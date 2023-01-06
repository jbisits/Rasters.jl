"""
    function convert_ocean_vars(raster::RasterStack, var_names::NamedTuple)
    function convert_ocean_vars(raster::Rasterseries, var_names::NamedTuple)
Convert ocean variables depth, practical salinity and potential temperature to pressure,
absolute salinity and conservative temperature. All conversions are done using the julia
implementation of TEOS-10 via [GibbsSeaWater](https://github.com/TEOS-10/GibbsSeaWater.jl).
The potential temperature and practical salinity variables in the `RasterStack` are replaced
by conservative temperature and absolute salinity. As pressure depends
on latitude and depth, it is added as a new variable --- that is, each longitude, latitude,
depth has a variable for pressure. A density variable is also computed which, by default, is
_in-situ_ density. Potential density at a reference pressure can be computed instead by
passing a the keyword argument `ref_pressure`.

The name of the variables for potential temperature and salinity
(either practical or absolute) must be passed in as a named tuple of the form
`(sp = :salt_name, pt = :potential_temp_name)` where `:potential_temp_name` and `:salt_name`
are the name of the potential temperature and salinity in the `Raster`.

**NOTE**
Currently this is only available for `RasterStacks`s or `RasterSeries`s with `X`, `Y`, `Z`
or `X`, `Y`, `Z`, `Ti` `dims`. If the `RasterStack/Series` has different dimensions e.g.
`X`, `Y`, `Ti` `convert_ocean_vars` will assume that the `Ti` dimension is depth (`Z`) and
compute incorrect variables. Though for pressure we require _at least_ `Y` (latitude) and
`Z` depth so it would not make sense to use this function if there is no depth `dim`.
If more methods are needed for different configurations of `dims` raise an issue on
GitHub.
"""
function convert_ocean_vars(raster::RasterStack, var_names::NamedTuple;
                            ref_pressure = nothing)

    Sₚ = raster[var_names.sp]
    θ = raster[var_names.pt]
    rs_dims = length(dims(Sₚ))==4 ? (dims(Sₚ, X), dims(Sₚ, Y), dims(Sₚ, Z), dims(Sₚ, Ti)) :
                                    (dims(Sₚ, X), dims(Sₚ, Y), dims(Sₚ, Z), nothing)
    p = depth_to_pressure(Sₚ, rs_dims)
    Sₐ = Sₚ_to_Sₐ(θ, Sₚ, p, rs_dims)
    Θ = θ_to_Θ(θ, Sₐ, rs_dims)
    converted_vars = (p = p, Sₐ = Sₐ, Θ = Θ)
    ρ = isnothing(ref_pressure) ? in_situ_density(Sₐ, Θ, p, rs_dims) :
                                  potential_density(Sₐ, Θ, ref_pressure, rs_dims)
    converted_vars = (p = p, Sₐ = Sₐ, Θ = Θ, ρ = ρ)

    return RasterStack(converted_vars, rs_dims)

end
function convert_ocean_vars(series::RasterSeries, var_names::NamedTuple; ref_pressure = nothing)

    rs_array = Array{RasterStack}(undef, length(series))
    for i ∈ eachindex(series)
        rs_array[i] = convert_ocean_vars(series[i], var_names; ref_pressure)
    end

    return RasterSeries(rs_array, dims(series, Ti))

end

"""
    function depth_to_pressure(raster::Raster)
Convert the depth dimension (`Z`) to pressure using `gsw_p_from_z`. Note that pressure
depends on depth and _latitude_ so the returned pressure is stored as a variable in the
resulting `Raster` rather than a vertical dimension.
"""
function depth_to_pressure(raster::Raster, rs_dims::Tuple)

    lons, lats, z, time = rs_dims
    p = similar(Array(raster))
    if isnothing(time)
        lats_array = repeat(Array(lats)'; outer = (length(lons), 1, length(z)))
        z_array = repeat(Array(z); outer = (1, length(lons), length(lats)))
        z_array = permutedims(z_array, (2, 3, 1))
        @. p = GibbsSeaWater.gsw_p_from_z(z_array, lats_array)
        rs_dims = (lons, lats, z)
    else
        lats_array = repeat(Array(lats)'; outer = (length(lons), 1, length(z), length(time)))
        z_array = repeat(Array(z); outer = (1, length(lons), length(lats), length(time)))
        z_array = permutedims(z_array, (2, 3, 1, 4))
        @. p = GibbsSeaWater.gsw_p_from_z(z_array, lats_array)
    end

    return Raster(p, rs_dims)

end

"""
    function Sₚ_to_Sₐ(raster::Raster, p::raster)
Convert a `Raster` of practical salinity (`Sₚ`) to absolute salinity (`Sₐ`) using
`gsw_sa_from_sp`.
"""
function Sₚ_to_Sₐ(θ::Raster, Sₚ::Raster, p::Raster, rs_dims::Tuple)

    lons, lats, z, time = rs_dims
    Sₐ = similar(Array(Sₚ), Union{Float64, Missing})
    if isnothing(time)
        for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

            Sₚ_profile = Sₚ[X(At(lon)), Y(At(lat))]
            θ_profile = θ[X(At(lon)), Y(At(lat))]
            find_nm = findall(.!ismissing.(θ_profile) .&& .!ismissing.(Sₚ_profile))
            p_profile = p[X(At(lon)), Y(At(lat)), Z(find_nm)]
            Sₐ[i, j, find_nm] = GibbsSeaWater.gsw_sa_from_sp.(Sₚ_profile[find_nm], p_profile, lon, lat)

        end
        rs_dims = (lons, lats, z)
    else
        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                Sₚ_profile = Sₚ[X(At(lon)), Y(At(lat)), Ti(t)]
                θ_profile = θ[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(θ_profile) .&& .!ismissing.(Sₚ_profile))
                p_profile = p[X(At(lon)), Y(At(lat)), Ti(t), Z(find_nm)]
                Sₐ[i, j, find_nm, t] = GibbsSeaWater.gsw_sa_from_sp.(Sₚ_profile[find_nm], p_profile, lon, lat)

            end
        end
    end

    return Raster(Sₐ, rs_dims)

end

"""
    function θ_to_Θ(raster::Raster, Sₐ::raster)
Convert a `Raster` of potential temperature (`θ`) to conservative temperature (`Θ`) using
`gsw_ct_from_pt`. This conversion depends on absolute salinity.
"""
function θ_to_Θ(θ::Raster, Sₐ::Raster, rs_dims::Tuple)

    lons, lats, z, time = rs_dims
    Θ = similar(Array(θ), Union{Float64, Missing})

    if isnothing(time)
        for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

            Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat))]
            θ_profile = θ[X(At(lon)), Y(At(lat))]
            find_nm = findall(.!ismissing.(θ_profile) .&& .!ismissing.(Sₐ_profile))
            Θ[i, j, find_nm] = GibbsSeaWater.gsw_ct_from_pt.(Sₐ_profile[find_nm], θ_profile[find_nm])

        end
        rs_dims = (lons, lats, z)
    else
        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat)), Ti(t)]
                θ_profile = θ[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(θ_profile) .&& .!ismissing.(Sₐ_profile))
                Θ[i, j, find_nm, t] = GibbsSeaWater.gsw_ct_from_pt.(Sₐ_profile[find_nm], θ_profile[find_nm])

            end
        end
    end

    return Raster(Θ, rs_dims)

end

"""
    function in_situ_density(Sₐ::Raster, Θ::Raster, p::Raster)
Compute in-situ density using `gsw_rho`.
"""
function in_situ_density(Sₐ::Raster, Θ::Raster, p::Raster, rs_dims::Tuple)

    lons, lats, z, time = rs_dims
    ρ = similar(Array(Θ))

    if isnothing(time)
        for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

            Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat))]
            Θ_profile = Θ[X(At(lon)), Y(At(lat))]
            find_nm = findall(.!ismissing.(Sₐ_profile) .&& .!ismissing.(Θ_profile))
            p_profile = p[X(At(lon)), Y(At(lat))]
            ρ[i, j, find_nm] = GibbsSeaWater.gsw_rho.(Sₐ_profile[find_nm], Θ_profile[find_nm], p_profile[find_nm])

        end
        rs_dims = (lons, lats, z)
    else
        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat)), Ti(t)]
                Θ_profile = Θ[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(Sₐ_profile) .&& .!ismissing.(Θ_profile))
                p_profile = p[X(At(lon)), Y(At(lat)), Z(find_nm), Ti(t)]
                ρ[i, j, find_nm, t] = GibbsSeaWater.gsw_rho.(Sₐ_profile[find_nm], Θ_profile[find_nm], p_profile)

            end
        end
    end

    return Raster(ρ, rs_dims)

end

"""
    function potential_density(Sₐ::Raster, Θ::Raster, p::Float64)
Compute potential density at reference pressure `p` using `gsw_rho`.
"""
function potential_density(Sₐ::Raster, Θ::Raster, p::Float64, rs_dims::Tuple)

    lons, lats, z, time = rs_dims
    ρ = similar(Array(Θ))

    if isnothing(time)
        for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

            Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat))]
            Θ_profile = Θ[X(At(lon)), Y(At(lat))]
            find_nm = findall(.!ismissing.(Sₐ_profile) .&& .!ismissing.(Θ_profile))
            ρ[i, j, find_nm] = GibbsSeaWater.gsw_rho.(Sₐ_profile[find_nm], Θ_profile[find_nm], p)

        end
        rs_dims = (lons, lats, z)
    else
        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat)), Ti(t)]
                Θ_profile = Θ[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(Sₐ_profile) .&& .!ismissing.(Θ_profile))
                ρ[i, j, find_nm, t] = GibbsSeaWater.gsw_rho.(Sₐ_profile[find_nm], Θ_profile[find_nm], p)

            end
        end
    end

    return Raster(ρ, rs_dims)

end
