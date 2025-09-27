using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

using Scanf
using GeometryBasics
using LinearAlgebra
using Interpolations
using Meshing

nbands = 18
eV_in_Ha = 27.211396641308;
joules_in_Hartree = 4.359748199E-18; # J/Ha
ħ = 1.054571817e-34; # J*s
planck = 6.62607015e-34;
m_e = 9.1093837139e-31; # Kg

function skiplines(io::IO, n::Int)
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io)
        read(io, Char) # actually skip the newline
    end
end

begin
    f = open("./DFT/data/fsurf72.bxsf", "r")
    skiplines(f, 9)
    _, Ef = @scanf(f, " Fermi Energy: %f", Float64)
    skiplines(f, 6)
    _, nbands, nx, ny, nz = @scanf(f, " %d %d %d %d", Int32, Int32, Int32, Int32)
    nkpts = nx * ny * nz
    bands = Vector{Array{Float64}}(undef, nbands)
    skiplines(f, 2)
    basis = Matrix{Float64}(undef, 3, 3)
    for i ∈ 1:3
        _, basis[1, i], basis[2, i], basis[3, i] = @scanf(f, " %f %f %f", Float64, Float64, Float64)
    end
    inv_basis = inv(basis)
    a = inv_basis[1, 2] * 2 # Ang
    sz = norm(basis[:, 1])
    for b ∈ 1:nbands
        skiplines(f, 2)
        bands[b] = data = Array{Float64}(undef, nx, ny, nz)
        for i ∈ 1:nx
            for j ∈ 1:ny
                for k ∈ 1:nz
                    _, data[i, j, k] = @scanf(f, " %f", Float64)
                end
            end
        end
    end
    close(f)
end

begin
    # XCrysDen grid are general, not periodic
    # => must exclude last k-points in each direction
    cnt = 0
    for b ∈ 11:14
        data = bands[b]
        @views cnt += @views count(E -> E > Ef, data[1:end-1, 1:end-1, 1:end-1])
    end
    periodic_nkpts = (nx - 1) * (ny - 1) * (nz - 1)
    println("n_h = $(cnt/periodic_nkpts)")
end


begin
    E = bands[11]
    X = LinRange(0, 1, nx)
    Y = LinRange(0, 1, ny)
    Z = LinRange(0, 1, nz)
    @views itp = interpolate(E, BSpline(Cubic(Periodic(OnGrid()))))
    itp = scale(itp, X, Y, Z)
    # itp = extrapolate(itp, Periodic(OnCell()));
end;

begin
    X = LinRange(0, 1, nx)
    Y = LinRange(0, 1, ny)
    Z = LinRange(0, 1, nz)

    min_v = fill(+Inf64, 2)
    max_v = fill(-Inf64, 2)

    weighted_v1 = zeros(2)
    weighted_v2 = zeros(2)

    weighted_im_cond = zeros(2)

    total_areas = zeros(2)

    GLp = [4 1 1; 1 4 1; 1 1 4] / 6
    GLw = [1, 1, 1] / 3
    nGLp = size(GLp)[2]

    for i ∈ 1:2
        E = bands[2*i+10]
        vertices, faces = isosurface(E, MarchingCubes(iso=Ef), X, Y, Z)
        vertices_as_mat = reshape(reinterpret(Float64, vertices), (3, :))
        @views itp = interpolate(E, BSpline(Cubic(Periodic(OnGrid()))))
        itp = scale(itp, X, Y, Z)

        for (a, b, c) ∈ faces
            v1 = basis * collect(vertices[b] .- vertices[a])
            v2 = basis * collect(vertices[c] .- vertices[a])
            triangle_area = 0.5 * norm(cross(v1, v2))

            @views integration_points = vertices_as_mat[:, [a, b, c]] * GLp

            vel1 = mapslices(integration_points, dims=1) do pos
                norm(inv_basis * Interpolations.gradient(itp, pos...))
            end

            vel2 = mapslices(integration_points, dims=1) do pos
                v = inv_basis * Interpolations.gradient(itp, pos...)
                return v ⋅ v
            end

            iM = map(eachcol(integration_points)) do pos
                inv_basis * inv_basis * Interpolations.hessian(itp, pos...)
            end

            im_cond = map(tr, iM)
            # m_therm = map(cbrt ∘ inv ∘ det, iM)

            min_v[i] = min(min_v[i], vel1...)
            max_v[i] = max(max_v[i], vel1...)

            total_areas[i] += triangle_area
            weighted_v1[i] += triangle_area * (GLw ⋅ vel1)
            weighted_v2[i] += triangle_area * (GLw ⋅ vel2)

            weighted_im_cond[i] += triangle_area * (GLw ⋅ im_cond)
        end
    end

    # We don't divide by hbar because our basis matrix
    # is already the reciprocal basis divided by 2π !!

    v_factor = 1e-10 * joules_in_Hartree / planck # Ha*Ang => m/s
    weighted_v1 .*= v_factor
    weighted_v2 .*= v_factor^2
    min_v .*= v_factor
    max_v .*= v_factor

    # division by 3 for the average of the traces
    m_factor = planck^2 * 1e20 / (joules_in_Hartree * m_e)     # 1/(Ha*Ang^2) => mₑ
    im_factor = (joules_in_Hartree * 1e-20 * m_e) / planck^2   # Ha*Ang^2 => 1/mₑ
    weighted_im_cond ./= m_factor * 3

    S = sum(total_areas)
    vel1 = sum(weighted_v1) / S
    vel2 = sum(weighted_v2) / S
    imeff = weighted_im_cond ./ total_areas
    m = 2 / sum(imeff)

    print("<|vf|> = $(vel1*1e-6) km/ms\n⎷(<|vf|²>) = $(sqrt(vel2)*1e-6) km/ms\n")
    print("m_c = 3/<Tr(inv(M))> = $(1 ./ imeff) => $m  (mₑ)\n")
end

begin
    X = LinRange(0, 1, nx)
    Y = LinRange(0, 1, ny)
    Z = LinRange(0, 1, nz)
    @views vertices, faces = isosurface(E, MarchingCubes(iso=Ef), X, Y, Z)
    vertices_as_mat = reshape(reinterpret(Float64, vertices), (3, :))'
    faces_as_mat = reshape(reinterpret(Int, faces), (3, :))'

    upper = sum(basis, dims=2)
    fig = Figure(size=(600, 600))
    ax = Axis3(
        fig[1, 1],
        aspect=(1, 1, 1),
        azimuth=1.2π,
        elevation=0.12π,
    )
    limits!(ax, 0, upper[1], 0, upper[2], 0, upper[3])
    path1 = [
        0 0 0;
        1 0 0;
        1 1 0;
        0 1 0;
        0 0 0;
        0 0 1;
        1 0 1;
        1 1 1;
        0 1 1;
        0 0 1
    ] * basis
    lines!(ax, path1, color=:black)

    path2 = [
        1 0 0;
        1 0 1;
        1 1 0;
        1 1 1;
        0 1 0;
        0 1 1
    ] * basis
    linesegments!(ax, path2, color=:black)
    mesh!(ax, vertices_as_mat, faces_as_mat)

    display(fig)
end


function band_speed(idx_band)
    X = LinRange(0, 1, nx)
    Y = LinRange(0, 1, ny)
    Z = LinRange(0, 1, nz)

    N = length(idx_band)
    weighted_v = zeros(N)

    for b ∈ 1:N
        E = bands[idx_band[b]]
        @views itp = interpolate(E, BSpline(Cubic(Periodic(OnGrid()))))
        itp = scale(itp, X, Y, Z)

        for i in 1:(nx-1)
            for j in 1:(ny-1)
                for k in 1:(nz-1)
                    v = inv_basis * Interpolations.gradient(itp, X[i], Y[j], Z[k])
                    weighted_v[b] += v ⋅ v
                end
            end
        end
    end

    # We don't divide by hbar because our basis matrix
    # is already the reciprocal basis divided by 2π !!

    v_factor = 1e-10 * joules_in_Hartree / planck # Ha*Ang => m/s
    weighted_v .*= v_factor^2
    weighted_v ./= (nx - 1) * (ny - 1) * (nz - 1)

    for b ∈ 1:N
        print("$(idx_band[b]): $(sqrt(weighted_v[b]))\n")
    end
end
