using Scanf
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

using GeometryBasics
using LinearAlgebra
using Interpolations
using Meshing

nbands = 18
nkpts = 7236
bands = Matrix{Float64}(undef, nkpts, nbands)

eV_in_Ha = 27.211396641308;
joules_in_Hartree = 4.359748199E-18; # J/Ha
ħ = 1.054571817e-34; # J*s
m_e = 9.1093837139e-31; # Kg

function skiplines(io::IO, n::Int) 
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io);
        read(io, Char); # actually skip the newline
    end
end

begin
    f = open("./DFT/data/fsurf.bxsf", "r")
    skiplines(f, 9);
    _, Ef = @scanf(f, " Fermi Energy: %f", Float64);
    skiplines(f, 6)
    _, nbands, nx, ny, nz = @scanf(f, " %d %d %d %d", Int32, Int32, Int32, Int32);
    nkpts = nx * ny * nz;
    bands = Vector{Array{Float64}}(undef, nbands);
    skiplines(f, 2);
    basis = Matrix{Float64}(undef, 3, 3);
    for i ∈ 1:3
        _, basis[1,i], basis[2,i], basis[3,i] = @scanf(f," %f %f %f", Float64, Float64, Float64);
    end
    sz = norm(basis[:,1]);
    for b ∈ 1:nbands
        skiplines(f,2);
        bands[b] = data = Array{Float64}(undef, nx, ny, nz);
        for i ∈ 1:nx
            for j ∈ 1:ny
                for k ∈ 1:nz
                    _, data[i,j,k] = @scanf(f, " %f", Float64)
                end
            end
        end
    end
    close(f); 
end

begin
    # XCrysDen grid are general, not periodic
    # => must exclude last k-points in each direction
    cnt = 0
    for b ∈ 11:14
        data = bands[b];
        @views cnt += @views count(E -> E > Ef, data[1:end-1, 1:end-1, 1:end-1]) 
    end
    periodic_nkpts = (nx-1)*(ny-1)*(nz-1)
    println("n_h = $(cnt/periodic_nkpts)")
end

E = bands[11];

begin
    dx = 1/(nx-1.0);
    dy = 1/(ny-1.0);
    dz = 1/(nz-1.0);
    X = LinRange(0,1,nx);
    Y = LinRange(0,1,ny);
    Z = LinRange(0,1,nz);
    @views itp = interpolate(E, BSpline(Cubic(Periodic(OnGrid()))));
    itp = scale(itp, X, Y, Z);
    # itp = extrapolate(itp, Periodic(OnCell()));
end;

begin
    inv_basis = inv(basis);

    X = LinRange(0,1,nx)
    Y = LinRange(0,1,ny)
    Z = LinRange(0,1,nz)

    weighted_iM = [zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3)];
    weighted_vf = zeros(4);
    total_areas = zeros(4);

    GLp = [4 1 1; 1 4 1; 1 1 4]/6;
    GLw = [1, 1, 1]/3;
    nGLp = size(GLp)[2];

    for i ∈ 1:4
        E = bands[i+10];
        vertices, faces = isosurface(E, MarchingCubes(iso=Ef), X, Y, Z);
        vertices_as_mat = reshape(reinterpret(Float64, vertices), (3,:));
        @views itp = interpolate(E, BSpline(Cubic(Periodic(OnCell()))));
        itp = scale(itp, X, Y, Z);
        for (a, b, c) ∈ faces
            v1 = basis*collect(vertices[b] .- vertices[a]);
            v2 = basis*collect(vertices[c] .- vertices[a]);
            triangle_area = 0.5 * norm(cross(v1, v2))

            @views integration_points = vertices_as_mat[:, [a,b,c]] * GLp;
            vel = mapslices(pos -> norm(inv_basis * Interpolations.gradient(itp, pos...)), integration_points, dims=1) * GLw;

            # baricenter = (vertices[a] .+ vertices[b] .+ vertices[c])./3;
            # vel = inv_basis * Interpolations.gradient(itp, baricenter...);
            # iM = inv_basis * inv_basis * Interpolations.hessian(itp, baricenter...);
            # weighted_iM[i] += triangle_area * iM;
            weighted_vf[i] += triangle_area * norm(vel);
            total_areas[i] += triangle_area;
        end
    end

    weighted_iM .*= (joules_in_Hartree * m_e) / (ħ^2 * 1e20); # 1/mₑ
    weighted_vf .*= 1e-10 * joules_in_Hartree / ħ; # m/s

    iM = sum(weighted_iM, dims=1)[1] / sum(total_areas);
    vf = sum(weighted_vf) / sum(total_areas);
end


begin
    X = LinRange(0,1,nx)
    Y = LinRange(0,1,ny)
    Z = LinRange(0,1,nz)
    @views vertices, faces = isosurface(E, MarchingCubes(iso=Ef), X, Y, Z);
    vertices_as_mat = reshape(reinterpret(Float64, vertices), (3,:))'
    faces_as_mat = reshape(reinterpret(Int, faces), (3,:))'

    upper = sum(basis, dims=2);
    fig = Figure(size=(600,600))
    ax = Axis3(
        fig[1,1], 
        aspect=(1,1,1),
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
        0 0 1;
    ] * basis;
    lines!(ax, path1, color=:black);

    path2 = [
        1 0 0;
        1 0 1;
        1 1 0;
        1 1 1;
        0 1 0;
        0 1 1;
    ]*basis;
    linesegments!(ax, path2, color=:black)
    mesh!(ax, vertices_as_mat, faces_as_mat)

    display(fig)
end