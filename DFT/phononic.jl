using Scanf
using JSON
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nkpts = 352
Ef = 0.527119

eV_in_Ha = 27.211396641308;
THz_in_meV = 0.2417990504024;
THz_in_Ha = THz_in_meV * 1e3 * eV_in_Ha;
K_in_meV = 11.60452500617; # K/meV
K_in_THz = K_in_meV / THz_in_meV;


kpts = 0:(nkpts-1);

function skiplines(io::IO, n::Int)
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io)
        read(io, Char) # actually skip the newline
    end
end

function load_bands(filepath::String)
    # Frequencies are in meV
    bands = Matrix{Float64}(undef, nkpts, 3)
    f = open(filepath, "r")
    skiplines(f, nkpts + 6)
    for k ∈ 1:nkpts
        c, bands[k, 1], bands[k, 2], bands[k, 3] =
            @scanf(f, "%*d %f %f %f \n", Float64, Float64, Float64)
        if c != 3
            println("Read only $c at line $k")
        end
    end
    close(f)
    bands *= THz_in_meV
    return bands
end


function load_dos(filepath::String)
    # energies are in Hartree
    f = open(filepath, "r")
    skiplines(f, 8)
    N = countlines(f)
    seekstart(f)
    dos = Matrix{Float64}(undef, N, 3)
    skiplines(f, 8)
    for (i, line) ∈ enumerate(eachline(f))
        _, E, d, dos[i, 3] = @scanf(line, " %f %f %f %*f %*f", Float64, Float64, Float64)
        dos[i, 1] = E * THz_in_Ha
        dos[i, 2] = d / eV_in_Ha * 1e-3
    end
    close(f)
    return dos
end

begin # load eph bands
    # Frequencies are in meV
    bands_ep = Matrix{Float64}(undef, nkpts, 3)
    f = open("./DFT/data/phbands_ep.agr", "r")
    skiplines(f, nkpts + 36)
    for b ∈ 1:3
        skiplines(f, 3)
        for k ∈ 1:nkpts
            _, bands_ep[k, b] = @scanf(f, "%*d %f \n", Float64)
        end
    end
    close(f)
end

bands = load_bands("./DFT/data/phbands.dat")
bands12 = load_bands("./DFT/data/phbands12.dat")

dos = load_dos("./DFT/data/phdos.dat")
dos12 = load_dos("./DFT/data/phdos12.dat")

begin
    yticks = 0:5:25
    fig = Figure(size=(700, 400), figure_padding=3)
    colsize!(fig.layout, 1, Relative(3 / 4))
    ax_bands = Axis(
        fig[1, 1],
        yticks=yticks,
        ylabel=L"E\quad\mathrm{(meV)}",
        xlabel=L"k\quad(2\pi/a)",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    )
    xlims!(ax_bands, 0, nkpts - 1)
    ax_dos = Axis(
        fig[1, 2],
        yticks=yticks,
        ylabel=L"E\quad\mathrm{(meV)}",
        yaxisposition=:right,
        ylabelrotation=-0.5π,
        xticks=(0:0.1:0.3, ["0", "0.1", "0.2", "0.3"]),
        xlabel=L"\mathrm{DOS}\ \mathrm{(meV^{-1})}"
    )
    xlims!(ax_dos, 0, nothing)
    ylims!(ax_dos, 0, 28)
    linkyaxes!(ax_bands, ax_dos)


    for b ∈ 1:3
        @views y = bands[:, b]
        lines!(ax_bands, kpts, y, color=:black)
    end
    @views lines!(ax_dos, dos[:, 2], dos[:, 1], color=:black)

    display(fig)
end

save("dft/fig/phbands.svg", fig, pt_per_unit=2)


begin
    @views d = Dict(
        "k18_q6" => Dict(
            "bands" => [
                bands[:, 1],
                bands[:, 2],
                bands[:, 3],
            ],
            "omega" => dos[:, 1],
            "dos" => dos[:, 2]
        ),
        "k12_q12" => Dict(
            "bands" => [
                bands12[:, 1],
                bands12[:, 2],
                bands12[:, 3],
            ],
            "omega" => dos12[:, 1],
            "dos" => dos12[:, 2]
        ),
    )
    open("DFT/export/phonons.json", "w") do f
        JSON.print(f, d)
    end
end


begin #experimental
    # Gamma-P-H
    zzz = [
        0.031 0.83 NaN;
        0.048 1.33 NaN;
        0.067 NaN 0.83;
        0.070 1.83 NaN;
        0.10 2.46 NaN;
        0.117 NaN 1.33;
        0.15 3.51 1.71;
        0.166 NaN 1.83;
        0.20 4.34 2.11;
        0.25 5.04 2.62;
        0.275 5.42 NaN;
        0.30 5.69 3.23;
        0.32 5.90 NaN;
        0.34 6.00 NaN;
        0.35 5.91 3.88;
        0.36 5.97 NaN;
        0.375 5.76 NaN;
        0.40 5.54 4.39;
        0.425 5.32 NaN;
        0.44 5.11 NaN;
        0.45 5.17 4.82;
        0.46 5.11 NaN;
        0.475 5.12 NaN;
        0.48 5.10 NaN;
        0.50 5.03 5.03;
        0.52 4.86 NaN;
        0.525 4.93 NaN;
        0.54 4.69 NaN;
        0.55 4.69 5.07;
        0.56 4.58 NaN;
        0.575 4.46 NaN;
        0.58 4.35 NaN;
        0.60 4.20 5.10;
        0.62 3.94 NaN;
        0.625 NaN 5.05;
        0.64 3.76 NaN;
        0.65 3.68 5.11;
        0.66 3.62 NaN;
        0.675 3.51 NaN;
        0.68 3.50 NaN;
        0.70 3.50 5.20;
        0.72 3.60 NaN;
        0.725 3.66 NaN;
        0.74 3.72 NaN;
        0.75 3.80 5.45;
        0.76 3.79 NaN;
        0.775 3.95 5.60;
        0.78 3.94 NaN;
        0.80 4.14 5.80;
        0.82 4.27 NaN;
        0.84 4.77 NaN;
        0.85 4.99 6.12;
        0.86 5.05 NaN;
        0.88 5.42 NaN;
        0.90 5.80 6.30;
        0.95 6.29 6.41;
        1.00 6.49 6.49
    ]

    # Gamma-N
    zz0 = [
        0.05 1.13 NaN NaN;
        0.075 1.62 NaN 0.78;
        0.10 2.12 0.73 1.04;
        0.125 2.55 0.93 1.31;
        0.15 3.00 1.19 1.57;
        0.20 3.83 1.76 2.15;
        0.25 4.32 2.39 2.71;
        0.30 4.86 3.17 3.24;
        0.338 5.05 NaN NaN;
        0.35 5.12 3.82 3.57;
        0.40 5.26 4.53 3.80;
        0.45 5.47 4.92 3.88;
        0.50 5.66 5.07 3.93
    ]

    # Gamma-H
    z00 = [
        0.031 0.50 NaN;
        0.044 0.75 NaN;
        0.059 1.00 NaN;
        0.074 1.25 NaN;
        0.088 NaN 0.50;
        0.091 1.50 NaN;
        0.10 1.63 NaN;
        0.143 NaN 0.73;
        0.145 NaN 0.75;
        0.15 2.39 NaN;
        0.193 NaN 0.93;
        0.20 3.10 NaN;
        0.205 NaN 1.00;
        0.213 NaN 1.05;
        0.233 NaN 1.14;
        0.243 NaN 1.19;
        0.25 3.73 NaN;
        0.25 NaN 1.25;
        0.253 NaN 1.25;
        0.273 NaN 1.37;
        0.291 NaN 1.50;
        0.293 NaN 1.51;
        0.30 4.35 1.54;
        0.313 NaN 1.63;
        0.333 NaN 1.78;
        0.343 NaN 1.87;
        0.35 4.82 NaN;
        0.353 NaN 1.97;
        0.375 5.05 NaN;
        0.393 NaN 2.28;
        0.40 5.29 2.27;
        0.425 5.46 NaN;
        0.443 NaN 2.82;
        0.45 5.60 NaN;
        0.455 NaN 3.00;
        0.46 5.63 NaN;
        0.48 5.72 NaN;
        0.493 NaN 3.38;
        0.50 5.71 3.55;
        0.505 NaN 3.50;
        0.52 5.81 NaN;
        0.54 5.79 NaN;
        0.545 NaN 4.00;
        0.55 5.82 NaN;
        0.58 5.75 NaN;
        0.60 5.71 4.57;
        0.65 5.55 5.04;
        0.70 5.50 NaN;
        0.725 5.50 NaN;
        0.75 5.59 5.90;
        0.775 5.68 NaN;
        0.80 5.70 6.23;
        0.84 5.92 NaN;
        0.85 6.03 NaN;
        0.86 6.10 NaN;
        0.88 6.15 NaN;
        0.90 6.29 6.42;
        0.92 6.32 NaN;
        0.94 6.37 NaN;
        0.95 6.40 6.47;
        0.97 6.43 NaN;
        1.00 6.49 6.49
    ]

    # 0  85  145  205  278  351
    # Γ  H   N    Γ    P    H
    zzz[:, 1] = (351 - 205) * zzz[:, 1] .+ 205
    zz0[:, 1] = 2*(145 - 205) * zz0[:, 1] .+ 205
    z00[:, 1] = 85 * z00[:, 1]

    pippo(A, col) = begin
        @views sel = (!).(isnan.(A[:, col]))
        return Dict(
            "k" => A[sel, 1],
            "f" => A[sel, col]
        )
    end

    d = [
        pippo(z00, 2),
        pippo(z00, 3),
        pippo(zz0, 2),
        pippo(zz0, 3),
        pippo(zz0, 4),
        pippo(zzz, 2),
        pippo(zzz, 3),
    ]
    open("DFT/export/phbands_experimental.json", "w") do f
        JSON.print(f, d)
    end
end

begin # Debye temperature
    @views E = dos12[:, 1]
    @views xx = E .^ 2
    a = findfirst(v -> v >= 0, E)
    G = map(50:500) do Δ
        sel = a:(a+Δ)
        @views y = dos12[sel, 2]
        @views A = (xx[sel]\y)[1]
        @views s = xx[sel] * A
        y_avg = sum(y) / length(y)
        R = sum((s .- y) .^ 2) / sum((y .- y_avg) .^ 2)
        return (A, R)
    end

    R, i = findmin(g -> g[2], G)
    println("Min residue ($(R)) with window of $(i + 49)")

    A = G[i][1] # 1/(meV * THz²)
    b = a + 49 + i
    sel = a:b

    # x^3 * A/3 = 3
    # x = cbrt(9/A)
    theta = cbrt(9 / A * THz_in_meV)
    sx = LinRange(0, theta, 50)
    sy = (sx .^ 2) * A


    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="ω/THz", ylabel="DOS (states/meV)")
    xlims!(ax, 0, nothing)
    @views lines!(ax, E, dos12[:, 2])
    lines!(ax, sx, sy)
    vlines!(ax, [E[b], theta], linestyle=:dash, color=:black)
    display(fig)

    println("Debye θ = $(theta)THz = $(theta*K_in_THz)K")
    println("A = $A /(meV * THz^2)")
end
