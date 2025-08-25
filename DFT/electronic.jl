using Scanf
using JSON
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nbands = 18
nkpts = 352
bands = Matrix{Float64}(undef, nkpts, nbands)
Ef = 0.527119

eV_in_Ha = 27.211396641308;
kpts = 0:(nkpts-1)

function skiplines(io::IO, n::Int)
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io)
        read(io, Char) # actually skip the newline
    end
end

begin # load bands
    # Energies are in eV. Zero set to efermi
    f = open("./DFT/data/ebands.dat", "r")
    skiplines(f, nkpts + 8)
    for k ∈ 1:nkpts
        skipchars(isdigit, f) # skip index of k-point
        for b ∈ 1:nbands
            _, bands[k, b] = @scanf(f, " %f", Float64)
        end
        skipchars(isspace, f) # skip final newline
    end
    close(f)
end

function smoothen_bands(swaps...; plot=false)
    # 1, 86, 146, 206, 279, 352
    smooth_bands = deepcopy(bands)
    for (a, b, start) ∈ swaps
        tmp = smooth_bands[start:end, a]
        smooth_bands[start:end, a] = smooth_bands[start:end, b]
        smooth_bands[start:end, b] = tmp
        tmp = smooth_bands[start:end, a+1]
        smooth_bands[start:end, a+1] = smooth_bands[start:end, b+1]
        smooth_bands[start:end, b+1] = tmp
    end
    if plot
        fig = Figure(figure_padding=3)
        colsize!(fig.layout, 1, Relative(3 / 4))
        ax = Axis(
            fig[1, 1],
            # ylabel=L"E - E_F\ \mathrm{(eV)}",
            xticks=(
                [1, 86, 146, 206, 279, 352],
                ["Γ", "H", "N", "Γ", "P", "H"]
            )
        )
        xlims!(ax, 1, nkpts)
        ylims!(ax, -4, nothing)
        vlines!(ax, collect(map(t -> t[3] - 0.5, swaps)), linewidth=1, linestyle=:dash, color=:gray80)
        for b ∈ 9:nbands
            s = isodd(b) ? :solid : :dash
            @views lines!(ax, 1:nkpts, smooth_bands[:, b], linewidth=1, linestyle=s, label="$b")
        end
        Legend(fig[1, 2], ax)
        display(fig)
    end
    return smooth_bands
end

# (9,11,56), (10, 12, 56),
# (9,13,126), (10, 14, 126)

begin # export bands as json
    sbands = smoothen_bands(
        (15, 17, 27), (13, 17, 32), (11, 17, 32), (9, 17, 56), (9, 11, 126),
        (13, 15, 176), (13, 15, 236), (9, 17, 256), (9, 11, 256), (11, 17, 314)
    )
    @views d = map(b -> sbands[:, b], 1:nbands)
    open("DFT/export/bands.json", "w") do f
        JSON.print(f, d)
    end
end
# ang = map(2:nkpts) do k; rad2deg(atan(bands[k,10] - bands[k-1,10])); end

begin # load dos
    # energies are in Hartree
    f = open("./DFT/data/edos.dat", "r")
    skiplines(f, 7)
    _, Ef = @scanf(readline(f), "# Fermi energy :       %f", Float64)
    skiplines(f, 6)
    N = countlines(f)
    seekstart(f)
    skiplines(f, 14)
    dos = Matrix{Float64}(undef, N, 3)
    for (i, line) ∈ enumerate(eachline(f))
        _, E, d, dos[i, 3] = @scanf(line, " %f %f %f", Float64, Float64, Float64)
        dos[i, 1] = (E - Ef) * eV_in_Ha
        dos[i, 2] = d / eV_in_Ha
    end
    close(f)
end

begin # load projected dos
    # energies are in Hartree
    f = open("./DFT/data/edos-orbitals.dat", "r")
    skiplines(f, 7)
    _, Ef = @scanf(f, "# Fermi energy :       %f", Float64)
    skiplines(f, 4)
    _, N = @scanf(f, "# at %d energies", Int32)
    skiplines(f, 8)
    energies = Vector{Float64}(undef, N)
    orb_dos = Matrix{Float64}(undef, N, 5)
    for (i, line) ∈ enumerate(eachline(f))
        _, energies[i], orb_dos[i, 1], orb_dos[i, 2],
        orb_dos[i, 3], orb_dos[i, 4], orb_dos[i, 5] = @scanf(line,
            " %f %f %f %f %f %f %*f %*f %*f %*f %*f",
            Float64, Float64, Float64, Float64, Float64, Float64
        )
    end
    map!(E -> (E - Ef) * eV_in_Ha, energies, energies)
    orb_dos ./= eV_in_Ha
    close(f)
end

begin # export dos as json
    @views a = 1999 + findfirst(v -> v > 0, dos[2000:end, 2])
    @views b = a - 1 + findlast(v -> v > 0, dos[a:end, 2])
    sel = a:b

    @views d = Dict(
        "N" => length(sel),
        "energy" => energies[sel],
        "full" => dos[sel, 2],
        "s" => orb_dos[sel, 1],
        "p" => orb_dos[sel, 2],
        "d" => orb_dos[sel, 3],
        "f" => orb_dos[sel, 4],
        "g" => orb_dos[sel, 5],
    )
    open("DFT/export/edos.json", "w") do f
        JSON.print(f, d)
    end
end


dos_sel = (N>>1):N;

begin
    fig = Figure(size=(700, 400), figure_padding=3)
    colsize!(fig.layout, 1, Relative(3 / 4))
    ax_bands = Axis(
        fig[1, 1],
        ylabel=L"E - E_F\ \mathrm{(eV)}",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    )
    xlims!(ax_bands, 0, nkpts - 1)
    ax_dos = Axis(
        fig[1, 2],
        ylabel=L"E - E_F\ \mathrm{(eV)}",
        yaxisposition=:right,
        ylabelrotation=-0.5π
    )
    xlims!(ax_dos, 0, nothing)
    ylims!(ax_dos, -7, 13)
    linkyaxes!(ax_bands, ax_dos)

    for b ∈ 9:nbands
        @views y = bands[:, b]
        scatter!(ax_bands, kpts, y, color=:black, markersize=3)
    end
    @views lines!(ax_dos, dos[dos_sel, 2], dos[dos_sel, 1])

    display(fig)
end

save("dft/fig/ebands.svg", fig, pt_per_unit=2)

begin # plot DOS only
    fig = Figure(figure_padding=3)
    ax = Axis(
        fig[1, 1],
        xlabel=L"(E - E_F)/\mathrm{eV}",
        ylabel="Density of states (electrons/eV/cell)"
    )
    xlims!(ax, -7, 27)
    ylims!(ax, 0, nothing)
    @views lines!(ax, dos[dos_sel, 1], dos[dos_sel, 2])
    display(fig)
end


begin # find electrons per cell at T=0K
    @views i = findfirst(E -> E > 0, dos[:, 1])
    E1 = dos[i-1, 1]
    E2 = dos[i, 1]

    D1 = dos[i-1, 2]
    D2 = dos[i, 2]

    m = (D2 - D1) / (E2 - E1)
    q = D2 - m * E2
    println("DOS(Ef) = $q")

    I1 = dos[i-1, 3]
    I2 = dos[i, 3]

    m = (I2 - I1) / (E2 - E1)
    # I2 = m*E2 + q
    q = I2 - m * E2
    println("N_inf = $q")

    # remove two previous peaks
    @views i = findfirst(E -> E > -10, dos[:, 1])
    println("N_top = $(q-dos[i,3])")
end
