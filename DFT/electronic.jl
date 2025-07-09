using Scanf
using CSV
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nbands = 30
nkpts = 352
bands = Matrix{Float64}(undef, nkpts, nbands)
Ef = 0.527119

eV_in_Ha = 27.211396641308;

begin # load data
    # Energies are in eV. Zero set to efermi
    tobeskipped = nkpts + nbands + 33;
    f = open("./DFT/data/bands.agr", "r");
    for _ ∈ 1:tobeskipped; readline(f); end
    for b ∈ 1:nbands
        readline(f) # &
        readline(f) # @target G0.S0
        readline(f) # @type xy
        for k ∈ 1:nkpts
            _, bands[k,b] = @scanf(readline(f), "%*d %f", Float64);
        end
    end
    close(f)
end

begin # plot bands
    x = 0:(nkpts-1)
    fig = Figure(figure_padding=3)
    ax = Axis(
        fig[1,1],
        ylabel=L"E(k) - E_F",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    );
    xlims!(ax, 0, nkpts-1);
    for b ∈ 9:nbands
        @views y = bands[:,b];
        scatter!(ax, x, y, color=:black, markersize=3)
    end
    display(fig)
end


begin
    Ef = 0.52650654;
    f = open("./DFT/data/dos.dat", "r");
    for _ ∈ 1:14; readline(f); end;
    N = countlines(f)
    seekstart(f);
    dos = Matrix{Float64}(undef, N, 3);
    for _ ∈ 1:14; readline(f); end;
    for (i,line) ∈ enumerate(eachline(f))
        _, E, d, dos[i,3] = @scanf(line, " %f %f %f", Float64, Float64, Float64);
        dos[i,1] = (E - Ef) * eV_in_Ha;
        dos[i,2] = d / eV_in_Ha;
    end
    close(f)
end

begin # plot DOS
    sel = (N>>1):N;
    fig = Figure(figure_padding=3);
    ax = Axis(
        fig[1,1],
        xlabel=L"(E - E_F)/\mathrm{eV}",
        ylabel="Density of states (electrons/eV/cell)"
    );
    xlims!(ax, -7, 27);
    ylims!(ax, 0, nothing)
    @views lines!(ax, dos[sel,1], dos[sel,2]);
    display(fig);
end


begin
    fig = Figure(size=(700, 400), figure_padding=3)
    colsize!(fig.layout, 1, Relative(3/4))
    ax_bands = Axis(
        fig[1,1],
        ylabel=L"E - E_F\ \mathrm{(eV)}",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    );
    xlims!(ax_bands,0, nkpts-1);
    ax_dos = Axis(
        fig[1,2],
        ylabel=L"E - E_F\ \mathrm{(eV)}",
        yaxisposition=:right,
        ylabelrotation=-0.5π
    );
    xlims!(ax_dos, 0, nothing);
    ylims!(ax_dos, -7, 27)
    linkyaxes!(ax_bands, ax_dos);

    for b ∈ 9:nbands
        @views y = bands[:,b];
        scatter!(ax_bands, x, y, color=:black, markersize=3)
    end
    @views lines!(ax_dos, dos[sel,2], dos[sel,1]);

    display(fig)
end

save("dft/fig/bands.svg", fig, pt_per_unit=2)
