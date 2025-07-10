using Scanf
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nkpts = 279
bands = Matrix{Float64}(undef, nkpts, 3)
Ef = 0.527119

eV_in_Ha = 27.211396641308;
kpts = 0:(nkpts-1)


begin # load bands
    # Energies are in meV
    tobeskipped = nkpts + 8;
    f = open("./DFT/data/phbands.dat", "r");
    for _ ∈ 1:tobeskipped; readline(f); end
    for k ∈ 1:nkpts
        _, bands[k,1], bands[k,2], bands[k,3] = 
        @scanf(f, "%*d %f %f %f \n", Float64, Float64, Float64);
    end
    close(f)
end


begin # load dos
    # energies are in Hartree
    tobeskipped = 8
    f = open("./DFT/data/phdos.dat", "r");
    for _ ∈ 1:tobeskipped; readline(f); end;
    N = countlines(f)
    seekstart(f);
    dos = Matrix{Float64}(undef, N, 3);
    for _ ∈ 1:tobeskipped; readline(f); end;
    for (i,line) ∈ enumerate(eachline(f))
        _, E, d, dos[i,3] = @scanf(line, " %f %f %f %*f %*f", Float64, Float64, Float64);
        dos[i,1] = E * eV_in_Ha * 1e3;
        dos[i,2] = d / eV_in_Ha;
    end
    close(f)
end

begin
    fig = Figure(size=(700, 400), figure_padding=3)
    colsize!(fig.layout, 1, Relative(3/4))
    ax_bands = Axis(
        fig[1,1],
        ylabel=L"E/\mathrm{meV}",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    );
    xlims!(ax_bands,0, nkpts-1);
    ax_dos = Axis(
        fig[1,2],
        ylabel=L"E/\mathrm{meV}",
        yaxisposition=:right,
        ylabelrotation=-0.5π
    );
    xlims!(ax_dos, 0, nothing);
    # ylims!(ax_dos, -7, 13)
    linkyaxes!(ax_bands, ax_dos);

    for b ∈ 1:3
        @views y = bands[:,b];
        scatter!(ax_bands, kpts, y, color=:black, markersize=3)
    end
    @views lines!(ax_dos, dos[:,2], dos[:,1]);

    display(fig)
end

save("dft/fig/phbands.svg", fig, pt_per_unit=2)

begin # plot DOS only
    fig = Figure(figure_padding=3);
    ax = Axis(
        fig[1,1],
        xlabel=L"(E - E_F)/\mathrm{eV}",
        ylabel="Density of states (electrons/eV/cell)"
    );
    xlims!(ax, -7, 27);
    ylims!(ax, 0, nothing)
    @views lines!(ax, dos[dos_sel,1], dos[dos_sel,2]);
    display(fig);
end
