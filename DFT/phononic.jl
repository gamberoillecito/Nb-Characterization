using Scanf
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nkpts = 279
bands = Matrix{Float64}(undef, nkpts, 3)
Ef = 0.527119

eV_in_Ha = 27.211396641308;
kpts = 0:(nkpts-1)

function skiplines(io::IO, n::Int) 
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io);
        read(io, Char); # actually skip the newline
    end
end

begin # load bands
    # Frequencies are in meV
    f = open("./DFT/data/phbands.dat", "r");
    skiplines(f, nkpts + 6);
    for k ∈ 1:nkpts
        _, bands[k,1], bands[k,2], bands[k,3] = 
        @scanf(f, "%*d %f %f %f \n", Float64, Float64, Float64);
    end
    close(f)
end


begin # load dos
    # energies are in Hartree
    f = open("./DFT/data/phdos.dat", "r");
    skiplines(f, 8)
    N = countlines(f)
    seekstart(f);
    dos = Matrix{Float64}(undef, N, 3);
    skiplines(f, 8)
    for (i,line) ∈ enumerate(eachline(f))
        _, E, d, dos[i,3] = @scanf(line, " %f %f %f %*f %*f", Float64, Float64, Float64);
        dos[i,1] = E * eV_in_Ha * 1e3;
        dos[i,2] = d / eV_in_Ha * 1e-3;
    end
    close(f)
end

begin
    yticks = 0:5:25
    fig = Figure(size=(700, 400), figure_padding=3)
    colsize!(fig.layout, 1, Relative(3/4))
    ax_bands = Axis(
        fig[1,1],
        yticks=yticks,
        ylabel=L"E\quad\mathrm{(meV)}",
        xlabel=L"k\quad(2\pi/a)",
        xticks=(
            [0, 85, 145, 205, 278, 351],
            ["Γ", "H", "N", "Γ", "P", "H"]
        )
    );
    xlims!(ax_bands,0, nkpts-1);
    ax_dos = Axis(
        fig[1,2],
        yticks=yticks,
        ylabel=L"E\quad\mathrm{(meV)}",
        yaxisposition=:right,
        ylabelrotation=-0.5π,
        xticks=(0:0.1:0.3, ["0","0.1","0.2", "0.3"]),
        xlabel=L"\mathrm{DOS}\ \mathrm{(meV^{-1})}"
    );
    xlims!(ax_dos, 0, nothing);
    ylims!(ax_dos, 0, 28)
    linkyaxes!(ax_bands, ax_dos);


    for b ∈ 1:3
        @views y = bands[:,b];
        lines!(ax_bands, kpts, y, color=:black)
    end
    @views lines!(ax_dos, dos[:,2], dos[:,1], color=:black);

    display(fig)
end

save("dft/fig/phbands.svg", fig, pt_per_unit=2)
