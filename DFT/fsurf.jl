using Scanf
using CairoMakie
update_theme!(theme_latexfonts());
update_theme!(Theme(fontsize=22));

nbands = 18
nkpts = 7236
bands = Matrix{Float64}(undef, nkpts, nbands)
Ef = 0.525983014

eV_in_Ha = 27.211396641308;

function skiplines(io::IO, n::Int) 
    for _ ∈ 1:n
        skipchars(c -> c != '\n', io);
        read(io, Char); # actually skip the newline
    end
end

begin
    bands = Vector{Matrix{Float64}}(undef, nbands);

    f = open("./DFT/data/fsurf.bxsf", "r")
    skiplines(f, 19);

    for b ∈ 1:nbands
        data = Matrix{Float64}(undef, nkpts, 7);
        bands[b] = data;
        skiplines(f, 3);
        for i ∈ 1:nkpts
            _, data[i,1], data[i,2], data[i,3], data[i,4], data[i,5], data[i,6], data[i,7] = @scanf(f, " %f %f %f %f %f %f %f", Float64, Float64, Float64, Float64, Float64, Float64, Float64);
        end
    end
    close(f); 
end

begin
    cnt = 0
    for b ∈ 11:14
        data = bands[b];
        cnt += @views count(E -> E < Ef, data[:,4]) 
    end
    println("n_h = $(cnt/nkpts)")
end