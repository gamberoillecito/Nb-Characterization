using Scanf
using CairoMakie
using LsqFit
using NonlinearSolve

begin
    f = open("data/heating_8_to_296.dat")
    N = countlines(f);
    seekstart(f);

    idxs = Vector{Int32}(undef, N);

    front_time = Vector{Float32}(undef, N);         # s
    front_temp = Vector{Float32}(undef, N);         # K
    front_temp_fluct = Vector{Float32}(undef, N);   # K
    front_2WR = Vector{Float32}(undef, N);          # Ohm

    back_time = Vector{Float32}(undef, N);         # s
    back_temp = Vector{Float32}(undef, N);         # K
    back_temp_fluct = Vector{Float32}(undef, N);   # K
    back_2WR = Vector{Float32}(undef, N);          # Ohm
    
    # parse each line of the file
    for i ∈ 1:N
        r, idxs[i], 
        front_time[i], front_temp[i], front_temp_fluct[i], front_2WR[i],
        back_time[i], back_temp[i], back_temp_fluct[i], back_2WR[i] =
        @scanf(f, " %*f %f %f %*f %*f %f %f %*f %*f %*f %*f %*f %f %f %*f %*f %f %f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f", 
        Int32, Float32, Float32, Float32, Float32, Float32, Float32, Float32, Float32);
        if r != 9; println("line $i, read $r/9"); end;
    end
    close(f)
end

cutoff = 35e-3; # 35 mK
# which are the points for which ΔTf < cutoff and ΔTb < cutoff?
# TODO: filter Tf and Tb separately
sel = [(
    front_temp[i] > 8 && back_temp[i] > 8
    && front_temp_fluct[i] <= cutoff && back_temp_fluct[i] <= cutoff
    # && abs(front_temp[i] - back_temp[i]) < 0.025
    ) for i ∈ 1:N
];

begin 
    # after selecting only the "good" points
    # we need to sort them by temperature
    # since they likely won't be sorted

    M = sum(sel);
    idx = Vector{Int64}(undef, M);
    
    Tf = front_temp[sel];
    Rf = front_2WR[sel];
    sortperm!(idx, Tf);
    Tf = Tf[idx];
    Rf = Rf[idx];
    

    Tb = back_temp[sel];
    Rb = back_2WR[sel];
    sortperm!(idx, Tb);
    Tb = Tb[idx];
    Rb = Rb[idx];
end;

begin
    # define equispaced sample of temperature and
    # average of points inside the bins centered around the sample T
    sample_T = 8:0.5:296;
    Ns = length(sample_T);
    clean_Rf = zeros(Ns);
    clean_Rb = zeros(Ns);
    jf = 1;
    jb = 1;
    for i ∈ 1:Ns
        lim = sample_T[i] + 0.25;
        count = 0;
        while jf < length(Rf) && Tf[jf] < lim
            count += 1;
            clean_Rf[i] += Rf[jf];
            jf += 1;
        end
        clean_Rf[i] /= count;
        count = 0;
        while jb < length(Rb) && Tb[jb] < lim
            count += 1;
            clean_Rb[i] += Rb[jb];
            jb += 1;
        end
        clean_Rb[i] /= count;
    end
end

begin
    # solve Van der Pau pairwise
    vdp(x, p) = exp(-π*p[1]/x) + exp(-π*p[2]/x) - 1;
    interp_R = map(zip(clean_Rf, clean_Rb)) do (f,b) 
        if isnan(f) || isnan(b)
            return NaN;
        else
            prob = IntervalNonlinearProblem(vdp, (1e-3, 160.), [f, b]);
            sol = solve(prob);
            return sol.u;
        end
    end
end


begin
    fig = Figure()
    ax = Axis(fig[1,1]);
    scatter!(ax, Tf, Rf, markersize=3, label="Front");
    lines!(ax, sample_T, clean_Rf, linestyle=:dash, color=:black);
    scatter!(ax, Tb, Rb, markersize=3, label="Back");
    lines!(ax, sample_T, clean_Rb, linestyle=:dash, color=:black);
    lines!(ax, sample_T, interp_R, linestyle=:dash, color=:black);
    axislegend(ax, position=:lt)
    display(fig)
end

begin
    # plot fluctuations
    fig = Figure(size=(400, 800))
    ax = Axis(fig[1,1], yticks=0:0.1:1);
    scatter!(ax, front_temp, front_temp_fluct, markersize=3, label="Front");
    scatter!(ax, back_temp, back_temp_fluct, markersize=3, label="Back");
    hlines!(ax, [cutoff], linestyle=:dash, color=:black);
    axislegend(ax, position=:lt)
    display(fig)
end