using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

include("load_files.jl")
include("reconstructTarget.jl")
include("runinitial.jl")
include("runtest.jl")
include("runsim.jl")
include("runperformance.jl")
include("funMovAvg.jl")
include("funRollingAvg.jl")
include("funSample.jl")
include("loglogfit.jl")
include("meanRates.jl")
include("psth.jl")

L_list = [2.0, 4.0, 6.0, 8.0, 10.0]
ll = 1 # L = 2.0

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
nn = 5

# dirData = "/data/kimchm/data/dale/janelia/trained/L$(ll)/taup$(nn)/"
dirdata = "/data/kimchm/data/dale/janelia/trained/basic/"
dirfig = "/data/kimchm/data/dale/janelia/figure/"
dirutarg = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
dirfig = "/data/kimchm/data/dale/janelia/figure/"

if ~ispath(dirfig)
    mkpath(dirfig)
end

println("Loading files")
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells = load_files(dirdata)

# run a trained network
rollingavg = true
xsum, times, ns, 
vsum, vbalE, vbalI, vplasE, vplasI = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,rollingavg,matchedCells)


figure(figsize=(6,5))
ncell = 250
axvspan(p.stim_on, p.stim_off, color="dodgerblue", alpha=0.3)
for ci = 1:ncell
    plot(times[ci,1:ns[ci]], ci*ones(ns[ci]), color="red", marker="o", ms=1, linestyle="")
end
for ci = p.Ne+1:p.Ne+ncell
    plot(times[ci,1:ns[ci]], (ci-p.Ne+ncell)*ones(ns[ci]), color="blue", marker="o", ms=1, linestyle="")
end
xlim([0,p.train_time])
ylim([0,2*ncell])
xlabel("time (ms)", fontsize=12)
ylabel("neuron", fontsize=12)
tight_layout()

savefig(dirfig * "trained_raster.png", dpi=300)


timev = p.dt*collect(1:p.Nsteps)

for nid = 1:10
    figure(figsize=(3.5,3))
    axvspan(p.stim_on, p.stim_off, color="dodgerblue", alpha=0.3)
    plot(timev, vsum[:,nid][:], c="black", linewidth=0.5, alpha=1)
    plot(timev, vbalE[:,nid][:], c="red", linewidth=0.5, alpha=0.3)
    plot(timev, vbalI[:,nid][:], c="blue", linewidth=0.5, alpha=0.3)
    plot(timev, vplasE[:,nid][:], c="red", linewidth=1.0, alpha=1)
    plot(timev, vplasI[:,nid][:], c="blue", linewidth=1.0, alpha=1)
    plot(timev, vplasE[:,nid][:] + vplasI[:,nid][:], c="darkorange", linewidth=1.0, alpha=1)
    plot(timev, ones(length(timev)), c="magenta", linewidth=1.5, linestyle="--", alpha=1)
    xlim([0,3000])
    ylim([-4,4])
    xlabel("time (ms)", fontsize=12)
    ylabel("synaptic input", fontsize=12)
    tight_layout()

    savefig(dirfig * "trained_syn$(nid).png", dpi=300)
    close()
end
