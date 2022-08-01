using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

include("load_files.jl")
include("runtest.jl")
# include("runperformance.jl")
# include("funMovAvg.jl")
# include("funRollingAvg.jl")
# include("funSample.jl")
# include("meanRates.jl")
# include("psth.jl")

L_list = [2.0, 4.0, 6.0, 8.0, 10.0]
ll = 1 # L = 2.0

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
nn = 5

# repi = parse(Int64,ARGS[1])
# repi = 1

# dirdata = "/data/kimchm/data/dale/janelia/trained/basic/"
# dirsim = "/data/kimchm/data/dale/janelia/trained/basic/sim/noavg/"
# dirutarg = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
# dirfig = "/data/kimchm/data/dale/janelia/figure/"

dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/L6/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/sim/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirutarg_pyr_lickright = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickright/"
dirutarg_pyr_lickleft = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickleft/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/"


if ~ispath(dirsim)
    mkpath(dirsim)
end

#----- lick right -----#
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells = load_files(dirdata, lickRL)
rollingavg = true
usum_R, times_R, ns_R = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,rollingavg)


#----- lick left -----#
println("Loading files")
lickRL = "left"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells = load_files(dirdata, lickRL)
rollingavg = true
usum_L, times_L, ns_L = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,rollingavg)

rate_R = zeros(p.Ncells)
rate_L = zeros(p.Ncells)
for ii = 1:p.Ncells
    rate_R[ii] = sum(times_R[ii, 1:ns_R[ii]] .> p.stim_off) / p.train_duration * 1000
    rate_L[ii] = sum(times_L[ii, 1:ns_L[ii]] .> p.stim_off) / p.train_duration * 1000
end

rdiff_pyr = rate_L[1:p.Ne] - rate_R[1:p.Ne]
rdiff_fs = rate_L[p.Ne+1:end] - rate_R[p.Ne+1:end]

pyr_sorted = rdiff_pyr[sortperm(rdiff_pyr)]
fs_sorted = rdiff_fs[sortperm(rdiff_fs)]


figure(figsize=(7,3.0))
subplot(121)
plot(pyr_sorted, color="red", marker="o", ms=2, linestyle="", label="Exc")
plot(collect(1:p.Ne), zeros(p.Ne), color="gray", linestyle="--")
legend(frameon=false, loc=2, fontsize=12)

subplot(122)
plot(fs_sorted, color="blue", marker="o", ms=2, linestyle="", label="Inh")
plot(collect(1:p.Ne), zeros(p.Ne), color="gray", linestyle="--")
legend(frameon=false, loc=2, fontsize=12)

tight_layout()

savefig(dirfig * "selectivity_pop.png", dpi=300)


