using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

include("load_files.jl")
include("runtestdrive.jl")
# include("runtest.jl")
# include("runperformance.jl")
# include("funMovAvg.jl")
# include("funRollingAvg.jl")
# include("funSample.jl")
# include("meanRates.jl")
# include("psth.jl")

L_list = [6.0, 12.0]
ll = 1 # L = 2.0

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
nn = 5

Lffwd_list = [100, 150, 200, 250, 300, 350, 400]
oo = 5

repi = parse(Int64,ARGS[1])
# repi = 1

# dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/saved/Lffwd3/lam2/"
dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/wffwd/Lffwd$(oo)/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirutarg_pyr_lickright = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickright/"
dirutarg_pyr_lickleft = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickleft/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


if ~ispath(dirsim)
    mkpath(dirsim)
    mkpath(dirfig)
end

#----- lick right -----#
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

# usum_R, times_R, ns_R = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
usum_R, times_R, ns_R = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1])

# fname_usum = dirsim * "usum_R_$(repi).jld"
fname_times = dirsim * "times_R_$(repi).jld"
fname_ns = dirsim * "ns_R_$(repi).jld"
# save(fname_usum,"usum", usum_R)
save(fname_times,"times", times_R)
save(fname_ns,"ns", ns_R)


#----- lick left -----#
println("Loading files")
lickRL = "left"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

# usum_L, times_L, ns_L = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
usum_L, times_L, ns_L = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2])

# fname_usum = dirsim * "usum_L_$(repi).jld"
fname_times = dirsim * "times_L_$(repi).jld"
fname_ns = dirsim * "ns_L_$(repi).jld"
# save(fname_usum,"usum", usum_L)
save(fname_times,"times", times_L)
save(fname_ns,"ns", ns_L)

