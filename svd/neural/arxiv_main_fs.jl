using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

include("load_files.jl")
include("matchfs_nrep.jl")
include("matchfs_noreplacement.jl")
include("explainedVar.jl")
include("runmatchfs.jl")
include("figure/plt_fs_pcor.jl")
include("figure/plt_fs_raster.jl")

temporalAvg = "noavg" # "avg", "noavg"

dirsim_bal = "/data/kimchm/data/dale/janelia/trained/basic/simffwd/balanced/" * temporalAvg * "/"
dirsim = "/data/kimchm/data/dale/janelia/trained/basic/sim/" * temporalAvg * "/"
dirtarget = "/data/kimchm/data/dale/janelia/s1alm/target/"
dirutargPyr = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
dirutargFS = "/data/kimchm/data/dale/janelia/s1alm/target/utargFS/"
dirfig = "/data/kimchm/data/dale/janelia/figure/lickright/" * temporalAvg * "/"

if ~ispath(dirfig)
    mkpath(dirfig)
end

# # load targets (Pyr)
# utargPyr = transpose(load(dirutargPyr * "utarg1Hz.jld", "utarg"))[:,:]

# load FS
rtargFS = transpose(load(dirtarget * "movingrate_FS.jld", "FS"))[1:100,:]
utargFS = load(dirutargFS * "utargFS.jld", "utargFS")[1:100,:]

# remove neurons with zero firing rate
idx = (sum(rtargFS, dims=1) .!= 0)[:]
rtargFS = rtargFS[:,idx]


# load usum (trained network)
nrep = 20
Ne = 2500
Ncells = 5000
navg_list = [50, 100, 200, 400]
navg = navg_list[4]

navg = 50
expvar_bal50, expvar_rate50, pcor_bal50, pcor_rate50, cells_rate50 = runmatchfs(dirsim, dirsim_bal, navg, utargFS, rtargFS)

navg = 400
expvar_bal400, expvar_rate400, pcor_bal400, pcor_rate400, cells_rate400 = runmatchfs(dirsim, dirsim_bal, navg, utargFS, rtargFS)

plt_fs_pcor(dirfig, expvar_bal400, expvar_rate400, pcor_bal400, pcor_rate400, 
expvar_bal50, expvar_rate50, pcor_bal50, pcor_rate50)

# plt_fs_raster(dirfig, dirsim, rtargFS, pcor_rate400, cells_rate400)

