using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD2

# include("load_files.jl")
# include("runinitial.jl")
# include("runtest.jl")
# # include("runtest_scaling.jl")
# include("runperformance.jl")
# include("runperformance_nrep.jl")
# include("funMovAvg.jl")
# include("funRollingAvg.jl")
# include("funSample.jl")
# include("loglogfit.jl")
# include("meanRates.jl")
include("psth.jl")

L_list = [2.0, 4.0, 6.0, 8.0, 10.0]
ll = 1 # L = 2.0

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
nn = 5

navg_list = [50, 100, 200, 400]
# oo = parse(Int, ARGS[1])
oo = 4


# dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirsim = "/data/grabelmz/trainedNetwork/sim/"
dirsimTrialAvg = "/data/grabelmz/trainedNetwork/sim/trialAveraged/"


# average trained network activities
nrep = 1
navg = navg_list[oo]
# usum_navg = zeros(100, 5000, nrep)
# rate_navg = zeros(100, 5000, nrep)

# #---------- generate usum ----------#
# for repi = 1:nrep
#     for avgi = 1:navg
#         ii = (repi-1)*navg + avgi
#         fname_usum = dirsim * "usum_rep$(ii).jld"
#         utmp = load(fname_usum, "usum")
#         usum_navg[:,:,repi] .+= utmp/navg
#     end
# end
# fname_usum_navg = dirsim * "usum_navg$(navg).jld"
# save(fname_usum_navg,"usum", usum_navg)


#---------- generate psth ----------#
lickRL = "right"
rate_navg = psth(dirsim, navg, nrep, lickRL)
fname_rate_navg = dirsimTrialAvg * "rate_R_navg$(navg).jld2"
save(fname_rate_navg,"rate", rate_navg)


lickRL = "left"
rate_navg = psth(dirsim, navg, nrep, lickRL)
fname_rate_navg = dirsimTrialAvg * "rate_L_navg$(navg).jld2"
save(fname_rate_navg,"rate", rate_navg)

