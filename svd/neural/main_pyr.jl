using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

include("load_files.jl")
include("runsim.jl")
include("runsimlong.jl")
include("runperformance_nrep.jl")
include("funMovAvg.jl")
include("funMovAvg_2d.jl")
include("funRollingAvg.jl")
include("funSample.jl")
include("psth.jl")
include("genSpikeRaster.jl")
include("genPyrRaster.jl")

# include("figure/plt_pyr_RL_raster.jl")
# include("figure/plt_exc_RL_raster.jl")
include("figure/plt_pyrexc_psth.jl")
include("figure/plt_pyrexc_psth_heatmap.jl")
include("figure/plt_pyr_raster.jl")
include("figure/plt_exc_raster.jl")
include("figure/plt_pyr_pcor.jl")
include("figure/plt_syn_dist.jl")
include("figure/plt_syn.jl")
include("figure/plt_fano.jl")
include("figure/plt_rate_dist_model.jl")
include("figure/plt_rate_dist_data.jl")

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm
gridspec = matplotlib.gridspec

navg_list = [50, 100, 200, 400]

dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/wffwd/Lffwd5/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


#----- load targets -----#
pyr_R = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[1:100,:]
pyr_L = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[1:100,:]
fs_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

#----- lick right -----#
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
ubal_avg_R, uebal_avg_R, uibal_avg_R, uplas_avg_R, uext_avg_R, 
ubal_R, uebal_R, uibal_R, uplas_R, uext_R, times_R, ns_R = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1])

#----- lick left -----#
println("Loading files")
lickRL = "left"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
ubal_avg_L, uebal_avg_L, uibal_avg_L, uplas_avg_L, uext_avg_L, 
ubal_L, uebal_L, uibal_L, uplas_L, uext_L, times_L, ns_L = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2])

# #################################################
# #----- Lick right / left population raster -----#
# #----- order neurons by selectivity -----#
# #################################################

# #----- pyramidal cells (data) -----#
# pyr_rate_R, pyr_rate_L, fs_rate_R, fs_rate_L, 
# pyr_cell_R, pyr_cell_L, fs_cell_R, fs_cell_L,
# pyr_units_R, pyr_units_L, fs_units_R, fs_units_L,
# pyr_spktime_R, pyr_eventtime_R, pyr_spktime_L, pyr_eventtime_L,
# fs_spktime_R, fs_eventtime_R, fs_spktime_L, fs_eventtime_L = genSpikeRaster()

# ns_R_pyr, times_R_pyr, 
# ns_L_pyr, times_L_pyr, cellsOrd_pyr, rate_diff_pyr = genPyrRaster(pyr_rate_R, pyr_rate_L, pyr_cell_R, pyr_cell_L, pyr_units_R, pyr_units_L, pyr_spktime_R, pyr_spktime_L, pyr_eventtime_R, pyr_eventtime_L)

# plt_pyr_RL_raster(dirfig, times_R_pyr, times_L_pyr, ns_R_pyr, ns_L_pyr, rate_diff_pyr, pyr_R, pyr_L, cellsOrd_pyr)

# #----- excitatory neuron (model) -----#
# fname_rate_R = dirsim * "rate_R_navg400.jld"
# fname_rate_L = dirsim * "rate_L_navg400.jld"
# rate_R = load(fname_rate_R, "rate")
# rate_L = load(fname_rate_L, "rate")
# rate_diff = mean(rate_R[:,:,1] - rate_L[:,:,1], dims=1)[:]
# Excrate_R = rate_R[:,1:p.Ne,1]
# Excrate_L = rate_L[:,1:p.Ne,1]

# cellsOrd = sortperm(rate_diff_pyr[almOrd]) ########## IMPORTANT!! Order exc neurons by selectivty of pyr cells ##########
# matchedCellsOrd = matchedCells[cellsOrd]

# plt_exc_RL_raster(dirfig, times_R, times_L, ns_R, ns_L, rate_diff, Excrate_R, Excrate_L, cellsOrd, matchedCellsOrd)

#############################################
#----- Lick right / left neuron raster -----#
#############################################

#----- plot pyr-exc performance ------#
pyr_pcor_R, pyr_pcor_L = plt_pyr_pcor(dirfig, dirsim, navg_list, pyr_R, pyr_L, almOrd, matchedCells)

#----- compare pyr and exc psth ------#
plt_pyrexc_psth(dirfig, dirsim, pyr_R, pyr_L, pyr_pcor_R, pyr_pcor_L, almOrd, matchedCells, "R")
plt_pyrexc_psth(dirfig, dirsim, pyr_R, pyr_L, pyr_pcor_R, pyr_pcor_L, almOrd, matchedCells, "L")
plt_pyrexc_psth_heatmap(dirfig, dirsim, pyr_R, pyr_L, pyr_pcor_R, pyr_pcor_L, almOrd, matchedCells)

#----- plot pyr, exc raster ------#
plt_exc_raster(dirfig, dirsim, pyr_R, pyr_L, pyr_pcor_R, pyr_pcor_L, almOrd, matchedCells)

plt_pyr_raster(dirfig, dirsim, pyr_cell_R, pyr_units_R, pyr_spktime_R, pyr_eventtime_R, pyr_rate_R, pyr_pcor_R,
pyr_cell_L, pyr_units_L, pyr_spktime_L, pyr_eventtime_L, pyr_rate_L, pyr_pcor_L, almOrd, matchedCells)

#########################
#----- fano factor -----#
#########################
ff_R, ff_L = fanofactor(dirsim, p)
plt_fano(dirfig, ff_R, ff_L)


#################################
#----- synaptic activities -----#
#################################

#----- distribution of synaptic inputs ------#
plt_syn_dist(dirfig, uplas_avg_R, ubal_avg_R, uebal_avg_R, uibal_avg_R, 
uplas_avg_L, ubal_avg_L, uebal_avg_L, uibal_avg_L)

#----- synaptic inputs to neurons ------#
plt_syn(dirfig, p, ubal_R, uebal_R, uibal_R, uplas_R)


######################################
#----- firing rate distribution -----#
######################################

#----- Trained model: lick right -----#
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

simtime = 200
times_long, ns_long = runsimlong(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1],simtime)
mrate_long = ns_long / simtime # simulation length: 20 sec

plt_rate_dist_model(dirfig, mrate_long)

#----- Neural data: Pyr/FS, lick right/left -----#
pyr_R_tmp = transpose(load(dirtarget_selectivity * "movingrate_Pyr_lickright.jld", "Pyr"))[1:100,:]
pyr_L_tmp = transpose(load(dirtarget_selectivity * "movingrate_Pyr_lickleft.jld", "Pyr"))[1:100,:]
fs_R_tmp = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L_tmp = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

plt_rate_dist_data(dirfig, pyr_R_tmp, pyr_L_tmp, fs_R_tmp, fs_L_tmp)





