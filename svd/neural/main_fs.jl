using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD
using HypothesisTests

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm
gridspec = matplotlib.gridspec

include("genSpikeRaster.jl")
include("genFSRaster.jl")
include("load_files.jl")
include("matchfs_nrep.jl")
include("matchfs_noreplacement.jl")
include("explainedVar.jl")
include("runmatchfs.jl")
include("runmatchfs_RL.jl")
include("runsim.jl")
include("load_1trials_model.jl")
include("load_1trials_data.jl")
include("figure/plt_fs_pcor.jl")
include("figure/plt_fsinh_psth.jl")
include("figure/plt_fsinh_psth_heatmap.jl")
include("figure/plt_inh_raster.jl")
# include("figure/plt_inh_RL_raster.jl")
include("figure/plt_RL_1trial.jl")
include("figure/plt_fs_raster.jl")
# include("figure/plt_fs_RL_raster.jl")
include("figure/plt_null_RL_raster.jl")

dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/wffwd/Lffwd5/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirsim_balanced = "/data/kimchm/data/dale/janelia/trained/balanced/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


# load FS
fs_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

####################################
#----- match fs and inh cells -----#
####################################
navg = 400 # 400
nrep = 1
# trained network
cells_matched, pcor_R, pcor_L, mse_R, mse_L = runmatchfs_RL(dirsim, navg, fs_R, fs_L, nrep)
# initial network
cells_matched_bal, pcor_R_bal, pcor_L_bal, mse_R_bal, mse_L_bal = runmatchfs_RL(dirsim_balanced, navg, fs_R, fs_L, nrep)

#----- fs performance -----#
plt_fs_pcor(dirfig, pcor_R, pcor_L, mse_R, mse_L, pcor_R_bal, pcor_L_bal, mse_R_bal, mse_L_bal)



###########################################
#----- plot fs and inh psth / raster -----#
###########################################
#----- compare fs and inh psth -----#
plt_fsinh_psth(dirfig, dirsim, fs_R, fs_L, pcor_R, pcor_L, cells_matched, "R")
plt_fsinh_psth(dirfig, dirsim, fs_R, fs_L, pcor_R, pcor_L, cells_matched, "L")
plt_fsinh_psth_heatmap(dirfig, dirsim, fs_R, fs_L, mse_R, mse_L, cells_matched)

#----- plot inh (model) raster -----#
plt_inh_raster(dirfig, dirsim, fs_R, fs_L, pcor_R, pcor_L, mse_R, mse_L, cells_matched)

# #----- plot fs (data) raster -----#
# pyr_rate_R, pyr_rate_L, fs_rate_R, fs_rate_L, 
# pyr_cell_R, pyr_cell_L, fs_cell_R, fs_cell_L,
# pyr_units_R, pyr_units_L, fs_units_R, fs_units_L,
# pyr_spktime_R, pyr_eventtime_R, pyr_spktime_L, pyr_eventtime_L,
# fs_spktime_R, fs_eventtime_R, fs_spktime_L, fs_eventtime_L = genSpikeRaster()

# plt_fs_raster(dirfig, dirsim, fs_cell_R, fs_units_R, fs_spktime_R, fs_eventtime_R, fs_rate_R, pcor_R, mse_R,
# fs_cell_L, fs_units_L, fs_spktime_L, fs_eventtime_L, fs_rate_L, pcor_L, mse_L, cells_matched)

# #----- selectivity: plot fs (data) RL raster -----#
# ns_R_fs, times_R_fs, 
# ns_L_fs, times_L_fs, cellsOrd_fs, rate_diff_fs = genFSRaster(fs_rate_R, fs_rate_L, fs_cell_R, fs_cell_L, fs_units_R, fs_units_L, fs_spktime_R, fs_spktime_L, fs_eventtime_R, fs_eventtime_L)

# plt_fs_RL_raster(dirfig, times_R_fs, times_L_fs, ns_R_fs, ns_L_fs, rate_diff_fs, cellsOrd_fs, fs_rate_R, fs_rate_L)


# #----- selectivity: plot inh (model) RL raster -----#
# # simulate: lick right 
# println("Loading files")
# lickRL = "right"
# p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
# wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
# ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
# ubal_avg_R, uebal_avg_R, uibal_avg_R, uplas_avg_R, uext_avg_R, 
# ubal_R, uebal_R, uibal_R, uplas_R, uext_R, times_R, ns_R = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1])

# # simulate: lick left 
# println("Loading files")
# lickRL = "left"
# p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
# wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
# ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
# ubal_avg_L, uebal_avg_L, uibal_avg_L, uplas_avg_L, uext_avg_L, 
# ubal_L, uebal_L, uibal_L, uplas_L, uext_L, times_L, ns_L = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2])

# # PSTH: inhibitory neuron
# Ncells = 5000
# Ne = 2500
# fname_rate_R = dirsim * "rate_R_navg400.jld"
# fname_rate_L = dirsim * "rate_L_navg400.jld"
# rate_R = load(fname_rate_R, "rate")
# rate_L = load(fname_rate_L, "rate")
# Excrate_R = rate_R[:,1:Ne,1]
# Excrate_L = rate_L[:,1:Ne,1]
# Excrate_diff = mean(Excrate_R - Excrate_L, dims=1)[:]
# Inhrate_R = rate_R[:,Ne+1:Ncells,1]
# Inhrate_L = rate_L[:,Ne+1:Ncells,1]
# Inhrate_diff = mean(Inhrate_R - Inhrate_L, dims=1)[:]

# plt_inh_RL_raster(dirfig, times_R, times_L, ns_R, ns_L, Inhrate_diff, Inhrate_R, Inhrate_L)

# # PSTH: initial balanced network
# Ncells = 5000
# Ne = 2500
# fname_rate_R_bal = dirsim_balanced * "rate_R_navg400.jld"
# fname_rate_L_bal = dirsim_balanced * "rate_L_navg400.jld"
# rate_R_bal = load(fname_rate_R_bal, "rate")
# rate_L_bal = load(fname_rate_L_bal, "rate")
# ns_R_bal = load(dirsim_balanced * "ns_R_1.jld")["ns"]
# ns_L_bal = load(dirsim_balanced * "ns_L_1.jld")["ns"]
# Inhrate_R_bal = rate_R_bal[:,Ne+1:Ncells,1]
# Inhrate_L_bal = rate_L_bal[:,Ne+1:Ncells,1]
# Inhrate_diff_bal = mean(Inhrate_R_bal - Inhrate_L_bal, dims=1)[:]

# plt_null_RL_raster(dirfig, ns_R_bal, ns_L_bal, Inhrate_diff_bal, Inhrate_R_bal, Inhrate_L_bal)

# # single-trial selectivity (model): correlation betwen single trials and PSTH
# pcor_RL_exc = load_exc_1trials(dirsim, Excrate_diff)
# pcor_RL_inh = load_inh_1trials(dirsim, Inhrate_diff)
# pcor_RL_null = load_inh_1trials(dirsim_balanced, Inhrate_diff_bal)

# # single-trial selectivity (data): 
# pcor_RL_fs = load_fs_1trials(fs_rate_R, fs_rate_L, fs_cell_R, fs_cell_L, fs_units_R, fs_units_L, fs_spktime_R, fs_spktime_L, fs_eventtime_R, fs_eventtime_L)
# pcor_RL_pyr = load_pyr_1trials(pyr_rate_R, pyr_rate_L, pyr_cell_R, pyr_cell_L, pyr_units_R, pyr_units_L, pyr_spktime_R, pyr_spktime_L, pyr_eventtime_R, pyr_eventtime_L)

# plt_RL_1trial(dirfig, pcor_RL_exc, pcor_RL_inh, pcor_RL_null, pcor_RL_pyr, pcor_RL_fs)
