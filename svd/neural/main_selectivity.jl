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
include("genFSRaster.jl")

include("figure/plt_pyr_selectivity.jl")
include("figure/plt_exc_selectivity.jl")
include("figure/plt_fs_selectivity.jl")
include("figure/plt_inh_selectivity.jl")
include("figure/plt_fsinh_selectivity.jl")
include("figure/plt_pyrexc_selectivity.jl")
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
dirsim_balanced = "/data/kimchm/data/dale/janelia/trained/balanced/"

#----- load targets -----#
pyr_R = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[1:100,:]
pyr_L = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[1:100,:]
fs_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

pyr_mean = ((mean(pyr_L, dims=1) + mean(pyr_R,dims=1))/2)[:]
pyr_diff = (mean(pyr_R,dims=1) - mean(pyr_L, dims=1))[:]

fs_mean = ((mean(fs_L, dims=1) + mean(fs_R,dims=1))/2)[:]
fs_diff = (mean(fs_R,dims=1) - mean(fs_L, dims=1))[:]
fs_diff_no0 = fs_diff[fs_diff .!= 0]
fs_mean_no0 = fs_mean[fs_diff .!= 0]

#----- load model neurons -----#
Ncells = 5000
Ne = 2500
fname_rate_R = dirsim * "rate_R_navg400.jld"
fname_rate_L = dirsim * "rate_L_navg400.jld"
rate_R = load(fname_rate_R, "rate")
rate_L = load(fname_rate_L, "rate")
rate_diff = mean(rate_R[:,:,1] - rate_L[:,:,1], dims=1)[:]

Excrate_R = rate_R[:,1:Ne,1]
Excrate_L = rate_L[:,1:Ne,1]
Excrate_mean = mean((Excrate_R + Excrate_L)/2, dims=1)[:]
Excrate_diff = mean(Excrate_R - Excrate_L, dims=1)[:]

Inhrate_R = rate_R[:,Ne+1:Ncells,1]
Inhrate_L = rate_L[:,Ne+1:Ncells,1]
Inhrate_mean = mean((Inhrate_R + Inhrate_L)/2, dims=1)[:]
Inhrate_diff = mean(Inhrate_R - Inhrate_L, dims=1)[:]

Excrate_diff_no0 = Excrate_diff[Excrate_diff .!= 0]
Excrate_mean_no0 = Excrate_mean[Excrate_diff .!= 0]
Inhrate_diff_no0 = Inhrate_diff[Inhrate_diff .!= 0]
Inhrate_mean_no0 = Inhrate_mean[Inhrate_diff .!= 0]


#----- load neural data -----#
pyr_rate_R, pyr_rate_L, fs_rate_R, fs_rate_L, 
pyr_cell_R, pyr_cell_L, fs_cell_R, fs_cell_L,
pyr_units_R, pyr_units_L, fs_units_R, fs_units_L,
pyr_spktime_R, pyr_eventtime_R, pyr_spktime_L, pyr_eventtime_L,
fs_spktime_R, fs_eventtime_R, fs_spktime_L, fs_eventtime_L = genSpikeRaster()

# simulate: lick right 
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
ubal_avg_R, uebal_avg_R, uibal_avg_R, uplas_avg_R, uext_avg_R, 
ubal_R, uebal_R, uibal_R, uplas_R, uext_R, times_R, ns_R = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1])

# simulate: lick left 
println("Loading files")
lickRL = "left"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
ubal_avg_L, uebal_avg_L, uibal_avg_L, uplas_avg_L, uext_avg_L, 
ubal_L, uebal_L, uibal_L, uplas_L, uext_L, times_L, ns_L = runsim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2])


####################################
#----- Pyramidal - Excitatory -----#
####################################

#----- pyramidal cells (data) -----#
ns_R_pyr, times_R_pyr, 
ns_L_pyr, times_L_pyr, cellsOrd_pyr, rate_diff_pyr = genPyrRaster(pyr_rate_R, pyr_rate_L, pyr_cell_R, pyr_cell_L, pyr_units_R, pyr_units_L, pyr_spktime_R, pyr_spktime_L, pyr_eventtime_R, pyr_eventtime_L)

plt_pyr_selectivity(dirfig, pyr_diff, pyr_R, pyr_L)

cellsOrd = sortperm(pyr_diff[almOrd]) ########## IMPORTANT!! Order exc neurons by selectivty of pyr cells ##########
matchedCellsOrd = matchedCells[cellsOrd]

plt_exc_selectivity(dirfig, Excrate_diff, Excrate_R, Excrate_L, matchedCellsOrd)

plt_pyrexc_selectivity(dirfig, pyr_diff, pyr_R, pyr_L, Excrate_diff, Excrate_R, Excrate_L, matchedCellsOrd)


#######################################
#----- Fast spiking - Inhibitory -----#
#######################################

plt_fs_selectivity(dirfig, fs_diff, fs_R, fs_L)

plt_inh_selectivity(dirfig, Inhrate_diff, Inhrate_R, Inhrate_L)

plt_fsinh_selectivity(dirfig, fs_diff, fs_R, fs_L, Inhrate_diff, Inhrate_R, Inhrate_L)

#----- plot fs (data) raster -----#
plt_fs_raster(dirfig, dirsim, fs_cell_R, fs_units_R, fs_spktime_R, fs_eventtime_R, fs_rate_R, pcor_R, mse_R,
fs_cell_L, fs_units_L, fs_spktime_L, fs_eventtime_L, fs_rate_L, pcor_L, mse_L, cells_matched)

#----- selectivity: plot fs (data) RL raster -----#
# ns_R_fs, times_R_fs, 
# ns_L_fs, times_L_fs, cellsOrd_fs, rate_diff_fs = genFSRaster(fs_rate_R, fs_rate_L, fs_cell_R, fs_cell_L, fs_units_R, fs_units_L, fs_spktime_R, fs_spktime_L, fs_eventtime_R, fs_eventtime_L)



#######################################
#----- Null model                -----#
#######################################

# PSTH: initial balanced network
Ncells = 5000
Ne = 2500
fname_rate_R_bal = dirsim_balanced * "rate_R_navg400.jld"
fname_rate_L_bal = dirsim_balanced * "rate_L_navg400.jld"
rate_R_bal = load(fname_rate_R_bal, "rate")
rate_L_bal = load(fname_rate_L_bal, "rate")
# ns_R_bal = load(dirsim_balanced * "ns_R_1.jld")["ns"]
# ns_L_bal = load(dirsim_balanced * "ns_L_1.jld")["ns"]
Inhrate_R_bal = rate_R_bal[:,Ne+1:Ncells,1]
Inhrate_L_bal = rate_L_bal[:,Ne+1:Ncells,1]
Inhrate_diff_bal = mean(Inhrate_R_bal - Inhrate_L_bal, dims=1)[:]

plt_null_selectivity(dirfig, Inhrate_diff_bal, Inhrate_R_bal, Inhrate_L_bal)



# #######################################
# #----- Single trial              -----#
# #######################################

# # single-trial selectivity (model): correlation betwen single trials and PSTH
# pcor_RL_exc = load_exc_1trials(dirsim, Excrate_diff)
# pcor_RL_inh = load_inh_1trials(dirsim, Inhrate_diff)
# pcor_RL_null = load_inh_1trials(dirsim_balanced, Inhrate_diff_bal)

# # single-trial selectivity (data): 
# pcor_RL_fs = load_fs_1trials(fs_rate_R, fs_rate_L, fs_cell_R, fs_cell_L, fs_units_R, fs_units_L, fs_spktime_R, fs_spktime_L, fs_eventtime_R, fs_eventtime_L)
# pcor_RL_pyr = load_pyr_1trials(pyr_rate_R, pyr_rate_L, pyr_cell_R, pyr_cell_L, pyr_units_R, pyr_units_L, pyr_spktime_R, pyr_spktime_L, pyr_eventtime_R, pyr_eventtime_L)

# plt_RL_1trial(dirfig, pcor_RL_exc, pcor_RL_inh, pcor_RL_null, pcor_RL_pyr, pcor_RL_fs)




#------------------------------------------------------------------------------------#





# figure()
# # hist(pyr_diff ./ pyr_mean, bins=100, histtype="step", color="black", density=true, label="Pyr")
# hist(pyr_diff, bins=100, histtype="step", color="black", density=true, label="Pyr")
# # hist(Excrate_diff_no0 ./ Excrate_mean_no0, bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# hist(Excrate_diff_no0, bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# xlabel("normalized selectivity")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_selectivity_normalized_ExcPyr0.png")


# figure()
# hist(fs_diff_no0 ./ fs_mean_no0, bins=30, histtype="step", color="black", density=true, label="FS")
# # hist(fs_diff_no0, bins=60, histtype="step", color="black", density=true, label="FS")
# hist(Inhrate_diff_no0 ./ Inhrate_mean_no0, bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# # hist(Inhrate_diff_no0, bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# xlim([-20, 15])
# xlabel("normalized selectivity")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_selectivity_normalized_InhFS0.png")








# figure()
# hist(log10.(abs.(pyr_diff)), bins=100, histtype="step", color="black", density=true, label="Pyr")
# # hist(log10.(abs.(fs_diff_no0)), bins=30, histtype="step", color="gray", density=true, label="FS")
# hist(log10.(abs.(Excrate_diff_no0)), bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# # hist(log10.(abs.(Inhrate_diff_no0)), bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# xlabel("log |selectivity|")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_selectivity_log_ExcPyr.png")



# figure()
# # hist(log10.(abs.(pyr_diff)), bins=100, histtype="step", color="black", density=true, label="Pyr")
# hist(log10.(abs.(fs_diff_no0)), bins=30, histtype="step", color="black", density=true, label="FS")
# # hist(log10.(abs.(Excrate_diff_no0)), bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# hist(log10.(abs.(Inhrate_diff_no0)), bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# xlabel("log |selectivity|")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_selectivity_log_InhFS.png")



# figure()
# hist(log10.(abs.(pyr_diff)), bins=100, histtype="step", color="black", density=true, label="Pyr")
# hist(log10.(abs.(fs_diff_no0)), bins=30, histtype="step", color="gray", density=true, label="FS")
# hist(log10.(abs.(Excrate_diff_no0)), bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# hist(log10.(abs.(Inhrate_diff_no0)), bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# xlabel("log |selectivity|")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_selectivity_log.png")





# figure()
# hist(log10.(abs.(pyr_diff ./ pyr_mean)), bins=100, histtype="step", color="black", density=true, label="Pyr")
# hist(log10.(abs.(pyr_mean)), bins=100, histtype="step", color="gray", density=true, label="Pyr")
# # hist(fs_diff_no0 ./ fs_mean_no0, bins=30, histtype="step", color="gray", density=true, label="FS")
# hist(log10.(abs.(Excrate_diff_no0 ./ Excrate_mean_no0)), bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# # hist(Inhrate_diff_no0 ./ Inhrate_mean_no0, bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# xlabel("normalized selectivity")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_tmp.png")


# figure()
# # hist(pyr_diff ./ pyr_mean, bins=100, histtype="step", color="black", density=true, label="Pyr")
# hist(log10.(abs.(fs_diff_no0 ./ fs_mean_no0)), bins=30, histtype="step", color="black", density=true, label="FS")
# hist(log10.(abs.(fs_mean_no0)), bins=30, histtype="step", color="gray", density=true, label="FS")
# # hist(Excrate_diff_no0 ./ Excrate_mean_no0, bins=100, histtype="step", color="limegreen", density=true, label="Exc")
# hist(log10.(abs.(Inhrate_diff_no0 ./ Inhrate_mean_no0)), bins=100, histtype="step", color="darkorange", density=true, label="Inh")
# hist(log10.(abs.(Inhrate_mean_no0)), bins=100, histtype="step", color="darkorange", density=true, alpha=0.5, label="Inh")
# # xlim([-1.5, 1.5])
# xlabel("normalized selectivity")
# ylabel("density")
# legend()
# tight_layout()
# savefig("fig_tmp1.png")




# #----------------------------------------------------------------------------------------------------#


# figure()
# prd = Erate_mean_R .* Erate_mean_L
# Erate_mean_R_no0 = Erate_mean_R[prd .!=0]
# Erate_mean_L_no0 = Erate_mean_L[prd .!=0]
# scatter(log10.(Erate_mean_R_no0), log10.(Erate_mean_L_no0), s=4)
# tight_layout()
# savefig("fig_E_rateRL.png")


# figure()
# scatter(log10.(Erate_mean), log10.(abs.(Excrate_diff)),s=9)
# xlim([0.2, 2.2])
# tight_layout()
# savefig("fig_selectivityE.png")

# figure()
# scatter(log10.(pyr_mean), log10.(abs.(pyr_diff)),s=9)
# tight_layout()
# savefig("fig_selectivityPYR.png")

# figure()
# scatter(log10.(Irate_mean), log10.(abs.(Inhrate_diff)),s=9)
# # xlim([0.2, 2.2])
# xlabel("log(inh rate)")
# ylabel("log|inh selectivity|")
# tight_layout()
# savefig("fig_selectivityI.png")

# figure()
# scatter(log10.(fs_mean), log10.(abs.(fs_diff)),s=9)
# tight_layout()
# savefig("fig_selectivityFS.png")



# alpha = 1.0
# exponent = 1.0 + alpha
# A = 1.0
# x = 0.1:0.01:50
# # y = log10.(A./x.^alpha)
# y = A./x.^exponent

# figure(); 
# cnt1, bin1 = hist(abs.(Excrate_diff), bins=200, histtype="step", density=true); 
# cnt2, bin2 = hist(abs.(rate_diff_pyr),bins=100, histtype="step",density=true); 
# close()

# figure(figsize=(2.5, 2.3))
# # plot(log10.(bin1[1:end-1]), log10.(cnt1), c="black", lw=0.8)
# # plot(log10.(bin2[1:end-1]), log10.(cnt2), c="limegreen", lw=0.8)
# # plot(log10.(x),log10.(y),c="red", linestyle="--", lw=1); 
# loglog(bin1[1:end-1], cnt1, c="black", lw=0.8, label="Exc")
# loglog(bin2[1:end-1], cnt2, c="limegreen", lw=0.8, label="Pyr")
# loglog(x,y,c="red", linestyle="--", lw=0.8); 
# legend(frameon=false, fontsize=7)
# xticks(fontsize=7)
# yticks(fontsize=7)
# xlabel("selectiviey", fontsize=7)
# ylabel("density", fontsize=7)
# # xlim([0, 1.5]); 
# xlim([0.1, 50])
# ylim([1e-3,1e2])
# tight_layout()
# savefig("fig_rdiff_pyr.png", dpi=600)




# alpha = 0.8
# exponent = 1.0 + alpha
# A = 0.8
# x = 0.1:0.01:50
# # y = log10.(A./x.^alpha)
# y = A./x.^exponent

# figure(); 
# cnt1, bin1 = hist(abs.(Inhrate_diff), bins=200, histtype="step", density=true); 
# cnt2, bin2 = hist(abs.(rate_diff_fs),bins=100, histtype="step",density=true); 
# close()

# figure(figsize=(2.5, 2.3))
# # plot(log10.(bin1[1:end-1]), log10.(cnt1), c="black", lw=0.8)
# # plot(log10.(bin2[1:end-1]), log10.(cnt2), c="darkorange", lw=0.8)
# # plot(log10.(x),log10.(y),c="red", linestyle="--", lw=1); 
# loglog(bin1[1:end-1], cnt1, c="black", lw=0.8, label="Inh")
# loglog(bin2[1:end-1], cnt2, c="darkorange", lw=0.8, label="FS")
# loglog(x,y,c="red", linestyle="--", lw=0.8); 
# legend(frameon=false, fontsize=7)
# xticks(fontsize=7)
# yticks(fontsize=7)
# xlabel("selectiviey (Hz)", fontsize=7)
# ylabel("density", fontsize=7)
# xlim([0.1, 50])
# ylim([1e-3,1e2])
# tight_layout()
# savefig("fig_rdiff_fs.png", dpi=600)


