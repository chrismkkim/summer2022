using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD
using MultivariateStats
using GLM
# using DataFrames


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
include("funMovAvg2D.jl")

include("figure/plt_fs_pcor.jl")
include("figure/plt_fsinh_psth.jl")
include("figure/plt_inh_raster.jl")
# include("figure/plt_inh_RL_raster.jl")
include("figure/plt_fs_raster.jl")
# include("figure/plt_fs_RL_raster.jl")
# include("figure/plt_null_RL_raster.jl")

include("figure/plt_pcomp.jl")
include("figure/plt_pcomp_pyrexc.jl")
include("figure/plt_pcomp_fsinh.jl")
include("figure/plt_pcomp_bal.jl")
include("figure/plt_pcor_pcomp.jl")


dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/wffwd/Lffwd5/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirsim_balanced = "/data/kimchm/data/dale/janelia/trained/balanced/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


# load data
pyr_R = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[1:100,:]
pyr_L = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[1:100,:]
fs_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

# load trained network
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells_excTrained,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)


####################################
#----- match fs and inh cells -----#
####################################
navg = 400 # 400
nrep = 1
# trained network
cells_matched, pcor_R, pcor_L, mse_R, mse_L = runmatchfs_RL(dirsim, navg, fs_R, fs_L, nrep)



#----- PCA of pyr rates -----#
# remove mean 
for ci = 1:size(pyr_R)[2]
    pyr_R[:,ci] .-= mean(pyr_R[:,ci])
    pyr_L[:,ci] .-= mean(pyr_L[:,ci])
end
Mpyr_R = fit(PCA, pyr_R; pratio=0.99)
Mpyr_L = fit(PCA, pyr_L; pratio=0.99)



#----- PCA of fs rates -----#
# remove mean 
for ci = 1:size(fs_R)[2]
    fs_R[:,ci] .-= mean(fs_R[:,ci])
    fs_L[:,ci] .-= mean(fs_L[:,ci])
end
Mfs_R = fit(PCA, fs_R; pratio=0.99)
Mfs_L = fit(PCA, fs_L; pratio=0.99)

fs_R_pcomp = transform(Mfs_R, fs_R)
fs_L_pcomp = transform(Mfs_L, fs_L)




################################################
# relationship between performance and PC1
################################################
plt_pcor_pcomp(dirfig, pcor_R, pcor_L, fs_R_pcomp, fs_L_pcomp)





########################
#----- Plot PCA's -----# 
########################

#----- PCA: PSTH of trained inhibitory neuron -----# 
startind = 7
wid = 0

Ncells = 5000
Ne = 2500
matchedCells_excTrained = sort(matchedCells_excTrained)
excNotTrained = deleteat!(collect(1:Ne), matchedCells_excTrained)

fname_rate_R = dirsim * "rate_R_navg400.jld"
fname_rate_L = dirsim * "rate_L_navg400.jld"
rate_R = load(fname_rate_R, "rate")[:,:,1]
rate_L = load(fname_rate_L, "rate")[:,:,1]
Excrate_R = funMovAvg2D(rate_R[startind:end, matchedCells_excTrained], wid)
Excrate_L = funMovAvg2D(rate_L[startind:end, matchedCells_excTrained], wid)
Inhrate_R = funMovAvg2D(rate_R[startind:end,Ne+1:Ncells], wid)
Inhrate_L = funMovAvg2D(rate_L[startind:end,Ne+1:Ncells], wid)
Excrate_R_NotTrained = funMovAvg2D(rate_R[startind:end, excNotTrained], wid)
Excrate_L_NotTrained = funMovAvg2D(rate_L[startind:end, excNotTrained], wid)

# remove mean 
for ci = 1:size(Excrate_R)[2]
    Excrate_R[:,ci] .-= mean(Excrate_R[:,ci])
    Excrate_L[:,ci] .-= mean(Excrate_L[:,ci])
end
for ci = 1:size(Inhrate_R)[2]
    Inhrate_R[:,ci] .-= mean(Inhrate_R[:,ci])
    Inhrate_L[:,ci] .-= mean(Inhrate_L[:,ci])
end
for ci = 1:size(Excrate_R_NotTrained)[2]
    Excrate_R_NotTrained[:,ci] .-= mean(Excrate_R_NotTrained[:,ci])
    Excrate_L_NotTrained[:,ci] .-= mean(Excrate_L_NotTrained[:,ci])
end
Mexc_R = fit(PCA, Excrate_R; pratio=1)
Mexc_L = fit(PCA, Excrate_L; pratio=1)
Minh_R = fit(PCA, Inhrate_R; pratio=1)
Minh_L = fit(PCA, Inhrate_L; pratio=1)
Mexc_R_NotTrained = fit(PCA, Excrate_R_NotTrained; pratio=1)
Mexc_L_NotTrained = fit(PCA, Excrate_L_NotTrained; pratio=1)

#----- PCA: PSTH of balanced inhibitory neuron -----# 
Ncells = 5000
Ne = 2500
fname_rate_R_bal = dirsim_balanced * "rate_R_navg400.jld"
fname_rate_L_bal = dirsim_balanced * "rate_L_navg400.jld"
rate_R_bal = load(fname_rate_R_bal, "rate")[:,:,1]
rate_L_bal = load(fname_rate_L_bal, "rate")[:,:,1]
Inhrate_R_bal = funMovAvg2D(rate_R_bal[startind:end,Ne+1:Ncells,1], wid)
Inhrate_L_bal = funMovAvg2D(rate_L_bal[startind:end,Ne+1:Ncells,1], wid)

# remove mean 
for ci = 1:size(Inhrate_R_bal)[2]
    Inhrate_R_bal[:,ci] .-= mean(Inhrate_R_bal[:,ci])
    Inhrate_L_bal[:,ci] .-= mean(Inhrate_L_bal[:,ci])
end
Mbal_R = fit(PCA, Inhrate_R_bal; pratio=0.99)
Mbal_L = fit(PCA, Inhrate_L_bal; pratio=0.99)


fname_pyr = "pyrexcR"
plt_pcomp_pyrexc(dirfig, fname_pyr, startind, Mpyr_R, Mexc_R)

fname_pyr = "pyrexcL"
plt_pcomp_pyrexc(dirfig, fname_pyr, startind, Mpyr_L, Mexc_L)


fname_pyr = "pyrexcR_NotTrained"
plt_pcomp_pyrexc(dirfig, fname_pyr, startind, Mpyr_R, Mexc_R_NotTrained)

fname_pyr = "pyrexcL_NotTrained"
plt_pcomp_pyrexc(dirfig, fname_pyr, startind, Mpyr_L, Mexc_L_NotTrained)


fname_fs = "fsinhR"
plt_pcomp_fsinh(dirfig, fname_fs, startind, Mfs_R, Minh_R)

fname_fs = "fsinhL"
plt_pcomp_fsinh(dirfig, fname_fs, startind, Mfs_L, Minh_L)


fname_bal = "bal"
plt_pcomp_bal(dirfig, fname_bal, startind, Mbal_L)




# fname_pyr = "pyr"
# plt_pcomp(dirfig, fname_pyr, 1, Mpyr_L, Mpyr_R)

# fname_fs = "fs"
# plt_pcomp(dirfig, fname_fs, 1, Mfs_L, Mfs_R)

# fname_exc = "exc"
# plt_pcomp(dirfig, fname_exc, startind, Mexc_L, Mexc_R)

# fname_inh = "inh"
# plt_pcomp(dirfig, fname_inh, startind, Minh_L, Minh_R)

# fname_bal = "bal"
# plt_pcomp(dirfig, fname_bal, startind, Mbal_L, Mbal_R)




# fs_R_tmp = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
# fs_L_tmp = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]
# fs_R_mean = mean(fs_R_tmp, dims=1)[:]
# fs_L_mean = mean(fs_L_tmp, dims=1)[:]
# fs_R_std = std(fs_R_tmp, dims=1)[:]
# fs_L_std = std(fs_L_tmp, dims=1)[:]

# fs_R_pcompsub = sum(abs.(fs_R_pcomp[2:end,:]),dims=1)[:]
# fs_L_pcompsub = sum(abs.(fs_L_pcomp[2:end,:]),dims=1)[:]

# RL_rate = hcat(fs_L_mean, fs_R_mean)
# RL_pcomp = hcat(fs_L_pcomp[1,:], fs_R_pcomp[1,:])
# RL_pcor = hcat(pcor_L, pcor_R)

# idx = argmin(abs.(RL_pcomp), dims=2)
# slope_pcor = RL_pcor[idx]

# idx = argmax(RL_rate, dims=2)
# rate_pcor = RL_pcor[idx]
# rate_diff = fs_R_mean - fs_L_mean





# figure()
# # plot(mean(fs_R, dims=2), c="red")
# plot(projection(Mfs_R)[:,1]*10, c="black")

# savefig("fs_pc1.png")


# figure()
# scatter(pcor_R, log.(fs_R_pc1), s=9, c="blue")
# scatter(pcor_L, log.(fs_L_pc1), s=9, c="red")
# xlabel("cor")
# ylabel("explained var by pc1")
# tight_layout()

# savefig("fs_pcor_expvar.png")




# figure(figsize=(3,3))
# for i = 1:4
#     subplot(2,2,i)
#     plot(pcor_R, fs_R_pcomp[i,:],marker="o", ms=1, mec="None", c="blue", linestyle="")
#     title("PC$(i)", fontsize=8)
#     xlim([-1,1])
#     ylim([-60,60])
#     xticks(fontsize=8)
#     yticks(fontsize=8)
#     xlabel("corr", fontsize=8)
#     ylabel("pcomp", fontsize=8)
# end

# tight_layout()

# savefig("fs_pcor_pcs_R.png", dpi=600)



# figure(figsize=(3,3))
# for i = 1:4
#     subplot(2,2,i)
#     plot(pcor_L, fs_L_pcomp[i,:],marker="o", ms=1, mec="None", c="red", linestyle="")
#     title("PC$(i)")
#     xlim([-1,1])
#     ylim([-60,60])
#     xticks(fontsize=8)
#     yticks(fontsize=8)
#     xlabel("corr", fontsize=8)
#     ylabel("pcomp", fontsize=8)
# end

# tight_layout()

# savefig("fs_pcor_pcs_L.png", dpi=600)




# figure(figsize=(1.8,1.5))

# plot(pcor_R, log.(abs.(fs_R_pcomp[1,:] .- 15)), marker="o", ms=1, mec="None", c="blue", alpha=1, linestyle="")
# plot(pcor_L, log.(abs.(fs_L_pcomp[1,:] .- 5)), marker="o", ms=1, mec="None", c="red", alpha=1, linestyle="")

# xlim([-1,1])
# ylim([-3,6])
# xticks(fontsize=8)
# yticks(fontsize=8)
# xlabel("corr", fontsize=8)
# ylabel(L"$\log \vert pc_1 - b_{R,L} \vert $", fontsize=8)
# tight_layout()

# savefig("fs_pcor_pcomp.png", dpi=600)




# figure(figsize=(1.7,1.5))

# plot(pcor_R, log.(abs.(fs_R_pcomp[3,:] .- 15)), marker="o", ms=1, mec="None", c="blue", alpha=1, linestyle="")
# plot(pcor_L, log.(abs.(fs_L_pcomp[3,:] .- 5)), marker="o", ms=1, mec="None", c="red", alpha=1, linestyle="")

# xlim([-0.75,1.2])
# ylim([-3,6])
# xticks(fontsize=8)
# yticks(fontsize=8)
# xlabel("cor", fontsize=8)
# ylabel("log |PC3 - b|", fontsize=8)
# tight_layout()

# savefig("fs_pcor_pcomp3.png", dpi=600)




