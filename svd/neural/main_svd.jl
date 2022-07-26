using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD
using MultivariateStats

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm
gridspec = matplotlib.gridspec


include("explainedVar.jl")
include("sharedVar.jl")

# include("figure/plt_sharedvar_pyrexc.jl")
# include("figure/plt_sharedvar_fsinh.jl")
include("figure/plt_sharedvar_fsinh_RL.jl")
# include("figure/plt_sharedvar_data_ei.jl")
# include("figure/plt_cormat.jl")
include("figure/plt_cormat_RL.jl")
include("funMovAvg2D.jl")

dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirsim_balanced = "/data/kimchm/data/dale/janelia/trained/balanced/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


startind = 7

#-------- load Pyr and FS --------#
rtargPyr_R = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[startind:100,:]
rtargFS_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[startind:100,:]
rtargPyr_L = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[startind:100,:]
rtargFS_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[startind:100,:]

#-------- remove mean --------#
#lick right
for ci = 1:size(rtargPyr_R)[2]
    rtargPyr_R[:,ci] .-= mean(rtargPyr_R[:,ci])
end
for ci = 1:size(rtargFS_R)[2]
    rtargFS_R[:,ci] .-= mean(rtargFS_R[:,ci])
end
# lick left
for ci = 1:size(rtargPyr_L)[2]
    rtargPyr_L[:,ci] .-= mean(rtargPyr_L[:,ci])
end
for ci = 1:size(rtargFS_L)[2]
    rtargFS_L[:,ci] .-= mean(rtargFS_L[:,ci])
end

# #-------- project to subspace --------#
# Mpyr_R = fit(PCA, rtargPyr_R; pratio=0.99)
# Mfs_R = fit(PCA, rtargFS_R; pratio=0.99)
# Mpyr_L = fit(PCA, rtargPyr_L; pratio=0.99)
# Mfs_L = fit(PCA, rtargFS_L; pratio=0.99)

###########################################
#-------- load network activities --------#
###########################################

#-------- trained network --------#
navg_list = [50, 100, 200, 400]
navg = navg_list[4]
Ne = 2500
Ncells = 5000
wid = 0

fname_rate_navg_R = dirsim * "rate_R_navg$(navg).jld"
fname_rate_navg_L = dirsim * "rate_L_navg$(navg).jld"
rate_navg_R = funMovAvg2D(load(fname_rate_navg_R, "rate")[startind:100,:,1], wid)
rate_navg_L = funMovAvg2D(load(fname_rate_navg_L, "rate")[startind:100,:,1], wid)

rateExc_R = rate_navg_R[:,1:Ne]
rateInh_R = rate_navg_R[:,Ne+1:Ncells]
rateExc_L = rate_navg_L[:,1:Ne]
rateInh_L = rate_navg_L[:,Ne+1:Ncells]
for ci = 1:Ne
    rateExc_R[:,ci] .-= mean(rateExc_R[:,ci])
    rateInh_R[:,ci] .-= mean(rateInh_R[:,ci])
    rateExc_L[:,ci] .-= mean(rateExc_L[:,ci])
    rateInh_L[:,ci] .-= mean(rateInh_L[:,ci])
end

#-------- initial balanced network --------#
fname_rate_navg_R_bal = dirsim_balanced * "rate_R_navg$(navg).jld"
fname_rate_navg_L_bal = dirsim_balanced * "rate_L_navg$(navg).jld"
rate_navg_R_bal = funMovAvg2D(load(fname_rate_navg_R_bal, "rate")[startind:100,:,1], wid)
rate_navg_L_bal = funMovAvg2D(load(fname_rate_navg_L_bal, "rate")[startind:100,:,1], wid)

rateExc_R_bal = rate_navg_R_bal[:,1:Ne]
rateInh_R_bal = rate_navg_R_bal[:,Ne+1:Ncells]
rateExc_L_bal = rate_navg_L_bal[:,1:Ne]
rateInh_L_bal = rate_navg_L_bal[:,Ne+1:Ncells]
for ci = 1:Ne
    rateExc_R_bal[:,ci] .-= mean(rateExc_R_bal[:,ci])
    rateInh_R_bal[:,ci] .-= mean(rateInh_R_bal[:,ci])
    rateExc_L_bal[:,ci] .-= mean(rateExc_L_bal[:,ci])
    rateInh_L_bal[:,ci] .-= mean(rateInh_L_bal[:,ci])
end

##########################################
#-------- SVD of trained network --------#
##########################################
#-------- SVD: corr of FS and INH --------#
idx_fs_R = .!isnan.(cor(rateInh_R, rtargFS_R, dims=1)[1,:])
idx_inh_R = .!isnan.(cor(rateInh_R, rtargFS_R, dims=1)[:,1])
idx_fs_L = .!isnan.(cor(rateInh_L, rtargFS_L, dims=1)[1,:])
idx_inh_L = .!isnan.(cor(rateInh_L, rtargFS_L, dims=1)[:,1])

rtargFS_trim_R = rtargFS_R[:,idx_fs_R]
rateInh_trim_R = rateInh_R[:,idx_inh_R]
rtargFS_trim_L = rtargFS_L[:,idx_fs_L]
rateInh_trim_L = rateInh_L[:,idx_inh_L]

# FS and trained inh
cor_fsinh_R = cor(rtargFS_trim_R, rateInh_trim_R, dims=1)
svd_fsinh_R = svd(cor_fsinh_R)
cor_fsinh_L = cor(rtargFS_trim_L, rateInh_trim_L, dims=1)
svd_fsinh_L = svd(cor_fsinh_L)

#-------- shared variance (FS and trained inh) --------#
svar_fsinh_cor_R, svar_fsinh_FS_R, svar_fsinh_INH_R,
svar_fsinh_FS_ksum_R, svar_fsinh_INH_ksum_R, 
svar_fsinh_FS_tsum_R, svar_fsinh_INH_tsum_R = sharedVar(svd_fsinh_R, rtargFS_trim_R, rateInh_trim_R)

svar_fsinh_cor_L, svar_fsinh_FS_L, svar_fsinh_INH_L,
svar_fsinh_FS_ksum_L, svar_fsinh_INH_ksum_L, 
svar_fsinh_FS_tsum_L, svar_fsinh_INH_tsum_L = sharedVar(svd_fsinh_L, rtargFS_trim_L, rateInh_trim_L)


##########################################
#-------- SVD of initial network --------#
##########################################
#-------- SVD: corr of FS and initial INH --------#
idx_fs_R_bal = .!isnan.(cor(rateInh_R_bal, rtargFS_R, dims=1)[1,:])
idx_inh_R_bal = .!isnan.(cor(rateInh_R_bal, rtargFS_R, dims=1)[:,1])
idx_fs_L_bal = .!isnan.(cor(rateInh_L_bal, rtargFS_L, dims=1)[1,:])
idx_inh_L_bal = .!isnan.(cor(rateInh_L_bal, rtargFS_L, dims=1)[:,1])

rtargFS_trim_R_bal = rtargFS_R[:,idx_fs_R_bal]
rateInh_trim_R_bal = rateInh_R_bal[:,idx_inh_R_bal]
rtargFS_trim_L_bal = rtargFS_L[:,idx_fs_L_bal]
rateInh_trim_L_bal = rateInh_L_bal[:,idx_inh_L_bal]

# FS and initial inh
cor_fsinh_R_bal = cor(rtargFS_trim_R_bal, rateInh_trim_R_bal, dims=1)
svd_fsinh_R_bal = svd(cor_fsinh_R_bal)
cor_fsinh_L_bal = cor(rtargFS_trim_L_bal, rateInh_trim_L_bal, dims=1)
svd_fsinh_L_bal = svd(cor_fsinh_L_bal)

#-------- shared variance (FS and trained inh) --------#
bal_fsinh_cor_R, bal_fsinh_FS_R, bal_fsinh_INH_R,
bal_fsinh_FS_ksum_R, bal_fsinh_INH_ksum_R, 
bal_fsinh_FS_tsum_R, bal_fsinh_INH_tsum_R = sharedVar(svd_fsinh_R_bal, rtargFS_trim_R_bal, rateInh_trim_R_bal)

bal_fsinh_cor_L, bal_fsinh_FS_L, bal_fsinh_INH_L,
bal_fsinh_FS_ksum_L, bal_fsinh_INH_ksum_L, 
bal_fsinh_FS_tsum_L, bal_fsinh_INH_tsum_L = sharedVar(svd_fsinh_L_bal, rtargFS_trim_L_bal, rateInh_trim_L_bal)




# plot - FS and trained inh
fname = "fsinh"
plt_sharedvar_fsinh_RL(dirfig, fname, 
svar_fsinh_cor_R, svar_fsinh_FS_R, svar_fsinh_INH_R, 
svar_fsinh_FS_ksum_R, svar_fsinh_FS_tsum_R, svar_fsinh_INH_ksum_R, svar_fsinh_INH_tsum_R,
svar_fsinh_cor_L, svar_fsinh_FS_L, svar_fsinh_INH_L, 
svar_fsinh_FS_ksum_L, svar_fsinh_FS_tsum_L, svar_fsinh_INH_ksum_L, svar_fsinh_INH_tsum_L,
bal_fsinh_cor_R, bal_fsinh_FS_R, bal_fsinh_INH_R, 
bal_fsinh_FS_ksum_R, bal_fsinh_FS_tsum_R, bal_fsinh_INH_ksum_R, bal_fsinh_INH_tsum_R,
bal_fsinh_cor_L, bal_fsinh_FS_L, bal_fsinh_INH_L, 
bal_fsinh_FS_ksum_L, bal_fsinh_FS_tsum_L, bal_fsinh_INH_ksum_L, bal_fsinh_INH_tsum_L)


fname = "fsinh"
plt_cormat_RL(dirfig, fname, cor_fsinh_R, svd_fsinh_R, cor_fsinh_L, svd_fsinh_L)


# fname = "dataei"
# plt_sharedvar_data_ei(fname, svar_FSPYR_cor, svar_Pyr, svar_FS, 
# svar_FS_ksum, svar_FS_tsum, svar_Pyr_ksum, svar_Pyr_tsum,
# svar_excinh_cor, svar_exc, svar_inh, 
# svar_exc_ksum, svar_exc_tsum, svar_inh_ksum, svar_inh_tsum)

# # plot - pca of PYR and FS
# plt_pca(dirfig, Mpyr, Mfs)


# # plot - FS and PYR
# fname = "data"
# plt_sharedvar_data(fname, svar_FSPYR_cor, svar_Pyr, svar_FS, 
# svar_FS_ksum, svar_FS_tsum, svar_Pyr_ksum, svar_Pyr_tsum)

# # plot - trained exc and inh
# fname = "excinh"
# plt_sharedvar_excinh(fname, svar_excinh_cor, svar_exc, svar_inh, 
# svar_exc_ksum, svar_exc_tsum, svar_inh_ksum, svar_inh_tsum)








# figure(figsize=(3.5,3))
# lensvd = length(expvar_trained)
# nsing = 20
# plot(collect(1:nsing), expvar_trained[1:nsing], color="red", marker="o", ms=2, linewidth=1)
# plot(collect(1:nsing), expvar_balanced[1:nsing], color="black", marker="o", ms=2, linewidth=1)
# xlabel("rank", fontsize=12)
# ylabel("singular value", fontsize=12)
# xticks(fontsize=12)
# yticks(fontsize=12)
# legend_elements = [Line2D([0], [0], color="red", lw=2, label="trained"),
#                    Line2D([0], [0], color="black", lw=2, label="not trained")]
# legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
# tight_layout()
# savefig(dirfig * "singval.png", dpi=300)



# figure(figsize=(3.5,3))
# imshow(cor_trained, cmap="bwr", vmin=-1,vmax=1, interpolation="None", aspect="auto")
# xlabel("inh neuron model", fontsize=12)
# ylabel("FS cells", fontsize=12)
# colorbar()
# tight_layout()
# savefig(dirfig * "cormat_trained.png", dpi=300)


# sortrow_trained = sortperm(svd_trained.U[:,1])
# figure(figsize=(3.5,3))
# imshow(cor_trained[sortrow_trained,:], cmap="bwr", vmin=-1,vmax=1, interpolation="None", aspect="auto")
# xlabel("inh neuron model", fontsize=12)
# ylabel("FS cells", fontsize=12)
# colorbar()
# tight_layout()
# savefig(dirfig * "cormat_trained_ordered.png", dpi=300)


# figure(figsize=(3.5,3))
# imshow(cor_balanced, cmap="bwr", vmin=-0.3, vmax=0.3, interpolation="None", aspect="auto")
# xlabel("inh neuron model", fontsize=12)
# ylabel("FS cells", fontsize=12)
# colorbar()
# tight_layout()
# savefig(dirfig * "cormat_balanced.png", dpi=300)


# sortrow_balanced = sortperm(svd_balanced.U[:,1])
# figure(figsize=(3.5,3))
# imshow(cor_balanced[sortrow_balanced,:], cmap="bwr", vmin=-0.3, vmax=0.3, interpolation="None", aspect="auto")
# xlabel("inh neuron model", fontsize=12)
# ylabel("FS cells", fontsize=12)
# colorbar()
# tight_layout()
# savefig(dirfig * "cormat_balanced_ordered.png", dpi=300)



