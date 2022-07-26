using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD
using MultivariateStats

# modifiable parameters
g_list = [1.0, 1.5]
ii = 1 # g = 1.0

# lam_list = collect(1.0:1.0:3.0) # 0.2
# lam_list = [0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
lam_list = [0.04, 0.08, 0.2, 0.6, 1.0]
jj = 1 

mu_list = [0.0, 2.0, 4.0, 8.0]
kk = 4 # mu = 8.0

L_list = [2.0, 4.0, 6.0, 8.0, 10.0]
# ll = parse(Int64,ARGS[1]) # L = 2.0
ll = 3

wpscale_list = [1.0, 2.0, 3.0]
mm = 2 # wpscale = 2.0. To improve dale's law try wpscale = 1.0 (Jan 6)

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
# nn = parse(Int64,ARGS[2])
nn = 5

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm
gridspec = matplotlib.gridspec


include("explainedVar.jl")
include("sharedVar.jl")

include("figure/plt_sharedvar_pyrexc.jl")
include("figure/plt_sharedvar_fsinh.jl")
# include("figure/plt_sharedvar_data_ei.jl")
include("figure/plt_pca.jl")
include("figure/plt_cormat.jl")

# temporalAvg = "noavg" # "avg", "noavg"

# dirsim = "/data/kimchm/data/dale/janelia/trained/basic/sim/" * temporalAvg * "/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
# dirsim_bal = "/data/kimchm/data/dale/janelia/trained/basic/simffwd/balanced/" * temporalAvg * "/"
# dirtarget = "/data/kimchm/data/dale/janelia/s1alm/target/lickright/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
# dirutargPyr = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
# dirutargFS = "/data/kimchm/data/dale/janelia/s1alm/target/utargFS/"
# dirfig = "/data/kimchm/data/dale/janelia/figure/ffwd/pca/" * temporalAvg * "/psth/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"

# if ~ispath(dirfig)
#     mkpath(dirfig)
# end

# lickRL = "right"
lickRL = "left"

startind = 7

#-------- load Pyr and FS --------#
if lickRL == "right"
    rtargPyr = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[startind:100,:]
    rtargFS = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[startind:100,:]
elseif lickRL == "left"
    rtargPyr = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[startind:100,:]
    rtargFS = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[startind:100,:]
end

for ci = 1:size(rtargPyr)[2]
    rtargPyr[:,ci] .-= mean(rtargPyr[:,ci])
end
for ci = 1:size(rtargFS)[2]
    rtargFS[:,ci] .-= mean(rtargFS[:,ci])
end

# project to subspace
Mpyr = fit(PCA, rtargPyr; pratio=0.99)
Mfs = fit(PCA, rtargFS; pratio=0.99)

#-------- load INH --------#
navg_list = [50, 100, 200, 400]
navg = navg_list[4]
Ne = 2500
Ncells = 5000

if lickRL == "right"
    fname_rate_navg = dirsim * "rate_R_navg$(navg).jld"
    rate_navg = load(fname_rate_navg, "rate")[startind:100,:,1]
elseif lickRL == "left"
    fname_rate_navg = dirsim * "rate_L_navg$(navg).jld"
    rate_navg = load(fname_rate_navg, "rate")[startind:100,:,1]
end
rateExc = rate_navg[:,1:Ne]
rateInh = rate_navg[:,Ne+1:Ncells]
for ci = 1:Ne
    rateExc[:,ci] .-= mean(rateExc[:,ci])
    rateInh[:,ci] .-= mean(rateInh[:,ci])
end

# #-------- load INH of balanced network --------#
# fname_bal_rate_navg = dirsim_bal * "bal_rate_navg$(navg).jld"
# rate_navg_bal = load(fname_bal_rate_navg, "rate")[3:100,:,1]
# rateExc_bal = rate_navg_bal[:,1:Ne]
# rateInh_bal = rate_navg_bal[:,Ne+1:Ncells]
# for ci = 1:Ne
#     rateExc_bal[:,ci] .-= mean(rateExc_bal[:,ci])
#     rateInh_bal[:,ci] .-= mean(rateInh_bal[:,ci])
# end


#-------- SVD: corr of FS and INH --------#
idx_pyr = .!isnan.(cor(rateExc, rtargPyr, dims=1)[1,:])
idx_exc = .!isnan.(cor(rateExc, rtargPyr, dims=1)[:,1])
idx_fs = .!isnan.(cor(rateInh, rtargFS, dims=1)[1,:])
idx_inh = .!isnan.(cor(rateInh, rtargFS, dims=1)[:,1])
# idx_bal = .!isnan.(cor(rtargFS, rateInh_bal, dims=1)[1,:])

rtargPyr_trim = rtargPyr[:,idx_pyr]
rateExc_trim = rateExc[:,idx_exc]
rtargFS_trim = rtargFS[:,idx_fs]
rateInh_trim = rateInh[:,idx_inh]

# rtargPyr_trim = rtargPyr[:,idx_pyr]
# rtargFS_trim = rtargFS[:,idx_fs]
# rateInh_bal = rateInh_bal[:,idx_bal]

# # FS and Pyr
# cor_fspyr = cor(rtargFS_trim, rtargPyr_trim, dims=1)
# svd_fspyr = svd(cor_fspyr)

# Pyr and trained exc
cor_pyrexc = cor(rtargPyr_trim, rateExc_trim, dims=1)
svd_pyrexc = svd(cor_pyrexc)

# FS and trained inh
cor_fsinh = cor(rtargFS_trim, rateInh_trim, dims=1)
svd_fsinh = svd(cor_fsinh)

# # Pyr and untrained balanced 
# cor_balanced_pyrexc = cor(rtargPyr_trim, rateExc_bal, dims=1)
# svd_balanced_pyrexc = svd(cor_balanced_pyrexc)

# # FS and untrained balanced 
# cor_balanced_fsinh = cor(rtargFS_trim, rateInh_bal, dims=1)
# svd_balanced_fsinh = svd(cor_balanced_fsinh)

# # trained exc and inh
# cor_excinh = cor(rateExc, rateInh, dims=1)
# svd_excinh = svd(cor_excinh)

#-------- shared variance (Pyr and trained exc) --------#
svar_pyrexc_cor, svar_pyrexc_Pyr, svar_pyrexc_Exc,
svar_pyrexc_Pyr_ksum, svar_pyrexc_Exc_ksum, 
svar_pyrexc_Pyr_tsum, svar_pyrexc_Exc_tsum = sharedVar(svd_pyrexc, rtargPyr_trim, rateExc_trim)

#-------- shared variance (FS and trained inh) --------#
svar_fsinh_cor, svar_fsinh_FS, svar_fsinh_INH,
svar_fsinh_FS_ksum, svar_fsinh_INH_ksum, 
svar_fsinh_FS_tsum, svar_fsinh_INH_tsum = sharedVar(svd_fsinh, rtargFS_trim, rateInh_trim)

# #-------- shared variance (Pyr and balanced) --------#
# svar_balanced_cor, svar_balanced_Pyr, svar_balanced_Exc,
# svar_balanced_Pyr_ksum, svar_balanced_Exc_ksum, 
# svar_balanced_Pyr_tsum, svar_balanced_Exc_tsum = sharedVar(svd_balanced_pyrexc, rtargPyr_trim, rateExc_bal)

# #-------- shared variance (FS and balanced) --------#
# svar_balanced_cor, svar_balanced_FS, svar_balanced_INH,
# svar_balanced_FS_ksum, svar_balanced_INH_ksum, 
# svar_balanced_FS_tsum, svar_balanced_INH_tsum = sharedVar(svd_balanced_fsinh, rtargFS_trim, rateInh_bal)

# #-------- shared variance (FS and Pyr) --------#
# svar_FSPYR_cor, svar_FS, svar_Pyr,
# svar_FS_ksum, svar_Pyr_ksum, 
# svar_FS_tsum, svar_Pyr_tsum = sharedVar(svd_fspyr, rtargFS_trim, rtargPyr_trim)

# #-------- shared variance (trained exc and inh) --------#
# svar_excinh_cor, svar_exc, svar_inh,
# svar_exc_ksum, svar_inh_ksum, 
# svar_exc_tsum, svar_inh_tsum = sharedVar(svd_excinh, usumExc, rateInh)


# plot - Pyr and trained exc
if lickRL == "right"
    fname = "pyrexc_R"
elseif lickRL == "left"
    fname = "pyrexc_L"
end
plt_sharedvar_pyrexc(dirfig, fname, svar_pyrexc_cor, svar_pyrexc_Pyr, svar_pyrexc_Exc, 
svar_pyrexc_Pyr_ksum, svar_pyrexc_Pyr_tsum, svar_pyrexc_Exc_ksum, svar_pyrexc_Exc_tsum)

# plot - FS and trained inh
if lickRL == "right"
    fname = "fsinh_R"
elseif lickRL == "left"
    fname = "fsinh_L"
end    
plt_sharedvar_fsinh(dirfig, fname, svar_fsinh_cor, svar_fsinh_FS, svar_fsinh_INH, 
svar_fsinh_FS_ksum, svar_fsinh_FS_tsum, svar_fsinh_INH_ksum, svar_fsinh_INH_tsum)


if lickRL == "right"
    fname = "pyrexc_R"
elseif lickRL == "left"
    fname = "pyrexc_L"
end
plt_cormat(dirfig, fname, cor_pyrexc, svd_pyrexc)


if lickRL == "right"
    fname = "fsinh_R"
elseif lickRL == "left"
    fname = "fsinh_L"
end
plt_cormat(dirfig, fname, cor_fsinh, svd_fsinh)


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



