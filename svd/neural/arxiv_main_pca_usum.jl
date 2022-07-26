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


# include("param.jl")
# include("genWeights.jl")
# include("genPlasticWeights.jl")
# include("convertWgtIn2Out.jl")
# include("genTarget.jl")
# include("genStim.jl")
# include("genCellsTrained.jl")
# include("runinitial.jl")
# include("runtrain.jl")
# # include("runtrain_tmp.jl")
# include("runtest.jl")
# include("funMovAvg.jl")
# include("funCorrTarg.jl")
# include("funCorrDecomp.jl")
# include("funSample.jl")
# include("funRollingAvg.jl")
# include("runperformance_train.jl")
# include("runperformance_test.jl")
# include("calcWeights.jl")
# include("replaceWp.jl")

# include("genWeightsFfwd.jl")
# include("genffwdRate.jl")
# include("genffwdRateSub.jl")
# include("runffwd.jl")
# include("psth.jl")
include("explainedVar.jl")
include("sharedVar.jl")

include("figure/plt_sharedvar_trained.jl")
include("figure/plt_sharedvar_data_ei.jl")
include("figure/plt_pca.jl")

temporalAvg = "noavg" # "avg", "noavg"

dirsim = "/data/kimchm/data/dale/janelia/trained/basic/sim/" * temporalAvg * "/"
dirsim_bal = "/data/kimchm/data/dale/janelia/trained/basic/simffwd/balanced/" * temporalAvg * "/"
dirtarget = "/data/kimchm/data/dale/janelia/s1alm/target/"
dirutargPyr = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
dirutargFS = "/data/kimchm/data/dale/janelia/s1alm/target/utargFS/"
dirfig = "/data/kimchm/data/dale/janelia/figure/ffwd/pca/" * temporalAvg * "/usum/"

if ~ispath(dirfig)
    mkpath(dirfig)
end

#-------- load Pyr and FS --------#
rtargPyr = transpose(load(dirtarget * "movingrate_Pyr1Hz.jld", "Pyr"))[3:100,:]
rtargFS = transpose(load(dirtarget * "movingrate_FS.jld", "FS"))[3:100,:]

#-------- load Pyr --------#
utargPyr = transpose(load(dirutargPyr * "utarg1Hz.jld", "utarg"))[3:100,:]
Npyr = size(utargPyr)[2]
for ci = 1:Npyr
    utargPyr[:,ci] = utargPyr[:,ci] .- mean(utargPyr[:,ci])
end
# project to subspace
Mpyr = fit(PCA, utargPyr; pratio=0.99)
# pcomp_pyr = transform(Mpyr, utargPyr)
# utargPyr_r = reconstruct(Mpyr,pcomp_pyr)


#-------- load FS --------#
utargFS = load(dirutargFS * "utargFS.jld", "utargFS")[3:100,:]
Nfs = size(utargFS)[2]
for ci = 1:Nfs
    utargFS[:,ci] = utargFS[:,ci] .- mean(utargFS[:,ci])
end
# project to subspace
Mfs = fit(PCA, utargFS; pratio=0.99)
# pcomp_fs = transform(Mfs, utargFS)
# utargFS_r = reconstruct(Mfs,pcomp_fs)


#-------- load INH --------#
navg_list = [50, 100, 200, 400]
navg = navg_list[4]
Ne = 2500
Ncells = 5000

fname_usum_navg = dirsim * "usum_navg$(navg).jld"
usum_navg = load(fname_usum_navg, "usum")[3:100,:,1]
usumExc = usum_navg[:,1:Ne]
usumInh = usum_navg[:,Ne+1:Ncells]
for ci = 1:Ne
    usumExc[:,ci] .-= mean(usumExc[:,ci])
    usumInh[:,ci] .-= mean(usumInh[:,ci])
end
# project to subspace
Minh = fit(PCA, usumInh; pratio=0.8)
# pcomp_inh = transform(Minh, usumInh)
# usumInh_r = reconstruct(Minh, pcomp_inh)

#-------- load INH of balanced network --------#
fname_bal_usum_navg = dirsim_bal * "bal_usum_navg$(navg).jld"
usum_navg_bal = load(fname_bal_usum_navg, "usum")[3:100,:,1]
usumInh_bal = usum_navg_bal[:,Ne+1:Ncells]
for ci = 1:Ne
    usumInh_bal[:,ci] .-= mean(usumInh_bal[:,ci])
end


#-------- SVD: corr of FS and INH --------#
idx_pyr = .!isnan.(cor(utargFS, utargPyr, dims=1)[1,:])
idx_fs = .!isnan.(cor(usumInh, utargFS, dims=1)[1,:])

utargPyr_trim = utargPyr[:,idx_pyr]
utargFS_trim = utargFS[:,idx_fs]

# FS and Pyr
cor_fspyr = cor(utargFS_trim, utargPyr_trim, dims=1)
svd_fspyr = svd(cor_fspyr)

# FS and trained inh
cor_fsinh = cor(utargFS_trim, usumInh, dims=1)
svd_fsinh = svd(cor_fsinh)

# FS and untrained balanced 
cor_balanced = cor(utargFS_trim, usumInh_bal, dims=1)
svd_balanced = svd(cor_balanced)

# trained exc and inh
cor_excinh = cor(usumExc, usumInh, dims=1)
svd_excinh = svd(cor_excinh)


#-------- shared variance (FS and trained inh) --------#
svar_fsinh_cor, svar_fsinh_FS, svar_fsinh_INH,
svar_fsinh_FS_ksum, svar_fsinh_INH_ksum, 
svar_fsinh_FS_tsum, svar_fsinh_INH_tsum = sharedVar(svd_fsinh, utargFS_trim, usumInh)

#-------- shared variance (FS and balanced) --------#
svar_balanced_cor, svar_balanced_FS, svar_balanced_INH,
svar_balanced_FS_ksum, svar_balanced_INH_ksum, 
svar_balanced_FS_tsum, svar_balanced_INH_tsum = sharedVar(svd_balanced, utargFS_trim, usumInh_bal)

# #-------- shared variance (FS and Pyr) --------#
# svar_FSPYR_cor, svar_FS, svar_Pyr,
# svar_FS_ksum, svar_Pyr_ksum, 
# svar_FS_tsum, svar_Pyr_tsum = sharedVar(svd_fspyr, utargFS_trim, utargPyr_trim)

# #-------- shared variance (trained exc and inh) --------#
# svar_excinh_cor, svar_exc, svar_inh,
# svar_exc_ksum, svar_inh_ksum, 
# svar_exc_tsum, svar_inh_tsum = sharedVar(svd_excinh, usumExc, usumInh)


# plot - FS and trained inh
fname = "navg$(navg)_fsinh"
plt_sharedvar_trained(dirfig, fname, svar_fsinh_cor, svar_fsinh_FS, svar_fsinh_INH, 
svar_fsinh_FS_ksum, svar_fsinh_FS_tsum, svar_fsinh_INH_ksum, svar_fsinh_INH_tsum, 
svar_balanced_INH_ksum, svar_balanced_INH_tsum)

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
