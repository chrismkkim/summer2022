using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD
using MAT

matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm
gridspec = matplotlib.gridspec


include("load_files.jl")
include("sim_selectivity.jl")
include("sim_perturbation.jl")
# include("sim_perturbation_subpop.jl")
# include("sim_perturbation_lickLeft.jl")
include("sim_perturbation_savedata.jl")
include("runtestdrive.jl")
include("runtestperturb.jl")
# include("runtestperturb_with_stim.jl")
# include("runtestperturb_lickLeft.jl")
include("runperformance_nrep.jl")
include("funMovAvg.jl")
include("funMovAvg2D.jl")
include("spk2rate.jl")
# include("runtest.jl")
# include("reconstructTarget.jl")
# include("runinitial.jl")
# include("runtest.jl")
# include("runtest_scaling.jl")
# include("runperformance.jl")
# include("funRollingAvg.jl")
# include("funSample.jl")
# include("loglogfit.jl")
# include("meanRates.jl")
# include("psth.jl")

include("figure/plt_pyr_rightleft.jl")
include("figure/plt_pyr_pcor.jl")
# include("figure/plt_pyr_raster.jl")
include("figure/plt_selectivity_pop.jl")
include("figure/plt_perturbation.jl")
# include("figure/plt_perturbation_lickRight.jl")
# include("figure/plt_perturbation_lickLeft.jl")

L_list = [2.0, 4.0, 6.0, 8.0, 10.0]
ll = 1 # L = 2.0

taup_list = [25.0, 50.0, 75.0, 100.0, 150.0]
nn = 5

navg_list = [50, 100, 200, 400]


# dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/saved/Lffwd3/lam2/"
dirdata = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/wffwd/Lffwd5/"
dirsim = "/data/kimchm/data/dale/janelia/trained/selectivity/drive/sim/"
dirperturb = "/data/kimchm/data/dale/janelia/trained/perturb/training_stimulus/"
dirtarget_selectivity = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/"
dirutarg_pyr_lickright = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickright/"
dirutarg_pyr_lickleft = "/data/kimchm/data/dale/janelia/s1alm/target/selectivity/utarg/pyr/lickleft/"
dirfig = "/data/kimchm/data/dale/janelia/figure/selectivity/drive/"


if ~ispath(dirfig)
    mkpath(dirfig)
end

# load files
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, "right")

# load targets: 
# --- trial-averaged firing rates 
# pyr_R: [time x neurons]
# pyr_L: [time x neurons]
# fs_R: [time x neurons]
# fs_L: [time x neurons]
pyr_R = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickright.jld", "Pyr"))[1:100,:]
pyr_L = transpose(load(dirtarget_selectivity * "movingrate_Pyr1Hz_lickleft.jld", "Pyr"))[1:100,:]
fs_R = transpose(load(dirtarget_selectivity * "movingrate_FS_lickright.jld", "FS"))[1:100,:]
fs_L = transpose(load(dirtarget_selectivity * "movingrate_FS_lickleft.jld", "FS"))[1:100,:]

fname_rate_R = dirsim * "rate_R_navg400.jld"
fname_rate_L = dirsim * "rate_L_navg400.jld"
rate_R = load(fname_rate_R, "rate")[:,:,1]
rate_L = load(fname_rate_L, "rate")[:,:,1]


stimtype = "cue" # cue, inh
ncell = 250 # 50, 250


#------------------------------------------------#
#-------- run perturbation experiments ----------#
#------------------------------------------------#
pert_on = 1500.0
pert_off = pert_on + 200
subpoptype = "exc"
codingdim_PR, codingdim_PL, codingdim_R, codingdim_L, 
homdim_PR, homdim_PL, homdim_R, homdim_L, 
codingdimInh_PR, codingdimInh_PL, codingdimInh_R, codingdimInh_L,
codingdim_PR_list, codingdim_PL_list, codingdim_R_list, codingdim_L_list, 
homdim_PR_list, homdim_PL_list, homdim_R_list, homdim_L_list = sim_perturbation(dirdata, rate_R, rate_L, pyr_R, pyr_L, pert_on, pert_off, subpoptype, ncell, stimtype)

#---------------------------------------------#
#-------- plot perturbation results ----------#
#---------------------------------------------#
perturbtype = "lickRight"
plt_perturbation(dirfig, 
codingdim_PR, codingdim_PL, codingdim_R, codingdim_L, 
homdim_PR, homdim_PL, homdim_R, homdim_L, 
codingdimInh_PR, codingdimInh_PL, codingdimInh_R, codingdimInh_L,
codingdim_PR_list, codingdim_PL_list, codingdim_R_list, codingdim_L_list, 
homdim_PR_list, homdim_PL_list, homdim_R_list, homdim_L_list,
pert_on, pert_off, perturbtype)


perturbtype = "lickLeft"
plt_perturbation(dirfig, 
codingdim_PR, codingdim_PL, codingdim_R, codingdim_L, 
homdim_PR, homdim_PL, homdim_R, homdim_L, 
codingdimInh_PR, codingdimInh_PL, codingdimInh_R, codingdimInh_L,
codingdim_PR_list, codingdim_PL_list, codingdim_R_list, codingdim_L_list, 
homdim_PR_list, homdim_PL_list, homdim_R_list, homdim_L_list,
pert_on, pert_off, perturbtype)

#---------------------------------------------------#
#-------- save perturbation data for Ran  ----------#
#---------------------------------------------------#
# ifsave = true
# subpoptype = "exc"
# rate_R_perturb, rate_R, rate_L = sim_perturbation_savedata(p, dirdata, dirperturb, rate_R, rate_L, pyr_R, pyr_L, pert_on, pert_off, subpoptype, ifsave, ncell, stimtype)

