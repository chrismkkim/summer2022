using Distributions
using PyCall
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Random
using SparseArrays
using JLD

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


include("param.jl")
include("genWeights.jl")
include("genPlasticWeights.jl")
include("convertWgtIn2Out.jl")
include("genTarget.jl")
include("genStim.jl")
include("genCellsTrained.jl")
include("runinitial.jl")
include("runtrain.jl")
# include("runtrain_tmp.jl")
include("runtest.jl")
include("funMovAvg.jl")
include("funCorrTarg.jl")
include("funCorrDecomp.jl")
include("funSample.jl")
include("funRollingAvg.jl")
include("runperformance_train.jl")
include("runperformance_test.jl")
include("calcWeights.jl")
include("replaceWp.jl")

include("genWeightsFfwd.jl")
include("genffwdRate.jl")
include("genffwdRateSub.jl")
include("runffwd.jl")
include("psth.jl")


dirdata = "/data/kimchm/data/dale/janelia/trained/basic/"
# dirdata = "/data/kimchm/data/dale/janelia/trained/rwr/"
dirtarget = "/data/kimchm/data/dale/janelia/s1alm/target/"
dirutarg = "/data/kimchm/data/dale/janelia/s1alm/target/utarg/"
dirutargFS = "/data/kimchm/data/dale/janelia/s1alm/target/utargFS/"
dirfig = "/data/kimchm/data/dale/janelia/figure/ffwd/"

if ~ispath(dirfig)
    mkpath(dirfig)
end

targetRate = load(dirtarget * "movingrate_Pyr1Hz.jld", "Pyr")
FSrate = load(dirtarget * "movingrate_FS.jld", "FS")

#-------- simulate inhibitory network with ffwd inputs from Pyr cells --------#
# inhibitory network network 
w0Index, w0Weights, nc0 = genWeightsFfwd(p)
# Pyr cells PSTH
ffwdRate = genffwdRate(p, targetRate)
# ffwdRate = genffwdRateSub(p, targetRate)

# simulate inhibitory network
nrep = 50
xsum = zeros(100, 5000, nrep)
maxTimes = round(Int,p.maxrate*p.train_time/1000)
times = zeros(p.Ncells,maxTimes, nrep)
ns = zeros(Int,p.Ncells, nrep)
for ii = 1:nrep
    println(ii)
    xsum[:,:,ii], times[:,:,ii], ns[:,ii] = runffwd(p,w0Index,w0Weights,nc0,ffwdRate)
end
usum = mean(xsum, dims=3)[:,:]

#-------- save usum (trial averaged) --------#
fname_usum = dirdata * "ffwd_usum.jld"
save(fname_usum,"usum", usum)

meanrate = psth(p, times, ns, nrep)
fname_meanrate = dirdata * "ffwd_meanrate.jld"
save(fname_meanrate,"meanrate", meanrate)


#-------- load FS --------#
utargFS = zeros(size(FSrate))
Nfs = size(utargFS)[1]
for ci = 1:Nfs
    utargFS[ci,:] = load(dirutargFS * "utarg$(ci).jld", "utarg")
end
utargFS = transpose(utargFS)[:,:]

#-------- load usum --------#
Ne = 2500
Ncells = 5000
# fname_usum = dirdata * "usum.jld"
# usum = load(fname_usum,"usum")
usumInh = usum[:,Ne+1:Ncells]

#-------- match FS and usum --------#
Ni = 2500
cells_matched = zeros(Nfs)
usumInh_matched = zeros(100, Nfs)
pcor_matched = zeros(Nfs)
for nid = 1:Nfs
    println(nid)
    ptmp = zeros(Ni)
    for ci = 1:Ni
        usumInh_slice = @view usumInh[1:100,ci]
        utargFS_slice = @view utargFS[1:100,nid]
        ptmp[ci] = cor(usumInh_slice, utargFS_slice)
    end
    maxcell = argmax(ptmp)
    cells_matched[nid] = maxcell
    usumInh_matched[:,nid] = usumInh[1:100,maxcell]
    pcor_matched[nid] = ptmp[maxcell]
end



figure(figsize=(3.5,3))
hist(pcor_matched, bins=50, range=(0,1), color="black", histtype="step")
xlabel("performance", fontsize=12)
ylabel("count", fontsize=12)
tight_layout()

savefig(dirfig * "ffwd_FSmatched_performance.png", dpi=300)




pcorSorted = reverse(sortperm(pcor_matched))

timev = 20*collect(1:100)


for ngrp = 1:12
    println(ngrp)
    figure(figsize=(10,10))
    for i = 1:25
        subplot(5,5,i)
        idx = (ngrp-1)*25 + i
        nid = pcorSorted[idx]
        plot(timev, utargFS[1:100,nid] .- mean(utargFS[1:100,nid]), c="red")
        plot(timev, usumInh_matched[1:100,nid] .- mean(usumInh_matched[1:100,nid]), c="black")
        title("cor: $(round(pcor_matched[nid], digits=2))", fontsize=12)
        xticks(fontsize=10)
        yticks(fontsize=10)
    end
    tight_layout()
    savefig(dirfig * "ffwd_FSmatched_shift-usum$(ngrp).png", dpi=300)
    close()
end




figure(figsize=(6,5))
ncell = 500
# axvspan(p.stim_on, p.stim_off, color="dodgerblue", alpha=0.3)
# for ci = 1:ncell
#     plot(times[ci,1:ns[ci]], ci*ones(ns[ci]), color="red", marker="o", ms=1, linestyle="")
# end
for ci = p.Ne+1:p.Ne+ncell
    plot(times[ci,1:ns[ci]], (ci-p.Ne)*ones(ns[ci]), color="blue", marker="o", ms=1, linestyle="")
end
xlim([0,p.train_time])
ylim([0,ncell])
xlabel("time (ms)", fontsize=12)
ylabel("neuron", fontsize=12)
tight_layout()

savefig(dirfig * "ffwd_raster.png", dpi=300)



# # select plastic weights to be trained
# wpWeightIn, wpWeightOut, wpIndexIn, wpIndexOut, wpIndexConvert, ncpIn, ncpOut = genPlasticWeights(p, w0Index, nc0, ns0, matchedCells)
# stim = genStim(p)


# # load targets
# xtarg = transpose(load(dirutarg * "utarg1Hz.jld", "utarg"))[:,:]


# #----------- save files --------------#
# fname_param = dirdata * "p.jld"
# fname_w0Index = dirdata * "w0Index.jld"
# fname_w0Weights = dirdata * "w0Weights.jld"
# fname_nc0 = dirdata * "nc0.jld"
# fname_wpIndexIn = dirdata * "wpIndexIn.jld"
# fname_wpIndexOut = dirdata * "wpIndexOut.jld"
# fname_wpIndexConvert = dirdata * "wpIndexConvert.jld"
# fname_ncpIn = dirdata * "ncpIn.jld"
# fname_ncpOut = dirdata * "ncpOut.jld"
# fname_uavg = dirdata * "uavg.jld"
# fname_stim = dirdata * "stim.jld"
# fname_almOrd = dirdata * "almOrd.jld"
# fname_matchedCells = dirdata * "matchedCells.jld"


# save(fname_param,"p", p)
# save(fname_w0Index,"w0Index", w0Index)
# save(fname_w0Weights,"w0Weights", w0Weights)
# save(fname_nc0,"nc0", nc0)
# save(fname_wpIndexIn,"wpIndexIn", wpIndexIn)
# save(fname_wpIndexOut,"wpIndexOut", wpIndexOut)
# save(fname_wpIndexConvert,"wpIndexConvert", wpIndexConvert)
# save(fname_ncpIn,"ncpIn", ncpIn)
# save(fname_ncpOut,"ncpOut", ncpOut)
# save(fname_uavg,"uavg", uavg)
# save(fname_stim,"stim", stim)
# save(fname_almOrd,"almOrd", almOrd)
# save(fname_matchedCells,"matchedCells", matchedCells)

# #----------- run train --------------#
# wpWeightIn, wpWeightOut, performance_test, frac_switch_exc, frac_switch_inh, 
# wp_rowsum_exc, wp_rowsum_inh, wp_std_exc, wp_std_inh,
# error = runtrain(dirdata,p,w0Index,w0Weights,nc0,stim,xtarg,wpIndexIn,wpIndexOut,wpIndexConvert,wpWeightIn,wpWeightOut,ncpIn,ncpOut,almOrd,matchedCells)


# #----------- save files --------------#
# fname_wpWeightIn = dirdata * "wpWeightIn.jld"
# fname_wpWeightOut = dirdata * "wpWeightOut.jld"
# fname_performance_test = dirdata * "_performance_test.jld"
# fname_frac_switch_exc = dirdata * "_frac_switch_exc.jld"
# fname_frac_switch_inh = dirdata * "_frac_switch_inh.jld"

# save(fname_wpWeightIn,"wpWeightIn", wpWeightIn)
# save(fname_wpWeightOut,"wpWeightOut", wpWeightOut)
# save(fname_performance_test,"performance_test", performance_test)
# save(fname_frac_switch_exc,"frac_switch_exc", frac_switch_exc)
# save(fname_frac_switch_inh,"frac_switch_inh", frac_switch_inh)
