function sim_selectivity(dirdata)


#----- lick right -----#
println("Loading files")
lickRL = "right"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
# usum_R, times_R, ns_R = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
ubal, uplas, uext, times_R, ns_R = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1])


#----- lick left -----#
println("Loading files")
lickRL = "left"
p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
# usum_L, times_L, ns_L = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
ubal, uplas, uext, times_L, ns_L = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2])

rate_R = zeros(p.Ncells)
rate_L = zeros(p.Ncells)
for ii = 1:p.Ncells
    rate_R[ii] = sum(times_R[ii, 1:ns_R[ii]] .> p.stim_off) / p.train_duration * 1000
    rate_L[ii] = sum(times_L[ii, 1:ns_L[ii]] .> p.stim_off) / p.train_duration * 1000
end

rdiff_exc = rate_L[1:p.Ne] - rate_R[1:p.Ne]
rdiff_inh = rate_L[p.Ne+1:end] - rate_R[p.Ne+1:end]

exc_sorted = rdiff_exc[sortperm(rdiff_exc)]
inh_sorted = rdiff_inh[sortperm(rdiff_inh)]

return rdiff_exc, rdiff_inh

end