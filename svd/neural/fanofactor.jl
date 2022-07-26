function fanofactor(dirsim, p)

    # load spikes
    navg = 50
    Ncells = 5000
    tmp = load(dirsim * "times_R_1.jld", "times")
    times_R = zeros(size(tmp)[1], size(tmp)[2], navg)
    times_L = zeros(size(tmp)[1], size(tmp)[2], navg)
    ns_R = zeros(Int,Ncells, navg)
    ns_L = zeros(Int,Ncells, navg)
    for avgi = 1:navg
        times_R[:,:,avgi] = load(dirsim * "times_R_$(avgi).jld", "times")
        times_L[:,:,avgi] = load(dirsim * "times_L_$(avgi).jld", "times")
        ns_R[:,avgi] = load(dirsim * "ns_R_$(avgi).jld", "ns")
        ns_L[:,avgi] = load(dirsim * "ns_L_$(avgi).jld", "ns")
    end

    # convert spike time to spike counts
    wbin = 1000 # ms
    tbin_R, spkcnt_R = spktime2cnt(p,wbin,navg,times_R,ns_R)
    tbin_L, spkcnt_L = spktime2cnt(p,wbin,navg,times_L,ns_L)

    # calculate Fano factor in t > stim_off
    #--- lick right ---#
    idx_R = tbin_R .> p.stim_off
    tbin_R_trim = tbin_R[idx_R]
    spkcnt_R_trim = spkcnt_R[:,idx_R,:]
    spkcnt_R_merge = sum(spkcnt_R_trim, dims=2)[:,1,:]
    ff_R = calcFano(p, spkcnt_R_merge)

    #--- lick left ---#
    idx_L = tbin_L .> p.stim_off
    tbin_L_trim = tbin_L[idx_L]
    spkcnt_L_trim = spkcnt_L[:,idx_L,:]
    spkcnt_L_merge = sum(spkcnt_L_trim, dims=2)[:,1,:]
    ff_L = calcFano(p, spkcnt_L_merge)

    return ff_R, ff_L

end


function spktime2cnt(p,wbin,ntest,times,ns)
    
    nbin = Int(p.train_time/wbin)
    tbin = collect(1:nbin)*wbin
    spkcnt = zeros(p.Ncells,nbin,ntest)
    for itest = 1:ntest
        for ci = 1:p.Ncells
            for ti = 1:Int(ns[ci,itest])
                idx = floor(Int,times[ci,ti,itest]/wbin)+1 # add 1 to start at index 1
                if idx <= nbin # idx = nbin+1 can happen. This is julia's numerical error.
                    spkcnt[ci,idx,itest] += 1 
                end
            end
        end
    end

    return tbin, spkcnt

end



function calcFano(p,spkcnt)

    spkcnt_ff    = zeros(p.Ncells)
    spkcnt_mean  = zeros(p.Ncells)
    spkcnt_var   = zeros(p.Ncells)

    # cross-trial mean, var and fano factor of spike counts
    for ci=1:p.Ncells # loop over cells
            spkcnt_mean[ci]  = mean(spkcnt[ci,:])
            spkcnt_var[ci]   = var(spkcnt[ci,:])
            # fano factor
            if spkcnt_mean[ci] > 0. 
                spkcnt_ff[ci] = spkcnt_var[ci]/spkcnt_mean[ci];
            else
                spkcnt_ff[ci] = 0.;
            end
    end    

    return spkcnt_ff    

end
    


# function calcFano(p,tbin,spkcnt)

#     spkcnt_ff    = zeros(p.Ncells,length(tbin))
#     spkcnt_mean  = zeros(p.Ncells,length(tbin))
#     spkcnt_var   = zeros(p.Ncells,length(tbin))

#     # cross-trial mean, var and fano factor of spike counts
#     for ci=1:p.Ncells # loop over cells
#         for j=1:length(tbin) # loop over bins
#             spkcnt_mean[ci,j]  = mean(spkcnt[ci,j,:]);
#             spkcnt_var[ci,j]   = var(spkcnt[ci,j,:]);
#             # fano factor
#             if spkcnt_mean[ci,j] > 0. 
#                 spkcnt_ff[ci,j] = spkcnt_var[ci,j]/spkcnt_mean[ci,j];
#             else
#                 spkcnt_ff[ci,j] = 0.;
#             end
#         end
#     end    

#     return spkcnt_ff    

# end
    
