function runmatchfs_RL(dirsim, navg, targ_R, targ_L, nrep)

    Ncells = 5000
    Ne = 2500
    
    fname_rate_navg_R = dirsim * "rate_R_navg$(navg).jld"
    fname_rate_navg_L = dirsim * "rate_L_navg$(navg).jld"
    rate_R = load(fname_rate_navg_R, "rate")
    rate_L = load(fname_rate_navg_L, "rate")
    rateInh_R = rate_R[:,Ne+1:Ncells,:]
    rateInh_L = rate_L[:,Ne+1:Ncells,:]
    

    # nrep = 20
    Nfs = size(targ_R)[2]
    Ni = size(rateInh_R)[2]
    cells_inh = collect(1:Ni)    
    cells_matched = zeros(Nfs)
    expvar_R = zeros(Nfs)
    expvar_L = zeros(Nfs)
    pcor_R = zeros(Nfs)
    pcor_L = zeros(Nfs)

    # loop FS cells
    for nid = 1:Nfs
        Ni = length(cells_inh)
        vtmp_R = zeros(Ni,nrep)
        vtmp_L = zeros(Ni,nrep)
        ptmp_R = zeros(Ni,nrep)        
        ptmp_L = zeros(Ni,nrep)
        # scan Inh neurons
        for (ind, ci) in enumerate(cells_inh)
            for repi = 1:nrep
                rateInh_R_slice = @view rateInh_R[:,ci,repi]
                rateInh_L_slice = @view rateInh_L[:,ci,repi]
                targ_R_slice = @view targ_R[:,nid]          
                targ_L_slice = @view targ_L[:,nid]                
                vtmp_R[ind,repi] = explainedVar(rateInh_R_slice, targ_R_slice)
                vtmp_L[ind,repi] = explainedVar(rateInh_L_slice, targ_L_slice)
                ptmp_R[ind,repi] = cor(rateInh_R_slice, targ_R_slice)
                ptmp_L[ind,repi] = cor(rateInh_L_slice, targ_L_slice)
            end
        end        
        vtmp_R_avg = mean(vtmp_R, dims=2)[:]
        vtmp_L_avg = mean(vtmp_L, dims=2)[:]
        ptmp_R_avg = mean(ptmp_R, dims=2)[:]
        ptmp_L_avg = mean(ptmp_L, dims=2)[:]

        # find the best matching Inh neuron
        maxidx = argmax(vtmp_R_avg + vtmp_L_avg)
        # maxidx = argmax(ptmp_R_avg + ptmp_L_avg)
        # best performance
        expvar_R[nid] = vtmp_R_avg[maxidx]
        expvar_L[nid] = vtmp_L_avg[maxidx]
        pcor_R[nid] = ptmp_R_avg[maxidx]
        pcor_L[nid] = ptmp_L_avg[maxidx]

        # best matching Inh neuron id
        maxcell = cells_inh[maxidx]
        cells_matched[nid] = maxcell

        # remove the matched cell
        filter!(x->x!=maxcell, cells_inh)
    end

    return cells_matched, pcor_R, pcor_L, expvar_R, expvar_L

    # # cells_matched_usum, usumInh_matched, expvar_matched_usum, pcor_matched_usum = matchfs_noreplacement(utargFS, usumInh)
    # cells_matched, rateInh_matched, expvar_matched, pcor_matched = matchfs_noreplacement(fs_rate, rateInh)
    # # cells_matched_bal, rateInh_matched_bal, expvar_matched_bal, pcor_matched_bal = matchfs_noreplacement(rtargFS, rateInh_bal)
    
    # return expvar_matched, pcor_matched, cells_matched

end