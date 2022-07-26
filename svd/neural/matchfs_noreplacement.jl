function matchfs_noreplacement(utarg, usum)

    # cnt = 1

    nrep = 20
    Nfs = size(utarg)[2]
    Ni = size(usum)[2]
    cells_inh = collect(1:Ni)    
    cells_matched = zeros(Nfs)
    usum_matched = zeros(100, Nfs)
    expvar_matched = zeros(Nfs)
    pcor_matched = zeros(Nfs)
    # loop FS cells
    for nid = 1:Nfs
        Ni = size(cells_inh)[1]
        vtmp = zeros(Ni,nrep)
        ptmp = zeros(Ni,nrep)
        # scan Inh neurons
        for (ind, ci) in enumerate(cells_inh)
            for repi = 1:nrep
                usum_slice = @view usum[:,ci,repi]
                utarg_slice = @view utarg[:,nid]                
                vtmp[ind,repi] = explainedVar(usum_slice, utarg_slice)
                ptmp[ind,repi] = cor(usum_slice, utarg_slice)
            end
        end        
        vtmp_avg = mean(vtmp, dims=2)[:]
        ptmp_avg = mean(ptmp, dims=2)[:]

        # find the best matching Inh neuron
        maxidx = argmax(vtmp_avg)
        # best performance
        expvar_matched[nid] = vtmp_avg[maxidx]
        pcor_matched[nid] = ptmp_avg[maxidx]

        # best matching Inh neuron id
        maxcell = cells_inh[maxidx]
        cells_matched[nid] = maxcell
        usum_matched[:,nid] = mean(usum[:,maxcell,:], dims=2)[:]

        # remove the matched cell
        filter!(x->x!=maxcell, cells_inh)
    end

    return cells_matched, usum_matched, expvar_matched, pcor_matched
end