function matchfs_nrep(utargFS, usumInh)

    cnt = 1

    nrep = 20
    Nfs = size(utargFS)[2]
    Ni = size(usumInh)[2]
    cells_matched = zeros(Nfs)
    usumInh_matched = zeros(100, Nfs)
    pcor_matched = zeros(Nfs)
    for nid = 1:Nfs
        ptmp = zeros(Ni,nrep)
        for ci = 1:Ni
            for repi = 1:nrep
                usumInh_slice = @view usumInh[:,ci,repi]
                utargFS_slice = @view utargFS[:,nid]
                ptmp[ci,repi] = cor(usumInh_slice, utargFS_slice)
                if isnan(ptmp[ci,repi])
                    println("ci ", ci)
                    cnt += 1
                end
            end
        end
        ptmp_avg = mean(ptmp, dims=2)[:]
        maxcell = argmax(ptmp_avg)
        cells_matched[nid] = maxcell
        usumInh_matched[:,nid] = mean(usumInh[:,maxcell,:], dims=2)[:]
        pcor_matched[nid] = ptmp_avg[maxcell]
    end

    return cells_matched, usumInh_matched, pcor_matched
end