function matchfs(utargFS, usumInh)

    Nfs = size(utargFS)[2]
    Ni = 2500
    cells_matched = zeros(Nfs)
    usumInh_matched = zeros(100, Nfs)
    pcor_matched = zeros(Nfs)
    for nid = 1:Nfs
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

    return cells_matched, usumInh_matched, pcor_matched
end