function runmatchfs(dirsim, navg, fs_rate, lickRL)

    Ncells = 5000
    Ne = 2500
    
    if lickRL == "right"
        fname_rate_navg = dirsim * "rate_R_navg$(navg).jld"
    elseif lickRL == "left"
        fname_rate_navg = dirsim * "rate_L_navg$(navg).jld"
    end
    rate_navg = load(fname_rate_navg, "rate")
    rateInh = rate_navg[:,Ne+1:Ncells,:]
    
    # load usum_bal (untrained balanced network)
    # fname_bal_usum_navg = dirsim_bal * "bal_usum_navg$(navg).jld"
    # fname_bal_rate_navg = dirsim_bal * "bal_rate_navg$(navg).jld"
    # usum_navg_bal = load(fname_bal_usum_navg, "usum")
    # rate_navg_bal = load(fname_bal_rate_navg, "rate")
    # usumInh_bal = usum_navg_bal[:,Ne+1:Ncells,:]
    # rateInh_bal = rate_navg_bal[:,Ne+1:Ncells,:]
    
    # rateInh = rateInh[:,1:end .!= 171,:] # remove inhibitory cell 171
    # rateInh_bal = rateInh_bal[:,1:end .!=1188,:]
    
    # cells_matched_usum, usumInh_matched, expvar_matched_usum, pcor_matched_usum = matchfs_noreplacement(utargFS, usumInh)
    cells_matched, rateInh_matched, expvar_matched, pcor_matched = matchfs_noreplacement(fs_rate, rateInh)
    # cells_matched_bal, rateInh_matched_bal, expvar_matched_bal, pcor_matched_bal = matchfs_noreplacement(rtargFS, rateInh_bal)
    
    return expvar_matched, pcor_matched, cells_matched

end