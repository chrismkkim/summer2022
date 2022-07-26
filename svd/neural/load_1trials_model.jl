function load_exc_1trials(dirsim, rate_diff)

    # inhibitory neurons only
    Ncells = 5000
    Ne = 2500
    ntrials = 100

    ns_size = size(load(dirsim * "ns_R_1.jld")["ns"])
    ns_R = zeros(ns_size[1],ntrials)
    ns_L = zeros(ns_size[1],ntrials)
    for ii = 1:ntrials
        ns_R[:,ii] = load(dirsim * "ns_R_$(ii).jld")["ns"]
        ns_L[:,ii] = load(dirsim * "ns_L_$(ii).jld")["ns"]
    end
    ns_R = ns_R[1:Ne,:]
    ns_L = ns_L[1:Ne,:]

    # order inhibitory neurons by mean rate difference
    cellsOrd = sortperm(rate_diff)

    # correlation: 1 trial vs. psth
    pearcorr = zeros(ntrials)
    rate_diff_1trial = (ns_R - ns_L) / 2
    for ii = 1:ntrials
        rate_diff_ord = rate_diff[cellsOrd]
        rate_diff_1trial_ord = rate_diff_1trial[cellsOrd,ii]
        pearcorr[ii] = cor(rate_diff_ord, rate_diff_1trial_ord)
    end

    return pearcorr

end


function load_inh_1trials(dirsim, rate_diff)

    # inhibitory neurons only
    Ncells = 5000
    Ne = 2500
    ntrials = 100

    ns_size = size(load(dirsim * "ns_R_1.jld")["ns"])
    ns_R = zeros(ns_size[1],ntrials)
    ns_L = zeros(ns_size[1],ntrials)
    for ii = 1:ntrials
        ns_R[:,ii] = load(dirsim * "ns_R_$(ii).jld")["ns"]
        ns_L[:,ii] = load(dirsim * "ns_L_$(ii).jld")["ns"]
    end
    ns_R = ns_R[Ne+1:Ncells,:]
    ns_L = ns_L[Ne+1:Ncells,:]

    # order inhibitory neurons by mean rate difference
    cellsOrd = sortperm(rate_diff)

    # correlation: 1 trial vs. psth
    pearcorr = zeros(ntrials)
    rate_diff_1trial = (ns_R - ns_L) / 2
    for ii = 1:ntrials
        rate_diff_ord = rate_diff[cellsOrd]
        rate_diff_1trial_ord = rate_diff_1trial[cellsOrd,ii]
        pearcorr[ii] = cor(rate_diff_ord, rate_diff_1trial_ord)
    end

    return pearcorr

end
