function sim_perturbation(dirdata, pyr_R, pyr_L, pert_on, pert_off)

    nrep = 10
    exc_Rp_list = zeros(nrep, 151)
    inh_Rp_list = zeros(nrep, 151)
    exc_R_list = zeros(nrep, 151)
    inh_R_list = zeros(nrep, 151)
    trained_Rp_list = zeros(nrep, 151)
    trained_R_list = zeros(nrep, 151)
    untrained_Rp_list = zeros(nrep, 151)
    untrained_R_list = zeros(nrep, 151)

    codingdim_Rp_list = zeros(nrep, 151)
    codingdim_R_list = zeros(nrep, 151)
    codingdim_L_list = zeros(nrep, 151)

    for repi = 1:nrep
        println(repi)

        #----- lick right -----#
        println("Loading files")
        lickRL = "right"
        p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
        wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
        ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
        # usum_R, times_R, ns_R = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
        isperturbed = false
        usum_R, times_R, ns_R = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off)

        isperturbed = true
        usum_R_perturb, times_R_perturb, ns_R_perturb = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off)

        #----- lick left -----#
        println("Loading files")
        lickRL = "left"
        p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
        wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
        ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)
        # usum_L, times_L, ns_L = runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)
        isperturbed = false
        usum_L, times_L, ns_L = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate[2],wpWeightFfwd[2],isperturbed, pert_on, pert_off)
        
        # coding dimension
        pyr_codingdim = mean(pyr_R - pyr_L, dims=1)[:][almOrd]
        pyr_codingdim = pyr_codingdim / sqrt(mean(pyr_codingdim.^2))
        
        rate_R_perturb = spk2rate(p, matchedCells, times_R_perturb, ns_R_perturb)
        rate_R = spk2rate(p, matchedCells, times_R, ns_R)
        rate_L = spk2rate(p, matchedCells, times_L, ns_L)
        rate_M = mean((rate_R + rate_L)/2, dims=1)[:]
    
        codingdim_R_perturb = mean((rate_R_perturb .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
        codingdim_R = mean((rate_R .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
        codingdim_L = mean((rate_L .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
        
        # mean rates
        exc_R_perturb = spk2rate(p, collect(1:p.Ne), times_R_perturb, ns_R_perturb)
        inh_R_perturb = spk2rate(p, collect(p.Ne+1:p.Ncells), times_R_perturb, ns_R_perturb)
        exc_R = spk2rate(p, collect(1:p.Ne), times_R, ns_R)
        inh_R = spk2rate(p, collect(p.Ne+1:p.Ncells), times_R, ns_R)
        untrained_R_perturb = spk2rate(p, deleteat!(collect(1:p.Ne), sort(matchedCells)), times_R_perturb, ns_R_perturb)
        untrained_R = spk2rate(p, deleteat!(collect(1:p.Ne), sort(matchedCells)), times_R, ns_R)

        exc_Rp_list[repi,:] = mean(exc_R_perturb, dims=2)[:]
        inh_Rp_list[repi,:] = mean(inh_R_perturb, dims=2)[:]
        exc_R_list[repi,:] = mean(exc_R, dims=2)[:]
        inh_R_list[repi,:] = mean(inh_R, dims=2)[:]
        trained_Rp_list[repi,:] = mean(rate_R_perturb, dims=2)[:]
        trained_R_list[repi,:] = mean(rate_R, dims=2)[:]
        untrained_Rp_list[repi,:] = mean(untrained_R_perturb, dims=2)[:]
        untrained_R_list[repi,:] = mean(untrained_R, dims=2)[:]

        codingdim_Rp_list[repi,:] = codingdim_R_perturb
        codingdim_R_list[repi,:] = codingdim_R
        codingdim_L_list[repi,:] = codingdim_L
    end

    codingdim_Rp_mean = mean(codingdim_Rp_list, dims=1)[:]
    codingdim_R_mean = mean(codingdim_R_list, dims=1)[:]
    codingdim_L_mean = mean(codingdim_L_list, dims=1)[:]

    exc_Rp_mean = mean(exc_Rp_list, dims=1)[:]
    inh_Rp_mean = mean(inh_Rp_list, dims=1)[:]
    exc_R_mean = mean(exc_R_list, dims=1)[:]
    inh_R_mean = mean(inh_R_list, dims=1)[:]
    trained_Rp_mean = mean(trained_Rp_list, dims=1)[:]
    trained_R_mean = mean(trained_R_list, dims=1)[:]
    untrained_Rp_mean = mean(untrained_Rp_list, dims=1)[:]
    untrained_R_mean = mean(untrained_R_list, dims=1)[:]


    return codingdim_R_mean, codingdim_L_mean, codingdim_Rp_mean, untrained_R_mean, untrained_Rp_mean, exc_R_mean, inh_R_mean, exc_Rp_mean, inh_Rp_mean, trained_R_mean, trained_Rp_mean
    
    end