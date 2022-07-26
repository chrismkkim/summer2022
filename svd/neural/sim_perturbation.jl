function sim_perturbation(dirdata, rate_R_avg, rate_L_avg, pyr_R, pyr_L, pert_on, pert_off, subpoptype, ncell, stimtype)

    npyr = size(pyr_R)[2]
    nexc = 2500
    ntrial = 15 # 15
    nsession = 10 # 25, 18
    # ncell = 50

    codingdim_PR_list = zeros(nsession, ntrial, 151)
    codingdim_PL_list = zeros(nsession, ntrial, 151)
    codingdim_R_list = zeros(nsession, ntrial, 151)
    codingdim_L_list = zeros(nsession, ntrial, 151)

    homdim_PR_list = zeros(nsession, ntrial, 151)
    homdim_PL_list = zeros(nsession, ntrial, 151)
    homdim_R_list = zeros(nsession, ntrial, 151)
    homdim_L_list = zeros(nsession, ntrial, 151)

    codingdimInh_PR_list = zeros(nsession, ntrial, 151)
    codingdimInh_PL_list = zeros(nsession, ntrial, 151)
    codingdimInh_R_list = zeros(nsession, ntrial, 151)
    codingdimInh_L_list = zeros(nsession, ntrial, 151)

    # There are nsession subpopulations. Select random ncells for each subpopulation. 
    # Random.seed!(333) # fixe the random seed to sample same set of neurons
    subpop = Vector{Array{Int,1}}()
    if subpoptype == "trained"
        cellshuffled = shuffle(collect(1:npyr))
    elseif subpoptype == "exc"
        cellshuffled = shuffle(collect(1:nexc))
    end
    for popi = 1:nsession
        start_idx = (popi-1)*ncell + 1
        end_idx = popi*ncell
        idx = collect(start_idx:end_idx)
        cellidx = cellshuffled[idx]
        push!(subpop, cellidx)
    end

    lickRL = "right"
    p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
    wpWeightOut, ncpIn, ncpOut, stim_R, almOrd, matchedCells,
    ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

    lickRL = "left"
    p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
    wpWeightOut, ncpIn, ncpOut, stim_L, almOrd, matchedCells,
    ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

    for repi = 1:ntrial
        println(repi)

        #----- lick right -----#
        println("lick right")
        isperturbed = false
        stim_training = copy(stim_R)
        stim_perturb = copy(stim_L)
        usum_R, times_R, ns_R = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off, stimtype)

        #----- lick left -----#
        println("lick left")
        isperturbed = false
        stim_training = copy(stim_L)
        stim_perturb = copy(stim_R)
        usum_L, times_L, ns_L = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[2],wpWeightFfwd[2],isperturbed, pert_on, pert_off, stimtype)
        
        #----- perturbation -----#
        println("perturbation - lick right")
        # trigger learned activity with stim_R, perturb with stim_L
        isperturbed = true            
        stim_training = copy(stim_R)
        stim_perturb = copy(stim_L)    
        usum_P, times_PR, ns_PR = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off, stimtype)

        println("perturbation - lick left")
        # trigger learned activity with stim_L, perturb with stim_R
        isperturbed = true            
        stim_training = copy(stim_L)
        stim_perturb = copy(stim_R)    
        usum_P, times_PL, ns_PL = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[2],wpWeightFfwd[2],isperturbed, pert_on, pert_off, stimtype)
  
        for popi = 1:nsession
            if subpoptype == "trained"
                matchedCells_subpop = matchedCells[subpop[popi]]
            elseif subpoptype == "exc"
                matchedCells_subpop = collect(1:nexc)[subpop[popi]]
            end
            inhCells = collect(nexc+1:2*nexc)

            rate_PR = spk2rate(p, matchedCells_subpop, times_PR, ns_PR)
            rate_PL = spk2rate(p, matchedCells_subpop, times_PL, ns_PL)
            rate_R = spk2rate(p, matchedCells_subpop, times_R, ns_R)
            rate_L = spk2rate(p, matchedCells_subpop, times_L, ns_L)
            rate_M = mean((rate_R + rate_L)/2, dims=1)[:]
            rate_M = zeros(size(rate_M)) # do not remove mean rate

            # coding dimension --- excitatory
            pyr_codingdim = mean(rate_R_avg - rate_L_avg, dims=1)[matchedCells_subpop]
            pyr_codingdim = pyr_codingdim / sqrt(mean(pyr_codingdim.^2))

            codingdim_PR_list[popi, repi,:] = mean((rate_PR .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
            codingdim_PL_list[popi, repi,:] = mean((rate_PL .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
            codingdim_R_list[popi, repi,:] = mean((rate_R .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
            codingdim_L_list[popi, repi,:] = mean((rate_L .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]

            # homogeneous dimension --- excitatory
            pyr_homdim = ones(size(rate_R)[2])
            pyr_homdim = pyr_homdim / sqrt(mean(pyr_homdim.^2))

            homdim_PR_list[popi, repi,:] = mean((rate_PR .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            homdim_PL_list[popi, repi,:] = mean((rate_PL .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            homdim_R_list[popi, repi,:] = mean((rate_R .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            homdim_L_list[popi, repi,:] = mean((rate_L .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            
            # coding dimension --- untrained inhibitory
            rateInh_PR = spk2rate(p, inhCells, times_PR, ns_PR)
            rateInh_PL = spk2rate(p, inhCells, times_PL, ns_PL)
            rateInh_R = spk2rate(p, inhCells, times_R, ns_R)
            rateInh_L = spk2rate(p, inhCells, times_L, ns_L)
            rateInh_M = mean((rateInh_R + rateInh_L)/2, dims=1)[:]

            fs_codingdim = mean(rate_R_avg - rate_L_avg, dims=1)[inhCells]
            fs_codingdim = fs_codingdim / sqrt(mean(fs_codingdim.^2))

            codingdimInh_PR_list[popi, repi, :] = mean((rateInh_PR .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]
            codingdimInh_PL_list[popi, repi, :] = mean((rateInh_PL .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]
            codingdimInh_R_list[popi, repi, :] = mean((rateInh_R .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]
            codingdimInh_L_list[popi, repi, :] = mean((rateInh_L .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]

        end

    end

    # excitatory
    codingdim_PR = dropdims(mean(codingdim_PR_list, dims=2), dims=2)
    codingdim_PL = dropdims(mean(codingdim_PL_list, dims=2), dims=2)
    codingdim_R = dropdims(mean(codingdim_R_list, dims=2), dims=2)
    codingdim_L = dropdims(mean(codingdim_L_list, dims=2), dims=2)
    homdim_PR = dropdims(mean(homdim_PR_list, dims=2), dims=2)
    homdim_PL = dropdims(mean(homdim_PL_list, dims=2), dims=2)
    homdim_R = dropdims(mean(homdim_R_list, dims=2), dims=2)
    homdim_L = dropdims(mean(homdim_L_list, dims=2), dims=2)

    # inhibitory
    codingdimInh_PR = dropdims(mean(codingdimInh_PR_list, dims=2), dims=2)
    codingdimInh_PL = dropdims(mean(codingdimInh_PL_list, dims=2), dims=2)
    codingdimInh_R = dropdims(mean(codingdimInh_R_list, dims=2), dims=2)
    codingdimInh_L = dropdims(mean(codingdimInh_L_list, dims=2), dims=2)


    return codingdim_PR, codingdim_PL, codingdim_R, codingdim_L, 
    homdim_PR, homdim_PL, homdim_R, homdim_L, 
    codingdimInh_PR, codingdimInh_PL, codingdimInh_R, codingdimInh_L,
    codingdim_PR_list, codingdim_PL_list, codingdim_R_list, codingdim_L_list, 
    homdim_PR_list, homdim_PL_list, homdim_R_list, homdim_L_list
    
    end