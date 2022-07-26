function sim_perturbation_subpop(dirdata, rate_R_avg, rate_L_avg, pyr_R, pyr_L, pert_on, pert_off, subpoptype)

    npyr = size(pyr_R)[2]
    nexc = 2500
    ntrial = 15 # 15
    nsession = 10 # 25, 18
    ncell = 50

    codingdim_Rp_list = zeros(nsession, ntrial, 151)
    codingdim_R_list = zeros(nsession, ntrial, 151)
    codingdim_L_list = zeros(nsession, ntrial, 151)

    homdim_Rp_list = zeros(nsession, ntrial, 151)
    homdim_R_list = zeros(nsession, ntrial, 151)

    codingdimInh_Rp_list = zeros(nsession, ntrial, 151)
    codingdimInh_R_list = zeros(nsession, ntrial, 151)
    codingdimInh_L_list = zeros(nsession, ntrial, 151)

    # There are nsession subpopulations. Select random ncells for each subpopulation. 
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
        usum_R, times_R, ns_R = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_R,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off)

        isperturbed = true
        # usum_R_perturb, times_R_perturb, ns_R_perturb = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_R,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off)
        usum_R_perturb, times_R_perturb, ns_R_perturb = runtestperturb_with_stim(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_R,stim_L,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off)

        #----- lick left -----#
        println("lick left")
        isperturbed = false
        usum_L, times_L, ns_L = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_L,ffwdRate[2],wpWeightFfwd[2],isperturbed, pert_on, pert_off)
        

        for popi = 1:nsession
            # almOrd_subpop = almOrd[subpop[popi]]
            if subpoptype == "trained"
                matchedCells_subpop = matchedCells[subpop[popi]]
            elseif subpoptype == "exc"
                matchedCells_subpop = collect(1:nexc)[subpop[popi]]
            end
            inhCells = collect(nexc+1:2*nexc)

            rate_R_perturb = spk2rate(p, matchedCells_subpop, times_R_perturb, ns_R_perturb)
            rate_R = spk2rate(p, matchedCells_subpop, times_R, ns_R)
            rate_L = spk2rate(p, matchedCells_subpop, times_L, ns_L)
            rate_M = mean((rate_R + rate_L)/2, dims=1)[:]

            # coding dimension --- excitatory
            pyr_codingdim = mean(rate_R_avg - rate_L_avg, dims=1)[matchedCells_subpop]
            pyr_codingdim = pyr_codingdim / sqrt(mean(pyr_codingdim.^2))

            codingdim_Rp_list[popi, repi,:] = mean((rate_R_perturb .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
            codingdim_R_list[popi, repi,:] = mean((rate_R .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]
            codingdim_L_list[popi, repi,:] = mean((rate_L .- reshape(rate_M,1,:)).*reshape(pyr_codingdim,1,:), dims=2)[:]

            # homogeneous dimension --- excitatory
            pyr_homdim = ones(size(rate_R)[2])
            pyr_homdim = pyr_homdim / sqrt(mean(pyr_homdim.^2))

            homdim_Rp_list[popi, repi,:] = mean((rate_R_perturb .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            homdim_R_list[popi, repi,:] = mean((rate_R .- reshape(rate_M,1,:)).*reshape(pyr_homdim,1,:), dims=2)[:]
            
            # coding dimension --- untrained inhibitory
            rateInh_R_perturb = spk2rate(p, inhCells, times_R_perturb, ns_R_perturb)
            rateInh_R = spk2rate(p, inhCells, times_R, ns_R)
            rateInh_L = spk2rate(p, inhCells, times_L, ns_L)
            rateInh_M = mean((rateInh_R + rateInh_L)/2, dims=1)[:]

            fs_codingdim = mean(rate_R_avg - rate_L_avg, dims=1)[inhCells]
            fs_codingdim = fs_codingdim / sqrt(mean(fs_codingdim.^2))

            codingdimInh_Rp_list[popi, repi, :] = mean((rateInh_R_perturb .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]
            codingdimInh_R_list[popi, repi, :] = mean((rateInh_R .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]
            codingdimInh_L_list[popi, repi, :] = mean((rateInh_L .- reshape(rateInh_M,1,:)).*reshape(fs_codingdim,1,:), dims=2)[:]

        end

    end

    # excitatory
    codingdim_Rp_mean = dropdims(mean(codingdim_Rp_list, dims=2), dims=2)
    codingdim_R_mean = dropdims(mean(codingdim_R_list, dims=2), dims=2)
    codingdim_L_mean = dropdims(mean(codingdim_L_list, dims=2), dims=2)
    homdim_Rp_mean = dropdims(mean(homdim_Rp_list, dims=2), dims=2)
    homdim_R_mean = dropdims(mean(homdim_R_list, dims=2), dims=2)

    # inhibitory
    codingdimInh_Rp_mean = dropdims(mean(codingdimInh_Rp_list, dims=2), dims=2)
    codingdimInh_R_mean = dropdims(mean(codingdimInh_R_list, dims=2), dims=2)
    codingdimInh_L_mean = dropdims(mean(codingdimInh_L_list, dims=2), dims=2)


    return codingdim_R_mean, codingdim_L_mean, codingdim_Rp_mean, homdim_R_mean, homdim_Rp_mean, codingdimInh_R_mean, codingdimInh_L_mean, codingdimInh_Rp_mean,
    codingdim_Rp_list, codingdim_R_list, homdim_Rp_list, homdim_R_list
    
    end