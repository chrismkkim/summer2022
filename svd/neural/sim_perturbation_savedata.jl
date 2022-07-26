function sim_perturbation_savedata(p, dirdata, dirperturb, rate_R_avg, rate_L_avg, pyr_R, pyr_L, pert_on, pert_off, subpoptype, ifsave, ncell, stimtype)

    npyr = size(pyr_R)[2]
    nexc = 2500
    ntrial = 15 # 15
    nsession = 10 # 25, 18
    # ncell = 50

    wid = 1
    nsteps = Int(p.train_time / wid)

    rate_R_perturb = zeros(nsteps, ncell, ntrial, nsession)
    rate_L_perturb = zeros(nsteps, ncell, ntrial, nsession)
    rate_R = zeros(nsteps, ncell, ntrial, nsession)
    rate_L = zeros(nsteps, ncell, ntrial, nsession)

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

    for triali = 1:ntrial
        println(triali)

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
        usum_PR, times_PR, ns_PR = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[1],wpWeightFfwd[1],isperturbed, pert_on, pert_off, stimtype)

        println("perturbation - lick left")
        # trigger learned activity with stim_L, perturb with stim_R
        isperturbed = true            
        stim_training = copy(stim_L)
        stim_perturb = copy(stim_R)    
        usum_PL, times_PL, ns_PL = runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate[2],wpWeightFfwd[2],isperturbed, pert_on, pert_off, stimtype)
  

        for popi = 1:nsession
            # almOrd_subpop = almOrd[subpop[popi]]
            if subpoptype == "trained"
                matchedCells_subpop = matchedCells[subpop[popi]]
            elseif subpoptype == "exc"
                matchedCells_subpop = collect(1:nexc)[subpop[popi]]
            end
            inhCells = collect(nexc+1:2*nexc)

            rate_R_perturb[:, :, triali, popi] = spk2rate_data(p, wid, matchedCells_subpop, times_PR, ns_PR)[2:end,:]
            rate_L_perturb[:, :, triali, popi] = spk2rate_data(p, wid, matchedCells_subpop, times_PL, ns_PL)[2:end,:]
            rate_R[:, :, triali, popi]         = spk2rate_data(p, wid, matchedCells_subpop, times_R, ns_R)[2:end,:]
            rate_L[:, :, triali, popi]         = spk2rate_data(p, wid, matchedCells_subpop, times_L, ns_L)[2:end,:]

        end

    end

    if ifsave

        # save in julia
        # rate_R : time x neuron x trial x session
        save(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_R_perturb.jld", "rate_R_perturb", rate_R_perturb)
        save(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_L_perturb.jld", "rate_L_perturb", rate_L_perturb)
        save(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_R.jld", "rate_R", rate_R)
        save(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_L.jld", "rate_L", rate_L)


        # save in matlab
        file_R_perturb = matopen(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_R_perturb.mat","w")
        file_L_perturb = matopen(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_L_perturb.mat","w")
        file_R = matopen(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_R.mat","w")
        file_L = matopen(dirperturb * "rate_" * stimtype * "_ncell$(ncell)_L.mat","w")

        write(file_R_perturb, "rate_R_perturb", rate_R_perturb)
        write(file_L_perturb, "rate_L_perturb", rate_L_perturb)
        write(file_R, "rate_R", rate_R)
        write(file_L, "rate_L", rate_L)

        close(file_R_perturb)
        close(file_L_perturb)
        close(file_R)
        close(file_L)


        # plot data
        figure()
        subplot(211)
        plot(mean(rate_R[:,:,:,1], dims=(2,3))[:], c="blue", lw=0.8, label="R")
        legend()
        ylim([0,30])

        subplot(212)
        plot(mean(rate_R_perturb[:,:,:,1], dims=(2,3))[:], c="gray", lw=0.8, label="R perturb")
        legend()
            # plot(mean(rate_L[:,:,:,i], dims=(2,3))[:], c="red", lw=0.8, label="L")
        ylim([0,30])
        tight_layout()

        savefig(dirperturb * "fig_session1_psth.png", dpi=300)
        savefig("fig_session1_psth.png", dpi=300)

    end


    return rate_R_perturb, rate_R, rate_L

end


function spk2rate_data(p, wid, matchedCells, times, ns)

    # wid = 1 # 20ms
    Nsteps = Int(p.train_time / wid) + 1
    Ntrained = length(matchedCells)
    rates = zeros(Nsteps, Ntrained)

    for ci = 1:Ntrained
        nid = matchedCells[ci]
        nspk = ns[nid]
        if nspk > 0
            for spk = 1:nspk
                spktime = times[nid,spk]
                kk = Int(floor(spktime/wid)) + 1
                rates[kk,ci] += 1 / (wid/1000)
            end
        end
    end

    return rates

end