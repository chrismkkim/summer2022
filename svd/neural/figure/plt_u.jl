function plt_u()

    println("Loading files")
    lickRL = "right"
    p, w0Index, w0Weights, nc0, wpIndexIn, wpIndexOut, wpIndexConvert, 
    wpWeightOut, ncpIn, ncpOut, stim, almOrd, matchedCells,
    ffwdRate, wpWeightFfwd = load_files(dirdata, lickRL)

    # ubal, uplas, uext, times, ns = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate,wpWeightFfwd[1])

    # ##############################################
    # # block plastic weights trained on lick left #
    # ##############################################
    # numExc = Int(p.Lexc / 2)
    # numInh = Int(p.Linh / 2)
    # rangeIn = [collect(numExc+1 : p.Lexc); collect(p.Lexc+numInh+1 : p.Lexc+p.Linh)]
    # wpWeightIn[:,rangeIn] .= 0
    # wpWeightOut = convertWgtIn2Out(p,ncpIn,wpIndexIn,wpIndexConvert,wpWeightIn,wpWeightOut)
    # ubal, uplas, uext, times, ns = runtestdrive(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim,ffwdRate,wpWeightFfwd[1])

    timev = p.dt*collect(1:p.Nsteps)

    for nid = 1:10
        # nid = 1
        # usum = ubal[:,nid][:] + uplas[:,nid][:] + uext[:,nid][:]

        figure(figsize=(3.5,3))
        axvspan(p.stim_on, p.stim_off, color="dodgerblue", alpha=0.3)
        # plot(timev, usum, c="black", linewidth=0.5, alpha=0.5)
        plot(timev, ubal[:,nid][:], c="black", linewidth=1.0, alpha=0.5, label="Ubal")
        plot(timev, uext[:,nid][:], c="red", linewidth=1.0, alpha=1, label="Uext")
        plot(timev, uplas[:,nid][:], c="darkorange", linewidth=1.0, alpha=1, label="Uplas")
        # plot(timev, ones(length(timev)), c="magenta", linewidth=1.5, linestyle="--", alpha=1)
        legend(loc=1, frameon=false, fontsize=12)
        xlim([0,3000])
        ylim([-1,2])
        xlabel("time (ms)", fontsize=12)
        ylabel("synaptic input", fontsize=12)
        tight_layout()

        savefig("trained_syn$(nid)_alt.png", dpi=300)


        close()
    end

end