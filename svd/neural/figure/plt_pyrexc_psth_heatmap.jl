function plt_pyrexc_psth_heatmap(dirfig, dirsim, rate_R, rate_L, pcor_R, pcor_L, almOrd, matchedCells)

    Ne = 2500
    Ncells = 5000

    start_idx = 5
    rate_R = rate_R[start_idx:end,:]
    rate_L = rate_L[start_idx:end,:]

    # load rate_navg
    rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[start_idx:end,1:Ne,:]
    rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[start_idx:end,1:Ne,:]
    
    # performance: Pyr - Exc
    # pcorMatched = (pcor_R + pcor_L)/2
    pcorMatched = pcor_R
    pcorSorted = reverse(sortperm(pcorMatched))

    ci = Int.(matchedCells[pcorSorted])
    ci_alm = Int.(almOrd[pcorSorted])
    rate_R_navg400_sorted = rate_R_navg400[:,ci,1]
    rate_L_navg400_sorted = rate_L_navg400[:,ci,1]
    rate_R_sorted = rate_R[:,ci_alm]
    rate_L_sorted = rate_L[:,ci_alm]

    # smooth the network model PSTH, then normalize by the max rate
    rate_R_navg400_sorted_smooth = funMovAvg_2d(copy(rate_R_navg400_sorted), 1)
    rate_L_navg400_sorted_smooth = funMovAvg_2d(copy(rate_L_navg400_sorted), 1)

    rate_R_navg400_sorted = rate_R_navg400_sorted ./ maximum(rate_R_navg400_sorted_smooth, dims=1)
    rate_L_navg400_sorted = rate_L_navg400_sorted ./ maximum(rate_L_navg400_sorted_smooth, dims=1)
    rate_R_sorted = rate_R_sorted ./ maximum(rate_R_sorted, dims=1)
    rate_L_sorted = rate_L_sorted ./ maximum(rate_L_sorted, dims=1)


    #----- lick right, inhibitory (model) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_R_navg400_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    # imshow(transpose(rate_R_navg400_sorted), vmin=0, vmax=20, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 1000], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Exc neuron", fontsize=7)
    tight_layout()
    
    # savefig("fig_lickR_exc.png", dpi=600)

    savefig(dirfig * "raster/pyrexc/fig_lickR_exc.pdf", dpi=600)

    #----- lick right, fast spiking (data) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_R_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    # imshow(transpose(rate_R_sorted), vmin=0, vmax=20, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 1000], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Pyr neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickR_pyr.png", dpi=600)

    savefig(dirfig * "raster/pyrexc/fig_lickR_pyr.pdf", dpi=600)


    

    #----- lick right, inhibitory (model) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_L_navg400_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 1000], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Exc neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickL_exc.png", dpi=600)

    savefig(dirfig * "raster/pyrexc/fig_lickL_exc.pdf", dpi=600)

    #----- lick right, fast spiking (data) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_L_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 1000], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Pyr neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickL_pyr.png", dpi=600)

    savefig(dirfig * "raster/pyrexc/fig_lickL_pyr.pdf", dpi=600)





end