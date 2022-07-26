function plt_fsinh_psth_heatmap(dirfig, dirsim, rate_R, rate_L, pcor_R, pcor_L, cells_matched)

    Ne = 2500
    Ncells = 5000

    start_idx = 5
    rate_R = rate_R[start_idx:end,:]
    rate_L = rate_L[start_idx:end,:]

    # load rate_navg
    rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[start_idx:end,Ne+1:Ncells,:]
    rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[start_idx:end,Ne+1:Ncells,:]
    
    # performance: Pyr - Exc
    # pcorMatched = (pcor_R + pcor_L)/2
    pcorMatched = pcor_R 
    pcorSorted = reverse(sortperm(pcorMatched))[1:295]


    ci = Int.(cells_matched[pcorSorted])
    ci_alm = pcorSorted
    rate_R_navg400_sorted = rate_R_navg400[:,ci,1]
    rate_L_navg400_sorted = rate_L_navg400[:,ci,1]
    rate_R_sorted = rate_R[:,ci_alm]
    rate_L_sorted = rate_L[:,ci_alm]

    rate_R_navg400_sorted = rate_R_navg400_sorted ./ maximum(rate_R_navg400_sorted, dims=1)
    rate_L_navg400_sorted = rate_L_navg400_sorted ./ maximum(rate_L_navg400_sorted, dims=1)
    rate_R_sorted = rate_R_sorted ./ maximum(rate_R_sorted, dims=1)
    rate_L_sorted = rate_L_sorted ./ maximum(rate_L_sorted, dims=1)

    #----- lick right, inhibitory (model) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_R_navg400_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 200], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Inh neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickR_inh.png", dpi=600)

    savefig(dirfig * "raster/fsinh/fig_lickR_inh.pdf", dpi=600)

    #----- lick right, fast spiking (data) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_R_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 200], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("FS neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickR_fs.png", dpi=600)

    savefig(dirfig * "raster/fsinh/fig_lickR_fs.pdf", dpi=600)


    

    #----- lick right, inhibitory (model) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_L_navg400_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 200], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("Inh neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickL_inh.png", dpi=600)

    savefig(dirfig * "raster/fsinh/fig_lickL_inh.pdf", dpi=600)

    #----- lick right, fast spiking (data) -----#
    figure(figsize=(2.0,1.5))
    imshow(transpose(rate_L_sorted), vmin=0, vmax=1, cmap="jet", interpolation="gaussian", aspect="auto")
    cbar = colorbar(ticks=[0,1])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.ax.tick_params(labelsize=7)
    # cbar.ticks([0,1],fontsize=7)
    xticks([0,50,96], [-2,-1,0], fontsize=7)
    # yticks([0,100,200],fontsize=7)
    yticks([0, 200], fontsize=7)
    xlabel("time (s)", fontsize=7)
    ylabel("FS neuron", fontsize=7)
    tight_layout()

    # savefig("fig_lickL_fs.png", dpi=600)

    savefig(dirfig * "raster/fsinh/fig_lickL_fs.pdf", dpi=600)




end