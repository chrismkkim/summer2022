function plt_perturbation(dirfig, codingdim_R_mean, codingdim_L_mean, codingdim_Rp_mean,
    untrained_R_mean, untrained_Rp_mean, exc_R_mean, inh_R_mean, exc_Rp_mean, inh_Rp_mean, trained_R_mean, trained_Rp_mean, pert_on, pert_off)


    timev = 20*collect(1:151)
    figure(figsize=(2.5,1.8))
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    plot(timev, codingdim_R_mean, color="blue", lw=0.8, label="right")
    plot(timev, codingdim_L_mean, color="red", lw=0.8, label="left")
    plot(timev, codingdim_Rp_mean, color="gray", lw=0.8, label="perturb")
    xlim([800, 3000])
    ylim([-2, 3])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)
    legend_elements = [Line2D([0], [0], color="blue", lw=1, label="right"),
                       Line2D([0], [0], color="red", lw=1, label="left"),
                       Line2D([0], [0], color="gray", lw=1, label="perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlabel("time (s)", fontsize=8)  
    ylabel("choice mode (a.u.)", fontsize=8)  
    tight_layout()

    savefig("fig_codingdim.png", dpi=600)
    # savefig(dirfig * "codingdim.png", dpi=600)


    figure(figsize=(2.0,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    plot(timev, untrained_R_mean, color="C2", lw=0.8, alpha=0.5)
    plot(timev, untrained_Rp_mean, color="C2", lw=0.8, label="untrained")
    plot(timev, trained_R_mean, color="C1", lw=0.8, alpha=0.5)
    plot(timev, trained_Rp_mean, color="C1", lw=0.8, label="trained")
    ylim([0, 15])
    xlabel("time (s)", fontsize=8)  
    ylabel("spiking rate (Hz)", fontsize=8)  
    text(1000,12,"not trained", fontsize=6, color="black")
    text(2200,3,"trained", fontsize=6, color="black")
    # legend(frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    savefig("fig_meanrate_trained.png", dpi=600)
    # savefig(dirfig * "meanrate_trained.png", dpi=600)


    figure(figsize=(2.0,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    plot(timev, exc_R_mean, color="blue", lw=0.8, alpha=0.5)
    plot(timev, inh_R_mean, color="red", lw=0.8, alpha=0.5)
    plot(timev, exc_Rp_mean, color="blue", lw=0.8, label="Exc")
    plot(timev, inh_Rp_mean, color="red", lw=0.8, label="Inh")
    ylim([0, 15])
    legend(loc=1, frameon=false, fontsize=8, handlelength=1)
    xlabel("time (s)", fontsize=8)  
    ylabel("spiking rate (Hz)", fontsize=8)  
    # text(0,10,"excitatory", fontsize=8, color="black")
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    savefig("fig_meanrate_exc.png", dpi=600)
    # savefig(dirfig * "meanrate_exc.png", dpi=600)

end