function plt_perturbation_lickRight(dirfig, codingdim_R_mean, codingdim_L_mean, codingdim_Rp_mean,
    homdim_R_mean, homdim_Rp_mean, 
    codingdimInh_R_mean, codingdimInh_L_mean, codingdimInh_Rp_mean, 
    codingdim_Rp_list, codingdim_R_list, homdim_Rp_list, homdim_R_list, pert_on, pert_off)


    timev = 20*collect(1:151)
    figure(figsize=(2.5,1.8))
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    plot(timev, mean(codingdim_R_mean, dims=1)[:], color="blue", lw=0.8, label="right")
    plot(timev, mean(codingdim_L_mean, dims=1)[:], color="red", lw=0.8, label="left")
    plot(timev, mean(codingdim_Rp_mean, dims=1)[:], color="gray", lw=0.8, label="perturb")
    xlim([800, 3000])
    ylim([-3, 3])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)
    legend_elements = [Line2D([0], [0], color="blue", lw=1, label="right"),
                       Line2D([0], [0], color="red", lw=1, label="left"),
                       Line2D([0], [0], color="gray", lw=1, label="perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlabel("time (s)", fontsize=8)  
    ylabel("choice mode (a.u.)", fontsize=8)  
    tight_layout()

    savefig("fig_codingdim_subpop.png", dpi=600)
    # savefig(dirfig * "codingdim.png", dpi=600)



    figure(figsize=(9,9))
    for i = 1:10
        subplot(3,4,i)
        axvspan(800,1000, color="dodgerblue", alpha=0.3)
        axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
        plot(timev, funMovAvg(codingdim_R_mean[i,:], 5), color="blue", lw=1.5, label="right")
        plot(timev, funMovAvg(codingdim_Rp_mean[i,:], 5), color="gray", lw=1.5, label="perturb")
        # error estimate
        rsem = std(funMovAvg2D(transpose(codingdim_R_list[i,:,:]), 5), dims=2)[:] /sqrt(10)
        rpsem = std(funMovAvg2D(transpose(codingdim_Rp_list[i,:,:]), 5), dims=2)[:] /sqrt(10)
        rmean = funMovAvg(codingdim_R_mean[i,:], 5)
        rpmean = funMovAvg(codingdim_Rp_mean[i,:], 5)
        fill_between(timev, rmean + rsem, rmean - rsem, color="blue", alpha=0.3)
        fill_between(timev, rpmean + rpsem, rpmean - rpsem, color="gray", alpha=0.3)
  
        xlim([800, 3000])
        ylim([-3, 3])
        xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
        yticks(fontsize=8)
    end
    tight_layout()

    savefig("fig_codingdim_sessions.png", dpi=600)
    # savefig(dirfig * "codingdim.png", dpi=600)
    
    


    figure(figsize=(2.5,1.8))
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    plot(timev, mean(codingdimInh_R_mean, dims=1)[:], color="darkorange", lw=0.8, label="right")
    plot(timev, mean(codingdimInh_Rp_mean, dims=1)[:], color="gray", lw=0.8, label="perturb")
    xlim([800, 3000])
    ylim([-3, 3])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)
    legend_elements = [Line2D([0], [0], color="darkorange", lw=1, label="right"),
                       Line2D([0], [0], color="gray", lw=1, label="perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlabel("time (s)", fontsize=8)  
    ylabel("choice mode (a.u.)", fontsize=8)  
    tight_layout()

    savefig("fig_codingdimInh_subpop.png", dpi=600)
    # savefig(dirfig * "codingdim.png", dpi=600)



    figure(figsize=(2.5,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    # plot(timev, untrained_R_mean, color="C2", lw=0.8, alpha=0.5)
    # plot(timev, untrained_Rp_mean, color="C2", lw=0.8, label="untrained")
    plot(timev, mean(homdim_R_mean, dims=1)[:], color="C2", lw=0.8, alpha=1, label="hom")
    plot(timev, mean(homdim_Rp_mean, dims=1)[:], color="gray", lw=0.8, label="perturb")
    ylim([-3, 3])
    xlabel("time (s)", fontsize=8)  
    ylabel("hom mode (a.u.)", fontsize=8)  
    legend_elements = [Line2D([0], [0], color="C2", lw=1, label="hom"),
                       Line2D([0], [0], color="gray", lw=1, label="perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   

    # text(1000,12,"not trained", fontsize=6, color="black")
    # text(2200,3,"trained", fontsize=6, color="black")
    # legend(frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    savefig("fig_homdim_subpop.png", dpi=600)
    # savefig(dirfig * "meanrate_trained.png", dpi=600)




    figure(figsize=(9,9))
    for i = 1:10
        subplot(3,4,i)
        axvspan(800,1000, color="dodgerblue", alpha=0.3)
        axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
        plot(timev, funMovAvg(homdim_R_mean[i,:], 5), color="C2", lw=1.5, label="right")
        plot(timev, funMovAvg(homdim_Rp_mean[i,:], 5), color="gray", lw=1.5, label="perturb")
        # error estimate
        rsem = std(funMovAvg2D(transpose(homdim_R_list[i,:,:]), 5), dims=2)[:] /sqrt(10)
        rpsem = std(funMovAvg2D(transpose(homdim_Rp_list[i,:,:]), 5), dims=2)[:] /sqrt(10)
        rmean = funMovAvg(homdim_R_mean[i,:], 5)
        rpmean = funMovAvg(homdim_Rp_mean[i,:], 5)
        fill_between(timev, rmean + rsem, rmean - rsem, color="C2", alpha=0.3)
        fill_between(timev, rpmean + rpsem, rpmean - rpsem, color="gray", alpha=0.3)
  
        xlim([800, 3000])
        ylim([-3, 3])
        xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
        yticks(fontsize=8)
    end
    tight_layout()

    savefig("fig_homdim_sessions.png", dpi=600)
    # savefig(dirfig * "codingdim.png", dpi=600)
    


    figure(figsize=(2.0,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, color="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    # plot(timev, untrained_R_mean, color="C2", lw=0.8, alpha=0.5)
    # plot(timev, untrained_Rp_mean, color="C2", lw=0.8, label="untrained")
    plot(timev, mean(homdim_Rp_mean, dims=1)[:] - mean(homdim_R_mean, dims=1)[:], color="C2", lw=0.8)
    plot(timev, mean(codingdim_Rp_mean, dims=1)[:] - mean(codingdim_R_mean, dims=1)[:], color="blue", lw=0.8)
    # plot(timev, mean(codingdimInh_Rp_mean, dims=1)[:] - mean(codingdimInh_R_mean, dims=1)[:], color="darkorange", lw=0.8)
    ylim([-2, 2])
    xlabel("time (s)", fontsize=8)  
    ylabel("diff (a.u.)", fontsize=8)  
    # text(1000,12,"not trained", fontsize=6, color="black")
    # text(2200,3,"trained", fontsize=6, color="black")
    # legend(frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    savefig("fig_diff_subpop.png", dpi=600)


    # figure(figsize=(2.0,1.8))
    # timev = 20*collect(1:151)
    # axvspan(800,1000, color="dodgerblue", alpha=0.3)
    # axvspan(pert_on,pert_off, color="dodgerblue", alpha=0.3)
    # plot(timev, exc_R_mean, color="blue", lw=0.8, alpha=0.5)
    # plot(timev, inh_R_mean, color="red", lw=0.8, alpha=0.5)
    # plot(timev, exc_Rp_mean, color="blue", lw=0.8, label="Exc")
    # plot(timev, inh_Rp_mean, color="red", lw=0.8, label="Inh")
    # ylim([0, 15])
    # legend(loc=1, frameon=false, fontsize=8, handlelength=1)
    # xlabel("time (s)", fontsize=8)  
    # ylabel("spiking rate (Hz)", fontsize=8)  
    # # text(0,10,"excitatory", fontsize=8, color="black")
    # xlim([800, 3000])
    # xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    # yticks(fontsize=8)

    # tight_layout()

    # savefig("fig_meanrate_exc.png", dpi=600)
    # savefig(dirfig * "meanrate_exc.png", dpi=600)

end