function plt_perturbation(dirfig, 
    codingdim_PR, codingdim_PL, codingdim_R, codingdim_L, 
    homdim_PR, homdim_PL, homdim_R, homdim_L, 
    codingdimInh_PR, codingdimInh_PL, codingdimInh_R, codingdimInh_L,
    codingdim_PR_list, codingdim_PL_list, codingdim_R_list, codingdim_L_list, 
    homdim_PR_list, homdim_PL_list, homdim_R_list, homdim_L_list,
    pert_on, pert_off, perturbtype)

    if perturbtype == "lickRight"

        codingdim_P = copy(codingdim_PR)
        codingdim_P_list = copy(codingdim_PR_list)
        codingdim_ref = copy(codingdim_R)
        codingdim_ref_list = copy(codingdim_R_list)
        codingdimInh_P = copy(codingdimInh_PR)
        codingdimInh_ref = copy(codingdimInh_R)
        homdim_P = copy(homdim_PR)
        homdim_P_list = copy(homdim_PR_list)
        homdim_ref = copy(homdim_R)
        homdim_ref_list = copy(homdim_R_list)
        label_ref = "Right"
        clr_ref = "blue"

    elseif perturbtype == "lickLeft"

        codingdim_P = copy(codingdim_PL)
        codingdim_P_list = copy(codingdim_PL_list)
        codingdim_ref = copy(codingdim_L)
        codingdim_ref_list = copy(codingdim_L_list)
        codingdimInh_P = copy(codingdimInh_PL)
        codingdimInh_ref = copy(codingdimInh_L)
        homdim_P = copy(homdim_PL)
        homdim_P_list = copy(homdim_PL_list)
        homdim_ref = copy(homdim_L)
        homdim_ref_list = copy(homdim_L_list)
        label_ref = "Left"
        clr_ref = "red"
    end

    timev = 20*collect(1:151)
    figure(figsize=(2.5,1.8))
    axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
    plot(timev, mean(codingdim_R, dims=1)[:], color="blue", lw=0.8, label="right")
    plot(timev, mean(codingdim_L, dims=1)[:], color="red", lw=0.8, label="left")
    plot(timev, mean(codingdim_P, dims=1)[:], color="gray", lw=0.8, label="perturb")
    xlim([800, 3000])
    ylim([-3, 3])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)
    legend_elements = [Line2D([0], [0], color="blue", lw=1, label="Right"),
                       Line2D([0], [0], color="red", lw=1, label="Left"),
                       Line2D([0], [0], color="gray", lw=1, label="Perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlabel("time (s)", fontsize=8)  
    ylabel("choice mode (a.u.)", fontsize=8)  
    tight_layout()

    if perturbtype == "lickRight"
        savefig("fig_R_codingdim.png", dpi=600)
    elseif perturbtype == "lickLeft"
        savefig("fig_L_codingdim.png", dpi=600)
    end
    # savefig(dirfig * "codingdim.png", dpi=600)



    figure(figsize=(2.5,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
    # plot(timev, untrained_R, color="C2", lw=0.8, alpha=0.5)
    # plot(timev, untrained_Rp, color="C2", lw=0.8, label="untrained")
    plot(timev, mean(homdim_ref, dims=1)[:], color="black", lw=0.8, alpha=1, label="Hom")
    plot(timev, mean(homdim_P, dims=1)[:], color="gray", lw=0.8, label="Perturb")
    ylim([-1, 10])
    xlabel("time (s)", fontsize=8)  
    ylabel("hom mode (a.u.)", fontsize=8)  
    legend_elements = [Line2D([0], [0], color="black", lw=1, label="Hom"),
                       Line2D([0], [0], color="gray", lw=1, label="Perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   

    # text(1000,12,"not trained", fontsize=6, color="black")
    # text(2200,3,"trained", fontsize=6, color="black")
    # legend(frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    if perturbtype == "lickRight"
        savefig("fig_R_homdim.png", dpi=600)
    elseif perturbtype == "lickLeft"
        savefig("fig_L_homdim.png", dpi=600)
    end
    # savefig(dirfig * "meanrate_trained.png", dpi=600)



    figure(figsize=(2.5,1.8))
    timev = 20*collect(1:151)
    axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
    plot(timev, mean(homdim_P, dims=1)[:] - mean(homdim_ref, dims=1)[:], color="black", lw=0.8)
    plot(timev, mean(codingdim_P, dims=1)[:] - mean(codingdim_ref, dims=1)[:], color=clr_ref, lw=0.8)
    ylim([-2, 2])
    xlabel("time (s)", fontsize=8)  
    ylabel("diff (a.u.)", fontsize=8)  
    legend_elements = [Line2D([0], [0], color=clr_ref, lw=1, label="Coding"),
                       Line2D([0], [0], color="black", lw=1, label="Hom")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlim([800, 3000])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)

    tight_layout()

    if perturbtype == "lickRight"
        savefig("fig_R_diff.png", dpi=600)
    elseif perturbtype == "lickLeft"
        savefig("fig_L_diff.png", dpi=600)
    end


    figure(figsize=(9,9))
    for i = 1:10
        subplot(3,4,i)
        axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
        axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
        plot(timev, funMovAvg(codingdim_ref[i,:], 5), color=clr_ref, lw=1.5, label=label_ref)
        plot(timev, funMovAvg(codingdim_P[i,:], 5), color="gray", lw=1.5, label="perturb")
        # error estimate
        rsem = std(funMovAvg2D(transpose(codingdim_ref_list[i,:,:]), 5), dims=2)[:] /sqrt(15) # 1/sqrt(ntrial)
        rpsem = std(funMovAvg2D(transpose(codingdim_P_list[i,:,:]), 5), dims=2)[:] /sqrt(15)
        rmean = funMovAvg(codingdim_ref[i,:], 5)
        rpmean = funMovAvg(codingdim_P[i,:], 5)
        fill_between(timev, rmean + rsem, rmean - rsem, color=clr_ref, alpha=0.3)
        fill_between(timev, rpmean + rpsem, rpmean - rpsem, color="gray", alpha=0.3)
  
        xlim([800, 3000])
        ylim([-4, 3])
        xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
        yticks(fontsize=8)
    end
    tight_layout()

    if perturbtype == "lickRight"
        savefig("fig_R_codingdim_sessions.png", dpi=600)
    elseif perturbtype == "lickLeft"        
        savefig("fig_L_codingdim_sessions.png", dpi=600)
    end        
    # savefig(dirfig * "codingdim.png", dpi=600)
    
    

    figure(figsize=(9,9))
    for i = 1:10
        subplot(3,4,i)
        axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
        axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
        plot(timev, funMovAvg(homdim_ref[i,:], 5), color="black", lw=1.5, label=label_ref)
        plot(timev, funMovAvg(homdim_P[i,:], 5), color="gray", lw=1.5, label="perturb")
        # error estimate
        rsem = std(funMovAvg2D(transpose(homdim_ref_list[i,:,:]), 5), dims=2)[:] /sqrt(15)
        rpsem = std(funMovAvg2D(transpose(homdim_P_list[i,:,:]), 5), dims=2)[:] /sqrt(15)
        rmean = funMovAvg(homdim_ref[i,:], 5)
        rpmean = funMovAvg(homdim_P[i,:], 5)
        fill_between(timev, rmean + rsem, rmean - rsem, color="black", alpha=0.3)
        fill_between(timev, rpmean + rpsem, rpmean - rpsem, color="gray", alpha=0.3)
  
        xlim([800, 3000])
        ylim([4, 10])
        xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
        yticks(fontsize=8)
    end
    tight_layout()

    if perturbtype == "lickRight"      
        savefig("fig_R_homdim_sessions.png", dpi=600)
    elseif perturbtype == "lickLeft"
        savefig("fig_L_homdim_sessions.png", dpi=600)
    end        
    # savefig(dirfig * "codingdim.png", dpi=600)
        


    figure(figsize=(2.5,1.8))
    axvspan(800,1000, facecolor="dodgerblue", alpha=0.3)
    axvspan(pert_on,pert_off, facecolor="dodgerblue", alpha=0.3)
    plot(timev, mean(codingdimInh_ref, dims=1)[:], color="darkorange", lw=0.8, label=label_ref)
    plot(timev, mean(codingdimInh_P, dims=1)[:], color="gray", lw=0.8, label="perturb")
    xlim([800, 3000])
    ylim([-3, 3])
    xticks([1000, 2000, 3000],[0, 1, 2], fontsize=8)
    yticks(fontsize=8)
    legend_elements = [Line2D([0], [0], color="darkorange", lw=1, label=label_ref),
                       Line2D([0], [0], color="gray", lw=1, label="Perturb")]
    legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")   
    xlabel("time (s)", fontsize=8)  
    ylabel("choice mode (a.u.)", fontsize=8)  
    tight_layout()

    if perturbtype == "lickRight"           
        savefig("fig_R_codingdimInh.png", dpi=600)
    elseif perturbtype == "lickLeft"      
        savefig("fig_L_codingdimInh.png", dpi=600)
    end
    # savefig(dirfig * "codingdim.png", dpi=600)





end