function plt_sharedvar_data_ei(fname, svar_cor1, svarX1, svarY1, svarX1_ksum, svarX1_tsum, svarY1_ksum, svarY1_tsum,
    svar_cor2, svarX2, svarY2, svarX2_ksum, svarX2_tsum, svarY2_ksum, svarY2_tsum)

    timev = 20*collect(2:100)

    fig=figure(figsize=(2.5,4.5))
    ncomp = 4
    sgn = [-1,1,1,1]
    for ii = 1:ncomp
        ax = fig.add_subplot(4,1,ii)
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)
        ylabel("SVC $(ii)", fontsize=12)
        xticks(color="None")
        yticks(color="None")
        xlim([0, 2000])
        # ylim([-4, 2.5])

        plot(timev, svarX1[2:end,ii], c="limegreen")
        plot(timev, svarY1[2:end,ii], c="darkorange")
        if ii == 1
            legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                               Line2D([0], [0], color="darkorange", lw=2, label="FS")]
            legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
        end

        # plot(timev, sgn[ii]*svarX2[2:end,ii], c="red", linestyle="-")
        # plot(timev, sgn[ii]*svarY2[2:end,ii], c="blue", linestyle="-")
        # if ii == 2
        #     legend_elements = [Line2D([0], [0], color="red", lw=2, label="Exc"),
        #                        Line2D([0], [0], color="blue", lw=2, label="Inh")]
        #     legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
        # end

        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -3.5*ones(length(tbar)), c="black", lw=2)
            text(1500, -3, "500ms", fontsize=8)
        end 
    end
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_pyrfs.png", dpi=300)
    
    


    fig=figure(figsize=(2.5,4.5))
    ncomp = 4
    sgn = [-1,1,1,1]
    for ii = 1:ncomp
        ax = fig.add_subplot(4,1,ii)
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)
        ylabel("SVC $(ii)", fontsize=12)
        xticks(color="None")
        yticks(color="None")
        xlim([0, 2000])
        # ylim([-4, 2.5])

        # plot(timev, svarX1[2:end,ii], c="limegreen")
        # plot(timev, svarY1[2:end,ii], c="darkorange")
        # if ii == 1
        #     legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
        #                        Line2D([0], [0], color="darkorange", lw=2, label="FS")]
        #     legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
        # end

        plot(timev, sgn[ii]*svarX2[2:end,ii], c="red", linestyle="-", lw=1)
        plot(timev, sgn[ii]*svarY2[2:end,ii], c="blue", linestyle="-", lw=1)
        if ii == 1
            legend_elements = [Line2D([0], [0], color="red", lw=2, label="Exc"),
                               Line2D([0], [0], color="blue", lw=2, label="Inh")]
            legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
        end

        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -3.5*ones(length(tbar)), c="black", lw=2)
            text(1500, -3, "500ms", fontsize=8)
        end
        
    end
    
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_ei.png", dpi=300)
    
    

    
    figure(figsize=(3.5,3))
    ncomp = 10
    plot(collect(1:ncomp), svar_cor1[1:ncomp], c="magenta", marker="o", ms=3)
    plot(collect(1:ncomp), svar_cor2[1:ncomp], c="dodgerblue", marker="o", ms=3, linestyle="--")
    # ylim([-0.1,1.1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("correlation " * L"($\rho$)", fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)

    legend_elements = [Line2D([0], [0], color="magenta", lw=2, label="Pyr-FS"),
                       Line2D([0], [0], color="dodgerblue", lw=2, label="E-I")]
    legend(loc=4, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     

    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_corall.png", dpi=300)
    
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 10
    plot(collect(1:ncomp), cumsum(svarX1_ksum[1:ncomp]/svarX1_tsum), c="limegreen", marker="o", ms=3)
    plot(collect(1:ncomp), cumsum(svarY1_ksum[1:ncomp]/svarY1_tsum), c="darkorange", marker="o", ms=3)
    # legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
    #                    Line2D([0], [0], color="darkorange", lw=2, label="FS")]
    # legend(loc=3, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     

    plot(collect(1:ncomp), cumsum(svarX2_ksum[1:ncomp]/svarX2_tsum), c="red", marker="o", ms=3, linestyle="--")
    plot(collect(1:ncomp), cumsum(svarY2_ksum[1:ncomp]/svarY2_tsum), c="blue", marker="o", ms=3, linestyle="--")
    legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                       Line2D([0], [0], color="darkorange", lw=2, label="FS"),
                       Line2D([0], [0], color="red", lw=2, label="Exc"),
                       Line2D([0], [0], color="blue", lw=2, label="Inh")]
    legend(loc=4, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     

    ylim([0,1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("cumulative var", fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_cumvar.png", dpi=300)
    
    

end

