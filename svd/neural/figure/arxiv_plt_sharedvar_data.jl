function plt_sharedvar_data(fname, svar_cor, svarX, svarY, svarX_ksum, svarX_tsum, svarY_ksum, svarY_tsum)

    timev = 20*collect(2:100)

    fig=figure(figsize=(5.0,4.5))
    ncomp = 4
    for ii = 1:ncomp
        subplot(2,2,ii)
        plot(timev, svarX[2:end,ii], c="limegreen")
        plot(timev, svarY[2:end,ii], c="darkorange")
        title("SVC $(ii)", fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        if ii > 1
            xticks(color="None")
            yticks(color="None")
        end
        if ii == 1
            legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                               Line2D([0], [0], color="darkorange", lw=2, label="FS")]
            legend(loc=4, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
        end
        xlim([0, 2000])
        ylim([-4, 2.5])
    end
    ax = fig.add_subplot(111)
    xlabel("time (ms)", fontsize=12, labelpad=30)
    ylabel("SVC", fontsize=12, labelpad=30)
    ax.set_frame_on(false)
    ax.tick_params(left=false, bottom=false)
    ax.tick_params(labelleft=false, labelbottom=false)
    
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * ".png", dpi=300)
    
    
    
    
    
    fig=figure(figsize=(5.0,4.5))
    ncomp = 4
    for ii = 1:ncomp
        subplot(2,2,ii)
        plot(svarX[2:end,ii], svarY[2:end,ii], c="magenta", marker="o", ms=1, linestyle="", label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")
        title("SVC $(ii)", fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        if ii > 1
            xticks(color="None")
            yticks(color="None")
        end
        xlim([-2, 1.5])
        ylim([-4, 2.5])
        legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")]
        legend(loc=2, handles=legend_elements, frameon=false, fontsize=12, handlelength=0)     
    end
    ax = fig.add_subplot(111)
    xlabel("Pyr cells", fontsize=12, labelpad=30)
    ylabel("FS cells", fontsize=12, labelpad=30)
    ax.set_frame_on(false)
    ax.tick_params(left=false, bottom=false)
    ax.tick_params(labelleft=false, labelbottom=false)
    
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_cor.png", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 10
    plot(collect(1:ncomp), svar_cor[1:ncomp], c="magenta", marker="o")
    ylim([-0.1,1.1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("correlation " * L"($\rho$)", fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_corall.png", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 10
    plot(collect(1:ncomp), svarX_ksum[1:ncomp]/svarX_tsum, c="limegreen", marker="o")
    plot(collect(1:ncomp), svarY_ksum[1:ncomp]/svarY_tsum, c="darkorange", marker="o")
    ylim([5e-4,1])
    yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("explained var", fontsize=12)
    legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                       Line2D([0], [0], color="darkorange", lw=2, label="FS")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_expvar.png", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 10
    plot(collect(1:ncomp), cumsum(svarX_ksum[1:ncomp]/svarX_tsum), c="limegreen", marker="o")
    plot(collect(1:ncomp), cumsum(svarY_ksum[1:ncomp]/svarY_tsum), c="darkorange", marker="o")
    legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                       Line2D([0], [0], color="darkorange", lw=2, label="FS")]
    legend(loc=2, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
    
    ylim([0,1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("cumulative var", fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_cumvar.png", dpi=300)
    
    

end

