matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

function plt_sharedvar_pyrexc(dirfig, fname, svar_cor, svarX, svarY, svarX_ksum, svarX_tsum, svarY_ksum, svarY_tsum)

    startind = 7
    timev = 20*collect(startind:100)
    
    fig=figure(figsize=(2.5,4.5))
    ncomp = 4
    sgn = [-1,1,1,-1]
    for ii = 1:ncomp
        ax = fig.add_subplot(4,1,ii)
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)
        ylabel("SVC $(ii)", fontsize=12)
        xlim([0, 2000])
        ylim([-140, 120])

        plot(timev, sgn[ii]*svarX[:,ii], c="limegreen", lw=2)
        plot(timev, sgn[ii]*svarY[:,ii], c="red", lw=2, linestyle="--")
        if ii == 1
            legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                               Line2D([0], [0], color="red", lw=2, label="Exc")]
            legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
        end
        if ii == 4
            tbar = 20*collect(75:100)
            plot(tbar, -110*ones(length(tbar)), c="black", lw=2)
            text(1500, -100, "500ms", fontsize=8)
        end        
    end
    
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * ".png", dpi=300)
    
    
    
    
    fig=figure(figsize=(1.7,4.5))
    ncomp = 4
    for ii = 1:ncomp
        ax = fig.add_subplot(4,1,ii)
        plot(svarX[:,ii], svarY[:,ii], c="magenta", marker="o", ms=0.5, linestyle="", label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")
        title("SVC $(ii)", fontsize=8)
        xticks(fontsize=8)
        yticks(fontsize=8)
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)

        if ii > 1
            xticks(color="None")
            yticks(color="None")
        end
        # xlim([-2, 1.5])
        # ylim([-4, 2.5])
        legend_elements = [Line2D([0], [0], color="black", lw=0, label=L"$\rho =$" * "$(round(svar_cor[ii], digits=2))")]
        legend(loc=4, handles=legend_elements, frameon=false, fontsize=8, handlelength=0)     
    end
    ax = fig.add_subplot(111)
    xlabel("FS cells", fontsize=12, labelpad=10)
    ylabel("Inh model", fontsize=12, labelpad=10)
    ax.set_frame_on(false)
    ax.tick_params(left=false, bottom=false)
    ax.tick_params(labelleft=false, labelbottom=false)
    
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_cor.png", bbox_inches="tight", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 50
    plot(collect(1:ncomp), svar_cor[1:ncomp], c="magenta", marker="o", ms=3)
    ylim([-0.1,1.1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("correlation " * L"($\rho$)", fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_corall.png", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 20
    plot(collect(1:ncomp), svarX_ksum[1:ncomp]/svarX_tsum, c="limegreen", marker="o", ms=3)
    plot(collect(1:ncomp), svarY_ksum[1:ncomp]/svarY_tsum, c="red", marker="o", ms=3)
    ylim([5e-4,1])
    yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("explained var", fontsize=12)
    legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="data (Pyr)"),
                       Line2D([0], [0], color="red", lw=2, label="model (Exc)")]
    legend(loc=1, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
    xticks(fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_expvar.png", dpi=300)
    
    
    
    figure(figsize=(3.5,3))
    ncomp = 20
    plot(collect(1:ncomp), cumsum(svarX_ksum[1:ncomp]/svarX_tsum), c="limegreen", marker="o", ms=3)
    plot(collect(1:ncomp), cumsum(svarY_ksum[1:ncomp]/svarY_tsum), c="red", marker="o", ms=3)
    # plot(collect(1:ncomp), cumsum(svarY_bal_ksum[1:ncomp]/svarY_bal_tsum), c="gray", linestyle="--")
    legend_elements = [Line2D([0], [0], color="limegreen", lw=2, label="Pyr"),
                       Line2D([0], [0], color="red", lw=2, label="Exc")]
    legend(loc=4, handles=legend_elements, frameon=false, fontsize=12, handlelength=1)     
    
    ylim([0,1])
    # yscale("log")
    xlabel("SVC", fontsize=12)
    ylabel("cumulative var", fontsize=12)
    xticks([0,5,10,15,20], fontsize=12)
    yticks(fontsize=12)
    tight_layout()
    
    savefig(dirfig * "SVC" * fname * "_cumvar.png", dpi=300)
    
    

end

