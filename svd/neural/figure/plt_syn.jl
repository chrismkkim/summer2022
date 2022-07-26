function plt_syn(dirfig, p, ubal_R, uebal_R, uibal_R, uplas_R)


    timev = p.dt*collect(1:p.Nsteps)

    for nid = 1:10


        # figure(figsize=(2.3,1.5))
        figure(figsize=(1.7,1.5))
        axvspan(p.stim_on, p.stim_off, color="dodgerblue", ec="None", alpha=0.3)
        plot(timev, ubal_R[:,nid][:], c="black", linewidth=0.1, alpha=1)
        plot(timev, uebal_R[:,nid][:], c="limegreen", linewidth=0.1, alpha=1)
        plot(timev, uibal_R[:,nid][:], c="darkorange", linewidth=0.1, alpha=1)
        plot(timev, uplas_R[:,nid][:], c="cyan", linewidth=0.8, alpha=1)
        # plot(timev, uext_R[:,nid][:], c="black", linewidth=0.5, alpha=1)
        plot(timev, ones(length(timev)), c="magenta", linewidth=0.8, linestyle="--", alpha=1)
        # legend_elements = [ Line2D([0], [0], color="blue", lw=1, label=L"$u_{bal}^E$"),    
        # Line2D([0], [0], color="black", lw=1, label=L"$u_{bal}^S$"),    
        # Line2D([0], [0], color="darkorange", lw=1, label=L"$u_{plas}$"),
        # Line2D([0], [0], color="red", lw=1, label=L"$u_{bal}^I$")]
        # legend(handles=legend_elements, frameon=false, fontsize=8, handlelength=1, bbox_to_anchor=(1.02,1.02), loc="upper left")     
    
        xticks([1000, 2000, 3000], [-2,-1,0],fontsize=8)
        yticks(fontsize=8)
        xlim([200,3000])
        ylim([-6,6])
        xlabel("time (s)", fontsize=8)
        ylabel("syn input", fontsize=8)
        tight_layout()

        savefig(dirfig * "u_R$(nid).png", dpi=600)
        savefig(dirfig * "u_R$(nid).pdf", dpi=600)
        close()
    end
    


end