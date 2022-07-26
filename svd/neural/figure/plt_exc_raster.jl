function plt_exc_raster(dirfig, dirsim, rate_R, rate_L, pcor_R, pcor_L, almOrd, matchedCells)

    Ne = 2500
    Ncells = 5000

    # load rate_navg
    rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,1:Ne,:]
    rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[:,1:Ne,:]
    
    # performance: Pyr - Exc
    pcorMatched = (pcor_R + pcor_L)/2
    pcorSorted = reverse(sortperm(pcorMatched))

    # load spikes
    navg = 50
    Ncells = 5000
    tmp = load(dirsim * "times_R_1.jld", "times")
    times_R = zeros(size(tmp)[1], size(tmp)[2], navg)
    times_L = zeros(size(tmp)[1], size(tmp)[2], navg)
    ns_R = zeros(Int,Ncells, navg)
    ns_L = zeros(Int,Ncells, navg)
    for avgi = 1:navg
        times_R[:,:,avgi] = load(dirsim * "times_R_$(avgi).jld", "times")
        times_L[:,:,avgi] = load(dirsim * "times_L_$(avgi).jld", "times")
        ns_R[:,avgi] = load(dirsim * "ns_R_$(avgi).jld", "ns")
        ns_L[:,avgi] = load(dirsim * "ns_L_$(avgi).jld", "ns")
    end
    times_R = times_R[1:Ne,:,:]
    times_L = times_L[1:Ne,:,:]
    ns_R = ns_R[1:Ne,:]
    ns_L = ns_L[1:Ne,:]

    #--------------------------------------------#
    #----- spike raster of selected neurons -----#
    #--------------------------------------------#
    # for kk = 1:length(pcorSorted)
    for kk = 1:300
        println(kk)

        nid = pcorSorted[kk]

        ci = Int(matchedCells[nid])
        ci_alm = Int(almOrd[nid])
        fig = figure(figsize=(0.8,1.3))
        gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0.001, left=0.3, bottom=0.2)

        ratemax = maximum([rate_R_navg400[:,ci,1][:]; rate_R[:,ci_alm][:]; rate_L_navg400[:,ci,1][:]; rate_L[:,ci_alm][:]])

        ax = fig.add_subplot(gs[1,1])
        for ii = 1:navg
            tidx_R = times_R[ci,:,ii] .> 1000
            plot(times_R[ci,tidx_R,ii], (ii+50)*ones(sum(tidx_R)), color="blue", marker="o", ms=0.5, mec="None", linestyle="")        
            xlim([1000,3000])
        end
        for ii = 1:navg
            tidx_L = times_L[ci,:,ii] .> 1000
            plot(times_L[ci,tidx_L,ii], ii*ones(sum(tidx_L)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
            xlim([1000,3000])
        end
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["bottom"].set_visible(false)
        tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
        # title("L:$(round(pcor_L[nid],digits=2)), R:$(round(pcor_R[nid],digits=2))", fontsize=8)


        ax = fig.add_subplot(gs[2,1])
        timev = 20*collect(1:100)
        # plot(timev, rate_R_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)
        plot(timev, rate_R_navg400[:,ci,1], color="blue", linewidth=0.8, alpha=1.0)  
        plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)    
        # plot(timev, rate_R[:,ci_alm], color="limegreen", linewidth=1.0, alpha=1.0, linestyle="-")    
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        xlim([0,2000])
        ylim([0, ratemax])
        xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
        yticks(fontsize=8)
        # ax.set_xlabel("time (s)", fontsize=8)
        # ax.set_ylabel("spk / s ", fontsize=8)


        # ax = fig.add_subplot(gs[1,1])
        # for ii = 1:navg
        #     tidx_L = times_L[ci,:,ii] .> 1000
        #     plot(times_L[ci,tidx_L,ii], ii*ones(sum(tidx_L)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
        #     xlim([1000,3000])
        # end
        # ax.spines["top"].set_visible(false)
        # ax.spines["right"].set_visible(false)
        # ax.spines["left"].set_visible(false)
        # ax.spines["bottom"].set_visible(false)
        # tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
        # title("L: $(round(pcor_L[nid],digits=2))", fontsize=8)


        # ax = fig.add_subplot(gs[2,1])
        # timev = 20*collect(1:100)
        # # plot(timev, rate_L_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)        
        # plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)  
        # # plot(timev, rate_L[:,ci_alm], color="limegreen", linewidth=1.0, alpha=1.0, linestyle="-")      
        # ax.spines["top"].set_visible(false)
        # ax.spines["right"].set_visible(false)
        # xlim([0,2000])
        # ylim([0, ratemax])
        # xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
        # yticks(fontsize=8)
        # ax.set_xlabel("time (s)", fontsize=8)
        # ax.set_ylabel("spk / s ", fontsize=8)

        # gs.update(left=0.3, bottom=0.3)
        # gs.tight_layout(fig)
        # savefig(dirfig * "fs_raster_pcor$(round(pcorMatched[nid], digits=2))_$(ci).png", dpi=300)
        savefig(dirfig * "raster/exc/exc_$(kk)" * "_L$(round(pcor_L[nid],digits=2))_R$(round(pcor_R[nid],digits=2))" * ".png", dpi=600)
        savefig(dirfig * "raster/exc/exc_$(kk)" * "_L$(round(pcor_L[nid],digits=2))_R$(round(pcor_R[nid],digits=2))" * ".pdf", dpi=600)
        close()
    end





end






    # #--------------------------------------------#
    # #----- spike raster of selected neurons -----#
    # #--------------------------------------------#
    # # for kk = 1:length(pcorSorted)
    # for kk = 1:5
    #     println(kk)

    #     nid = pcorSorted[kk]

    #     ci = Int(matchedCells[nid])
    #     ci_alm = Int(almOrd[nid])
    #     fig = figure(figsize=(1.5,1.5))
    #     gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.25, hspace=0.001, left=0.3, bottom=0.3)

    #     ratemax = maximum([rate_R_navg400[:,ci,1][:]; rate_R[:,ci_alm][:]; rate_L_navg400[:,ci,1][:]; rate_L[:,ci_alm][:]])

    #     ax = fig.add_subplot(gs[1,2])
    #     for ii = 1:navg
    #         tidx_R = times_R[ci,:,ii] .> 1000
    #         plot(times_R[ci,tidx_R,ii], ii*ones(sum(tidx_R)), color="blue", marker="o", ms=0.5, mec="None", linestyle="")        
    #         xlim([1000,3000])
    #     end
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     ax.spines["left"].set_visible(false)
    #     ax.spines["bottom"].set_visible(false)
    #     tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    #     title("R: $(round(pcor_R[nid],digits=2))", fontsize=8)


    #     ax = fig.add_subplot(gs[2,2])
    #     timev = 20*collect(1:100)
    #     # plot(timev, rate_R_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)
    #     plot(timev, rate_R_navg400[:,ci,1], color="blue", linewidth=0.8, alpha=1.0)    
    #     plot(timev, rate_R[:,ci_alm], color="limegreen", linewidth=1.0, alpha=1.0, linestyle="-")    
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     xlim([0,2000])
    #     ylim([0, ratemax])
    #     xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
    #     tick_params(labelleft=false)
    #     ax.set_xlabel("time (s)", fontsize=8)


    #     ax = fig.add_subplot(gs[1,1])
    #     for ii = 1:navg
    #         tidx_L = times_L[ci,:,ii] .> 1000
    #         plot(times_L[ci,tidx_L,ii], ii*ones(sum(tidx_L)), color="red", marker="o", ms=0.5, mec="None", linestyle="")        
    #         xlim([1000,3000])
    #     end
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     ax.spines["left"].set_visible(false)
    #     ax.spines["bottom"].set_visible(false)
    #     tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
    #     title("L: $(round(pcor_L[nid],digits=2))", fontsize=8)


    #     ax = fig.add_subplot(gs[2,1])
    #     timev = 20*collect(1:100)
    #     # plot(timev, rate_L_navg50[:,ci,1], color="gray", linewidth=0.8, alpha=0.5)        
    #     plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)  
    #     plot(timev, rate_L[:,ci_alm], color="limegreen", linewidth=1.0, alpha=1.0, linestyle="-")      
    #     ax.spines["top"].set_visible(false)
    #     ax.spines["right"].set_visible(false)
    #     xlim([0,2000])
    #     ylim([0, ratemax])
    #     xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
    #     yticks(fontsize=8)
    #     ax.set_xlabel("time (s)", fontsize=8)
    #     ax.set_ylabel("spk / s ", fontsize=8)

    #     # gs.update(left=0.3, bottom=0.3)
    #     # gs.tight_layout(fig)
    #     # savefig(dirfig * "fs_raster_pcor$(round(pcorMatched[nid], digits=2))_$(ci).png", dpi=300)
    #     savefig(dirfig * "raster/exc/exc_$(kk).png", dpi=600)
    #     savefig(dirfig * "raster/exc/exc_$(kk).pdf", dpi=600)
    #     close()
    # end

