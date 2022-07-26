function plt_pyr_raster_arxiv(dirfig, dirsim, rtargPyr, almOrd, matchedCells)

    # load rate_navg
    rate_navg50 = load(dirsim * "rate_navg50.jld", "rate")
    rate_navg400 = load(dirsim * "rate_navg400.jld", "rate")
    
    # performance: Pyr - Exc
    nrep = 20
    pcor_navg400 = runperformance_nrep(rtargPyr, rate_navg400, almOrd, matchedCells, nrep)
    pcorSorted = reverse(sortperm(pcor_navg400))
    pcor_high = pcorSorted[pcor_navg400[pcorSorted] .> 0.9]
    pcor_med = pcorSorted[pcor_navg400[pcorSorted] .> 0.5]
    pcor_low = pcorSorted[pcor_navg400[pcorSorted] .> 0.2]

    # load spikes
    navg = 50
    Ncells = 5000
    tmp = load(dirsim * "times_rep1.jld", "times")
    times = zeros(size(tmp)[1], size(tmp)[2], navg)
    ns = zeros(Int,Ncells, navg)
    for avgi = 1:navg
        times[:,:,avgi] = load(dirsim * "times_rep$(avgi).jld", "times")
        ns[:,avgi] = load(dirsim * "ns_rep$(avgi).jld", "ns")
    end


    #--------------------------------------------#
    #----- spike raster of selected neurons -----#
    #--------------------------------------------#
    for kk = 1:10
        # nid = pcor_high[end-kk+1]
        # nid = pcor_med[end-kk+1]
        nid = pcor_low[end-kk+1]

        ci = matchedCells[nid]
        ci_alm = almOrd[nid]
        figure(figsize=(3.5,3))
        ax = subplot(211)
        for ii = 1:navg
            tidx = times[ci,:,ii] .> 1000
            plot(times[ci,tidx,ii], ii*ones(sum(tidx)), color="black", marker="o", ms=1, linestyle="")        
            xlim([1000,3000])
        end
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.spines["left"].set_visible(false)
        ax.spines["bottom"].set_visible(false)
        tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)

        ax = subplot(212)
        timev = 20*collect(1:100)
        plot(timev, rtargPyr[:,ci_alm], color="red", linewidth=1, alpha=1.0)    
        plot(timev, rate_navg50[:,ci,1], color="gray", linewidth=1, alpha=0.5)
        plot(timev, rate_navg400[:,ci,1], color="black", linewidth=1, alpha=1.0)    
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        xlim([0,2000])
        xticks(fontsize=12)
        yticks(fontsize=12)
        xlabel("time (ms)", fontsize=12)
        ylabel("spk rate (Hz)", fontsize=12)

        tight_layout()

        savefig(dirfig * "raster_pcor$(round(pcor_navg400[nid], digits=2))_$(ci).png", dpi=300)
        close()
    end



    # # print all Pyr / Exc 
    # for ngrp = 1:20

    #     println(ngrp)

    #     figure(figsize=(20,20))
    #     for i = 1:100
    #         subplot(10,10,i)
    #         idx = (ngrp-1)*100 + i
    #         nid = pcorSorted[idx]
    #         ci_alm = almOrd[nid]
    #         ci = matchedCells[nid]
    #         plot(timev, rtargPyr[:,ci_alm], color="red", linewidth=1, alpha=1.0)    
    #         plot(timev, rate_navg50[:,ci,1], color="gray", linewidth=1, alpha=0.5)
    #         plot(timev, rate_navg400[:,ci,1], color="black", linewidth=1, alpha=1.0)    
    #         title("cor: $(round(pcor_navg400[nid], digits=2))", fontsize=12)
    #         xticks(fontsize=10)
    #         yticks(fontsize=10)
    #     end
    #     tight_layout()

    #     savefig(dirfig * "pyr_usum$(ngrp).png", dpi=300)
    #     close()

    # end




end