function plt_inh_raster(dirfig, dirsim, fs_rate_R, fs_rate_L, pcor_R, pcor_L, mse_R, mse_L, cellsMatched)

    Ne = 2500
    Ncells = 5000

    # load rate_navg
    rate_R_navg50 = load(dirsim * "rate_R_navg50.jld", "rate")[:,Ne+1:Ncells,:]
    rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[:,Ne+1:Ncells,:]
    rate_L_navg50 = load(dirsim * "rate_L_navg50.jld", "rate")[:,Ne+1:Ncells,:]
    rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[:,Ne+1:Ncells,:]
    
    # performance: Pyr - Exc
    pcorMatched = (pcor_R + pcor_L)/2
    pcorSorted = reverse(sortperm(pcorMatched))

    mseMatched = (mse_R + mse_L)/2
    mseSorted = reverse(sortperm(mseMatched))

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
    times_R = times_R[Ne+1:Ncells,:,:]
    times_L = times_L[Ne+1:Ncells,:,:]
    ns_R = ns_R[Ne+1:Ncells,:]
    ns_L = ns_L[Ne+1:Ncells,:]

    #--------------------------------------------#
    #----- spike raster of selected neurons -----#
    #--------------------------------------------#
    for kk = 1:length(mseSorted)
    # for kk = 1:10
        println(kk)

        nid = mseSorted[kk]

        ci = Int(cellsMatched[nid])
        ci_alm = nid
        fig = figure(figsize=(0.8,1.3))
        gs = gridspec.GridSpec(2, 1, figure=fig, wspace=0, hspace=0.001, left=0.3, bottom=0.2)

        ratemax = maximum([rate_R_navg400[:,ci,1][:]; fs_rate_R[:,ci_alm][:]; rate_L_navg400[:,ci,1][:]; fs_rate_L[:,ci_alm][:]])

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
        # title("R: $(round(pcor_R[nid],digits=2))", fontsize=8)


        ax = fig.add_subplot(gs[2,1])
        timev = 20*collect(1:100)
        plot(timev, rate_R_navg400[:,ci,1], color="blue", linewidth=0.8, alpha=1.0)    
        plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)   
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        xlim([0,2000])
        ylim([0, ratemax])
        xticks([0,1000,2000],["-2","-1","0"],fontsize=8)
        yticks(fontsize=8)
        # tick_params(labelleft=false)
        # ax.set_xlabel("time (s)", fontsize=8)

        savefig(dirfig * "raster/inh/inh_$(kk)" * "_R$(round(pcor_R[nid],digits=2))" * "_L$(round(pcor_L[nid],digits=2))" * ".png", dpi=600)
        savefig(dirfig * "raster/inh/inh_$(kk)" * "_R$(round(pcor_R[nid],digits=2))" * "_L$(round(pcor_L[nid],digits=2))" * ".pdf", dpi=600)
        close()
    end





end