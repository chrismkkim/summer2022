function plt_pyrexc_psth(dirfig, dirsim, rate_R, rate_L, pcor_R, pcor_L, almOrd, matchedCells, lickRL)

    Ne = 2500
    Ncells = 5000

    start_idx = 5
    rate_R = rate_R[start_idx:end,:]
    rate_L = rate_L[start_idx:end,:]

    # load rate_navg
    rate_R_navg400 = load(dirsim * "rate_R_navg400.jld", "rate")[start_idx:end,1:Ne,end]
    rate_L_navg400 = load(dirsim * "rate_L_navg400.jld", "rate")[start_idx:end,1:Ne,end]
    
    # smooth the network model PSTH, then normalize by the max rate
    rate_R_navg400_smooth = funMovAvg_2d(copy(rate_R_navg400), 1)
    rate_L_navg400_smooth = funMovAvg_2d(copy(rate_L_navg400), 1)
    
    # # normalize by max rate
    # rate_R_navg400 = rate_R_navg400 ./ maximum(rate_R_navg400_smooth, dims=1)
    # rate_L_navg400 = rate_L_navg400 ./ maximum(rate_L_navg400_smooth, dims=1)
    # rate_R = rate_R ./ maximum(rate_R, dims=1)
    # rate_L = rate_L ./ maximum(rate_L, dims=1)

    # performance: Pyr - Exc
    # pcorMatched = (pcor_R + pcor_L)/2
    pcorMatched = copy(pcor_R)
    pcorSorted = reverse(sortperm(pcorMatched))
    
    # shuffle pcorSorted
    # Random.seed!(1234)
    # pcorSorted = shuffle(pcorSorted)
    ngrp = 10
    pcorSorted_array = zeros(Int,ngrp,16)
    for grpi = 1:ngrp
        for coli = 1:4
            numCell = length(pcorSorted)
            numInc = 200 #Int(floor(0.15*numCell))
            start_ind = (coli-1)*numInc + 1
            end_ind = coli*numInc 
            rndCells = sort(shuffle(collect(1:numInc))[1:4])
            pcorSorted_array[grpi, (coli-1)*4+1:coli*4] = pcorSorted[start_ind : end_ind][rndCells]
        end
    end

    ncol = 3
    nrow = 3
    ncells = ncol * nrow
    # ngrp = floor(Int,length(pcor_R) / ncells) + 1
        
    #----------- plot lick right ------------#
    for grpi = 1:ngrp #6
        println(grpi)

        fig = figure(figsize=(2.2,1.5))
        gs = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.7, hspace=0.3, left=0.2, right = 0.95, bottom=0.25, top=0.95)

        for kk = 1:ncells

            # celli = (grpi-1)*ncells + kk
            # nid = pcorSorted[celli]
            nid = pcorSorted_array[grpi, kk]

            ci = Int(matchedCells[nid])
            ci_alm = Int(almOrd[nid])

            ratemax = maximum([7; rate_R_navg400[:,ci,1][:]; rate_R[:,ci_alm][:]; rate_L_navg400[:,ci,1][:]; rate_L[:,ci_alm][:]])

            coli = kk % ncol
            if coli == 0 coli = ncol end
            rowi = floor(Int,(kk-1)/ncol)+1
            ax = fig.add_subplot(gs[rowi,coli])

            timev = 20*collect(start_idx:100)
            # if lickRL == "R"                
                plot(timev, rate_R_navg400[:,ci,1], color="dodgerblue", linewidth=0.8, alpha=1.0)                    
                plot(timev, rate_R[:,ci_alm], color="black", linewidth=0.8, alpha=1.0, linestyle="-")    
                title("$(round(pcor_R[nid], digits=2))", fontsize=6, pad=-0.5)
            ymax = 5*(floor(maximum([maximum(rate_L_navg400[:,ci,1]), maximum(rate_L[:,ci_alm]), maximum(rate_R_navg400[:,ci,1]), maximum(rate_R[:,ci_alm])])/5) + 1)
            # elseif lickRL == "L"
            #     plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)  
            #     plot(timev, rate_L[:,ci_alm], color="black", linewidth=0.8, alpha=1.0, linestyle="-")   
            #     title("$(round(pcor_L[nid], digits=2))", fontsize=6, pad=-0.5)   
            # end            
            ax.spines["top"].set_visible(false)
            ax.spines["right"].set_visible(false)
            ax.spines["left"].set_linewidth(0.5)
            ax.spines["bottom"].set_linewidth(0.5)
            ax.xaxis.set_tick_params(width=0.5, length=2)
            ax.yaxis.set_tick_params(width=0.5, length=2)
            xlim([0,2000])
            ylim([0, ymax+1])
            if kk == 7
                xticks([0, 1000,2000], [-2,-1,0], fontsize=7)
                yticks([ymax], fontsize=7)
                xlabel("time (s)", fontsize=7)
                ylabel("spk/s", fontsize=7)
            else
                xticks([1000,2000], fontsize=0)
                yticks([ymax], fontsize=7)
                tick_params(labelbottom=false)
            end            

        end

        # savefig("fig_pyrexc$(grpi)_R.png", dpi=600)
        savefig(dirfig * "raster/exc/fig_pyrexc$(grpi)_R.pdf", dpi=600)
        close()

    end



    #----------- plot lick left ------------#
    for grpi = 1:ngrp #6
        println(grpi)

        fig = figure(figsize=(2.2,1.5))
        gs = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.7, hspace=0.3, left=0.2, right = 0.95, bottom=0.25, top=0.95)

        for kk = 1:ncells

            # celli = (grpi-1)*ncells + kk
            # nid = pcorSorted[celli]
            nid = pcorSorted_array[grpi, kk]

            ci = Int(matchedCells[nid])
            ci_alm = Int(almOrd[nid])

            ratemax = maximum([7; rate_R_navg400[:,ci,1][:]; rate_R[:,ci_alm][:]; rate_L_navg400[:,ci,1][:]; rate_L[:,ci_alm][:]])

            coli = kk % ncol
            if coli == 0 coli = ncol end
            rowi = floor(Int,(kk-1)/ncol)+1
            ax = fig.add_subplot(gs[rowi,coli])

            timev = 20*collect(start_idx:100)
            # if lickRL == "R"                
            #     plot(timev, rate_R_navg400[:,ci,1], color="dodgerblue", linewidth=0.8, alpha=1.0)                    
            #     plot(timev, rate_R[:,ci_alm], color="black", linewidth=0.8, alpha=1.0, linestyle="-")    
            #     title("$(round(pcor_R[nid], digits=2))", fontsize=6, pad=-0.5)
            # elseif lickRL == "L"
                plot(timev, rate_L_navg400[:,ci,1], color="red", linewidth=0.8, alpha=1.0)  
                plot(timev, rate_L[:,ci_alm], color="black", linewidth=0.8, alpha=1.0, linestyle="-")   
                title("$(round(pcor_L[nid], digits=2))", fontsize=6, pad=-0.5)   
            ymax = 5*(floor(maximum([maximum(rate_L_navg400[:,ci,1]), maximum(rate_L[:,ci_alm]), maximum(rate_R_navg400[:,ci,1]), maximum(rate_R[:,ci_alm])])/5) + 1) 
            # end
            
            ax.spines["top"].set_visible(false)
            ax.spines["right"].set_visible(false)
            ax.spines["left"].set_linewidth(0.5)
            ax.spines["bottom"].set_linewidth(0.5)
            ax.xaxis.set_tick_params(width=0.5, length=2)
            ax.yaxis.set_tick_params(width=0.5, length=2)
            xlim([0,2000])
            ylim([0, ymax+1])
            if kk == 7
                xticks([0, 1000,2000], [-2,-1,0], fontsize=7)
                yticks([ymax], fontsize=7)
                xlabel("time (s)", fontsize=7)
                ylabel("spk/s", fontsize=7)
                # title("$(round(pcor_L[nid], digits=2))", fontsize=7, pad=-0.9)   
            else
                xticks([1000,2000], fontsize=0)
                yticks([ymax], fontsize=7)
                tick_params(labelbottom=false)
            end            

        end

        # savefig("fig_pyrexc$(grpi)_L.png", dpi=600)
        savefig(dirfig * "raster/exc/fig_pyrexc$(grpi)_L.pdf", dpi=600)
        close()

    end    





end