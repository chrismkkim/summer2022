function plt_cormat_RL(dirfig, fname, cor_fsinh_R, svd_fsinh_R, cor_fsinh_L, svd_fsinh_L)

    svnorm_R = svd_fsinh_R.S / sum(svd_fsinh_R.S)
    svnorm_L = svd_fsinh_L.S / sum(svd_fsinh_L.S)
    
    ###################################################
    #----- singular values of correlation matrix -----#
    ###################################################
    fig = figure(figsize=(1.5,1.4))
    ax = fig.add_subplot(111)
    nsing = 20
    plot(collect(1:nsing), svnorm_R[1:nsing], color="blue", marker="o", ms=2, lw=0.8)
    plot(collect(1:nsing), svnorm_L[1:nsing], color="red", marker="o", ms=2, lw=0.8)
    ax.set_yscale("log")
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    xlabel("rank", fontsize=7)
    ylabel("singular value", fontsize=7)
    xticks(fontsize=7)
    yticks(fontsize=7)
    ylim([3e-3,1])
    legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="right"),
                       Line2D([0], [0], color="red", lw=1.5, label="left")]
    legend(handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
    tight_layout()
    
    savefig(dirfig * "cormat_sv_" * fname * ".png", dpi=600)
    savefig(dirfig * "cormat_sv_" * fname * ".pdf", dpi=600)
    
    
    ################################
    #----- correlation matrix -----#
    ################################
    # lick right
    sortrow_R = sortperm(svd_fsinh_R.U[:,1])
    fig = figure(figsize=(1.6,1.4))
    ax = fig.add_subplot(111)
    imshow(cor_fsinh_R[sortrow_R,:], cmap="bwr", vmin=-1,vmax=1, interpolation="None", aspect="auto")
    xticks(fontsize=7)
    yticks(fontsize=7)
    xlabel("Inh neuron", fontsize=7)
    ylabel("FS cell", fontsize=7)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.tick_params(labelsize=8)
    tight_layout()
    
    savefig(dirfig * "cormat_" * fname * "_R.png", dpi=600)
    savefig(dirfig * "cormat_" * fname * "_R.pdf", dpi=600)
    
    # lick left
    sortrow_L = sortperm(svd_fsinh_L.U[:,1])
    fig = figure(figsize=(1.6,1.4))
    ax = fig.add_subplot(111)
    imshow(cor_fsinh_L[sortrow_L,:], cmap="bwr", vmin=-1,vmax=1, interpolation="None", aspect="auto")
    xticks(fontsize=7)
    yticks(fontsize=7)
    xlabel("Inh neuron", fontsize=7)
    ylabel("FS cell", fontsize=7)
    cbar = colorbar(ticks=[-1,0,1])
    cbar.ax.tick_params(labelsize=8)
    tight_layout()
    
    savefig(dirfig * "cormat_" * fname * "_L.png", dpi=600)
    savefig(dirfig * "cormat_" * fname * "_L.pdf", dpi=600)
    
    
    end
    
    