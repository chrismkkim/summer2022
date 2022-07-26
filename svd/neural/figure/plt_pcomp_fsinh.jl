function plt_pcomp_fsinh(dirfig, fname, startind, Mfs, Minh)

    timev_pyr = 20*collect(1:100)
    timev_exc = 20*collect(startind:100)
    
    fig=figure(figsize=(1.0,0.7))
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.2, hspace=0.1, left=0, bottom=0, right=1, top=0.9)
    ncomp = 6
    # sgn_fs = [-1,-1,-1,-1]
    if fname == "fsinhR"
        sgn_fs = [-1, 1, -1, 1, -1, -1]
        sgn_inh = [-1, 1, -1, 1, 1, 1]
    elseif fname == "fsinhL"
        sgn_fs = [-1, 1, -1, -1, 1, 1]
        sgn_inh = [-1, 1, 1, -1, 1, 1]
    end
    
    for ii = 1:ncomp
        ax = fig.add_subplot(gs[ii])
        # ax = fig.add_subplot(2,3,ii)
        ax.set_frame_on(false)
        ax.tick_params(left=false, bottom=false)
        ax.tick_params(labelleft=false, labelbottom=false)
        title("PC $(ii)", fontsize=5, pad=-1)
        xticks(color="None", fontsize=0)
        yticks(color="None", fontsize=0)
        xlim([0, 2000])
        ylim([-0.25, 0.25])
    
        plot(timev_exc, sgn_inh[ii]*Minh.proj[:,ii], c="darkorange", lw=1.1)
        plot(timev_pyr, sgn_fs[ii]*Mfs.proj[:,ii], c="black", lw=1.1, linestyle="--")
        
    end
    
    # tight_layout()
    
    # savefig("pca_" * fname * ".png", dpi=600)
    savefig(dirfig * "pca_" * fname * ".pdf", dpi=600)
    
    
    
    
    
    ncomp = 9
    expvar_exc = cumsum(principalvars(Minh))[1:ncomp] / tvar(Minh)
    expvar_pyr = cumsum(principalvars(Mfs))[1:ncomp] / tvar(Mfs)
    
    fig = figure(figsize=(1.4,1.3))
    ax = fig.add_subplot(111)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    plot(collect(1:ncomp), expvar_exc, color="darkorange", marker="o", lw=1.1, ms=3, mec="None", clip_on=false)
    plot(collect(1:ncomp), expvar_pyr, color="black", marker="o", lw=1.1, ms=3, mec="None", clip_on=false, linestyle="-")
    # plot(collect(1:ncomp), expvar_exc, color="limegreen", lw=1.1, clip_on=false)
    # plot(collect(1:ncomp), expvar_pyr, color="black", lw=1.1, clip_on=false, linestyle="--")
    legend_elements = [Line2D([0], [0], color="darkorange", lw=0, label="Inhibitory"),
                       Line2D([0], [0], color="black", lw=0, label="Fast spiking")]
    legend(handles=legend_elements, frameon=false, fontsize=7, handlelength=0, loc=4, labelcolor="linecolor")     
    
    xlabel("PC", fontsize=7)
    ylabel("cumulative var", fontsize=7)
    xticks([1, 3, 5, 7, 9], fontsize=7)
    yticks(fontsize=7)
    ylim([0,1])
    tight_layout()
    # savefig(dirfig * fname * "_cumvar.png", dpi=600)
    # savefig(dirfig * "pca_cumvar_" * fname * ".png", dpi=600)
    # savefig("pca_cumvar_" * fname * ".png", dpi=600)
    savefig(dirfig * "pca_cumvar_" * fname * ".pdf", dpi=600)
    
    
    end
    