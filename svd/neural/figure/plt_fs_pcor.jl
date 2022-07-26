matplotlib = pyimport("matplotlib")
Line2D = matplotlib.lines.Line2D
cmap = matplotlib.cm

function plt_fs_pcor(dirfig, pcor_R, pcor_L, mse_R, mse_L, pcor_R_bal, pcor_L_bal, mse_R_bal, mse_L_bal)

############################
#----- trained nework -----#
############################
figure(figsize=(1.5,1.4))
ax = subplot(111)
# hist([pcor_R_bal; pcor_L_bal], bins=50, range=(-1,1), color="gray", alpha=0.5)
hist(pcor_R_bal, bins=50, range=(-1,1), color="gray", alpha=0.5)
cnt_R, bins_R = hist(pcor_R, bins=50, range=(-1,1), lw=0.8, color="blue", histtype="step")
cnt_L, bins_L = hist(pcor_L, bins=50, range=(-1,1), lw=0.8, color="red", histtype="step")
pmax_R = maximum(cnt_R)
pmax_L = maximum(cnt_L)
ptmpR = pcor_R
ptmpL = pcor_L
medR = mean(ptmpR[.!isnan.(ptmpR)])
medL = mean(ptmpL[.!isnan.(ptmpL)])
plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", lw=0.8, linestyle="--")
plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", lw=0.8, linestyle="--")

ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="Right"),
                   Line2D([0], [0], color="red", lw=1.5, label="Left"),
                   Line2D([0], [0], color="gray", lw=1.5, label="Null")]
legend(loc="upper left", bbox_to_anchor=[0,1], borderaxespad=0.0, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
xticks([-1,-0.5,0,0.5,1], ["-1", "-0.5", "0", "0.5", "1"], fontsize=7)
yticks(fontsize=7)

xlabel("correlation", fontsize=7)
ylabel("neuron count", fontsize=7)
tight_layout()

savefig(dirfig * "fs_pcor.png", dpi=600)
savefig(dirfig * "fs_pcor.pdf", dpi=600)




figure(figsize=(1.5,1.4))
ax = subplot(111)
hist(mse_R_bal, bins=50, range=(0,1), color="gray", alpha=0.5)
cnt_R, bins_R = hist(mse_R, bins=50, range=(0,1), lw=0.8, color="blue", histtype="step")
cnt_L, bins_L = hist(mse_L, bins=50, range=(0,1), lw=0.8, color="red", histtype="step")
vmax_R = maximum(cnt_R)
vmax_L = maximum(cnt_L)
vtmpR = mse_R
vtmpL = mse_L
medR = mean(vtmpR[.!isnan.(vtmpR)])
medL = mean(vtmpL[.!isnan.(vtmpL)])
plot(medR*ones(11), collect(0:vmax_R/10:vmax_R), color="blue", lw=0.8, linestyle="--")
plot(medL*ones(11), collect(0:vmax_L/10:vmax_L), color="red", lw=0.8, linestyle="--")

ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="Right"),
                   Line2D([0], [0], color="red", lw=1.5, label="Left"),
                   Line2D([0], [0], color="gray", lw=1.5, label="Null")]
legend(loc="upper left", bbox_to_anchor=[0,1], borderaxespad=0.0, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
xticks([0,0.5,1], fontsize=7)
yticks(fontsize=7)

xlabel("1 - nMSE", fontsize=7)
ylabel("count", fontsize=7)
tight_layout()

savefig(dirfig * "fs_mse.png", dpi=600)
savefig(dirfig * "fs_mse.pdf", dpi=600)



# ###########################################
# #----- balanced nework (not trained) -----#
# ###########################################
# figure(figsize=(1.7,1.5))
# ax = subplot(111)
# cnt_R, bins_R = hist(pcor_R_bal, bins=50, range=(-1,1), lw=0.8, color="blue", histtype="step")
# cnt_L, bins_L = hist(pcor_L_bal, bins=50, range=(-1,1), lw=0.8, color="red", histtype="step")
# pmax_R = maximum(cnt_R)
# pmax_L = maximum(cnt_L)
# ptmpR = pcor_R_bal
# ptmpL = pcor_L_bal
# medR = mean(ptmpR[.!isnan.(ptmpR)])
# medL = mean(ptmpL[.!isnan.(ptmpL)])
# plot(medR*ones(11), collect(0:pmax_R/10:pmax_R), color="blue", lw=0.8, linestyle="--")
# plot(medL*ones(11), collect(0:pmax_L/10:pmax_L), color="red", lw=0.8, linestyle="--")

# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="right"),
#                     Line2D([0], [0], color="red", lw=1.5, label="left")]
# legend(loc=2, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
# xticks([-1,-0.5,0,0.5,1], ["-1", "-0.5", "0", "0.5", "1"], fontsize=7)
# yticks(fontsize=7)

# xlabel("performance", fontsize=7)
# ylabel("count", fontsize=7)
# tight_layout()

# savefig(dirfig * "fs_pcor_bal.png", dpi=600)
# savefig(dirfig * "fs_pcor_bal.pdf", dpi=600)




# figure(figsize=(1.7,1.5))
# ax = subplot(111)
# cnt_R, bins_R = hist(mse_R_bal, bins=50, range=(0,1), lw=0.8, color="blue", histtype="step")
# cnt_L, bins_L = hist(mse_L_bal, bins=50, range=(0,1), lw=0.8, color="red", histtype="step")
# vmax_R = maximum(cnt_R)
# vmax_L = maximum(cnt_L)
# vtmpR = mse_R
# vtmpL = mse_L
# medR = mean(vtmpR[.!isnan.(vtmpR)])
# medL = mean(vtmpL[.!isnan.(vtmpL)])
# plot(medR*ones(11), collect(0:vmax_R/10:vmax_R), color="blue", lw=0.8, linestyle="--")
# plot(medL*ones(11), collect(0:vmax_L/10:vmax_L), color="red", lw=0.8, linestyle="--")

# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# legend_elements = [Line2D([0], [0], color="blue", lw=1.5, label="right"),
#                     Line2D([0], [0], color="red", lw=1.5, label="left")]
# legend(loc=2, handles=legend_elements, frameon=false, fontsize=7, handlelength=1)     
# xticks([0,0.5,1], fontsize=7)
# yticks(fontsize=7)

# xlabel("1 - nMSE", fontsize=7)
# ylabel("count", fontsize=7)
# tight_layout()

# savefig(dirfig * "fs_mse_bal.png", dpi=600)
# savefig(dirfig * "fs_mse_bal.pdf", dpi=600)




# figure(figsize=(3.5,3))
# ax = subplot(211)
# cntbal, binsbal = hist(pcor_bal50, bins=100, range=(-1,1), color="darkgray", log=true)
# cnt, bins = hist(pcor_rate50, bins=60, range=(-0.2,1), color="black", histtype="step", log=true)
# pmax = maximum(cntbal)
# plot(median(pcor_rate50)*ones(11), collect(0:pmax/10:pmax), color="tab:red", linestyle="--")

# legend_elements = [Line2D([0], [0], color="darkgray", lw=2, label="untrained"),
#                    Line2D([0], [0], color="black", lw=2, label="trained")]
# legend(handles=legend_elements, frameon=false, fontsize=10, handlelength=1, loc=2)     
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# ax.tick_params(labelbottom=false)
# yticks(fontsize=12)
# title(L"$N_{trial}=50$", fontsize=12)


# ax = subplot(212)
# cntbal, binsbal = hist(pcor_bal400, bins=100, range=(-1,1), color="darkgray", log=true)
# cnt, bins = hist(pcor_rate400, bins=60, range=(-0.2,1), color="black", histtype="step", log=true)
# pmax = maximum(cntbal)
# plot(median(pcor_rate400)*ones(11), collect(0:pmax/10:pmax), color="tab:red", linestyle="--")

# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# # ax.tick_params(labelbottom=false)
# xticks(fontsize=12)
# yticks(fontsize=12)
# xlabel("correlation", fontsize=12)
# ylabel("count", fontsize=12)
# title(L"$N_{trial}=400$", fontsize=12)

# tight_layout()

# savefig(dirfig * "fs_pcor.png", dpi=300)





# figure(figsize=(3.5,3))
# ax = subplot(211)
# cntbal, binsbal = hist(expvar_bal50, bins=50, range=(-1,1), color="darkgray", log=true)
# cnt, bins = hist(expvar_rate50, bins=50, range=(-1,1), color="black", histtype="step", log=true)
# pmax = maximum(cntbal)
# plot(median(expvar_rate50)*ones(11), collect(0:pmax/10:pmax), color="tab:red", linestyle="--")

# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# ax.tick_params(labelbottom=false)
# yticks(fontsize=12)
# title(L"$N_{trial}=50$", fontsize=12)


# ax = subplot(212)
# cntbal, binsbal = hist(expvar_bal400, bins=50, range=(-1,1), color="darkgray", log=true)
# cnt, bins = hist(expvar_rate400, bins=50, range=(-1,1), color="black", histtype="step", log=true)
# pmax = maximum(cntbal)
# plot(median(expvar_rate400)*ones(11), collect(0:pmax/10:pmax), color="tab:red", linestyle="--")

# legend_elements = [Line2D([0], [0], color="darkgray", lw=2, label="untrained"),
#                    Line2D([0], [0], color="black", lw=2, label="trained")]
# legend(handles=legend_elements, frameon=false, fontsize=10, handlelength=1, loc=2)     
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# # ax.tick_params(labelbottom=false)
# xticks(fontsize=12)
# yticks(fontsize=12)
# xlabel(L"$1 - nMSE$", fontsize=12)
# ylabel("count", fontsize=12)
# title(L"$N_{trial}=400$", fontsize=12)

# tight_layout()

# savefig(dirfig * "fs_expvar.png", dpi=300)






# figure(figsize=(3.5,3))
# ax = subplot(111)
# hist(pcor_bal400, bins=60, range=(-0.2,1), color="darkgray", log=true)
# hist(pcor_rate50, bins=60, range=(-0.2,1), color="black", histtype="step", log=true)
# hist(pcor_rate400, bins=60, range=(-0.2,1), color="black", histtype="step", log=true)
# legend_elements = [Line2D([0], [0], color="darkgray", lw=2, label="untrained"),
#                    Line2D([0], [0], color="black", lw=2, label="trained (spk)")]
# legend(handles=legend_elements, frameon=false, fontsize=12, handlelength=1, loc=1)     
# ax.spines["top"].set_visible(false)
# ax.spines["right"].set_visible(false)
# xticks(fontsize=12)
# yticks(fontsize=12)
# xlabel("correlation", fontsize=12)
# ylabel("count", fontsize=12)
# tight_layout()

# savefig(dirfig * "fs_pcor.png", dpi=300)


end