using JLD2
using PyPlot
using Statistics

dirNetwork = "/home/grabelmz/code/distributedActivity-main/data_network/"
dirALM = "/home/grabelmz/code/distributedActivity-main/data_alm/"
dirsimTrialAvg = "/data/grabelmz/trainedNetwork/sim/trialAveraged/"

almOrd = load(dirNetwork * "almOrd.jld", "almOrd")
matchedCells = load(dirNetwork * "matchedCells.jld", "matchedCells")

rtarg_lickright = load(dirALM * "target16_er.jld", "targetShifted_16_exc_right")
rtarg_lickright_inh = load(dirALM * "target16_ir.jld", "targetShifted_16_inh_right")
rtarg_lickleft_inh = load(dirALM * "target16_il.jld", "targetShifted_16_inh_left")
rtarg_lickleft = load(dirALM * "target16_el.jld", "targetShifted_16_exc_left")
selectCellsRight = load(dirNetwork * "selectCellsRight.jld", "selectCellsRight")
selectCellsLeft = load(dirNetwork * "selectCellsLeft.jld", "selectCellsLeft")

# neural data
rtarg_lickright = rtarg_lickright[:,selectCellsRight]
rtarg_lickleft = rtarg_lickleft[:,selectCellsLeft]

# rnn model
trained_lickright = dropdims(load(dirsimTrialAvg * "rate_R_navg400.jld2", "rate"), dims=3)
trained_lickleft = dropdims(load(dirsimTrialAvg * "rate_L_navg400.jld2", "rate"), dims=3)

# sort the real and model neurons
rtarg_lickright_sorted = rtarg_lickright[:,almOrd]
rtarg_lickleft_sorted = rtarg_lickleft[:,almOrd]

trained_lickright_sorted = trained_lickright[:,matchedCells]
trained_lickleft_sorted = trained_lickleft[:,matchedCells]

dt_data = 33.0
timev_data = dt_data * collect(1:59)

dt_model = 20.0
timev_model = dt_model * collect(1:100)


nid_list = rand(1:2500,16)
figure(figsize=(10,10))
title("lick right", fontsize=20)
for ci = 1:16
    subplot(4,4,ci)
    nid = nid_list[ci]
    plot(timev_data, rtarg_lickright_sorted[:,nid], color="black", label="data")
    plot(timev_model, trained_lickright_sorted[:,nid], color="blue", label="model")
    if ci == 1
        legend(frameon=false)
    end
    xlabel("time (ms)")
    ylabel("spike rate (Hz)")
end
tight_layout()

savefig("fig_trainedmodel_exc_right.png")


figure(figsize=(10,10))
title("lick left", fontsize=20)
for ci = 1:16
    subplot(4,4,ci)
    nid = nid_list[ci]
    plot(timev_data, rtarg_lickleft_sorted[:,nid], color="black", label="data")
    plot(timev_model, trained_lickleft_sorted[:,nid], color="red", label="model")
    if ci == 1
        legend(frameon=false)
    end
    xlabel("time (ms)")
    ylabel("spike rate (Hz)")
end
tight_layout()

savefig("fig_trainedmodel_exc_left.png")



figure(figsize=(10,10))
for ci = 1:16
    subplot(4,4,ci)
    mid = rand(1:500)
    # plot(timev_data, rtarg_lickright_inh[:,mid], color="red")

    shiftIndexToInhibitory = 2500
    nid = rand(1:2500) + shiftIndexToInhibitory # add 2500 to select inhibitory neurons
    plot(timev_model, trained_lickright[:,nid,1], color="blue", label="lick right")
    plot(timev_model, trained_lickleft[:,nid,1], color="red", label="lick left")
    if ci == 1
        legend(frameon=false)
    end
    xlabel("time (ms)")
    ylabel("spike rate (Hz)") 
end
tight_layout()

savefig("fig_trainedmodel_inh.png")


# PCA of 
# (1) exc data - lick right / lick left
# (2) inh data - lick right / lick left
# (3) exc model - lick right / lick left
# (4) inh model - lick right / lick left <===== interesting result

# Does the PCA of inh model match with PCA of inh data???




# center the data and model activity
# exc model
train_lickleft_tmp = trained_lickleft[:, 1:2500]
train_lickright_tmp = trained_lickright[:, 1:2500]
train_lickleft_meanzero = train_lickleft_tmp .- mean(train_lickleft_tmp, dims=1)
train_lickright_meanzero = train_lickright_tmp .- mean(train_lickright_tmp, dims=1)

# data
rtarg_lickleft_exc_meanzero = rtarg_lickleft_sorted .- mean(rtarg_lickleft_sorted, dims=1)
rtarg_lickright_exc_meanzero = rtarg_lickright_sorted .- mean(rtarg_lickright_sorted, dims=1)


dt_data = 33.0
dt_model = 20.0

using MultivariateStats, Plots
#---- exc model -----
npca = 10
# Left
pca_model_exc_left = fit(PCA, train_lickleft_meanzero, maxoutdim = npca)
pc_model_exc_left = projection(pca_model_exc_left::PCA)
expvar_model_exc_left = principalvars(pca_model_exc_left) / var(pca_model_exc_left)

# Right
pca_model_exc_right = fit(PCA, train_lickright_meanzero, maxoutdim = npca)
pc_model_exc_right = projection(pca_model_exc_right::PCA)
expvar_model_exc_right = principalvars(pca_model_exc_right) / var(pca_model_exc_right)

#---- exc data -----
# Left 
pca_data_exc_left = fit(PCA, rtarg_lickleft_exc_meanzero, maxoutdim = npca)
pc_data_exc_left = projection(pca_data_exc_left::PCA)
expvar_data_exc_left = principalvars(pca_data_exc_left) / var(pca_data_exc_left)

# Right
pca_data_exc_right = fit(PCA, rtarg_lickright_exc_meanzero, maxoutdim = npca)
pc_data_exc_right = projection(pca_data_exc_right::PCA)
expvar_data_exc_right = principalvars(pca_data_exc_right) / var(pca_data_exc_right)

# Exc
#PCA 1
sgn = -1
figure(); 
plot(collect(1:100) * dt_model, sgn * pc_model_exc_left[:,1], color = "red", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_exc_left[:,1], color = "black", label = "Data");
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC1: Lick Left Exc"); legend(frameon=false)
savefig("PCA1ExcL")

figure();
plot(collect(1:100) * dt_model, pc_model_exc_right[:,1], color = "blue", label = "Model"); 
plot(collect(1:59) * dt_data , pc_data_exc_right[:,1], color = "black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC1: Lick Right Exc"); legend(frameon=false)
savefig("PCA1ExcR")

#PCA 2
figure(); 
plot(collect(1:100) * dt_model, pc_model_exc_left[:,2], color = "red", label = "Model"); 
plot(collect(1:59) * dt_data , pc_data_exc_left[:,2], color = "black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC2: Lick Left Exc"); legend(frameon=false)
savefig("PCA2ExcL")

figure(); 
plot(collect(1:100) * dt_model, sgn * pc_model_exc_right[:,2], color = "blue", label = "Model"); 
plot(collect(1:59) * dt_data , pc_data_exc_right[:,2], color = "black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC2: Lick Right Exc"); legend(frameon=false)
savefig("PCA2ExcR")

#PCA 3
figure(); 
plot(collect(1:100) * dt_model , sgn * pc_model_exc_left[:,3], color = "red", label = "Model"); 
plot(collect(1:59) * dt_data , pc_data_exc_left[:,3], color = "black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC3: Lick Left Exc"); legend(frameon=false)
savefig("PCA3ExcR")

figure(); 
plot(collect(1:100) * dt_model, sgn * pc_model_exc_right[:,3], color = "blue", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_exc_right[:,3], color = "black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC3: Lick Right Exc"); legend(frameon=false)
savefig("PCA3ExcL")

# explained variance
figure()
plot(collect(1:npca), cumsum(expvar_model_exc_left), marker="o", color="red", label="Model")
plot(collect(1:npca), cumsum(expvar_data_exc_left), marker="o", color="black", label="Data")
xlabel("PC"); ylabel("Cumulative Variance"); title("Explained Variance: Exc Left"); legend(frameon=false)
savefig("ExpvarExcL")

figure()
plot(collect(1:npca), cumsum(expvar_model_exc_right), marker="o", color="blue", label="Model")
plot(collect(1:npca), cumsum(expvar_data_exc_right), marker="o", color="black", label="Data")
xlabel("PC"); ylabel("Cumulative Variance"); title("Explained Variance: Exc Right"); legend(frameon=false)
savefig("ExpvarExcR")



# center the data and model activity
# model
train_lickleft_tmp = trained_lickleft[:, 2501:5000]
train_lickright_tmp = trained_lickright[:, 2501:5000]
train_lickleft_meanzero = train_lickleft_tmp .- mean(train_lickleft_tmp, dims=1)
train_lickright_meanzero = train_lickright_tmp .- mean(train_lickright_tmp, dims=1)

# data
rtarg_lickleft_inh_meanzero = rtarg_lickleft_inh .- mean(rtarg_lickleft_inh, dims=1)
rtarg_lickright_inh_meanzero = rtarg_lickright_inh .- mean(rtarg_lickright_inh, dims=1)


#---- inh model -----#
npca = 10
# Left 
pca_model_inh_left = fit(PCA, train_lickleft_meanzero, maxoutdim = npca)
pc_model_inh_left = projection(pca_model_inh_left::PCA)
expvar_model_inh_left = principalvars(pca_model_inh_left) / var(pca_model_inh_left)

# Right
pca_model_inh_right = fit(PCA, train_lickright_meanzero, maxoutdim = npca)
pc_model_inh_right = projection(pca_model_inh_right::PCA)
expvar_model_inh_right = principalvars(pca_model_inh_right) / var(pca_model_inh_right)

#----- inh data -----#
# Left
pca_data_inh_left = fit(PCA, rtarg_lickleft_inh_meanzero, maxoutdim = npca)
pc_data_inh_left = projection(pca_data_inh_left::PCA)
expvar_data_inh_left = principalvars(pca_data_inh_left) / var(pca_data_inh_left)

# Right
pca_data_inh_right = fit(PCA, rtarg_lickright_inh_meanzero, maxoutdim = npca)
pc_data_inh_right = projection(pca_data_inh_right::PCA)
expvar_data_inh_right = principalvars(pca_data_inh_right) / var(pca_data_inh_right)


# Inh
#PCA 1
sgn = -1
figure(); 
plot(collect(1:100) * dt_model, pc_model_inh_left[:,1], color="red", label = "Model"); 
plot(collect(1:59) * dt_data, sgn * pc_data_inh_left[:,1], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC1: Lick Left Inh"); legend(frameon=false)
savefig("PCA1InhL")

figure(); 
plot(collect(1:100) * dt_model, pc_model_inh_right[:,1], color="blue", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_inh_right[:,1], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC1: Lick Right Inh"); legend(frameon=false)
savefig("PCA1InhR")

#PCA 2
figure(); 
plot(collect(1:100) * dt_model, pc_model_inh_left[:,2], color="red", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_inh_left[:,2], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC2: Lick Left Inh"); legend(frameon=false)
savefig("PCA2InhL")

figure(); 
plot(collect(1:100) * dt_model, sgn * pc_model_inh_right[:,2], color="blue", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_inh_right[:,2], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC2: Lick Right Inh"); legend(frameon=false)
savefig("PCA2InhR")

#PCA 3
figure(); 
plot(collect(1:100) * dt_model, pc_model_inh_left[:,3], color="red", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_inh_left[:,3], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC3: Lick Left Inh"); legend(frameon=false)
savefig("PCA3InhL")

figure(); 
plot(collect(1:100) * dt_model, sgn * pc_model_inh_right[:,3], color="blue", label = "Model"); 
plot(collect(1:59) * dt_data, pc_data_inh_right[:,3], color="black", label = "Data"); 
xlabel("Time (ms)"); ylabel("Neural Activity"); title("PC3: Lick Right Inh"); legend(frameon=false)
savefig("PCA3InhR")


# explained variance
figure()
plot(collect(1:npca), cumsum(expvar_model_inh_left), marker="o", color="red", label="Model")
plot(collect(1:npca), cumsum(expvar_data_inh_left), marker="o", color="black", label="Data")
xlabel("PC"); ylabel("Cumulative Variance"); title("Explained Variance: Inh Left"); legend(frameon=false)
savefig("ExpvarInhL")

figure()
plot(collect(1:npca), cumsum(expvar_model_inh_right), marker="o", color="blue", label="Model")
plot(collect(1:npca), cumsum(expvar_data_inh_right), marker="o", color="black", label="Data")
xlabel("PC"); ylabel("Cumulative Variance"); title("Explained Variance: Inh Right"); legend(frameon=false)
savefig("ExpvarInhR")
