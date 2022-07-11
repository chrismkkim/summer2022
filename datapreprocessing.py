import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def directory(num): # locate the data
    dir = ['150817', '150818', '150819', '150820', '150821','150825', '150826', '150827', '150828', '150831', '150901']
    postfile = ['post_150817_001_ch2-PnevPanResults-170808-190057.mat', 'post_150818_001_ch2-PnevPanResults-170808-180842.mat', 'post_150819_001_ch2-PnevPanResults-170815-163235.mat', 'post_150820_001_ch2-PnevPanResults-170808-185044.mat', 'post_150821_001-002_ch2-PnevPanResults-170808-184141.mat', 'post_150825_001-002-003_ch2-PnevPanResults-170814-191401.mat', 'post_150826_001_ch2-PnevPanResults-170808-002053.mat', 'post_150827_001_ch2-PnevPanResults-170807-171156.mat', 'post_150828_001-002_ch2-PnevPanResults-170807-204746.mat', 'post_150831_001-002_ch2-PnevPanResults-170807-193348.mat', 'post_150901_001_ch2-PnevPanResults-170807-072732.mat']
    morefile = ['more_150817_001_ch2-PnevPanResults-170808-190057.mat', 'more_150818_001_ch2-PnevPanResults-170808-180842.mat', 'more_150819_001_ch2-PnevPanResults-170815-163235.mat', 'more_150820_001_ch2-PnevPanResults-170808-185044.mat', 'more_150821_001-002_ch2-PnevPanResults-170808-184141.mat', 'more_150825_001-002-003_ch2-PnevPanResults-170814-191401.mat', 'more_150826_001_ch2-PnevPanResults-170808-002053.mat', 'more_150827_001_ch2-PnevPanResults-170807-171156.mat', 'more_150828_001-002_ch2-PnevPanResults-170807-204746.mat', 'more_150831_001-002_ch2-PnevPanResults-170807-193348.mat', 'more_150901_001_ch2-PnevPanResults-170807-072732.mat']
    os.chdir('C:\\Users\\mzgra\\Documents\\NIH\\data\\' + dir[num])
    post = scipy.io.loadmat(postfile[num])
    more = scipy.io.loadmat(morefile[num])
    return post, more

post17, more17 = directory(0) 
post18, more18 = directory(1)
post19, more19 = directory(2)
post20, more20 = directory(3)
post21, more21 = directory(4)
post25, more25 = directory(5)
post26, more26 = directory(6)
post31, more31 = directory(9)


def vars(post, more): # define variables needed
    outcomes = post['outcomes']
    outcomes = outcomes[0]
    coms = post['firstSideTryAl_COM']
    time = coms['time']
    traces = coms['traces']
    inhibit_cells = more['inhibitRois_pix'][0]
    allResp = post['allResp_HR_LR'][0]
    time = time[0][0][0]
    traces1 = traces[0][0] 
    return outcomes, time, traces1, inhibit_cells, allResp

outcomes17, time17, traces17, inhibit_cells17, allResp17 = vars(post17, more17)
outcomes18, time18, traces18, inhibit_cells18, allResp18 = vars(post18, more18)
outcomes19, time19, traces19, inhibit_cells19, allResp19 = vars(post19, more19)
outcomes20, time20, traces20, inhibit_cells20, allResp20 = vars(post20, more20)
outcomes21, time21, traces21, inhibit_cells21, allResp21 = vars(post21, more21)
outcomes25, time25, traces25, inhibit_cells25, allResp25 = vars(post25, more25)
outcomes26, time26, traces26, inhibit_cells26, allResp26 = vars(post26, more26)
outcomes31, time31, traces31, inhibit_cells31, allResp31 = vars(post31, more31)

def keepnums(outcomes, traces, allResp): # find which indicies to keep in each dataset
    trials = np.shape(traces)[2]
    nums = [] # gather number of trial indicies that have data
    for i in range(trials): # find all the indices of trials that have no NaNs
       if all(np.isnan(traces[:,0,i])) == False:
           nums.append(i)    
     
    keep = [] # keep all trial indices that were successful outcome
    for i in nums:
        if outcomes[i] == 1:
            keep.append(i)
            
    allResp = allResp[keep] # find the values of decisions made
    
    L = []
    R = []
    for i in range(len(keep)): # create new indicies of decisions 0 = L 1 = R
       if allResp[i] == 0:
           L.append(i)
       else:
           R.append(i)
           
    print('Trials with Data:', len(nums))
    print('Sucessful Trials:', len(keep))
    print('L:', len(L))
    print('R:', len(R))
    
    return keep, L, R


keep17, L17, R17 = keepnums(outcomes17, traces17, allResp17)
keep18, L18, R18 = keepnums(outcomes18, traces18, allResp18)
keep19, L19, R19 = keepnums(outcomes19, traces19, allResp19)
keep20, L20, R20 = keepnums(outcomes20, traces20, allResp20)
keep21, L21, R21 = keepnums(outcomes21, traces21, allResp21)
keep25, L25, R25 = keepnums(outcomes25, traces25, allResp25)
keep26, L26, R26 = keepnums(outcomes26, traces26, allResp26)
keep31, L31, R31 = keepnums(outcomes31, traces31, allResp31)


def avetrials(keep, L, R, traces): # average each dataset over the trials
    ntraces = traces[:,:,keep] # set new traces to be all the traces trials that have data
    Ltraces = ntraces[:,:,L]
    Rtraces = ntraces[:,:,R]
    Lavetrial = Ltraces.mean(axis=2)
    Ravetrial = Rtraces.mean(axis=2)
    return Lavetrial, Ravetrial
    
Lave17, Rave17 = avetrials(keep17, L17, R17, traces17)
Lave18, Rave18 = avetrials(keep18, L18, R18, traces18)
Lave19, Rave19 = avetrials(keep19, L19, R19, traces19)
Lave20, Rave20 = avetrials(keep20, L20, R20, traces20)
Lave21, Rave21 = avetrials(keep21, L21, R21, traces21)
Lave25, Rave25 = avetrials(keep25, L25, R25, traces25)
Lave26, Rave26 = avetrials(keep26, L26, R26, traces26)
Lave31, Rave31 = avetrials(keep31, L31, R31, traces31)

def heatmap(Lave, Rave): # heatmaps of individual datasets
    Ladata = pd.DataFrame(Lave)
    Radata = pd.DataFrame(Rave)
    Lad = Ladata.transpose()
    Rad = Radata.transpose()
    
    plt.figure()
    ax = plt.axes()
    sns.heatmap(Lad, ax=ax, cmap = "Blues", vmax = 0.1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Neuron')
    ax.set_title('L averaged data')
      
    plt.clf
    plt.figure()
    ax = plt.axes()
    sns.heatmap(Rad, ax=ax, cmap = "Blues", vmax = 0.1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Neuron')
    ax.set_title('R averaged data')
    
heatmap(Lave17, Rave17)
heatmap(Lave18, Rave18)
heatmap(Lave19, Rave19)
heatmap(Lave20, Rave20)
heatmap(Lave21, Rave21)
heatmap(Lave25, Rave25)
heatmap(Lave26, Rave26)
heatmap(Lave31, Rave31)

    
def PCA(Lave, Rave, time, inhibit_cells): #PCA analysis of individual datasets
    from sklearn.decomposition import PCA
    L = pd.DataFrame(Lave)
    R = pd.DataFrame(Rave)
    Lt = L.transpose()
    Rt = R.transpose()
    # R
    pca = PCA()
    pf = pca.fit(Rt)
    comp = pca.components_

    # excitatory / inhibitory
    excidx = (inhibit_cells == 0)
    inhidx = (inhibit_cells == 1)
    
    tnp = Rt.to_numpy()
    texc = tnp[excidx,:]
    tinh = tnp[inhidx,:]
    
    pca_exc = PCA()
    pf_exc = pca_exc.fit(texc)
    comp_exc = pca_exc.components_
    
    pca_inh = PCA()
    pf_inh = pca_inh.fit(tinh)
    comp_inh = pca_inh.components_
    
    plt.figure()
    plt.plot(time, comp_exc[0,:], label = 'Excitatory')
    plt.plot(time, comp_inh[0,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[0]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[0]))
    plt.xlabel("PCA 1 R")
     
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 R")
    
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 R")
    
    # L
    pca = PCA()
    pf = pca.fit(Lt)
    comp = pca.components_

    # excitatory / inhibitory
    excidx = (inhibit_cells == 0)
    inhidx = (inhibit_cells == 1)
    
    tnp = Lt.to_numpy()
    texc = tnp[excidx,:]
    tinh = tnp[inhidx,:]
    
    pca_exc = PCA()
    pf_exc = pca_exc.fit(texc)
    comp_exc = pca_exc.components_
    
    pca_inh = PCA()
    pf_inh = pca_inh.fit(tinh)
    comp_inh = pca_inh.components_
    
    plt.figure()
    plt.plot(time, comp_exc[0,:], label = 'Excitatory')
    plt.plot(time, comp_inh[0,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[0]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[0]))
    plt.xlabel("PCA 1 L")
    
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 L")
    
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 L")


PCA(Lave17, Rave17, time17, inhibit_cells17)
PCA(Lave18, Rave18, time18, inhibit_cells18)
PCA(Lave19, Rave19, time19, inhibit_cells19)
PCA(Lave20, Rave20, time20, inhibit_cells20)
PCA(Lave21, Rave21, time21, inhibit_cells21)
PCA(Lave25, Rave25, time25, inhibit_cells25)
PCA(Lave26, Rave26, time26, inhibit_cells26)
PCA(Lave31, Rave31, time31, inhibit_cells31)

    
def smoothdata(Lave, Rave, traces): # smooth each dataset neuron
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        #gaussian = np.exp(-(time/sigma)**2/2)
        y_smooth = np.convolve(y, box, mode = 'same')
        return y_smooth
    
    t = np.shape(traces)[0]
    n = np.shape(traces)[1]
    smoothdataL = np.zeros((t, n)) 
    for i in range(n):
            smoothdataL[:,i] = smooth(Lave[:,i], 5)
            
    smoothdataR = np.zeros((t, n)) 
    for i in range(n):
            smoothdataR[:,i] = smooth(Rave[:,i], 5)
    return smoothdataL, smoothdataR

# figure out how to throw out neruons that dont look like anything

smoothdataL17, smoothdataR17 = smoothdata(Lave17, Rave17, traces17)
smoothdataL18, smoothdataR18 = smoothdata(Lave18, Rave18, traces18)
smoothdataL19, smoothdataR19 = smoothdata(Lave19, Rave19, traces19)
smoothdataL20, smoothdataR20 = smoothdata(Lave20, Rave20, traces20)
smoothdataL21, smoothdataR21 = smoothdata(Lave21, Rave21, traces21)
smoothdataL25, smoothdataR25 = smoothdata(Lave25, Rave25, traces25)
smoothdataL26, smoothdataR26 = smoothdata(Lave26, Rave26, traces26)
smoothdataL31, smoothdataR31 = smoothdata(Lave31, Rave31, traces31)


def reducetime(time, smoothdataL, smoothdataR): # retrieve the reduced time steps for spikes
    b = (time >= -938.52) & (time <= 938.52)
    dataL = smoothdataL[b,:]
    dataR = smoothdataR[b,:]
    time = time[b] # reduced time steps for combination plots
    return dataL, dataR, time

dataL17, dataR17, time = reducetime(time17, smoothdataL17, smoothdataR17)
dataL18, dataR18, time = reducetime(time18, smoothdataL18, smoothdataR18)
dataL19, dataR19, time = reducetime(time19, smoothdataL19, smoothdataR19)
dataL20, dataR20, time = reducetime(time20, smoothdataL20, smoothdataR20)
dataL21, dataR21, time = reducetime(time21, smoothdataL21, smoothdataR21)
dataL25, dataR25, time = reducetime(time25, smoothdataL25, smoothdataR25)
dataL26, dataR26, time = reducetime(time26, smoothdataL26, smoothdataR26)
dataL31, dataR31, time = reducetime(time31, smoothdataL31, smoothdataR31)



finalL = np.hstack((dataL17, dataL18, dataL19, dataL20, dataL21, dataL25, dataL26, dataL31)) # combine all the neurons on reduced time steps
finalR = np.hstack((dataR17, dataR18, dataR19, dataR20, dataR21, dataR25, dataR26, dataR31))
finalinhibit = np.hstack((inhibit_cells17, inhibit_cells18, inhibit_cells19, inhibit_cells20, inhibit_cells21, inhibit_cells25, inhibit_cells26, inhibit_cells31))

def combinedPCA(finalL, finalR, finalinhibit, time): #PCA analysis of the combined data 
    from sklearn.decomposition import PCA
    L = pd.DataFrame(finalL)
    R = pd.DataFrame(finalR)
    Lt = L.transpose()
    Rt = R.transpose()
    # R
    pca = PCA()
    pf = pca.fit(Rt)
    comp = pca.components_

    # excitatory / inhibitory
    excidx = (finalinhibit == 0)
    inhidx = (finalinhibit == 1)
    
    tnp = Rt.to_numpy()
    texc = tnp[excidx,:]
    tinh = tnp[inhidx,:]
    
    pca_exc = PCA()
    pf_exc = pca_exc.fit(texc)
    comp_exc = pca_exc.components_
    
    pca_inh = PCA()
    pf_inh = pca_inh.fit(tinh)
    comp_inh = pca_inh.components_
    
    plt.figure()
    plt.plot(time, comp_exc[0,:], label = 'Excitatory')
    plt.plot(time, comp_inh[0,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[0]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[0]))
    plt.xlabel("PCA 1 R")
     
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 R")
    
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 R")
    
    # L
    pca = PCA()
    pf = pca.fit(Lt)
    comp = pca.components_

    # excitatory / inhibitory
    excidx = (finalinhibit == 0)
    inhidx = (finalinhibit == 1)
    
    tnp = Lt.to_numpy()
    texc = tnp[excidx,:]
    tinh = tnp[inhidx,:]
    
    pca_exc = PCA()
    pf_exc = pca_exc.fit(texc)
    comp_exc = pca_exc.components_
    
    pca_inh = PCA()
    pf_inh = pca_inh.fit(tinh)
    comp_inh = pca_inh.components_
    
    plt.figure()
    plt.plot(time, comp_exc[0,:], label = 'Excitatory')
    plt.plot(time, comp_inh[0,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[0]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[0]))
    plt.xlabel("PCA 1 L")
    
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 L")
    
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 L")

    
combinedPCA(finalL, finalR, finalinhibit, time)
