import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pdb
 
os.chdir('/home/grabelmz/FN_dataSharing/data/mouse1_fni16/150818')
os.getcwd()
post = scipy.io.loadmat('post_150818_001_ch2-PnevPanResults-170808-180842.mat')
more = scipy.io.loadmat('more_150818_001_ch2-PnevPanResults-170808-180842.mat')
#date = scipy.io.loadmat('150817_001.mat')
#data.keys()
outcomes = post['outcomes']
outcomes = outcomes[0]

def PCA(post, more):
    coms = post['firstSideTryAl_COM']
    time = coms['time']
    traces = coms['traces']
    inhibit_cells = more['inhibitRois_pix']
    time = time[0][0][0]
    traces1 = traces[0][0] # extract the traces proper format
    np.shape(traces1) # 68x420x84
    timedim = np.shape(traces1)[0]
    neurons = np.shape(traces1)[1]
    trials = np.shape(traces1)[2]
    print(trials)
    nums = [] # gather number of trial indicies. it is same regardless of neuron
    for i in range(trials): # find all the indices of trials that have no NaNs
        if all(np.isnan(traces1[:,0,i])) == False:
            nums.append(i)    
    print(nums)
  
    keep = []
    for i in nums:
        if outcomes[i] == 1:
            keep.append(i)
        
    allResp = post['allResp_HR_LR'][0]
    print(allResp[keep])
    print(keep)
    print(outcomes)
    print(outcomes[keep])
    
    ntraces = traces1[:,:,keep] # set new traces to be all the traces trials that have data
    np.shape(ntraces) 
        
    # Average each neuron over the 14 trials at each timestep
    
    avetrial = ntraces.mean(axis=2) #average over the trial axis with 14 trials
    np.shape(avetrial) # (68x420)!
    avetrial[:,4] # neuron 5 for all time
    np.shape(avetrial[:,4]) # 68 time steps
    plt.plot(time, avetrial[:,1]) # plot 1 averaged neuron
    plt.plot(time, ntraces[:,1,1])
    plt.xlabel('Time')
    plt.ylabel('Activity')
    
    plt.plot(time, avetrial[:,0:5]) # or a few
    plt.xlabel('Time')
    plt.ylabel('Activity')
    
    # convert final data to dataframe to export to excel
    averageddata = pd.DataFrame(avetrial)
    # # every row is a time step
    # averageddata.index = np.round(time,2)
    averageddata.iloc[:, 0:20]
    # every column is a neuron
    
    
    # averageddata.to_excel('traces.xlsx') # export data to excel
    
    
    plt.figure(figsize=(6,10))
    
    for i in range(len(keep)):
        tracesdf = pd.DataFrame(ntraces[:,:,i])
        tracesdf.index = np.round(time,2)
        tt = tracesdf.transpose()
        plt.subplot(10,6,i+1)
        sns.heatmap(tt, xticklabels = False, yticklabels = False, cmap = "Blues", vmax = 0.1)

    plt.savefig('heatmap_all.png')

    # transpose for heatmap
    t = averageddata.transpose()
    np.shape(t)
    
    ax = plt.axes()
    sns.heatmap(t, ax=ax, cmap = "Blues", vmax = 0.1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Neuron')
    
    plt.savefig('heatmap_avg.png')
    
    
    
    #PCA output time trajectories first dim neuron second time
    from sklearn.decomposition import PCA
    import plotly.express as px
    
    pca = PCA()
    pf = pca.fit(t)
    comp = pca.components_
    plt.plot(time, comp[2,:])
    
    pca.explained_variance_ratio_

    # excitatory / inhibitory
    inhcells_copy = np.copy(inhibit_cells)
    excidx = (inhibit_cells == 0)[0,:]
    inhidx = (inhibit_cells == 1 )[0,:]

    # pdb.set_trace()

    tnp = t.to_numpy()
    texc = tnp[excidx,:]
    tinh = tnp[inhidx,:]

    pca_exc = PCA()
    pf_exc = pca_exc.fit(texc)
    comp_exc = pca_exc.components_

    pca_inh = PCA()
    pf_inh = pca_inh.fit(tinh)
    comp_inh = pca_inh.components_

    plt.figure()
    plt.plot(time, comp_exc[0,:])
    plt.plot(time, comp_inh[0,:])
    plt.savefig('pc_excinh0.png')
    


    plt.figure()
    plt.plot(time, comp_exc[1,:])
    plt.plot(time, comp_inh[1,:])
    plt.savefig('pc_excinh1.png')
    


    plt.figure()
    plt.plot(time, comp_exc[2,:])
    plt.plot(time, comp_inh[2,:])
    plt.savefig('pc_excinh2.png')
    


PCA(post, more)