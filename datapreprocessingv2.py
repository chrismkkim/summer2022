###########################################################################################################################################
###########################################################################################################################################
#
# Data Preprocessing for http://repository.cshl.edu/id/eprint/36980/ the post and more data.
#
# Functions: dates, vars, keepnums, avetrials, heatmap, PCA, smoothdata, reducetime, combinedPCA, EI
# dates:
# takes two numbered dates as input e.g. '150817'. Links to directory where these numbered files are located. Will have to change directory
# based on where you have the data. retrieves the post and more file dicionaries and links them together for output along with the folder names.
# vars:
# takes the data dictionary from dates and the day as input. Outputs the outcomes, time, traces, inhibitory neurons, and decision made.
# puts each days data in another dicrionary keyed by the day. Must start iterating over the folder names to input and output the functions 
# keepnums:
# takes the outcomes, traces, and decision individual days as input. Outputs keep, L, and R. keep denotes all the indices of sucessful trials
# in that dataset that have data not NaN. L and R are the indices of neurons where the L or R decision was made by the mouse. It also prints
# out the numbers for a visual reference.
# avetrials:
# averages the the trials of the neurons over time. First gathers all the traces data with trial data from keep and splits into L and R traces
# averages these traces individually over the trail axis. However if either L or R is empty it just returns the empty traces data to output.
# heatmap:
# creates heatmaps for the averaged data split by the decision L and R. First converts data to dataframe and transposes it for plotting on 
# the correct axes. Uncomment the code below the function to plot the heatmap. Select days individually to see plots saved
# PCA:
# performs PCA analysis on the Lave and Rave data. Also takes time and inhibit_cells as input to determine the excitatory and inhibitory neurons 
# and the time axis for plotting. To plot all days PCA plots uncomment the code below the function. You will have to again change the names 
# automatically to get all plots otherwise it will override each iteration. Or just run one day you want to see in particular 
# smoothdata:
# smooths the data by convolving the averaged data by a set number of boxes. outputs the smoothed data in both L and R
# reducetime:
# reduces the timesteps to the minimum of all the data which was 59. selects the min range from the time variable and subsets the data to
# include only the neurons in that time frame to combine all days together later.
# combinedPCA:
# performs the same PCA analysis as the PCA function except once on the combined data.
# EI: 
# splits the combined L and R data into excitatory and inhibitory neuron datasets for the training data for the NN. Outputs to csv
#
# post and more datafiles are dictioaries holding data on specific variables. 
# when changing mouse change the directory to the mouse# and fni# the days in that directory and the timec variable day below.
# use pdb.set_trace() anywhere in the code to debug or manually run code in the linux bash console in biowulf
###########################################################################################################################################
###########################################################################################################################################

from bdb import set_trace
import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pdb

def dates(date1, date2):
    folders = os.listdir('/data/grabelmz/FN_dataSharing/data/mouse1_fni16') # have to change directory to mouse2_fni17 or mouse1_fni16 
    ifolders = []
    for i in range(len(folders)): # append folder names as an int
        n = int(folders[i])
        ifolders.append(n)
        
    sfolders = sorted(ifolders) # make sure they are sorted in acending order
    i1 = sfolders.index(date1)
    i2 = sfolders.index(date2)
    rfolders = sfolders[i1:i2+1] # get the subset of folders inputted
    rfolders = [str(x) for x in rfolders] # convert back to string
    datadic = {}
    for i in rfolders: # for each date folder retireve and store post and more dicionaries into another dictionary.
        dir = '/data/grabelmz/FN_dataSharing/data/mouse1_fni16/' + i # have to change directory to mouse2_fni17 or mouse1_fni16 
        print(i)
        files = os.listdir(dir)
        post = scipy.io.loadmat(dir + '/' + files[0])
        more = scipy.io.loadmat(dir + '/' + files[1])
        post.update(more)
        datadic[i] = post
    return datadic, rfolders

datadic, folders = dates(150817, 151029)  # all days fni16 150817-151029. all days fni17 150814-151001

def vars(datadic, day): # define variables needed
    outcomes = datadic[day]['outcomes'] # data on if a trial was sucessful 
    outcomes = outcomes[0]
    coms = datadic[day]['firstSideTryAl_COM'] # data on time and traces of neurons
    time = coms['time']
    traces = coms['traces']
    inhibit_cells = datadic[day]['inhibitRois_pix'][0] # data on inbitory vs excitatory neurons
    allResp = datadic[day]['allResp_HR_LR'][0] # data on decision mouse made. Left or Right. High or Low
    time = time[0][0][0]
    traces1 = traces[0][0] 
    return outcomes, time, traces1, inhibit_cells, allResp

outcomes = {}
time = {}
traces = {}
inhibit_cells = {}
allResp = {}
for i in folders: # store the vars needed in dictionaries based on the day as a key.
    outcomes[i], time[i], traces[i], inhibit_cells[i], allResp[i] = vars(datadic, i)


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
           L.append(i) #low freq
       else:
           R.append(i) #high freq
           
    print('Trials with Data:', len(nums))
    print('Sucessful Trials:', len(keep))
    print('L:', len(L))
    print('R:', len(R))
    
    return keep, L, R

keep = {}
L = {}
R = {}
for i in folders:
    keep[i], L[i], R[i] = keepnums(outcomes[i], traces[i], allResp[i])



def avetrials(keep, L, R, traces): # average each dataset over the trials
    ntraces = traces[:,:,keep] # set new traces to be all the traces trials that have data
    Ltraces = ntraces[:,:,L]
    Rtraces = ntraces[:,:,R]
    if len(L) == 0:
        Lavetrial = Ltraces
    else:
        Lavetrial = Ltraces.mean(axis=2)
    if len(R) == 0:
        Ravetrial = Rtraces
    else:
        Ravetrial = Rtraces.mean(axis=2)
    
    return Lavetrial, Ravetrial

Lave = {}
Rave = {}
for i in folders:
    Lave[i], Rave[i] = avetrials(keep[i], L[i], R[i], traces[i])


#pdb.set_trace()
def heatmap(Lave, Rave): # heatmaps of individual datasets

    if np.size(Lave) != 0:
        Ladata = pd.DataFrame(Lave)
        Lad = Ladata.transpose()
        plt.figure()
        ax = plt.axes()
        sns.heatmap(Lad, ax=ax, cmap = "Blues", vmax = 0.1)
        ax.set_xlabel('Time')
        ax.set_ylabel('Neuron')
        ax.set_title('L averaged data')
        #plt.savefig('heatmapLave.png')
        
    plt.clf
    if np.size(Rave) != 0:    
        Radata = pd.DataFrame(Rave)
        Rad = Radata.transpose()       
        plt.figure()
        ax = plt.axes()
        sns.heatmap(Rad, ax=ax, cmap = "Blues", vmax = 0.1)
        ax.set_xlabel('Time')
        ax.set_ylabel('Neuron')
        ax.set_title('R averaged data')
        #plt.savefig('heatmapRave.png') # will override save files unless names are changed with each iteration
        
# turn on for heatmaps   
#for i in folders:
#    heatmap(Lave[i], Rave[i])


#pdb.set_trace()

    
def PCA(Lave, Rave, time, inhibit_cells): #PCA analysis of individual datasets
    from sklearn.decomposition import PCA

    # R
    if np.size(Rave) != 0: # if R is not empty run PCA
        
        R = pd.DataFrame(Rave)
        Rt = R.transpose() 
        
        pca = PCA()
        pf = pca.fit(Rt)
        comp = pca.components_

        # excitatory / inhibitory neurons index
        excidx = (inhibit_cells == 0) 
        inhidx = (inhibit_cells == 1)
        
        tnp = Rt.to_numpy() # back to np
        texc = tnp[excidx,:] # index the neruons by exh and inh
        tinh = tnp[inhidx,:]
        
        pca_exc = PCA() # perform PCA on exc and inh
        pf_exc = pca_exc.fit(texc)
        comp_exc = pca_exc.components_ # get the components 
        
        pca_inh = PCA()
        pf_inh = pca_inh.fit(tinh)
        comp_inh = pca_inh.components_
        
        plt.figure() # plot the first PC by time
        plt.plot(time, comp_exc[0,:], label = 'Excitatory')
        plt.plot(time, comp_inh[0,:], label = 'Inhibitory')
        plt.legend()
        plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[0]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[0]))
        plt.xlabel("PCA 1 R")
        
        plt.figure() # plot the second PC by time
        plt.plot(time, comp_exc[1,:], label = 'Excitatory')
        plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
        plt.legend()
        plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
        plt.xlabel("PCA 2 R")
        
        plt.figure() # plot the third PC by time
        plt.plot(time, comp_exc[2,:], label = 'Excitatory')
        plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
        plt.legend()
        plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
        plt.xlabel("PCA 3 R")
        
    # L
    if np.size(Lave) != 0: # if L is not empty do the same thing as above.
        
        L = pd.DataFrame(Lave)
        Lt = L.transpose()

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

# turn on for PCA plots
#for i in folders:
#    PCA(Lave[i], Rave[i], time[i], inhibit_cells[i])
    
def smoothdata(Lave, Rave, traces): # smooth each dataset neuron
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode = 'same')
        return y_smooth
    
    t = np.shape(traces)[0] # time dim
    n = np.shape(traces)[1] # neuron dim
    if np.size(Lave) != 0:    # if non empty smooth the neurons
        smoothdataL = np.zeros((t, n)) 
        for i in range(n):
            smoothdataL[:,i] = smooth(Lave[:,i], 5)
    elif np.size(Lave) == 0: # else return nothing
        smoothdataL = np.zeros((t, n))
    if np.size(Rave) != 0:    
        smoothdataR = np.zeros((t, n)) 
        for i in range(n):
            smoothdataR[:,i] = smooth(Rave[:,i], 5)
    elif np.size(Rave) == 0:
        smoothdataR = np.zeros((t, n))
    return smoothdataL, smoothdataR

smoothdataL = {}
smoothdataR = {}
for i in folders: # iterate by day in the dicionary. 
    smoothdataL[i], smoothdataR[i] = smoothdata(Lave[i], Rave[i], traces[i])
 
    
def reducetime(time, smoothdataL, smoothdataR): # retrieve the reduced time steps for spikes
    b = (time >= -938.52) & (time <= 938.52)
    dataL = smoothdataL[b,:]
    dataR = smoothdataR[b,:]
    time = time[b] # reduced time steps for combination plots
    return dataL, dataR, time

dataL = {}
dataR = {}
timec = {} # time variable is the same regardless of day. We just need one day for plotting 
for i in folders:
    dataL[i], dataR[i], timec[i] = reducetime(time[i], smoothdataL[i], smoothdataR[i])


for i in folders: # delete the days where either L or R is empty or both are. We need the data to be symmetric  
    if sum(sum(dataL[i])) == 0 or sum(sum(dataR[i])) == 0:
        del dataL[i]
        del dataR[i]
    
finalL = np.hstack(dataL[i] for i in sorted(dataL)) # combine all nonempty days into one dataframe
finalR = np.hstack(dataR[i] for i in sorted(dataR))
finalinhibit = np.hstack(inhibit_cells[i] for i in sorted(dataL)) # similarly combine all the exc and inh data with the same days

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
    #plt.savefig('combinedpcaR1m2l.png')
     
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 R")
    #plt.savefig('combinedpcaR2m2l.png')
        
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 R")
    #plt.savefig('combinedpcaR3m2l.png')
        
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
    #plt.savefig('combinedpcaL1m2l.png')
    
    plt.figure()
    plt.plot(time, comp_exc[1,:], label = 'Excitatory')
    plt.plot(time, comp_inh[1,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[1]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[1]))
    plt.xlabel("PCA 2 L")
    #plt.savefig('combinedpcaL2m2l.png')
        
    plt.figure()
    plt.plot(time, comp_exc[2,:], label = 'Excitatory')
    plt.plot(time, comp_inh[2,:], label = 'Inhibitory')
    plt.legend()
    plt.title('Explained Var Exc: {}'.format(pca_exc.explained_variance_ratio_[2]) + '\nExplained Var Inh: {}'.format(pca_inh.explained_variance_ratio_[2]))
    plt.xlabel("PCA 3 L")    
    #plt.savefig('combinedpcaL3m2l.png')
    

combinedPCA(finalL, finalR, finalinhibit, timec['150817'])


def EI(finalL, finalR, finalinhibit): # split combined data into exc and inh datasets
    excidx = (finalinhibit == 0)
    inhidx = (finalinhibit == 1)
    EL = finalL[:, excidx]
    IL = finalL[:, inhidx]
    ER = finalR[:, excidx]
    IR = finalR[:, inhidx]
    # ELm = EL.mean(axis=0)
    # ILm = IL.mean(axis=0)
    # ERm = ER.mean(axis=0)
    # IRm = IR.mean(axis=0)
    return EL, IL, ER, IR

EL, IL, ER, IR = EI(finalL, finalR, finalinhibit) 
EL = pd.DataFrame(EL)
IL = pd.DataFrame(IL)
ER = pd.DataFrame(ER)
IR = pd.DataFrame(IR)

EL.to_csv('target_fni16_exc_left.csv')
IL.to_csv('target_fni16_inh_left.csv')
ER.to_csv('target_fni16_exc_right.csv')
IR.to_csv('target_fni16_inh_right.csv')


pdb.set_trace()

# other analysis

#identify neurons with activity > 0.1
n = np.shape(finalL)[1]
idxt = {}
ind = [] 
for i in range(n):
    idxt[i] = np.where(finalL[:,i] > 0.1)
    if len(idxt[i][0]) != 0:
        ind.append(i)

finalL_avg = np.mean(finalL, axis=0)
finalL_avg_log = np.log10(finalL_avg)
finalL_avg_log = finalL_avg_log[~np.isnan(finalL_avg_log)]

finalR_avg = np.mean(finalR, axis=0)
finalR_avg_log = np.log10(finalR_avg)
finalR_avg_log = finalR_avg_log[~np.isnan(finalR_avg_log)]

plt.figure()
plt.hist(finalL_avg_log, bins=100, histtype='step')
plt.hist(finalR_avg_log, bins=100, histtype='step')
plt.tight_layout()
plt.savefig('fig_hist1_m1.png')

pdb.set_trace()
