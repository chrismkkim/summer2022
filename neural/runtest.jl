using Distributions

function runtest(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim)

# copy param
nloop = copy(p.nloop) # train param
penlambda = copy(p.penlambda)
penmu = copy(p.penmu)
frac = copy(p.frac)
learn_every = copy(p.learn_every)
stim_on = copy(p.stim_on)
stim_off = copy(p.stim_off)
train_time = copy(p.train_time)
dt = copy(p.dt) # time param
Nsteps = copy(p.Nsteps) 
Ncells = copy(p.Ncells) # network param
Ne = copy(p.Ne)
Ni = copy(p.Ni)
taue = copy(p.taue) # neuron param
taui = copy(p.taui)
sqrtK = copy(p.sqrtK)
threshe = copy(p.threshe)
threshi = copy(p.threshi)
refrac = copy(p.refrac)
vre = copy(p.vre)
muemin = copy(p.muemin) # external input
muemax = copy(p.muemax)
muimin = copy(p.muimin)
muimax = copy(p.muimax)
tauedecay = copy(p.tauedecay) # synaptic time
tauidecay = copy(p.tauidecay)
taudecay_plastic = copy(p.taudecay_plastic)
maxrate = copy(p.maxrate)

# set up variables
mu = zeros(Ncells)
mu[1:Ne] = (muemax-muemin)*rand(Ne) .+ muemin
mu[(Ne+1):Ncells] = (muimax-muimin)*rand(Ni) .+ muimin

thresh = zeros(Ncells)
thresh[1:Ne] .= threshe
thresh[(1+Ne):Ncells] .= threshi

tau = zeros(Ncells)
tau[1:Ne] .= taue
tau[(1+Ne):Ncells] .= taui

maxTimes = round(Int,maxrate*train_time/1000)
times = zeros(Ncells,maxTimes)
ns = zeros(Int,Ncells)

forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
forwardInputsI = zeros(Ncells)
forwardInputsPE = zeros(Ncells)
forwardInputsPI = zeros(Ncells)
forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
forwardInputsIPrev = zeros(Ncells)
forwardInputsPEPrev = zeros(Ncells)
forwardInputsPIPrev = zeros(Ncells)
forwardSpike = zeros(Ncells)
forwardSpikePrev = zeros(Ncells)

xedecay = zeros(Ncells)
xidecay = zeros(Ncells)
xpdecay = zeros(Ncells)
xpedecay = zeros(Ncells)
xpidecay = zeros(Ncells)
synInputBalanced = zeros(Ncells)
synInputPlastic = zeros(Ncells)

v = rand(Ncells) #membrane voltage 

lastSpike = -100.0*ones(Ncells) #time of last spike
  
t = 0.0
bias = zeros(Ncells)

learn_nsteps = Int(p.train_duration / p.learn_every)
learn_seq = 1
example_neurons = 25
wid = 50
widInc = Int(2*wid/p.learn_every - 1)

xsum = zeros(learn_nsteps,Ncells)
xsumcnt = zeros(learn_nsteps)
tseq = 1


for ti=1:Nsteps
    if mod(ti,Nsteps/100) == 1  #print percent complete
        print("\r",round(Int,100*ti/Nsteps))
    end

    t = dt*ti;
    forwardInputsE .= 0.0;
    forwardInputsI .= 0.0;
    # forwardInputsP .= 0.0;
    forwardInputsPE .= 0.0;
    forwardInputsPI .= 0.0;
    forwardSpike .= 0.0;
    
    for ci = 1:Ncells
        xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]/tauedecay
        xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]/tauidecay
        xpedecay[ci] += -dt*xpedecay[ci]/taudecay_plastic + forwardInputsPEPrev[ci]/taudecay_plastic
        xpidecay[ci] += -dt*xpidecay[ci]/taudecay_plastic + forwardInputsPIPrev[ci]/taudecay_plastic
        synInputBalanced[ci] = xedecay[ci] + xidecay[ci]
        synInputPlastic[ci] = xpedecay[ci] + xpidecay[ci]
        synInput = synInputBalanced[ci] + synInputPlastic[ci]

        # # save rolling average for analysis
        # if rollingavg 
        #     if t > Int(p.stim_off) && t <= Int(p.train_time) && mod(t,1.0) == 0
        #         xsum[:,ci], xsumcnt = funRollingAvg(p,t,wid,widInc,learn_nsteps,xsum[:,ci],xsumcnt,synInput+mu[ci],ci)
        #     end
        # end

        # if rollingavg 
        #     if t > Int(stim_off) && t <= Int(train_time) && mod(t,learn_every) == 0
        #         xsum[tseq,ci] = synInput + mu[ci]
        #         if ci == Ncells
        #             tseq += 1
        #         end
        #     end
        # end
        
        # external input
        if t > stim_on && t < stim_off
            bias[ci] = mu[ci] + stim[ti-Int(stim_on/dt),ci]
        else
            bias[ci] = mu[ci]
        end

        #not in refractory period
        if t > (lastSpike[ci] + refrac)  
            v[ci] += dt*((1/tau[ci])*(bias[ci]-v[ci] + synInput))
            if v[ci] > thresh[ci]  #spike occurred
                v[ci] = vre
                forwardSpike[ci] = 1.
                lastSpike[ci] = t
                ns[ci] = ns[ci]+1
                if ns[ci] <= maxTimes
                    times[ci,ns[ci]] = t
                end
                for j = 1:nc0[ci]
                    wgt = w0Weights[j,ci]
                    cell = w0Index[j,ci]
                    if wgt > 0  #E synapse
                        forwardInputsE[cell] += wgt
                    elseif wgt < 0  #I synapse
                        forwardInputsI[cell] += wgt
                    end
                end #end loop over synaptic projections
                for j = 1:ncpOut[ci]
                    wgtp = wpWeightOut[j,ci]
                    cell = Int(wpIndexOut[j,ci])
                    if ci <= Ne # E synapse
                        forwardInputsPE[cell] += wgtp
                    elseif ci > Ne # I synapse
                        forwardInputsPI[cell] += wgtp
                    end
                end
            end #end if(spike occurred)
        end #end not in refractory period
    end #end loop over neurons

    forwardInputsEPrev = copy(forwardInputsE)
    forwardInputsIPrev = copy(forwardInputsI)
    forwardInputsPEPrev = copy(forwardInputsPE)
    forwardInputsPIPrev = copy(forwardInputsPI)
    forwardSpikePrev = copy(forwardSpike) # if training, compute spike trains

end #end loop over time
print("\r")

# for k = 1:learn_nsteps
#     xsum[k,:] = xsum[k,:]/xsumcnt[k]
# end

return xsum, times, ns 


end
