using Distributions

function runtestperturb(p,w0Index,w0Weights,nc0,wpIndexOut,wpWeightOut,ncpOut,stim_training,stim_perturb,ffwdRate,wpWeightFfwd, isperturbed, pert_on, pert_off, stimtype)

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
times_ffwd = zeros(p.Lffwd, maxTimes)
ns_ffwd = zeros(Int, p.Lffwd)

forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
forwardInputsI = zeros(Ncells)
forwardInputsP = zeros(Ncells)
forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
forwardInputsIPrev = zeros(Ncells)
forwardInputsPPrev = zeros(Ncells)
forwardSpike = zeros(Ncells)
forwardSpikePrev = zeros(Ncells)

xedecay = zeros(Ncells)
xidecay = zeros(Ncells)
xpdecay = zeros(Ncells)
synInputBalanced = zeros(Ncells)

v = zeros(Ncells) #membrane voltage 

lastSpike = -100.0*ones(Ncells) #time of last spike
  
t = 0.0
bias = zeros(Ncells)

learn_nsteps = Int((p.train_time - p.stim_off)/p.learn_every)
learn_seq = 1
example_neurons = 25
wid = 50
widInc = Int(2*wid/p.learn_every - 1)

xtotal = zeros(learn_nsteps,Ncells)

rndphase = rand(Ncells)
rndperturb = 6.0*(rand(Ncells) .- 0.5) # perturb inhibitory neurons
# perturb = 6.0*(rand(Ncells) .- 0.5) # perturb exc and inh neurons


for ti=1:Nsteps
    if mod(ti,Nsteps/100) == 1  #print percent complete
        print("\r",round(Int,100*ti/Nsteps))
    end

    t = dt*ti;
    forwardInputsE .= 0.0;
    forwardInputsI .= 0.0;
    forwardInputsP .= 0.0;
    forwardSpike .= 0.0;
    rndFfwd = rand(p.Lffwd)

    for ci = 1:Ncells
        xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]/tauedecay
        xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]/tauidecay
        xpdecay[ci] += -dt*xpdecay[ci]/taudecay_plastic + forwardInputsPPrev[ci]/taudecay_plastic
        synInputBalanced[ci] = xedecay[ci] + xidecay[ci]
        synInput = synInputBalanced[ci] + xpdecay[ci]
        
        # external input
        if t > Int(stim_on) && t < Int(stim_off)
            bias[ci] = mu[ci] + stim_training[ti-Int(stim_on/dt),ci]
        else
            bias[ci] = mu[ci]
        end

        # perturbation 
        if isperturbed

            # perturb using stimulus that triggers the opposite learned activity
            if stimtype == "cue"
                if t > pert_on && t < pert_on + Int(stim_off - stim_on)
                    bias[ci] = mu[ci] + stim_perturb[ti-Int(pert_on/dt),ci]
                end
            end
            
            # perturb inhibitory neurons with random stimulus
            if stimtype == "inh"
                if ci > Ne && t > pert_on && t < pert_off
                    bias[ci] += rndperturb[ci] * sin(2*pi*(t/100 - rndphase[ci]))
                end
            end

            # # perturb exc / inh neurons with random stimulus
            # if t > pert_on && t < pert_off
            #     bias[ci] += perturb[ci] * sin(2*pi*(t/100 - rndphase[ci]))
            # end

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
                    cell = Int(wpIndexOut[j,ci])
                    forwardInputsP[cell] += wpWeightOut[j,ci]
                end
            end #end if(spike occurred)
        end #end not in refractory period
    end #end loop over neurons

    # External input to trained excitatory neurons
    if ti > Int(stim_off/dt)
        tidx = ti - Int(stim_off/dt)
        for ci = 1:p.Lffwd
            # # if training, filter the spikes
            # s[ci] += -dt*s[ci]/taudecay_plastic + ffwdSpikePrev[ci]/taudecay_plastic

            # if Poisson neuron spiked
            if rndFfwd[ci] < ffwdRate[tidx,ci]/(1000/p.dt)
                # ffwdSpike[ci] = 1.
                ns_ffwd[ci] = ns_ffwd[ci]+1
                if ns_ffwd[ci] <= maxTimes
                    times_ffwd[ci,ns_ffwd[ci]] = t
                end
                for j = 1:Ne
                    forwardInputsP[j] += wpWeightFfwd[j,ci]
                end #end loop over synaptic projections
            end #end if spiked
        end #end loop over ffwd neurons
    end #end ffwd input

    forwardInputsEPrev = copy(forwardInputsE)
    forwardInputsIPrev = copy(forwardInputsI)
    forwardInputsPPrev = copy(forwardInputsP)
    forwardSpikePrev = copy(forwardSpike) # if training, compute spike trains

end #end loop over time
print("\r")


return xtotal, times, ns


end
