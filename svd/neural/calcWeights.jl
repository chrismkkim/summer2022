function calc_rowsum_wp(wpWeightIn)
    #-------------------------------------------------#
    # plastic weights: wpWeightIn = Ncells x indegree #
    #-------------------------------------------------#

    # number of plastic synapses per neuron
    Lexc = Int(size(wpWeightIn)[2]/2)
    Linh = Lexc

    # plastic weights: exc and inh rowsum
    wp_rowsumE = mean(wpWeightIn[:,1:Lexc], dims=2) 
    wp_rowsumI = mean(wpWeightIn[:,Lexc+1:Lexc+Linh], dims=2)

    return wp_rowsumE, wp_rowsumI
end


function calc_std_wp(wpWeightIn)
    #-------------------------------------------------#
    # plastic weights: wpWeightIn = Ncells x indegree #
    #-------------------------------------------------#

    # number of plastic synapses per neuron
    Lexc = Int(size(wpWeightIn)[2]/2)
    Linh = Lexc

    # plastic weights: exc and inh rowsum
    wp_stdE = std(wpWeightIn[:,1:Lexc], dims=2) 
    wp_stdI = std(wpWeightIn[:,Lexc+1:Lexc+Linh], dims=2)

    return wp_stdE, wp_stdI
end


function calc_frac_switch(p, wpWeightOut)

    # fraction of synapses that switched
    we_frac = wpWeightOut[:,1:p.Ne]
    Nesyn = sum(abs.(we_frac) .> 0)
    Nesyn_switch = sum(we_frac .< 0)

    wi_frac = wpWeightOut[:,p.Ne+1:p.Ncells]
    Nisyn = sum(abs.(wi_frac) .> 0)
    Nisyn_switch = sum(wi_frac .> 0)

    frac_switch_exc = Nesyn_switch/Nesyn
    frac_switch_inh = Nisyn_switch/Nisyn
    
    return frac_switch_exc, frac_switch_inh

end