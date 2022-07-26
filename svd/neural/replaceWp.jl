function replaceWp(p,ncpIn,wpIndexIn,wpIndexConvert,wpWeightIn,wpWeightOut, dirData, iloop)

    # start_time = time()

    idz = falses(p.Lexc)
    for ci = 1:p.Ncells
        # println(ci)
        # find incoming synapses that violate Dale's principle
        exc_syn = @view wpWeightIn[ci,1:p.Lexc]
        inh_syn = @view wpWeightIn[ci,p.Lexc+1:end]
        exc_idx = findall(x->x<0, exc_syn)
        inh_idx = findall(x->x>0, inh_syn)
        Nesyn = length(exc_idx)
        Nisyn = length(inh_idx)

        # replace exc syn that violate Dale's principle
        if Nesyn > 0
            # delete syn that violate Dale's principle
            dale_idx = deleteat!(collect(1:p.Lexc), exc_idx)
            # select random syn respecting Dale's principle
            rand_idx = sort(shuffle(dale_idx)[1:Nesyn])
            # replaced by the mean of random and nondale synapses
            replace_syn = (exc_syn[exc_idx] + exc_syn[rand_idx])/2
            wpWeightIn[ci,exc_idx] = replace_syn
            wpWeightIn[ci,rand_idx] = replace_syn
        end

        # replace inh syn that violate Dale's principle
        if Nisyn > 0
            # delete syn that violate Dale's principle
            dale_idx = deleteat!(collect(1:p.Linh), inh_idx)
            # select random syn respecting Dale's principle
            rand_idx = sort(shuffle(dale_idx)[1:Nisyn])
            # replace by random synapse
            replace_syn = (inh_syn[inh_idx] + inh_syn[rand_idx])/2
            wpWeightIn[ci,inh_idx .+ p.Lexc] = replace_syn
            wpWeightIn[ci,rand_idx .+ p.Lexc] = replace_syn
            ## wpWeightIn[ci,inh_idx .+ p.Lexc] = (inh_syn[dale_syn])[rand(1:length(dale_syn), Nisyn)]
        end
    end

    wpWeightOut = convertWgtIn2Out(p,ncpIn,wpIndexIn,wpIndexConvert,wpWeightIn,wpWeightOut)

    # elapsed_time = time()-start_time
    # writedlm(dirData * "_replaceWp_loop$(iloop).txt", elapsed_time)    

    return wpWeightIn, wpWeightOut
end