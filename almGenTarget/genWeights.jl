function genWeights(p)

nc0Max = Int(p.Ncells*p.pree) # outdegree
nc0 = Int.(nc0Max*ones(p.Ncells))
w0Index = zeros(Int,nc0Max,p.Ncells)
w0Weights = zeros(nc0Max,p.Ncells)
for i = 1:p.Ncells
    postcells = filter(x->x!=i, collect(1:p.Ncells)) # remove autapse
    w0Index[1:nc0Max,i] = sort(shuffle(postcells)[1:nc0Max]) # fixed outdegree nc0Max
    nexc = sum(w0Index[1:nc0Max,i] .<= p.Ne) # number of exc synapses
    if i <= p.Ne
        w0Weights[1:nexc,i] .= p.jee  ## EE weights
        w0Weights[nexc+1:nc0Max,i] .= p.jie  ## IE weights
    else
        w0Weights[1:nexc,i] .= p.jei  ## EI weights
        w0Weights[nexc+1:nc0Max,i] .= p.jii  ## II weights
    end
end


# nc0Max = Int(p.Ncells*p.pree) + floor(Int, 6*sqrt(p.Ncells*p.pree*(1-p.pree))) # outdegree = mean + 6*sd
# nc0 = ones(Int, p.Ncells)
# w0Index = zeros(Int,nc0Max,p.Ncells)
# w0Weights = zeros(nc0Max,p.Ncells)
# for i = 1:p.Ncells
#     # generate postsynaptic cells
#     v = rand(p.Ncells) .< p.pree
#     postcells = findall(v.==1)
#     postcells = filter(x->x!=i, postcells) # remove autapse
#     Npostcells = length(postcells)

#     # define nc0, w0Index
#     nc0[i] = Npostcells    
#     w0Index[1:Npostcells,i] = postcells

#     # define w0Weights
#     nexc = sum(postcells .<= p.Ne) # number of exc synapses
#     if i <= p.Ne
#         w0Weights[1:nexc,i] .= p.jee  ## EE weights
#         w0Weights[nexc+1:Npostcells,i] .= p.jie  ## IE weights
#     else
#         w0Weights[1:nexc,i] .= p.jei  ## EI weights
#         w0Weights[nexc+1:Npostcells,i] .= p.jii  ## II weights
#     end
# end



return w0Index, w0Weights, nc0


end
