####################################################################
# test #############################################################
####################################################################
function simple_bn()
    #=
    This will generate b_n for small system sizes
    =#
    # setting parameters, jzs is a list of J_z
    jzs = [.2,.3] # [.2,.3,.4,.5,.6]
    g = .3; gamma = .9; j = 1.0
    threshold = 1e-13
    # number of b_n to calculate
    steps = 600
    Ls = [6] # list of system size lengths
    Bs = zeros(steps,length(Ls),length(jzs))
    for (jzi,jz) in enumerate(jzs)
	for (Li,L) in enumerate(Ls)
	    println("jzi = $(jzi), Li = $(Li)")
            bsx = ED.get_bs(L,j,jz,g,gamma,threshold,steps) # calculates b_n
	    Bs[:,Li,jzi] = bsx
	end
    end
    # plots b_n
    fig,ax = plt.subplots(length(jzs),1)
    for (jzi,jz) in enumerate(jzs)
	for (Li,L) in enumerate(Ls)
	    ax[jzi].plot(Bs[:,Li,jzi])
	end
    end
    fig.savefig("simple_bn.png")
end
