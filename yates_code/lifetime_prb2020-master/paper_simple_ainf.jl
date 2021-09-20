function simple_ainf()
    # setting parameters, jzs is a list of J_z
    jzs = [.2,.3] # [.2,.3,.4,.5,.6]
    g = .3; gamma = .9; j = 1.0
    # setting system sizes (list)
    Ls = [4,6] #[6,8]
    # times we want to evaluate Ainf
    times = 10 .^ range(-1,4,length = 1000)
    Ainfs = zeros(length(Ls),length(jzs),length(times))
    for (jzi,jz) in enumerate(jzs)
	for (Li,L) in enumerate(Ls)
	    println("jzi = $(jzi), Li = $(Li)")
	    Ainfs[Li,jzi,:] = ED.get_Ainf(L,g,j,jz,gamma,times)
	end
    end
    fig,ax = plt.subplots(length(jzs),1)
    for (jzi,jz) in enumerate(jzs)
	for (Li,L) in enumerate(Ls)
	    ax[jzi].plot(times,Ainfs[Li,jzi,:])
	end
        ax[jzi].set_xscale("log")
    end 
    plt.show()
end
