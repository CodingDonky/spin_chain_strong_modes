function fig3()
    #=
    This will create Fig3 in the paper. 
    =#

    # setting up plot #################################################
    fig,ax = plt.subplots(2,2,figsize = [3.4,2.],
			  sharey = "row",
			  sharex = true)
    fig.subplots_adjust(wspace = 0.0,hspace = 0.0)
    ax[2,1].set_xlabel("n",labelpad = 0,fontsize = "small")
    ax[1,1].set_ylabel("\$b_n \$",labelpad = 0,fontsize = "small")
    ax[2,1].set_ylabel("\$|\\phi|^2,|\\eta|^2 \$",
		       labelpad = 0,fontsize = "small")
    ax[2,1].set_yscale("log")

    # setting up toy models and plotting ###############################
    L = 20
    delta = 2
    alpha = .3
    rho = .7
    bs,Nstar = gen_bs(L,alpha,delta,rho) # first toy model. gen_bs() is located in paper_util.jl
    phi_o,phi_e = zero_modes(bs) # get zero modes for toy model. zero_modes() is located in paper_util.jl

    ax[1,1].plot(1:length(bs),bs,"-o",ms = 1.5,lw = .5)
    ax[2,1].plot(1:2:length(phi_o), phi_o[1:2:end].^2,"-o",ms = 1.5,lw = .5,
		 label = "\$ \\phi\$")
    ax[2,1].plot(2:2:length(phi_e), phi_e[2:2:end].^2,"-s",ms = 1.5,lw = .5,
		 label = "\$ \\eta\$")

    # second toy model
    A = 1.0
    beta = 4
    bs = gen_bs3(10,2*L+1,alpha,delta,A,beta) #gen_bs3() located in paper_util.jl
    phi_o,phi_e = zero_modes(bs) #zero_modes() is located in paper_util.jl
    ax[1,2].plot(1:length(bs),bs,"-o",ms = 1.5,lw = .5)
    ax[2,2].plot(1:2:length(phi_o), phi_o[1:2:end].^2,"-o",ms = 1.5,lw = .5,
		 label = "\$ \\phi\$")
    ax[2,2].plot(2:2:length(phi_e), phi_e[2:2:end].^2,"-s",ms = 1.5,lw = .5,
		 label = "\$ \\eta\$")
    ax[2,2].set_xlabel("n",labelpad = 0,fontsize = "small")

    ax[2,2].legend(loc = "lower right",
		   ncol = 2,
		   fontsize = "small",
		   handletextpad = .5,
		   handlelength = 1.3,
		   columnspacing = .5
		   )

    # best fit line to check decay is suff
    # x1 = length(phi_o)-1
    # x2 = x1 - 10
    # y1 = phi_o[end-1].^2
    # y2 = phi_o[end-11].^2
    # a1 = log(y2/y1)/log(x1/x2)
    # b1 = log(y1) + a1*log(x1)
    # println(a1)
    # ax[2,2].plot(1:2:length(phi_o),(1:2:length(phi_o)).^(-a1).*exp(b1),
    #            "--")

    ax[1,1].tick_params(axis = "both",which = "both",labelsize = "small",
			direction = "in")
    ax[2,1].tick_params(axis = "both",which = "both",labelsize = "small",
			direction = "in")
    ax[1,2].tick_params(axis = "both",which = "both",labelsize = "small",
			direction = "in")
    ax[2,2].tick_params(axis = "both",which = "both",labelsize = "small",
			direction = "in")

    fig.savefig("fig3.png",bbox_inches = "tight",dpi = 800)    
    fig.show()
end
