function supp2_run()
    #= 
    This will make the krylov time-evolve data for supp2, 
    better to run it and save data.
    =#
    times = 10 .^ range(-1,2.5,length = 40)
    Ls = [20,40]
    L2 = 10000

    A = .5
    beta = 10
    alpha = .3
    delta = 1
    x0 = 10

    Bs = []
    Ainfs = []
    for (Li,L) in enumerate(Ls)
	println("Li = $(Li) ###############################")
	bs = gen_bs3(x0,L,alpha,delta,A,beta) #gen_bs3() is located in paper_util.jl
	bs = extend_bs(L2,bs)
	# krylov-time evolve, takes some time
	signal = time_evolve(bs,times)  #time_evolve() is located in paper_util.jl
	push!(Bs,bs)
	push!(Ainfs,signal)
    end
    # saving data
    save("data/supp2_data.jld",
	 "Bs",Bs,
	 "times",times,
	 "Ainfs",Ainfs,
	 "Ls",Ls,
	 "L2",L2,
	 "A",A,
	 "beta",beta,
	 "alpha",alpha,
	 "delta",delta,
	 "x0",x0
	 )
end

function supp2()
    #=
    This will create supp2 in the paper.
    =#
    # loading data from supp2_run() ###############################
    data = load("data/supp2_data.jld")
    Bs = data["Bs"]
    times = data["times"]
    Ainfs = data["Ainfs"]
    Ls = data["Ls"]
    L2 = data["L2"]
    alpha = data["alpha"]
    delta = data["delta"]
    x0 = data["x0"]
    beta = data["beta"]
    A = data["A"]

    # setting up plot ##############################################
    fig,ax = plt.subplots(3,1,figsize = [3.4,4.4])
    fig.subplots_adjust(hspace = .3)
    # running through data sets ###################################
    ms_tmp = [4,2]
    for (Li, L) in enumerate(Ls)
	bs = Bs[Li]
	ainfs = Ainfs[Li]
	phi_o,phi_e = zero_modes(bs[1:17])

	engs = range(-.1,.1,length = 200)
	# generating numerical greens functions
	data = num_greens_func(bs[1:L],bs[L+1],engs)
	# bs = gen_bs3(x0,L,alpha,delta,A,beta)
	engI = find_half_width(bs[1:L-1],bs[L])
	# plotting bs
	ax[1].plot((1:50), bs[1:50],
		   "-o",ms = ms_tmp[Li],lw = .5,color = "C$(Li-1)",
		   label = "L = $(L)")
	# plotting lorentizan for greens functions
	ax[2].plot(engs,abs.(imag.(data)),
		   lw = ms_tmp[Li]/1.25,color = "C$(Li-1)")
	ax[2].axvline(x = -engI,lw = ms_tmp[Li]/1.25,color = "C$(Li-1)")
	ax[2].axvline(x =  engI,lw = ms_tmp[Li]/1.25,color = "C$(Li-1)")
	ax[2].axhline(y = phi_o[1]^2/engI,lw = ms_tmp[Li]/1.25,color = "C$(Li-1)")
	# ax[1].axhline(y = t0)

	# plotting krylov-ainf from supp2_run
	ax[3].plot(times,ainfs,
		   "-o",
		   lw = .5,ms = ms_tmp[Li]/2.0,color = "C$(Li-1)",
		   label = "\$ L = $(L)\$")
	if Li == 2
	    ax[3].plot(times,phi_o[1]^2 .* exp.(-abs.(engI).*times),
		       lw = 5,alpha = .2,color = "black")
	end

    end
    # rest is plotting ########################################################
    ax[3].set_xscale("log")
    ax[1].set_xlabel("n",labelpad = 0.0,fontsize = "small")
    ax[1].set_ylabel("\$ b_n \$",labelpad = 0.0,fontsize = "small")
    ax[2].set_xlabel("\$ E \$",labelpad = 0.0,fontsize = "small")
    ax[3].set_xlabel("time",labelpad = 0.0,fontsize = "small")
    ax[2].set_ylabel("\$|\\Im(G)| \$",labelpad = 0.0,fontsize = "small")
    ax[3].set_ylabel("\$ A_{\\infty} \$",labelpad = 0.0,fontsize = "small")
    ax[1].tick_params(axis = "both",which = "both",labelsize = "small",
		      direction = "in")
    ax[2].tick_params(axis = "both",which = "both",labelsize = "small",
		      direction = "in")
    ax[3].tick_params(axis = "both",which = "both",labelsize = "small",
		      direction = "in")
    fig.savefig("supp2.png",bbox_inches = "tight",dpi = 800)
    fig.show()
end
