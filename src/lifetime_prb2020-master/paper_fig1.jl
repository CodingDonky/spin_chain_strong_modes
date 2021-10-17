function fig1()
    #=
    This will create fig1 in the paper.
    =#

    # Loading data ############################################################
    (datas_bs,datas_ainfs,Ls,
     plot_label,other_labels,D,
     ds,j,threshold,steps) = load_data4(2) # 1,2,3

    # data = load("data/bs_june25_quick_all.jld") # no gram-schmidt L=6-14
    data2 = load("data/bs_june25_slow_all.jld") # gram-schmidt L=6-12
    Bs = data2["Bs"]

    # setting up plotting #####################################################
    fig0,ax0 = plt.subplots(length(ds),1,figsize = [3.4,5])
    fig0.subplots_adjust(hspace = .2)
    times2 = 10 .^ range(-2,5,length = 100)

    # running through data sets ###############################################
    for (di,d) in enumerate(ds)
	for l in 1:length(Ls)
	    # getting ED ainf data
	    ainf = datas_ainfs[d]["ainfs"][l]
	    plt_times = datas_ainfs[d]["times"]
	    # plotting ED ainf data
	    ax0[di].plot(plt_times,ainf,"-o",lw = .5,ms = 1.5,
			 color = "C$(l-1)",label = "\$ L = $(Ls[l])\$")

	    # only plots estimate for last gamma
	    if l == length(Ls)-1
		bs_tmp = Bs[di][l][1:400]
		# analytical gamma, no R
		gammaA,prefactor = analytic_gamma_noR(bs_tmp)
		# numerical gamma
		engI = find_half_width(bs_tmp[1:end-1],bs_tmp[end])
		println("L = $(l): gammaA=$(gammaA), gammaB=$(engI), rel err=$(abs(gammaA-engI)/engI)")
		# plotting estimate
		ax0[di].plot(times2,prefactor .* exp.(-abs(gammaA) .* times2),
			     lw = 2,alpha = .5,color = "black")
	    end
	end
	# plotting settings
	ax0[di].set_xscale("log")
	ax0[di].text(.05,.05,
		     "\$ $(plot_label[1]) = $(plot_label[2][di])\$",
		     transform = ax0[di].transAxes,
		     fontsize = "small")
	ax0[di].set_ylabel("\$ A_{\\infty}\$",labelpad = 0.0,fontsize = "small")
	ax0[di].tick_params(axis = "both",which = "both",labelsize = "small",
			    direction = "in")
    end

    # the rest is plotting ############################################################
    ax0[end].set_xlabel("time",labelpad = 0.0,fontsize = "small")
    ax0[1].text(.07,.24,
		"\\begin{align*} $(other_labels[1]) &= $(other_labels[2])\\\\[-5pt] $(other_labels[3]) &= $(other_labels[4]) \\end{align*}",
		transform = ax0[1].transAxes,
		fontsize = "small")
    ax0[1].set_xlim(10^(-1),10^5)
    ax0[2].set_xlim(10^(-1),10^4)
    ax0[3].set_xlim(10^(-1),10^3)
    ax0[4].set_xlim(10^(-1),10^(2.4))
    ax0[5].set_xlim(10^(-1),10^(2.4))

    for l in 1:length(Ls)
	ax0[end].plot([],[],lw = 4,color = "C$(l-1)",
		      label = "\$L=$(Ls[l])\$")
    end
    lines = ax0[end].get_lines()
    ax0[1].legend(loc = "upper center",
		  handles = [lines[end-i] for i in 4:-1:0],
		  fontsize = "small",
		  framealpha = 1.0,
		  bbox_to_anchor = (.500,1.4),
		  handlelength = .5,
		  handletextpad = .15,
		  ncol = 5,
		  columnspacing = .5,
		  labelspacing = .2                    
		  )
    fig0.savefig("fig1.png",bbox_inches = "tight",dpi = 800)
    fig0.show()
end
