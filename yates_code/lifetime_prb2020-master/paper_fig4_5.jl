function fig4_5()
    #=
    This will create figures 4 and 5 in paper.
    =#
    # IMPORTING DATA ###############################################################
    (datas_bs,datas_ainfs,Ls,
    plot_label,other_labels,D,
    ds,j,threshold,steps) = load_data4(2) # 1,2,3
    data = load("data/bs_june25_quick_all.jld")
    Ls = data["Ls"]
    Bs = data["Bs"]
    jzs = data["jzs"]
    data2 = load("data/bs_june25_slow_all.jld")
    Bs2 = data2["Bs"]
    # SETTING UP PLOTS ###############################################################
    fig2,ax2 = plt.subplots(2,length(ds),
			    sharey = "row",
			    sharex = "row",
			    figsize = [6.8,2.4],
			    gridspec_kw = Dict("height_ratios" => (1,2))
			    )
    fig2.subplots_adjust(wspace = 0.0,hspace = 0.3)
    linewidths = [5,4,3,2]
    # SETTING UP PARAMS ###############################################################
    engsLims = [.005,.025,.1,.2,.5]
    nL = 40 # number of sites!! NOT number of b_ns (nL-1)
    nRs = [40,100,200,400] # same as above

    data3_ed = zeros(length(ds))
    data3_approx = zeros(length(nRs),length(ds))
    data3_toy = zeros(length(ds))
    # LOOPING THROUGH Jzs ###############################################################
    for (di,d) in enumerate(ds)
	l = 4 # working only with L = 12
	# last bn considered coupling to metal
	# and value of metallic hopping
	bs = Bs2[di][l][1:400]
	Nb = length(bs)
	phi_o,phi_e = zero_modes(bs[1:29])

	ainf = datas_ainfs[d]["ainfs"][l]
	plt_times = datas_ainfs[d]["times"]
	# determining gamma from ed ainf
	for dati in length(ainf):-1:1
	    if ainf[dati] > phi_o[1]^2*exp(-1)
		data3_ed[di] = 1/plt_times[dati]
		break
	    end
	end

	# for fig4 #####################################################
	# getting h_n and \tilde{h}_n
	bs_sign = (-1).^(1:Nb-1)
	bs_avg = (bs[2:end] .+ bs[1:end-1])./2
	bs_dif = (bs[1:end-1].- bs[2:end])./2 .* bs_sign

	# getting alpha,delta with best fit
	x_tmp = 1:20
	y_tmp = bs_avg[1:20]
	alpha,delta = simple_best_fit(x_tmp,y_tmp)
	alpha = round(alpha,digits = 2)
	delta = round(delta,digits = 2)
	#plotting
	ax2[1,di].plot(1:70,delta .+ alpha.*(1:70),"--",color = "black",
		       label = "\$ $(alpha)n + $(delta) \$")
	ax2[1,di].text(.025,.6,
		       "\\begin{align*} \\alpha &= $(alpha)\\\\[-5pt]\\delta &= $(delta)\\end{align*}",
		       transform = ax2[1,di].transAxes,
		       fontsize = "small")
	# getting moving average
	navg = 7
	bs_dif_avg = copy(bs_dif[1:end-(navg-1)]) ./navg
	for i in 1:navg-1
	    bs_dif_avg += bs_dif[1+i:end-(navg-1)+i] ./navg
	end
	# plotting
	pltEnd = 300
	ax2[1,di].plot(1:pltEnd,bs_avg[1:pltEnd],"-o",lw = .5,ms = 2.5)
	ax2[2,di].plot(1:pltEnd,bs_dif[1:pltEnd],"-o",lw = .5,ms = 1.5,alpha = .7,
		       label = "\$\\tilde{h}(n)\$")
	ax2[2,di].plot(1:pltEnd,bs_dif_avg[1:pltEnd],
		       label = "\$\\langle \\tilde{h}(n)\\rangle_{$(navg)}\$")
	# finding x_0 for square wave
	A = round(bs_dif_avg[1],digits = 2)
	x0 = 1
	for k in length(bs_dif_avg):-1:1
	    if (bs_dif_avg[k]) > A/2
		x0 = k
		break
	    end
	end
	ax2[2,di].text(.4,.8,
		       "\\begin{align*} M_0 &= $(A*2)\\\\[-5pt] n_0 &= $(x0)\\end{align*}",
		       transform = ax2[2,di].transAxes,
		       fontsize = "small",
		       bbox = Dict("boxstyle"=>"square,pad=0.3",
				   "facecolor"=>"white",
				   "alpha"=>.75))
	square = zeros(pltEnd)
	square[1:x0] .= A
	ax2[2,di].plot(1:pltEnd,square,label = "\$\\frac{M_0}{2} \\theta(n_0 -n)\$")
	# moving on to fig 5 ########################################################
	# toy model
	bs_tmp = gen_bs3(x0,50,alpha,delta,A,100)
	# analytic gamma
        # "analytic_gamma" is in paper_util.jl
        # there I fixed an error with "R" in the original plots in the first
        # version of the paper.
	gammaA_tmp,prefactor = analytic_gamma(bs_tmp)
	data3_toy[di] = gammaA_tmp
	for nRi in 1:length(nRs)
	    nR = nRs[nRi]
	    gammaA,prefactor = analytic_gamma(bs[1:nR])
	    data3_approx[nRi,di] = gammaA
	end
	# rest is plotting ############################################################
	ax2[1,di].set_xlim(-1,40); ax2[1,di].set_ylim(-1,15)
	ax2[2,di].set_xscale("log")
	ax2[2,di].set_xlabel("\$ n \$",labelpad = 0,fontsize = "small")
	ax2[1,di].set_xlabel("\$ n \$",labelpad = 0,fontsize = "small")
	ax2[1,di].tick_params(axis = "both",which = "both",labelsize = "small",
			      direction = "in")
	ax2[2,di].tick_params(axis = "both",which = "both",labelsize = "small",
			      direction = "in")
	ax2[1,di].text(.55,.05,
		       "\$ $(plot_label[1]) = $(plot_label[2][di])\$",
		       transform = ax2[1,di].transAxes,
		       fontsize = "small")
    end

    # plotting fig 5 ###############################################################33
    fig3,ax3 = plt.subplots(2,1,figsize = [3.4,5])
    fig3.subplots_adjust(hspace = 0.3)
    ax3[1].plot(1. ./jzs, data3_ed,"-o",ms = 7,lw = .5,color = "black",
	     label = "\$ \\Gamma^{\\text{ED}}\$")
    ax3[1].tick_params(axis = "both",which = "both",labelsize = "small",
		    direction = "in")
    ax3[1].plot(1. ./jzs, data3_toy,"-o",ms = 8,fillstyle = "none",lw = .5,
	     color = "black",
	     label = "\$ \\Gamma_A^{\\text{toy}}\$")
    for k in 1:length(nRs)
	ax3[1].plot(1. ./jzs, data3_approx[k,:],"-s",ms = 5,lw = .5,
		 color = "C$(k-1)",
		 label = "\$ \\Gamma_A^{($(nRs[k]))}\$")
    end

    x2 = [40,100,200,400]
    for k in 1:length(jzs)
        ax3[2].plot(x2,data3_approx[:,k],"-o",color = "C$(k-1)",
                    label = "\$\\Gamma_A^{(N)}\$")

        ax3[2].axhline(y = data3_ed[k],linestyle = "--",color = "C$(k-1)",
                       label = "\$\\Gamma^{ED}\$")
    end

    for k in 1:2
        ax3[k].set_yscale("log")
        ax3[k].set_ylabel("\$\\Gamma \$",
		          fontsize = "small",
		          labelpad = 0.0)
    end
    ax3[2].set_ylim(3e-4,2e0)
    ax3[2].set_xscale("log")
    ax3[1].set_xlabel("\$ 1/J_z \$",
		      fontsize = "small",
		      labelpad = 0.0)
    ax3[2].set_xlabel("\$ N \$",
		      fontsize = "small",
		      labelpad = 0.0)
    ax3[1].legend(loc = "upper right",
	          fontsize = "small",
	          handlelength = 1.5,
	          handletextpad = .2,
	          framealpha = 0.0,
	          ncol = 2,
	          columnspacing = .5)
    lines = ax3[2].get_lines()
    legs = []
    for k in 1:length(jzs)
        push!(legs,ax3[2].legend(
            handles = lines[2*k-1:2*k],
            title =  "\$ J_z = $(jzs[k])\$",
            loc = "upper right",
            bbox_to_anchor = (.03 + k/length(jzs),1.05),
	    fontsize = "small",
            title_fontsize = "small",
	    handlelength = 1.5,
	    handletextpad = .2,
	    framealpha = 0.0,
	    ncol = 1,
	    columnspacing = .5))
    end
    for l in legs[1:end-1]
        ax3[2].add_artist(l)
    end
    ax2[2,end].legend(loc = "lower right",
		      fontsize = "small",
		      bbox_to_anchor = (1.05,-.05),
		      handlelength = 1.25,
		      handletextpad = .2,
		      framealpha = .75,
		      ncol = 3,
		      columnspacing = .5)
    ax2[1,1].text(.58,.23,
		  "\\begin{align*} $(other_labels[1]) &= $(other_labels[2])\\\\[-5pt] $(other_labels[3]) &= $(other_labels[4]) \\end{align*}",
		  transform = ax2[1,1].transAxes,
		  fontsize = "small")
    ax2[1,1].set_ylabel("\$h(n)\$",labelpad = 0,fontsize = "small")
    ax2[2,1].set_ylabel("\$\\tilde{h}(n)\$",labelpad = 0,fontsize = "small")
    fig2.savefig("fig4.png",bbox_inches = "tight",dpi = 800)
    fig3.savefig("fig5.png",bbox_inches = "tight",dpi = 800)
end
