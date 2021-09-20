function fig2()
    #=
    This will create fig2 in the paper.
    =#

    # Loading data ###############################################################
    (datas_bs,datas_ainfs,Ls,
    plot_label,other_labels,D,
    ds,j,threshold,steps) = load_data4(2) # 1,2,3

    data = load("data/bs_june25_quick_all.jld")# without gram-schmidt L = 6-14
    Ls = data["Ls"]
    Bs = data["Bs"]

    data2 = load("data/bs_june25_slow_all.jld")# with gram-schmidt L=6-12
    Bs2 = data2["Bs"]

    # Setting up plots ##########################################################
    fig,ax = plt.subplots(2,length(ds),figsize = [6.8,3.],sharex = true,sharey = "row")
    fig.subplots_adjust(hspace = 0.0,wspace = 0.0)

    # running through data sets #################################################
    for (di,d) in enumerate(ds)
	for l in 1:length(Ls)
	    bs = Bs[di][l][1:400-1]
	    Nb = length(bs)
	    phi_o,phi_e = zero_modes(bs)
	    phi_o = phi_o ./ phi_o[1] # undo normalization
	    # plotting phi_o
	    ax[2,di].plot(1:2:Nb+1,phi_o[1:2:end].^2,
			  lw = .5,
			  color = "C$(l-1)",
			  label = "\$ L = $(Ls[l])\$")
	    # plotting bs
	    ax[1,di].plot(1:length(bs),bs,
			  lw = .5,
			  color = "C$(l-1)",
			  label = "\$ L = $(Ls[l])\$")
	    if l<length(Ls)
		# if L != 14, plot gram-schmidt data set
		bs_slow = Bs2[di][l][1:400-1]
		# plot bs
		ax[1,di].plot(1:length(bs_slow),bs_slow,lw = 2,
			      alpha = .2,
			      color = "C$(l-1)",
			      label = "\$ L = $(Ls[l])\$")
		Nb = length(bs_slow)
		phi_o,phi_e = zero_modes(bs_slow)
		phi_o = phi_o ./ phi_o[1] # undo normalization
		# plot phi_o
		ax[2,di].plot(1:2:Nb+1,phi_o[1:2:end].^2,
			      lw = 2,
			      alpha = .2,
			      color = "C$(l-1)")
	    end
	end
	# rest is plotting ######################################################
	ax[2,di].set_xlabel("\$ n \$",labelpad = 0.0,fontsize = "small")
	ax[1,di].set_xscale("log")
	ax[2,di].set_xscale("log")
	ax[2,di].set_yscale("log")
	ax[1,di].text(.05,.8,
		      "\$ $(plot_label[1]) = $(plot_label[2][di])\$",
		      transform = ax[1,di].transAxes,
		      fontsize = "small")

	ax[1,di].tick_params(axis = "both",which = "both",labelsize = "small",
			     direction = "in")
	ax[2,di].tick_params(axis = "both",which = "both",labelsize = "small",
			     direction = "in")
    end
    ax[1,1].set_ylabel("\$ b_n \$",labelpad = 0.0,fontsize = "small")
    ax[2,1].set_ylabel("\$|\\phi|^2 \$",labelpad = 0.0,fontsize = "small")
    ax[1,1].text(.083,.6,
		 "\\begin{align*} $(other_labels[1]) &= $(other_labels[2])\\\\[-5pt] $(other_labels[3]) &= $(other_labels[4]) \\end{align*}",
		 transform = ax[1,1].transAxes,
		 fontsize = "small")
    for l in 1:length(Ls)
	ax[1,1].plot([],[],lw = 4,color = "C$(l-1)",
		     label = "\$ L = $(Ls[l])\$")
    end
    lines = ax[1,1].get_lines()
    ax[2,end].legend(loc = "lower right",
		     handles = [lines[end-i] for i in 4:-1:0],
		     fontsize = "small",
		     framealpha = 0.0,
		     bbox_to_anchor = (1.07,-.05),
		     handlelength = 1.0,
		     handletextpad = .5,
		     ncol = 1,
		     columnspacing = .5,
		     labelspacing = .1                    
		     )
    fig.savefig("fig2.png",bbox_inches = "tight",dpi = 800)
    fig.show()
end
