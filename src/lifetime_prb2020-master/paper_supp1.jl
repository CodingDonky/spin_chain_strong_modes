function supp1()
    #=
    This will create supp1 in the paper.
    =#
    # Loading data ###########################################################
    (datas_bs,datas_ainfs,Ls,
    plot_label,other_labels,D,
    ds,j,threshold,steps) = load_data4(2) # 1,2,3
    data = load("data/bs_june25_quick_all.jld")
    Ls = data["Ls"]
    Bs = data["Bs"]
    data2 = load("data/bs_june25_slow_all.jld")
    Bs2 = data2["Bs"]

    # setting up figure ####################################################
    fig3,ax3 = plt.subplots(length(ds),figsize = [3.4,3.4],sharex = true)
    fig3.subplots_adjust(hspace = 0.0)
    # running through data sets ###########################################
    for (di,d) in enumerate(ds)
	for l in 1:length(Ls)-1
	    bs = Bs[di][l][1:400-1]
	    bs_slow = Bs2[di][l][1:400-1]
	    Nb = length(bs)
	    # plotting error between bs with and without gram-schmidt
	    ax3[di].plot(1:400-1,
			 abs.(bs[1:400-1] - bs_slow[1:400-1]),
			 "-o",lw = .5,ms = 1.5,color  = "C$(l-1)")
	end
	# rest is plotting #####################################################
	ax3[di].text(.025,.8,
		     "\$ $(plot_label[1]) = $(plot_label[2][di])\$",
		     transform = ax3[di].transAxes,
		     fontsize = "small")
	ax3[di].set_ylabel("\$\\text{err}(b_n)\$",labelpad = 0.0,
			   fontsize = "small")
	ax3[di].set_yscale("log")
	ax3[di].tick_params(axis = "both",which = "both",labelsize = "small",
			    direction = "in")
    end
    ax3[1].text(.04,.35,
		"\\begin{align*} $(other_labels[1]) &= $(other_labels[2])\\\\[-5pt] $(other_labels[3]) &= $(other_labels[4]) \\end{align*}",
		transform = ax3[1].transAxes,
		fontsize = "small")
    ax3[end].set_xlabel("n",labelpad = 0.0)
    for l in 1:length(Ls)
	ax3[1].plot([],[],lw = 4,color = "C$(l-1)",
		      label = "\$ L = $(Ls[l])\$")
    end
    lines = ax3[1].get_lines()
    ax3[1].legend(loc = "lower right",
		  handles = [lines[end-i] for i in 4:-1:1],
		  fontsize = "small",
		  framealpha = 0.0,
		  # bbox_to_anchor = (1.0,-.05),
		  handlelength = 1.0,
		  handletextpad = .5,
		  ncol = 2,
		  columnspacing = .5,
		  labelspacing = .1                    
		  )
    fig3.savefig("supp1.png",bbox_inches = "tight",dpi = 800)
    fig3.show()
end
