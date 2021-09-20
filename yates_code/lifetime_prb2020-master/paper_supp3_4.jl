function supp3()
    #=
    Comparing analytical estimate of gamma to 
    numerical solution of gamma, for toy models. 
    Similar to supp4(), but varying different parameters.
    =#
    L = 300
    x0s = 11:8:251
    As = [1.,1.,.5,1.]
    betas = [10,10,10,10]
    alphas = [.2,1.2,.2,.2]
    deltas = [1.5,1.5,1.5,3]
    gammas = []
    gammasA = []
    # running through parameters ###########################################33
    for i in 1:4
	A = As[i]
	M = 2*A
	beta = betas[i]
	alpha = alphas[i]
	delta = deltas[i]
	gamma_tmp =zeros(length(x0s))
	gammaA_tmp = zeros(length(x0s))
	for (x0i,x0) in enumerate(x0s)
	    bs = gen_bs3(x0,L,alpha,delta,A,beta)
	    engI = find_half_width(bs[1:end-1],bs[end])
	    gamma_tmp[x0i] = engI
	    g_tmp,prefactor = analytic_gamma(bs)
	    println("x0 = $(x0): approx=$(g_tmp), exact=$(engI), rel err=$(abs(g_tmp-engI)/engI)")
	    gammaA_tmp[x0i] = g_tmp
	end
	push!(gammas,gamma_tmp)
	push!(gammasA,gammaA_tmp)
    end
    # plotting ####################################################################
    fig,ax = plt.subplots(1,1,figsize = [6.8,1.8])
    for i in 1:4
	beta = betas[i]
	alpha = alphas[i]
	delta = deltas[i]
	A = As[i]; M0 = 2*A
	ax.plot(x0s,gammas[i],"-x",
		label = "\$ \\alpha = $(alpha), \\delta = $(delta), M_0 = $(M0)\$",
		lw = .5,
		ms= 3.0)
	ax.plot(x0s,gammasA[i],"-o",
		label = "\$ \\alpha = $(alpha), \\delta = $(delta), M_0 = $(M0)\$",
		lw = .5,
		ms= 2.0)
    end
    ax.tick_params(axis = "both",which = "both",labelsize = "small",
		   direction = "in")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("\$ n_0 \$")
    ax.set_ylabel("\$ \\Gamma^{($(L))}, \\Gamma_A \$")
    lines = ax.get_lines()
    leg0 = ax.legend(loc = "lower center",
		     handles = [lines[2*j-1] for j in 1:4],
		     title = "\$\\Gamma^{($(L))}\$",
		     bbox_to_anchor = [.4725,-.06],
		     fontsize = "small",
		     handlelength = 1.,
		     labelspacing = .1,
		     framealpha = 0.0,
		     handletextpad = .4)
    ax.legend(loc = "lower left",
	      handles = [lines[2*j] for j in 1:4],
	      title = "\$\\Gamma_A\$",
	      bbox_to_anchor = [-.01,-.06],
	      fontsize = "small",
	      handlelength = 1.,
	      labelspacing = .1,
	      framealpha = 0.0,
	      handletextpad = .4)
    ax.add_artist(leg0)
    fig.show()
end

function supp4()
    #=
    Similar to supp3(), compares analytical estimate of 
    gamma to numerical solution for toy model. 
    Varying different parameters compared to supp3().
    =#
    L = 150
    As = range(.1,2.0,length = 50)
    x0s = [75,75,75,130]
    betas = [10,10,10,10]
    alphas = [.2,1.2,.2,.2]
    deltas = [4.5,4.5,9,9]
    gammas = []
    gammasA = []
    # running through parameters ##################################################
    for i in 1:4
	x0 = x0s[i]
	beta = betas[i]
	alpha = alphas[i]
	delta = deltas[i]
	gamma_tmp =zeros(length(As))
	gammaA_tmp = zeros(length(As))
	for (Ai,A) in enumerate(As)
	    bs = gen_bs3(x0,L,alpha,delta,A,beta)
	    engI = find_half_width(bs[1:end-1],bs[end])
	    gamma_tmp[Ai] = engI
	    g_tmp,prefactor = analytic_gamma(bs)
	    gammaA_tmp[Ai] = g_tmp
	end
	push!(gammas,gamma_tmp)
	push!(gammasA,gammaA_tmp)
    end
    # plotting ######################################################################
    fig,ax = plt.subplots(1,1,figsize = [6.8,1.8])
    for i in 1:4
	beta = betas[i]
	alpha = alphas[i]
	delta = deltas[i]
	x0 = x0s[i]
	ax.plot(As .* 2,gammas[i],"-x",
		label = "\$ \\alpha = $(alpha), \\delta = $(delta), n_0 = $(x0)\$",
		lw = .5,
		ms= 3.0)
	ax.plot(As .* 2,gammasA[i],"-o",
		label = "\$ \\alpha = $(alpha), \\delta = $(delta), n_0 = $(x0)\$",
		lw = .5,
		ms= 2.0)
    end
    ax.tick_params(axis = "both",which = "both",labelsize = "small",
		   direction = "in")
    ax.set_yscale("log")
    ax.set_xlabel("\$ M_0 \$")
    ax.set_ylabel("\$ \\Gamma^{($(L))}, \\Gamma_A \$")
    lines = ax.get_lines()
    leg0 = ax.legend(loc = "lower center",
		     handles = [lines[2*j-1] for j in 1:4],
		     title = "\$\\Gamma^{($(L))}\$",
		     bbox_to_anchor = [.4725,-.06],
		     fontsize = "small",
		     handlelength = 1.,
		     labelspacing = .1,
		     framealpha = 0.0,
		     handletextpad = .4)
    ax.legend(loc = "lower left",
	      handles = [lines[2*j] for j in 1:4],
	      title = "\$\\Gamma_A\$",
	      bbox_to_anchor = [-.01,-.06],
	      fontsize = "small",
	      handlelength = 1.,
	      labelspacing = .1,
	      framealpha = 0.0,
	      handletextpad = .4)
    ax.add_artist(leg0)
    fig.show()
end
