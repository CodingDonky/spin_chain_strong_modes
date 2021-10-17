# Utility functions ###############################
function gen_bs(L,alpha,delta,rho)
    #=
    Generate toy bn for two slope model. 
    Only used for Fig3. 
    =#
    bs = [alpha]
    Nstar = 1
    switch = false
    for k in 1:L
	odd = alpha*(2*k+1)
	even = delta + rho*alpha*(2*k)
	if even < odd && !switch
	    Nstar = 2*k+1
	    switch = true
	end
	push!(bs,even);push!(bs,odd)
    end
    return bs,Nstar
end

function extend_bs(L2,bs)
    #=
    Extends bns by copying the last 
    hopping L2 times. 
    =#
    bs2 = copy(bs)
    tmpend = bs2[end]
    for k in 1:L2
	push!(bs2,tmpend)
    end
    return bs2
end

function gen_bs3(x0,L2,alpha,delta,A,beta)
    #=
    Second toy model in paper. A = M_0/2.
    Generates ramp, flat part can be added later on.
    =#
    bs = zeros(L2)
    for Li in 1:L2
	bs[Li] = alpha*Li + delta + (-1)^Li*A/((Li/x0)^(beta) + 1)
    end
    return bs
end

function time_evolve(bs,times)
    signal = zeros(length(times))
    L = length(bs)
    M = spdiagm(1 => bs,-1=>bs)
    op = zeros(Complex{Float64},L+1); op[1] = 1.0
    for (ti,time) in enumerate(times)
        tmp4,info = exponentiate(M,im*time,op,tol = (1e-2)/times[end],
                                 maxiter = 4000)
        tmp5  = real.(op'*tmp4)
        println("time step $(ti): info $(info), data point = $(tmp5)")
        if info.converged == 0
            break
        else
            signal[ti] = tmp5
        end
    end
    return signal
end

function load_data1()
    # this data will generate nice lifetime plots like Nayak17.pdf
    # loading up bns, from lanczos algo
    d1 = load("data/bs_jan4_test5.jld") # .25, .35, .9 = (g,jz,gamma)
    d2 = load("data/bs_jan_test.jld")   # .35, .35, .9 
    d3 = load("data/bs_jan4_test2.jld") # .45, .35, .9 
    d4 = load("data/bs_jan7_test6.jld") # .55, .35, .9 
    d5 = load("data/bs_jan7_test7.jld") # .65, .35, .9 

    # LOADING UP Ainfs, from ed
    e1 = load("data/ainf_jan4_test5.jld") # .25,.35,.9 = (g,jz,gamma)
    e2 = load("data/ainf_jan12_2.jld")    # .35,.35,.9 
    e3 = load("data/ainf_jan2.jld")       # .45,.35,.9 
    e4 = load("data/ainf_jan7_test6.jld") # .55,.35,.9 
    e5 = load("data/ainf_jan7_test7.jld") # .65,.35,.9 

    datas_bs = [d1,d2,d3,d4,d5]
    datas_ainfs = [e1,e2,e3,e4,e5]
    Ls = d1["Ls"] # all the same 
    jz = d1["jz"] # all the same
    gamma = d1["gamma"] # all the same
    j = d1["j"] # all the same
    threshold = d1["threshold"] # all the same
    steps = d1["steps"] # all the same
    gs = []
    for (di,d) in enumerate(datas_bs)
        push!(gs,d["g"])
    end
    return [datas_bs,
            datas_ainfs,
            Ls, gamma, j, threshold, steps, gs, jz]
end

function load_data2()
    # this data will generate nice lifetime plots like Nayak17.pdf
    # loading up bns, from lanczos algo
    d1 = load("data/bs_2_jan7.jld")       # .3, .2, .9 = (g,jz,gamma)
    d2 = load("data/bs_jan19_test8.jld")  # .3, .3, .9               
    d3 = load("data/bs_jan19_test9.jld")  # .3, .4, .9               
    d4 = load("data/bs_jan19_test10.jld") # .3, .5, .9               
    d5 = load("data/bs_jan19_test11.jld") # .3, .6, .9               

    # LOADING UP Ainfs, from ed
    e1 = load("data/ainf_l14_dec27_jzp2_gammap9_gp3.jld") # .3,.2,.9 = (g,jz,gamma)
    e2 = load("data/ainf_jan19_test8.jld")                # .3,.3,.9 
    e3 = load("data/ainf_jan19_test9.jld")                # .3,.4,.9 
    e4 = load("data/ainf_jan19_test10.jld")               # .3,.5,.9 
    e5 = load("data/ainf_jan19_test11.jld")               # .3,.6,.9 

    datas_bs = [d1,d2,d3,d4,d5]
    datas_ainfs = [e1,e2,e3,e4,e5]
    Ls = d1["Ls"] # all the same 
    g = d1["g"] # all the same
    gamma = d1["gamma"] # all the same
    j = d1["j"] # all the same
    threshold = d1["threshold"] # all the same
    steps = d1["steps"] # all the same
    jzs = []
    for (di,d) in enumerate(datas_bs)
        push!(jzs,d["jz"])
    end
    return [datas_bs,
            datas_ainfs,
            Ls, gamma, j, threshold, steps, g, jzs]
end

function load_data3()
    # loading up bns, from lanczos algo
    d1 = load("data/bs_jan22_test12.jld") # .2 ,.2 , 1.0 = (g,jz,gamma)
    d2 = load("data/bs_jan22_test13.jld") # .3 ,.3 , 1.0
    d3 = load("data/bs_jan22_test14.jld") # .4 ,.4 , 1.0
    d4 = load("data/bs_jan22_test15.jld") # .5 ,.5 , 1.0
    d5 = load("data/bs_jan22_test16.jld") # .6 ,.6 , 1.0
    d6 = load("data/bs_jan22_test17.jld") # .7 ,.7 , 1.0
    d7 = load("data/bs_jan22_test18.jld") # .1 ,.1 , 1.0
    d8 = load("data/bs_jan22_test19.jld") # .25,.25, 1.0
    d9 = load("data/bs_jan22_test20.jld") # .35,.35, 1.0
    # LOADING UP Ainfs, from ed
    e1 = load("data/ainf_jan21_test12.jld") # .2 ,.2 , 1.0 = (g,jz,gamma)
    e2 = load("data/ainf_jan21_test13.jld") # .3 ,.3 , 1.0
    e3 = load("data/ainf_jan21_test14.jld") # .4 ,.4 , 1.0
    e4 = load("data/ainf_jan21_test15.jld") # .5 ,.5 , 1.0
    e5 = load("data/ainf_jan21_test16.jld") # .6 ,.6 , 1.0
    e6 = load("data/ainf_jan21_test17.jld") # .7 ,.7 , 1.0
    e7 = load("data/ainf_jan21_test18.jld") # .1 ,.1 , 1.0
    e8 = load("data/ainf_jan21_test19.jld") # .25,.25, 1.0
    e9 = load("data/ainf_jan21_test20.jld") # .35,.35, 1.0
    datas_bs = [d7,d1,d8,d2,d9,d3,d4,d5,d6]
    datas_ainfs = [e7,e1,e8,e2,e9,e3,e4,e5,e6]
    Ls = d1["Ls"] # all the same 
    gamma = d1["gamma"] # all the same
    j = d1["j"] # all the same
    threshold = d1["threshold"] # all the same
    steps = d1["steps"] # all the same
    gs = []
    jzs = []
    for (di,d) in enumerate(datas_bs)
        push!(gs,d["g"])
        push!(jzs,d["jz"])
    end
    return [datas_bs,
            datas_ainfs,
            Ls, gamma, j, threshold, steps, gs, jzs]
end

function load_data4(dataset)
    # we want dataset 2 for a_inf data, this is from Jan 2020
    # this is from Jan 2020    
    if dataset == 1
	datas_bs,datas_ainfs,Ls,gamma,j,threshold,steps,gs,jzs = load_data1()
	plot_label = ["g",gs]
	other_labels = ["\\gamma",gamma,"J_z",jzs]
	D = length(datas_bs)
	ds = 1:D
    elseif dataset ==2
	datas_bs,datas_ainfs,Ls,gamma,j,threshold,steps,gs,jzs = load_data2()
	plot_label = ["J_z",jzs]
	other_labels = ["\\gamma",gamma,"g",gs]
	D = length(datas_bs)
	ds = 1:D
    else
	datas_bs,datas_ainfs,Ls,gamma,j,threshold,steps,gs,jzs = load_data3()
	plot_label = ["J_z",jzs] # and gzs
	other_labels = ["\\gamma",gamma,"g","J_z"]
	D = length(datas_bs)
	ds = 2:D-1 # 3:D-1
    end

    return (datas_bs,datas_ainfs,Ls,
	    plot_label,other_labels,D,
	    ds,j,threshold,steps)
end

function zero_modes(bs)
    N = length(bs)
    if N%2 == 0
        error("N has to be odd")
    end
    psi_even = zeros(N+1)
    psi_odd = zeros(N+1)
    psi_odd[1] = 1.0
    psi_even[2] = 1.0
    i = 1
    while 2*i-1 < N 
        psi_odd[2*i+1]  = -bs[2*i-1]/bs[2*i]*psi_odd[2*i-1]
        psi_even[2*i+2] = -bs[2*i]/bs[2*i+1]*psi_even[2*i]
        i +=1
    end
    psi_even = psi_even ./ sqrt(psi_even'*psi_even)
    psi_odd = psi_odd ./ sqrt(psi_odd'*psi_odd)
    return (psi_odd,psi_even)
end

function simple_best_fit(x,y)
    # least squares
    N = length(x)
    xavg = sum(x)/N
    yavg = sum(y)/N
    m = (y'*x - N*xavg*yavg)/(x'*x - xavg*xavg*N)
    b = yavg - m*xavg
    return m,b
end

function num_greens_func(bs,t0,engs)
    H = ED.make_M(bs) .* (1 + 0.0*im)
    data = zeros(Complex{Float64},length(engs))
    vec = zeros(length(bs)+1);vec[1] = 1.0
    for (engi,eng) in enumerate(engs)
	H[end,end] = 1/2*(eng - im *sqrt(4*t0^2-eng^2))
	data[engi] =((eng*I - H) \ vec)[1]
    end
    return data
end

function find_half_width(bs,t0,epsilon = 1e-16,max_iter= 100,engL = 0.0,engR = 2.0)
    #=
    Numerical solution for Gamma, searching (binary search) for half-width.
    bs -> hoppings for H_N
    t0 -> coupling and metal hopping
    epsilon -> precision of search
    max_iter -> max number of iterations for search
    engL -> starting left end of search
    engR -> starting right end of search 

    engR, if set too high can give errors for complex values, because
    it will be beyond the  band width of the system. 
    Setting it too low, will cause the binary search to miss the correct
    value of Gamma, the default value should work well for most cases.
    =#

    # making H_N
    Nb = length(bs)
    H = ED.make_M(bs) .*(1+0.0*im)
    vec = zeros(Nb+1);vec[1] = 1.0
    HB = copy(H)
    eng = 0.0
    # setting self-eng
    HB[end,end] = 1/2*(eng - im *sqrt(4t0^2-eng^2))
    # finding half-max -> val0
    val0 = abs(imag(((eng*I - HB) \ vec)[1]))/2 #1/2 val at e=0
    # setting up binary search
    iter = 0
    engI = (engR + engL)/2.0	
    res = 1.0
    # running binary search
    while iter <max_iter
	HB[end,end] = 1/2*(engI - im *sqrt(4t0^2-engI^2))
	val = abs(imag(((engI*I - HB) \ vec)[1]))
	res = val - val0
	if abs(res) < epsilon
	    break
	elseif res > 0
	    engL = engI
	else # if res < 0
	    engR = engI
	end
	iter += 1
	engI = (engR + engL)/2.0
    end
    println(iter," ",engI," ",abs(res)," ",val0," ",abs(res/val0))
    return engI
end

function find_corrected_half_width(bs,sigma,epsilon=1e-16,max_iter=100,engL=0.0,engR=2.0)
    eng = 0.0
    iter = 0
    vec0 = zeros(length(bs)+1); vec0[1] = 1
    H2 = ED.make_M(bs) .*(1+0.0*im)
    H2[end,end] = sigma 
    val0 = abs(imag(((eng*I - H2) \ vec0)[1]))/2 #1/2 val at e=0
    engI = (engR + engL)/2.0	
    res = 1.0
    while iter <max_iter
	val = abs(imag(((engI*I - H2) \ vec0)[1]))
	res = val - val0
	if abs(res) < epsilon
	    break
	elseif res > 0
	    engL = engI
	else # if res < 0
	    engR = engI
	end
	iter += 1
	engI = (engR + engL)/2.0
    end
    println(iter," ",engI," ",abs(res)," ",val0," ",abs(res/val0))
    return engI
end

function analytic_gamma_noR(bs)
    #=
    analytic estimate for gamma, 
    requires bs, from 1 to N (N must be even)
    bs[1:N-1] for H_N
    bs[end] for H_N-metal coupling
    =#
    gr = -im/bs[end] # bulk greens function
    phi_o,phi_e = zero_modes(bs[1:end-1]) # zero modes (normalized)
    Np = 1/phi_o[1]^2 
    phi_o = phi_o.^2 .* Np # undoing normalization
    gammaA = abs(imag(phi_o[end-1]/(gr*Np)))
    return gammaA,1/Np
end

function analytic_gamma(bs)
    sigma = -im*bs[end]
    phi_o,phi_e = zero_modes(bs[1:end-1])
    Np = 1/phi_o[1]^2
    Ne = 1/phi_e[2]^2
    phi_o = phi_o.^2 .* Np
    R = phi_o[end-1]^2*Ne*bs[end]^4/((sigma^2)*bs[1]^2)
    gammaA = abs(imag(phi_o[end-1]*bs[end]^2/(sigma*(Np + R))))
    return gammaA,1/(Np + R)

    # these are the old implementations
    # they are missing a factor of 1/(b_1^2)
    # in the R
    
    # gr = -im/bs_tmp[end]
    # phi_o,phi_e = zero_modes(bs_tmp[1:49])
    # Np = 1/phi_o[1]^2
    # Ne = 1/phi_e[2]^2
    # phi_o = phi_o.^2 .* Np
    # gammaA_tmp = abs(imag(phi_o[49]/(gr*(Np + phi_o[49]^2*Ne/gr^2))))

    # gr = -im/bs[nR]
    # phi_o,phi_e = zero_modes(bs[1:nR-1])
    # Np = 1/phi_o[1]^2
    # Ne = 1/phi_e[2]^2
    # phi_o = phi_o.^2 .* Np
    # gammaA = abs(imag(phi_o[nR-1]/(gr*(Np + phi_o[nR-1]^2*Ne/gr^2))))

    # gr = -im/bs[end]
    # phi_o,phi_e = zero_modes(bs[1:end-1])
    # Np = 1/phi_o[1]^2
    # phi_o = phi_o.^2 .* Np
    # g_tmp = abs(imag(phi_o[end-1]/(gr*Np)))
end
