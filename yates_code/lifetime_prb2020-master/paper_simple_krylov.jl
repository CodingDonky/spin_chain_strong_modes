function simple_krylov()
    #=
    Creates plot of kyrlov time evolution 
    for both toy model and for real bn
    =#
    
    # chose times to evaluate, cannot be too
    # far into the future,
    times = 10 .^ range(-1,2,length = 40)

    L = 40 # ramp length
    L2 = 1000 # plateau length, increase for longer times, 10-20k for t~1k
    # toy model params
    A = .5
    beta = 10
    alpha = .3
    delta = 1
    x0 = 10
    # generating toy model
    bs_toy = gen_bs3(x0,L,alpha,delta,A,beta)
    # adding plateau
    bs_toy = extend_bs(L2,bs_toy)
    # krylov-time evolve, takes some time
    signal_toy = time_evolve(bs_toy,times)

    # plotting toy results
    fig,ax = plt.subplots(2,2)
    ax[1,1].plot(1:length(bs_toy),bs_toy)
    ax[2,1].plot(times,signal_toy)
    ax[2,1].set_xscale("log")
    ax[1,1].set_xscale("log")

    data2 = load("data/bs_june25_slow_all.jld")# with gram-schmidt L=6-12
    Bs = data2["Bs"]
    jzs = data2["jzs"]
    jzi = 1 # 1,2,3,4,5,6 => J_z = .2,.3,.4,.5,.6 # instead of looping over all of them
    Li = 1 # 1,2,3,4 => L = 6,8,10,12 # instead of looping over all of them
    
    bs = Bs[jzi][Li][1:400-1] # get real bn
    bs = extend_bs(L2,bs) # extend real bn
    signal = time_evolve(bs,times) # get krylov-time evolve

    # plotting "real" results
    ax[1,2].plot(1:length(bs),bs)
    ax[2,2].plot(times,signal)
    ax[1,2].set_xscale("log")
    ax[2,2].set_xscale("log")
    fig.savefig("simple_krylov.png")
end
