push!(LOAD_PATH,pwd()*"/ed_main/src/")
module email_response
using LinearAlgebra
using PyPlot
using ED

function test()
    steps = 20 # number of lanczos steps
    times = 10 .^ range(-1,4,length = 100) # times for ainf
    # vary Jz
    jzs = range(0,.3,length = 100)
    L = 6
    # want negative coeffs
    j = -1.0
    gamma = -1.0 # J_y = 0 
    # storage for energies
    us = zeros(2^L,length(jzs))
    Ainfs = [] # storage for Ainfs
    Bs = []    # storage for Bn's 
    jz2 = [] # store particular jz values for plotting
    for (jzi,jz) in enumerate(jzs)
        # their relation (negated)
        g =  2*sqrt(jz^2 + jz)
        print("g = $(g)\r")
        # get ham
        H = ED.H_Delta(L,j,jz,g,gamma)
        # fix edges
        H = H - g/2 .*(ED.sigZ(L,1) + ED.sigZ(L,L))
        # get eigenvalues
        u = eigvals(H)
        us[:,jzi] = u
        # for some jz values, find Bn's and ainfs
        if jzi%15==0 || jzi == 2
            push!(jz2,jz) # store for plotting
            
            # do lanczos algorithm
            function L_op(op) # defining L
                return H*op - op*H
            end
            # starting operator is sigma_1^y
            O1 = -im.*ED.sigZ(L,1)*ED.sigX(L,1)
            phi = copy(O1) # build phi in ED
            coeff = 1.0 # coeff for phi
            O2 = L_op(O1)
            b1 = sqrt(tr(O2'*O2)/(2^L))
            O2 = O2 ./ b1
            bs = [b1]
            Os = [O2,O1]
            for n in 1:steps
                O_tmp = L_op(Os[1]) - bs[end] .* Os[2]
                b_tmp = sqrt(tr(O_tmp'*O_tmp)/(2^L))
                O_tmp = O_tmp ./ b_tmp
                push!(bs,b_tmp)
                Os = [O_tmp,Os[1]]
                if n%2==1
                    coeff = coeff*(-bs[end-1]/bs[end])
                    phi = phi + coeff.*O_tmp
                end
            end
            # save resulting Bn's
            push!(Bs,bs)
            # normalize phi
            phi = phi ./ sqrt(tr(phi'*phi)/(2^L))
            
            
            # calculate Ainf.. doing it by hand to
            # avoid sign issues in parameters
            # get full spectrum, not just eigenvalues
            u,v = eigen(H)
            # swap out operators to see different
            # time evolutions, phi or sigma_1^x
            # currently, its sigma_1^y
            # (phi will send the sholder of the exp to 1.0)
            sig1 =  copy(O1) # copy(phi)# ED.sigX(L,
            ainf = zeros(length(times))
            for (ti,time) in enumerate(times)
                U = v*Diagonal(exp.(-im*(time) .* u))*(v')
                ainf[ti] = real(tr(U'*sig1*U*sig1)/(2^L))
            end
            push!(Ainfs,ainf)
            println("\n ground state energies:\n",u[1])
            println(u[2])
            println("-(L-1)(jz + 1):\n", -(L-1)*(jz + 1))
            gs1 = v[:,1]
            gs2 = v[:,2]
            println(  "tr(y_1 phi)/(2^L): ",abs(tr(phi'*O1)/(2^L))^2)
            println(  "|<1|phi|1>|^2: ",abs(gs1'*phi*gs1)^2)
            println(  "|<2|phi|2>|^2: ",abs(gs2'*phi*gs2)^2)
            println(  "|<2|phi|1>|^2: ",abs(gs2'*phi*gs1)^2)
            println(  "|<1|phi|2>|^2: ",abs(gs1'*phi*gs2)^2)
        end
    end


    
    # rest is plotting
    fig,ax = plt.subplots(3,1)
    # plotting spectrum vs jz
    ax[1].plot(jzs,us',color = "black",
               lw = .5)
    # plotting their prediction of ground state energy
    ax[1].plot(jzs,-(L-1).*(jzs .+ 1) ,color = "red")
    # for the particular jz values, show Bn's calculated
    for k in 1:length(jz2)
        jz = round(jz2[k],sigdigits = 2)
        ax[1].axvline(x = jz2[k],color = "C$(k-1)",zorder = 0,
                      label = "jz = $(jz)")
        ax[2].plot(1:length(Bs[k]),real.(Bs[k]),"-o",color = "C$(k-1)",
                   label = "jz = $(jz)")
        ax[3].plot(times,Ainfs[k],"-",color = "C$(k-1)",
                   label = "jz = $(jz)")
    end
    ax[1].set_xlabel("jz")
    ax[2].set_xlabel("n")
    ax[3].set_xlabel("time")
    ax[1].set_ylabel("energy")
    ax[2].set_ylabel("bn")
    ax[3].set_ylabel("ainf")
    ax[3].set_xscale("log")
    ax[1].legend(fontsize = "small",
                 ncol = 2,
                 loc = "upper right",
                 labelspacing = .5,
                 columnspacing = .5,
                 framealpha = 1.0)
    ax[2].legend(fontsize = "small",
                 ncol = 2,
                 loc = "lower right",
                 labelspacing = .5,
                 columnspacing = .5,
                 framealpha = 1.0)
    ax[3].legend(fontsize = "small",
                 ncol = 2,
                 loc = "lower left",
                 labelspacing = .5,
                 columnspacing = .5,
                 framealpha = 1.0)
end
@timev test()
end
