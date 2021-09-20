# THIS FILE NEEDS ED CODE TO RUN
# MODIFY LOAD_PATH ACCORDINGLY
push!(LOAD_PATH,pwd()*"/ed_main/src/")
module tmp
using LinearAlgebra
using PyPlot
using ED
using Roots


# # check with get_bs.jl and gen_bs.jl...
# push!(LOAD_PATH,pwd()*"/../ED/src/")
using KrylovKit
using SparseArrays
using JLD


#=
This produces the time evolution via Kyrlov methods
(KrylovKit) for the bn of the SYK model
=#

function krylov_time_evolveSYK()
    # system parameters
    L = 8
    cs=[0,1,2,3,4,5,6,7,8,9,10]
   
    # times = 10 .^ range(-1,3,length = 40) # log scale time points
    times = range(0,1,length=10) # linear scale time points
    signal = zeros(length(times),length(cs))

    for (ci, c) in enumerate(cs)

     bs=[]
       
    for i in 0:L
    	bi = sqrt((i+1)*(i+c))
        push!(bs,bi) 
	end
	bs = Float64.(bs) 	
       
	println("c = ", "  ",c)

	println("bn= ", bs)
		
        # make a HK -> M
        M = spdiagm(1 => bs,-1=>bs)
	#println("M=   ", M," " ,"size of M=  ", size(M))
	#println("length of bs+1=  ", length(bs)+1)
        # start with wavefunction on site 1

	op = zeros(length(bs)+1); op[1] = 1.0

	
        # time evolve (KrylovKit)

	
	for (ti,time) in enumerate(times)
	    
	    tmp,info = exponentiate(M,im*time,op, tol = (1e-2)/times[end],maxiter=1000)

            println(info)


            tmp1  = real.(op'*tmp)
            # if tmp1 <.05 || info.converged == 0
            if info.converged == 0
                break
            else 
                signal[ti,ci] = tmp1
        	end

		if abs.(signal[ti,ci])<=exp(-1)*abs.(signal[1,ci])
		   println("Ainf has decayed: "," ","ti= ",ti)
		    end

		    println("time =", time, "  ","signal = ", signal[ti,ci])		    

		end

	
	end

    # rest is plotting
    
    markers = ["-o","-s","-^","-*","->"]
    
    fig,ax = plt.subplots(1,1)
    	   for (ci, c) in enumerate(cs)
        ax.plot(times,signal[:,ci],"-o",ms= 4.,lw = 1.0, label = "\$c = $(c)\$")
	 end
	ax.axhline(y = exp(-1),color = "red",alpha = .5)
    ax.set_xscale("log")
    ax.set_xlabel("\$ t \$")
    ax.set_ylabel("\$ A_\\infty(t) \$")
    ax.legend()
    plt.show()
    
end

#= This performs time-evolution with uniform b's=#


function krylov_time_metal()


 # system parameters
    L = 40
    cs=[1,2,3] 
   
    # times = 10 .^ range(-1,3,length = 40) # log scale time points
    times = range(0,10,length=10) # linear scale time points
    signal = zeros(length(times),length(cs))

    for (ci, c) in enumerate(cs)

     bs=[]
       
    for i in 0:L
    	bi = c #Uniform hopping
        push!(bs,bi) 
	end
	bs = Float64.(bs) 	
       
	println("c = ", "  ",c)

	println("bn= ", bs)
		
        # make a HK -> M
        M = spdiagm(1 => bs,-1=>bs)
	#println("M=   ", M," " ,"size of M=  ", size(M))
	#println("length of bs+1=  ", length(bs)+1)
        # start with wavefunction on site 1

	op = zeros(length(bs)+1); op[1] = 1.0

	
        # time evolve (KrylovKit)

	
	for (ti,time) in enumerate(times)
	    
	    tmp,info = exponentiate(M,im*time,op, tol = (1e-2)/times[end],maxiter=1000)

            println(info)

	   

            tmp1  = real.(op'*tmp)
            # if tmp1 <.05 || info.converged == 0
            if info.converged == 0
                break
            else 
                signal[ti,ci] = tmp1
        	end

		if abs.(signal[ti,ci])<=exp(-1)*abs.(signal[1,ci])
		   println("Ainf has decayed: "," ","ti= ",ti)
		    end

		 println("time =", time, "  ","signal = ", signal[ti,ci])   		    

		end


 		

	end

    
    # rest is plotting
    
    markers = ["-o","-s","-^","-*","->"]
    
    fig,ax = plt.subplots(1,1)
    	   for (ci, c) in enumerate(cs)
        ax.plot(times,signal[:,ci],"-o",ms= 4.,lw = 1.0, label = "\$c = $(c)\$")
	 end
	ax.axhline(y = exp(-1),color = "red",alpha = .5)
    #ax.set_xscale("log")
    ax.set_xlabel("\$ t \$")
    ax.set_ylabel("\$ A_\\infty(t) \$")
    ax.legend()
    plt.show()
    
end



#@timev krylov_time_evolveSYK()

@timev krylov_time_metal()

end