#=
This module is where common functions that are
closer to this current project will go.
=#

using LinearAlgebra
using Printf

# Expects Hermitian matrix, orders from low to high
function diagonalize(mat)
    u,v = eigen(mat)
    u = real.(u)
    order = sortperm(u)
    u = u[order]
    v = v[:,order]
    return u,v
end

# finds value of Z_2 for a list of states
function pairities(L,vs)
    pairities = zeros(2^L)
    for n in 1:2^L
        pairities[n] = Z2(L,Array{Complex{Float64},1}(vs[:,n]))
    end
    return pairities
end

# Expects Hermitian matrix, orders from low to high
function diagonalize_vals(mat)
    u = eigvals(mat)
    u = real.(u)
    order = sortperm(u)
    u = u[order]
    return u
end


# creates propagator from eigenbasis
function propagator(u,v,phase)
    return v*Diagonal(exp.(phase.*u))*v'
end


# If you have data for different Ls, Jzs,gs
# then print it out with this function.
# its assuming a single float for each entry.
function printout(data,Ls,jzs,gs,message;scientific = false)
    println(message)
    for (Li,L) in enumerate(Ls)
        @printf("L=% 5d|",L)
        for k in 1:length(jzs)
            @printf("jz=% 7.2f|",jzs[k])
        end
        println("")
        for (gi,g) in enumerate(gs)
            @printf("g=% 5.2f|",g)
            for (jzi,jz) in enumerate(jzs)
                if scientific
                    @printf("% 10.2e|",data[Li,gi,jzi])
                else
                    @printf("% 10.2f|",data[Li,gi,jzi])
                end
            end
            println("")
        end
        println("")
    end
end


# If you have data for different Ls, gs, no jzs
# its assuming a single float for each entry.
function printout(data,Ls,gs,message;scientific = false)
    println(message)
    @printf("%7s|","")
    for k in 1:length(gs)
        @printf("g=% 8.2f|",gs[k])
    end
    println("")
    for (Li,L) in enumerate(Ls)
        @printf("L=% 5d|",L)
        # @printf("%8s","")
        for (gi,g) in enumerate(gs)
            if scientific
                @printf("% 10.2e|",data[Li,gi])
            else
                @printf("% 10.2f|",data[Li,gi])
            end
        end
        println("")
    end
end


#=
The wavefunction will be in the block-diagonal basis.
desired time using sparse matrices and KrylovKit.
[(p = 0,m = 0 sector), (p = 0, m = 1), (p = 1, m = 0),(p = 1, m = 1)]
=#
function TimeEvolve(L::Int64,Jx::Float64,Jz::Float64,mu::Float64,time::Float64,overlaps)
    newstate = []
    for p in 0:1, m in 0:1
        # get wavefunction in this block
        index = 2*p + m + 1
        state = overlaps[index]
        # if state has no weight in this sector, just
        # return zero (Krylov stuff will fail otherwise)
        if norm(state) > 0.0
            # create the sparse matrices
            hzz = Hzz(L,Jz,p,m)
            hxx = Hxx(L,Jx,p,m)
            hz = Hz(L,mu,p,m)
            h = hxx + hz + hzz
            tmp1 = exponentiate(h,-im*time,state)
            state[:] = tmp1[1]
        end
        push!(newstate,state)
    end
    return newstate
end


function TimeEvolve_simple(h,time,state)
    tmp = exponentiate(h,-im*time,state,verbosity = 1)
    return tmp[1]
end

