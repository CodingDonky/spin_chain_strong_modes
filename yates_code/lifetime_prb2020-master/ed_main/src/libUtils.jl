#=
This is where assorted functions go..
=#

using LinearAlgebra
using Dates

function get_uvs(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64)
    #=
    This will find the eigenspectrum the smart way and 
    create the eigenvalue matrix U and the 
    eigenvector matrix V. Note that set L>=  L0 if you want to
    use binary basis for the computation for L<L0.
    =#
    us = []; vs = []
    if L >= 14 && L%2 == 0
        for (p1i, p) in enumerate(0:1), (mi,m) in enumerate(0:1)
            index = 2*p + m + 1
            u,v = eigen(H_Delta(L,jx,jz,g,gamma,p,m))
            us = vcat(us,u)
            vtmp = zeros(Complex{Float64},2^L,length(v[1,:]))
            for i in 1:length(v[1,:])
                tmp2 = block_to_binary(L,v[:,i],p,m)
                vtmp[:,i] = tmp2
            end
            if length(vs) == 0
                vs = vtmp
            else
                vs = hcat(vs,vtmp)
            end
        end
    else
        us,vs = eigen(H_Delta(L,jx,jz,g,gamma))
    end
    return us, vs
end


function get_Ainf(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64,times::AbstractArray)
    #=
    This is partially optimized.
    =#
    ainfs = zeros(Complex{Float64}, length(times))
    sigx = Complex{Float64}.(sigX(L,1))
    us,vs = get_uvs(L,g,jx,jz,gamma)
    sigxV = vs'*sigx*vs
    for ti in 1:length(times)
        print("ed: ti = $ti\r")
        U = Diagonal(exp.(-im*times[ti].*us))
        ainfs[ti] = 1/(2^L)*tr(U'*sigxV*U*sigxV)
    end
    return real.(ainfs)
end

function get_Ainf_Edge_mode(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64,times::AbstractArray,
                coeffs::Array{Float64}, Ops)
    #=
    same as get_Ainf, but now try to use edge mode as starting point
    =#
    ainfs = zeros(Complex{Float64}, length(times))
    #this will hold the overlaps with the krylov basis
    ainfs2 = zeros(Complex{Float64},length(Ops),length(times))
    # this will hold the autocorrelations of the krylov basis
    ainfs3 = zeros(Complex{Float64},length(Ops),length(times))

    edgemode = zeros(2^L,2^L)
    for (ci,coeff) in enumerate(coeffs)
        edgemode += Ops[ci].*coeff
    end

    sigx = sigX(L,1)
    us,vs = get_uvs(L,g,jx,jz,gamma)
    edge_time = vs'*edgemode*vs

    Ops2 = []
    for i in 1:length(Ops)
        push!(Ops2,vs'*Ops[i]*vs)
        # just add a dagger in trace below..
        # if i%2 == 0
        #     # add an imaginary factor to make these 
        #     # operators hermitian..
        #     push!(Ops2,im .* vs'*Ops[i]*vs)
        # else
        #     push!(Ops2,vs'*Ops[i]*vs)
        # end
    end

    for ti in 1:length(times)
        print("ed: ti = $ti\r")
        U = Diagonal(exp.(-im*times[ti].*us))
        ainfs[ti] = 1/(2^L)*tr(U'*edge_time*U*edge_time)
        for i in 1:length(Ops)
            ainfs2[i,ti] = 1/(2^L)*tr(U'*edge_time'*U*Ops2[i])
            ainfs3[i,ti] = 1/(2^L)*tr(U'*Ops2[i]'*U*Ops2[i])
        end
    end
    return real.(ainfs), real.(ainfs2), real.(ainfs3)
end



function get_L_ainf(M,times)
    u,v = eigen(M)
    v1 = v[1,:]
    v1c = conj.(v1)
    ainf = zeros(length(times))
    for (ti,time) in enumerate(times)
        print("time = $ti\r")
        ainf[ti] = real(transpose(v1)*Diagonal(exp.(im .* u .* time))*v1c)
    end 
    return ainf
end 

function make_M(bs)
    len = length(bs)+1
    M = zeros(len,len)
    for (bi,b) in enumerate(bs)
        M[bi,bi+1] = b
        M[bi+1,bi] = b
    end 
    return M
end


function make_M2(bs)
    len = length(bs)+1
    M = zeros(len,len)
    for (bi,b) in enumerate(bs)
        M[bi,bi+1] = -b
        M[bi+1,bi] =  b
    end 
    return M
end


function make_test_bs(L,alpha,delta,offset)
    bs = [offset + n*alpha + (1 + (-1)^n)/2*delta for n in 1:L]
    return bs
end


function make_test_bs2(L,alpha,delta,offset)
    bs = [offset + n*alpha + (1 + (-1)^n)/2*(delta/n) for n in 1:L]
    return bs
end



function make_my_edge_mode(bs)
    coeffs = []
    errors = []
    norms = []
    if length(bs)%2 == 1
        ratios = bs[1:2:end-1]./bs[2:2:end-1]
    else
        ratios = bs[1:2:end]./bs[2:2:end]
    end
    for n in 1:length(ratios)
        tmp = (-1)^(n+1)
        for m in 1:n-1
            tmp = tmp*ratios[m]
        end
        push!(coeffs,tmp)
        push!(errors,tmp*bs[2*n-1])
        push!(norms,sqrt(coeffs'*coeffs))
    end 
    return coeffs,errors,norms
end


function make_my_edge_mode2(bs)
    coeffs = []
    errors = []
    norms = []
    if length(bs)%2 == 1
        ratios = bs[1:2:end-1]./bs[2:2:end-1]
    else
        ratios = bs[1:2:end]./bs[2:2:end]
    end
    for n in 1:length(ratios)
        tmp = 1
        for m in 1:n-1
            tmp = tmp*ratios[m]
        end
        push!(coeffs,tmp)
        push!(errors,tmp*bs[2*n-1])
        push!(norms,sqrt(coeffs'*coeffs))
    end 
    return coeffs,errors,norms
end


function get_lifetime(times,data,cutoff = 0.1,window = 10)
    lifetime = 0
    for (ti,time) in enumerate(times[1:end-window])
        if 1/window*sum((data[ti:ti+window]))<cutoff # && time > 20
            lifetime = time
            break
        end
    end
    return lifetime
end


function get_fit(xs,ys)
    N1 = length(xs)
    avgX = sum(xs)/N1
    avgY = sum(ys)/N1
    slope = (xs'*ys - avgX*avgY*N1)/(xs'*xs - avgX^2*N1)
    inter = (xs'*ys - slope*(xs'*xs))/(avgX*N1)
    return slope, inter

end


function get_bn_fit(bs)
    xs = 1:length(bs)
    y1 = bs[1:2:end]
    x1 = xs[1:2:end]
    y2 = bs[2:2:end]
    x2 = xs[2:2:end]
    alpha_1,alpha0_1 = get_fit(x1,y1)
    alpha_2,alpha0_2 = get_fit(x2,y2)

    alpha =  (alpha_1 + alpha_2)/2
    delta = alpha0_2 - alpha0_1
    offset = alpha0_1
    return alpha, delta, offset
end



function find_domains(bs)
    #=
    This will find the positions where 
    staggering runs into an error.
    assumes the correct pattern is..
    down up down up down up...
    =#
    grow_shrink = []
    # values = []
    for bi in 1:2:length(bs)-1
        tmp = bs[bi]/bs[bi+1]
        if tmp>1.0
            push!(grow_shrink,1.0)
            push!(grow_shrink,1.0)
        else 
            push!(grow_shrink,-1.0)
            push!(grow_shrink,-1.0)
        end
    end
    return grow_shrink
end


function get_Ayy_inf(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64,times::AbstractArray)
    #=
    Unoptimized autocorrelation for y1-y1
    =#
    ainfs = zeros(Complex{Float64}, length(times))
    sigx = sigX(L,1)
    sigz = sigZ(L,1)
    sigy = im.*(sigz*sigx)
    us,vs = eigen(H_Delta(L,jx,jz,g,gamma))
    sigyV = vs'*sigy*vs
    for ti in 1:length(times)
        print("ed: ti = $ti\r")
        U = Diagonal(exp.(-im*times[ti].*us))
        ainfs[ti] = 1/(2^L)*tr(U'*sigyV*U*sigyV)
    end
    return real.(ainfs)
end


function get_Azz_inf(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64,times::AbstractArray)
    #=
    Unoptimized autocorrelation for z1-z1
    =#
    ainfs = zeros(Complex{Float64}, length(times))
    sigz = sigZ(L,1)
    us,vs = eigen(H_Delta(L,jx,jz,g,gamma))
    sigzV = vs'*sigz*vs
    for ti in 1:length(times)
        print("ed: ti = $ti\r")
        U = Diagonal(exp.(-im*times[ti].*us))
        ainfs[ti] = 1/(2^L)*tr(U'*sigzV*U*sigzV)
    end
    return real.(ainfs)
end


function get_Ainf_phase_floquet(L::Int64,
                                g::Float64,
                                jx::Float64,
                                jz::Float64,
                                gamma::Float64,
                                times::AbstractArray)
    #=
    =#
    ainfs = zeros(Complex{Float64}, length(times))
    sigz = sigZ(L,1)
    us,vs = eigen(H_Delta(L,jx,jz,g,gamma))
    sigzV = vs'*sigz*vs
    for ti in 1:length(times)
        print("ed: ti = $ti\r")
        U = Diagonal(exp.(-im*times[ti].*us))
        ainfs[ti] = 1/(2^L)*tr(U'*sigzV*U*sigzV)
    end
    return real.(ainfs)
end



function get_spectrum(L::Int64,g::Float64,jx::Float64,jz::Float64, gamma::Float64;outfile = "tmp_out.txt")
    #=
    This is partially optimized.
    THIS FUNCTION APPEARS TO BE NOT NEEDED...
    =#
    us = []; vs = []
    open(outfile,"a") do file
        write(file,"starting ED calculation\n")
    end

    if L >= 4 && L%2 == 0
        for (p1i, p) in enumerate(0:1), (mi,m) in enumerate(0:1)
            open(outfile,"a") do file
                write(file,"L = $(L), g = $(g), jz = $(jz), gamma = $(gamma), p = $(p), m = $(m), time = $(Dates.now())\n")
            end
            index = 2*p + m + 1
            u,v = eigen(H_Delta(L,jx,jz,g,gamma,p,m))
            us = vcat(us,u)
            vtmp = zeros(Complex{Float64},2^L,length(v[1,:]))
            for i in 1:length(v[1,:])
                tmp2 = block_to_binary(L,v[:,i],p,m)
                vtmp[:,i] = tmp2
            end
            if length(vs) == 0
                vs = vtmp
            else
                vs = hcat(vs,vtmp)
            end
        end
    else
        us,vs = eigen(H_Delta(L,jx,jz,g,gamma))
    end
    return (us,vs)
end



















