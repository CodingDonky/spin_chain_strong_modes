#=
Quantum spin chain functions for ED.
This is for open boundary conditions.
Block diagonalization takes advantage of spin parity and
inversion symmetry.

This code was based off of the helpful notes by Prof. Sandvik
https://physics.bu.edu/~sandvik
https://physics.bu.edu/~sandvik/perimeter/

"binary" is the "simple" basis (bit basis), not the
block-diagonalized basis (p,m basis)
p = 0,1 and m = 0,1
=#

# binary basis: \sigma_i^x \sigma_{i+1}^x
function Hxx(L::Int64, Jx::Float64)
    Hxx = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            j = i+1
            b = bitflip(a,i,j)
            Hxx[b+1,a+1] += Jx
        end
    end
    return Hxx
end

# sparse binary basis: \sigma_i^x \sigma_{i+1}^x
function HxxSP(L::Int64, Jx::Float64)
    Hxx = spzeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            j = i+1
            b = bitflip(a,i,j)
            Hxx[b+1,a+1] += Jx
        end
    end
    return Hxx
end

# binary basis \sigma_i^y \sigma_{i+1}^y
function Hyy(L::Int64,Jy::Float64)
    Hyy = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            j = i+1
            b = bitflip(a,i,j)
            coeff = Jy*bitsign(a,i)*bitsign(a,j)
            # Hyy[b+1,a+1] += coeff
            Hyy[b+1,a+1] -= coeff
        end
    end
    return Hyy
end

# sparse binary basis \sigma_i^y \sigma_{i+1}^y
function HyySP(L::Int64,Jy::Float64)
    Hyy = spzeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            j = i+1
            b = bitflip(a,i,j)
            coeff = Jy*bitsign(a,i)*bitsign(a,j)
            # Hyy[b+1,a+1] += coeff
            Hyy[b+1,a+1] -= coeff
        end
    end
    return Hyy
end


# block basis: \sigma_i^x \sigma_{i+1}^x
function Hxx(L::Int64,Jx::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hxx = zeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        Ra = findPeriod(a,L)
        for i in 0:(L-2)
            b = bitflip(a,i,i+1)
            Rb = findPeriod(b,L)
            c,l = findRep(b,L)
            k = locateState(c,ref_states)
            if k>0
                Hxx[k,n] += Jx*(-1)^(l*m)*sqrt(Ra/Rb)
            end
        end
    end
    return Hxx
end

# block basis: \sigma_i^x \sigma_{i+1}^x, SPARSE
function HxxSP(L::Int64,Jx::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hxx = spzeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        Ra = findPeriod(a,L)
        for i in 0:(L-2)
            b = bitflip(a,i,i+1)
            Rb = findPeriod(b,L)
            c,l = findRep(b,L)
            k = locateState(c,ref_states)
            if k>0
                Hxx[k,n] += Jx*(-1)^(l*m)*sqrt(Ra/Rb)
            end
        end
    end
    return Hxx
end


# block basis \sigma_i^y \sigma_{i+1}^y
function Hyy(L::Int64,Jy::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hyy = zeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        Ra = findPeriod(a,L)
        for i in 0:(L-2)
            b = bitflip(a,i,i+1)
            Rb = findPeriod(b,L)
            c,l = findRep(b,L)
            k = locateState(c,ref_states)
            coeff = Jy*bitsign(a,i)*bitsign(a,i+1)
            if k>0
                # Hyy[k,n] += (-1)^(l*m)*sqrt(Ra/Rb)*coeff
                Hyy[k,n] -= (-1)^(l*m)*sqrt(Ra/Rb)*coeff
            end
        end
    end
    return Hyy
end

# block basis \sigma_i^y \sigma_{i+1}^y
function HyySP(L::Int64,Jy::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hyy = spzeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        Ra = findPeriod(a,L)
        for i in 0:(L-2)
            b = bitflip(a,i,i+1)
            Rb = findPeriod(b,L)
            c,l = findRep(b,L)
            k = locateState(c,ref_states)
            coeff = Jy*bitsign(a,i)*bitsign(a,i+1)
            if k>0
                # Hyy[k,n] += (-1)^(l*m)*sqrt(Ra/Rb)*coeff
                Hyy[k,n] -= (-1)^(l*m)*sqrt(Ra/Rb)*coeff
            end
        end
    end
    return Hyy
end



# periodic boundary conditions, binary basis
function HxxP(L::Int64, Jx::Float64)
    pbc_xx = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        b = bitflip(a,L-1,0)
        pbc_xx[b+1,a+1] += Jx
    end
    return Hxx(L,Jx) + pbc_xx
end


# binary basis: \sigma_i^z \sigma_{i+1}^z
function Hzz(L::Int64, Jz::Float64)
    Hzz = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            Hzz[a+1,a+1] += Jz*bitsign(a,i)*bitsign(a,i+1)
        end
    end
    return Hzz
end

# binary basis: \sigma_i^z \sigma_{i+1}^z
function HzzSP(L::Int64, Jz::Float64)
    Hzz = spzeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-2)
            Hzz[a+1,a+1] += Jz*bitsign(a,i)*bitsign(a,i+1)
        end
    end
    return Hzz
end


# block basis: \sigma_i^z \sigma_{i+1}^z
function Hzz(L::Int64,Jz::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hzz = zeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        for i in 0:(L-2)
            Hzz[n,n] += Jz*bitsign(a,i)*bitsign(a,i+1)
        end
    end
    return Hzz
end

# block basis: \sigma_i^z \sigma_{i+1}^z
function HzzSP(L::Int64,Jz::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hzz = spzeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        for i in 0:(L-2)
            Hzz[n,n] += Jz*bitsign(a,i)*bitsign(a,i+1)
        end
    end
    return Hzz
end


# periodic boundary conditions, binary basis
function HzzP(L::Int64, Jz::Float64)
    pbc_zz = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        pbc_zz[a+1,a+1] += Jz*bitsign(a,L-1)*bitsign(a,0)
    end
    return Hzz(L,Jz) + pbc_zz
end


# binary basis: \sigma_i^z
function Hz(L::Int64, g::Float64)
    Hz = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-1)
            Hz[a+1,a+1] += g*bitsign(a,i)
        end
    end
    return Hz
end


# binary basis: \sigma_i^z
function sigZ(L::Int64, i::Int64)
    sigZ = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        sigZ[a+1,a+1] += bitsign(a,i-1)
    end
    return sigZ
end


# binary basis: \sigma_i^z
function HzSP(L::Int64, g::Float64)
    Hz = spzeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-1)
            Hz[a+1,a+1] += g*bitsign(a,i)
        end
    end
    return Hz
end

# block basis: \sigma_i^z
function Hz(L::Int64,mu::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hz = zeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        for i in 0:(L-1)
            Hz[n,n] += mu*bitsign(a,i)
        end
    end
    return Hz
end

# block basis: \sigma_i^z
function HzSP(L::Int64,mu::Float64,p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    Hz = spzeros(Complex{Float64},length(ref_states),length(ref_states))
    for (n,a) in enumerate(ref_states)
        for i in 0:(L-1)
            Hz[n,n] += mu*bitsign(a,i)
        end
    end
    return Hz
end


# binary basis: \sigma_i^z
function Hz_non_const(L::Int64, gs::Array{Float64})
    Hz = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        for i in 0:(L-1)
            Hz[a+1,a+1] += gs[i+1]*bitsign(a,i)
        end
    end
    return Hz
end


# binary basis: Hz + Hxx + Hzz
function H(L::Int64,Jx::Float64,Jz::Float64,g::Float64)
    hzz = Hzz(L,Jz)
    hxx = Hxx(L,Jx)
    hz = Hz(L,g)
    return hxx + hz + hzz
end

# binary basis: Hz + Hxx + Hzz
function H_SP(L::Int64,Jx::Float64,Jz::Float64,g::Float64)
    hzz = HzzSP(L,Jz)
    hxx = HxxSP(L,Jx)
    hz = HzSP(L,g)
    return hxx + hz + hzz
end


# binary basis: Hz + Hxx + Hzz
function H_Delta(L::Int64, Jx::Float64, Jz::Float64, g::Float64, gamma::Float64)
    hzz = Hzz(L,Jz)
    hxx = Hxx(L,Jx*(1+gamma)/2)
    hyy = Hyy(L,Jx*(1-gamma)/2)
    hz = Hz(L,g)
    return hxx + hyy + hz + hzz
end


# binary basis: Hz + Hxx + Hzz
function H_Delta_non_const(L::Int64, Jx::Float64, Jz::Float64, gs::Array{Float64}, gamma::Float64)
    hzz = Hzz(L,Jz)
    hxx = Hxx(L,Jx*(1+gamma)/2)
    hyy = Hyy(L,Jx*(1-gamma)/2)
    hz = Hz_non_const(L,gs)
    return hxx + hyy + hz + hzz
end


# binary basis: Hz + Hxx + Hzz
function H_Delta_SP(L::Int64, Jx::Float64, Jz::Float64, g::Float64, gamma::Float64)
    hzz = HzzSP(L,Jz)
    hxx = HxxSP(L,Jx*(1+gamma)/2)
    hyy = HyySP(L,Jx*(1-gamma)/2)
    hz = HzSP(L,g)
    return hxx + hyy + hz + hzz
end


# block basis: Hz + Hxx + Hzz
function H_Delta(L::Int64, Jx::Float64, Jz::Float64, g::Float64, gamma::Float64,
                 p::Int64, m::Int64)
    hzz = Hzz(L,Jz,p,m)
    hxx = Hxx(L,Jx*(1+gamma)/2,p,m)
    hyy = Hyy(L,Jx*(1-gamma)/2,p,m)
    hz = Hz(L,g,p,m)
    return hxx + hyy + hz + hzz
end

# block basis: Hz + Hxx + Hzz
function H_Delta_SP(L::Int64, Jx::Float64, Jz::Float64, g::Float64, gamma::Float64,
                 p::Int64, m::Int64)
    hzz = HzzSP(L,Jz,p,m)
    hxx = HxxSP(L,Jx*(1+gamma)/2,p,m)
    hyy = HyySP(L,Jx*(1-gamma)/2,p,m)
    hz = HzSP(L,g,p,m)
    return hxx + hyy + hz + hzz
end

# block basis: Hz + Hzz + Hxx
function H(L::Int64,Jx::Float64,Jz::Float64,g::Float64,p::Int64,m::Int64)
    hzz = Hzz(L,Jz,p,m)
    hxx = Hxx(L,Jx,p,m)
    hz = Hz(L,g,p,m)
    return hzz + hxx + hz
end


# periodic boundary conditions binary basis: Hz + Hxx + Hzz
function HP(L::Int64,Jx::Float64,Jz::Float64,g::Float64)
    hzz = HzzP(L,Jz)
    hxx = HxxP(L,Jx)
    hz = Hz(L,g)
    return hxx + hz + hzz
end

# binary basis, create \sigma_i^x
# Can only be in the binary basis as it
# does not commute with parity or mirror.
function sigX(L::Int64,i::Int64)
    sigX = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        b = xor(a,2^(i-1))
        sigX[a+1,b+1] += 1.0
    end
    return sigX
end

# binary basis, will apply sigX on a vector..
# should be faster than matrix multiplication
function sigX_V(L::Int64,site::Int64,vec::Array{Complex{Float64},1})
    ret = zeros(Complex{Float64},2^L)
    for j in 0:(2^(site -1)-1)
        # for k in 0:2^site:2^L-1
        #     println("j = $(j), k = $(k), $(bitstring(Int8(j + k))) gets $(bitstring(Int8(j + k + 2^(site-1))))")
        # end
        # display(collect(j:2^site:2^L) .|> Int8 .|> bitstring)
        # display(collect(j+2^(site-1):2^site:2^L) .|> Int8 .|> bitstring)
        ret[j+1            : 2^site : end] = vec[j+2^(site -1)+1 : 2^site : end]
        ret[j+2^(site-1)+1 : 2^site : end] = vec[j+1             : 2^site : end]
    end
    return ret
end


# binary basis, creates majorana at site i
# like sigX, can only be written in the binary
# basis.
function majorana(L::Int64,site::Int64)
    majorana = zeros(2^L,2^L)
    for a in 0:(2^L-1)
        tmp = 1
        for k in 0:(site-2)
            tmp *= bitsign(a,k)
        end
        b = xor(a,2^(site-1))
        majorana[a+1,b+1] += tmp
    end
    return majorana
end

# Binary basis: strong zero mode
function Psi0(L::Int64,g::Float64)
    psi0 = zeros(2^L,2^L)
    for j in 1:L
        psi0 += g^(j-1) .* majorana(L,j)
    end
    inv_norm = sqrt((1-g^2)/(1-g^(2*L)))
    psi0 = inv_norm .* psi0
    return psi0
end

# binary basis, calculate < Z_2 >
function Z2(L::Int64,state::Array{Complex{Float64},1})
    ans = 0.0
    for k in 0:2^L-1
        ans += (-1)^(bitsum(k,L))*abs(state[k+1])^2
    end
    return ans
end

# binary basis, calculate < I >
function mirror(L::Int64,state::Array{Complex{Float64},1})
    ans = 0.0
    for k in 0:2^L-1
        q = bitreverse(k,L)
        ans += conj(state[k+1])*state[q+1]
    end
    return ans
end


# takes block basis state into a binary basis vector
function block_to_binary(L::Int64,state::Array{Array{Complex{Float64},1},1})
    state_ans = zeros(Complex{Float64},2^L)
    for (p1i,p) in enumerate(0:1), (mi,m) in enumerate(0:1)
        ref_states = findBasis(L,p,m)
        index = 2*p + m + 1
        state_tmp = state[index]
        for (ni,n) in enumerate(ref_states)
            R = findPeriod(n,L)
            state_ans[n + 1] += state_tmp[ni]*sqrt(R)/2
            nr = bitreverse(n,L)
            state_ans[nr + 1] += (-1)^m*state_tmp[ni]*sqrt(R)/2
        end
    end
    return state_ans
end


# takes block basis state into a binary basis vector
# now takes specific block sectors
function block_to_binary(L::Int64,state::Array{Complex{Float64},1},p::Int64,m::Int64)
    state_ans = zeros(Complex{Float64},2^L)
    ref_states = findBasis(L,p,m)
    for (ni,n) in enumerate(ref_states)
        R = findPeriod(n,L)
        state_ans[n + 1] += state[ni]*sqrt(R)/2
        nr = bitreverse(n,L)
        state_ans[nr + 1] += (-1)^m*state[ni]*sqrt(R)/2
    end
    return state_ans
end



# takes binary basis vector into block basis
function binary_to_block(L::Int64,state::Array{Complex{Float64},1})
    overlaps :: Array{Array{Complex{Float64},1},1} = []
    for (p1i,p) in enumerate(0:1), (mi,m) in enumerate(0:1)
        ref_states = findBasis(L,p,m)
        overlaps_tmp = zeros(Complex{Float64},length(ref_states))
        for l in 0:(2^L-1)
            k1 = locateState(l,ref_states)
            k2 = locateState(bitreverse(l,L),ref_states)
            if k1>0
                N1 = 2/sqrt(findPeriod(ref_states[k1],L))
                overlaps_tmp[k1] += (1/N1*state[l+1])
            end
            if k2>0
                N2 = 2/sqrt(findPeriod(ref_states[k2],L))
                overlaps_tmp[k2] += ((-1)^m/N2*state[l+1])
            end
        end
        push!(overlaps,overlaps_tmp)
    end
    return overlaps
end

# takes binary basis vector into block basis, works for specific sector
function binary_to_block(L::Int64,state::Array{Complex{Float64},1},p::Int64,m::Int64)
    ref_states = findBasis(L,p,m)
    overlaps = zeros(Complex{Float64},length(ref_states))
    for l in 0:(2^L-1)
        k1 = locateState(l,ref_states)
        k2 = locateState(bitreverse(l,L),ref_states)
        if k1>0
            N1 = 2/sqrt(findPeriod(ref_states[k1],L))
            overlaps[k1] += (1/N1*state[l+1])
        end
        if k2>0
            N2 = 2/sqrt(findPeriod(ref_states[k2],L))
            overlaps[k2] += ((-1)^m/N2*state[l+1])
        end
    end
    return overlaps
end



#=
This will find eigenstate pairs via sigma_x
THIS IS SUPER SLOW... NOT USING BLAS,
NOT SURE WHY..

I NOW KNOW WHY!! IF YOU MULTIPLY A COMPLEX
MATRIX AGAINST A REAL MATRIX THE WHOLE THING
DEFAULTS TO SOME AWFUL SINGLE THREADED ROUTINE,
(MUCH WORSE PERFORMANCE THAN SINGLE THREADED
BLAS...)

YOU WOULD THINK IT WOULD PROMOTE THE FLOAT TO
COMPLEX, IT DOES NOT
=#
function find_pairs(L,v,site)
    pairs = zeros(Int64,2^L)
    sigx = sigX(L,site)
    println("before matrix mult")
    coeffs = abs.(v'*sigx*v).^2
    println("after matrix mult")
    for n in 1:2^L
        println("n = $n")
        pairs[n] = argmax(coeffs[:,n])
    end
    return pairs
end
