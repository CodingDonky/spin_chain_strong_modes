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

#To verify a command press "?" in the julia repl. Then it goes into Help mode.
#Type the command and julia gives all information about it

#=
Notes... 
There are functions called makeM which generate either 
the M matrix from the Fabian paper, or they take bns and produce 
the resulting Krylov Hamiltonian. Also, I often denote W matrices 
or HK or L matrices as M...

Much of the interacting data is either taken from the set below 
(Lanczos on HF) or its generated on the fly for small system sizes
so nothing in this file will take too long to produce.
=#



##################################################################################################
##################################################################################################
##################################################################################################

# toy data for interacting Floquet problem
# b_ns from Lanczos on Floquet Hamiltonian.
# system size is only 8 so not huge, but nice to have. Used in a couple spots.
# L = 8; jx = 1.0; g = .3; jz = 0.05; 
# Ts = 3/2 .*[.5,1.5,3.6,5.0,8.25]
Bs_5 = [[0.153668, 0.502557, 0.186668, 0.719168, 0.571393, 0.76389, 0.476687,
         1.25847, 0.769723, 1.15104, 1.13946, 1.37072, 1.04943, 1.40107,
         1.64001, 1.81812, 1.35673, 2.16312, 1.54627, 2.04627, 1.73094,
         2.22801, 1.69252, 1.65182, 1.58381, 1.72861, 1.68271, 1.81182,
         1.41862, 1.70555, 1.62641, 1.6823, 1.44449, 1.93471, 1.52069, 1.68394,
         1.57737, 1.6513, 1.88154, 1.86481, 2.02362, 1.87306, 1.84631, 1.77601,
         1.63427, 1.68665, 1.65522, 1.49915, 1.54126, 1.69039, 1.63431,
         1.58785, 1.70147, 1.53124, 1.66036, 1.57741, 1.59625, 1.67516, 1.599],
        [0.825275, 4.03712, 2.40067, 1.99139, 2.04133, 4.19992, 1.20796,
         2.89188, 3.24449, 2.94345, 2.01691, 4.04408, 2.43163, 3.49843,
         2.23785, 3.53926, 2.84441, 2.87791, 2.66057, 3.84921, 2.15255,
         3.30876, 2.65454, 3.24853, 2.7015, 3.0083, 2.97535, 3.30304, 2.5134,
         3.14043, 3.028, 3.02587, 3.03898, 3.09654, 2.99348, 3.03498, 2.90518,
         3.25604, 3.04237, 3.14561, 2.85129, 3.14447, 2.88055, 3.21971,
         2.88583, 3.10514, 2.76236, 3.06633, 2.81844, 3.12007, 2.84779,
         2.89981, 2.83939, 2.91201, 2.91757, 2.97014, 2.72452, 2.99674,
         2.98215],
        [1.65196, 2.79529, 1.17607, 4.02754, 2.02995, 3.92737, 3.09954, 3.3617,
         2.01674, 4.21213, 2.20123, 3.60392, 2.97296, 3.09578, 2.88263,
         3.65624, 2.78952, 3.1277, 3.22115, 2.85406, 3.04938, 3.47661, 2.72165,
         3.05168, 3.25492, 3.00395, 2.87637, 3.50783, 2.75012, 3.1458, 3.25741,
         2.94533, 2.99949, 3.24913, 2.96769, 2.9676, 3.31032, 3.01729, 2.95757,
         3.30517, 2.93738, 3.03126, 3.11529, 3.05235, 2.86891, 3.26741,
         2.97969, 3.07794, 3.13071, 3.04701, 2.96173, 3.14313, 2.93334,
         2.97225, 3.11706, 3.11186, 2.96609, 3.23135, 3.15586],
        [2.58818, 3.20078, 2.87247, 2.6967, 3.4334, 2.86531, 2.93879, 3.20209,
         3.03736, 2.88638, 3.31663, 2.99104, 2.94974, 3.18511, 3.08727,
         2.91986, 3.30999, 3.12533, 3.1162, 3.26369, 3.10212, 2.95698, 3.20072,
         3.06834, 3.07565, 3.20586, 3.06128, 2.87293, 3.03928, 3.09881,
         3.14008, 3.1332, 3.00609, 2.97769, 3.25052, 3.15786, 3.01757, 3.1181,
         3.19346, 3.08667, 3.14063, 3.10249, 3.01823, 3.10078, 3.17703,
         2.94533, 3.13593, 3.17021, 3.05802, 3.07668, 3.08843, 3.07137,
         3.17476, 3.01251, 3.13075, 3.04946, 3.10602, 3.13025, 3.21137],
        [2.99623, 1.48492, 3.35252, 3.54601, 2.12053, 2.76806, 4.22026,
         2.06471, 3.38667, 3.64851, 2.55608, 3.0502, 3.66245, 2.54851, 3.30759,
         3.3733, 2.95884, 3.13133, 3.2659, 2.88824, 3.11182, 3.17299, 3.06892,
         3.12376, 3.10426, 3.04621, 3.11319, 3.08492, 3.10696, 3.11056,
         3.09775, 3.10165, 3.08663, 3.05611, 3.08955, 3.10213, 3.11528, 3.1382,
         3.19319, 3.09454, 3.17018, 3.08758, 3.14018, 3.09987, 3.09252,
         3.06454, 3.05832, 3.11759, 3.11992, 3.12002, 3.0727, 3.11838, 3.1136,
         3.10883, 3.15252, 3.05341, 3.11574, 3.05304, 3.08169]]
##################################################################################################
##################################################################################################
#=
dot product for Arnoldi in operator space
L = system size
A = operator 1
B = operator 2
=#
function op_dot(L,A,B)
    return 1/(2^L) * tr(A'*B)
end

##################################################################################################
##################################################################################################
#=
commutation for Arnoldi in operator space
H = Hamiltonian spin space
O = operator spin space
=#
function op_H(H,O)
    return H*O - O*H
end

##################################################################################################
##################################################################################################
#=
conjugation for Arnoldi in operator space
U = propagator, spin space
O = operator, spin space
=#
function op_conj(U,O)
    return U'*O*U
end


##################################################################################################
##################################################################################################
#=
This shifts the energies into the FBZ..
engs = list of energies
=#
function shift_energies(engs)
    new_engs = mod.(engs .- pi,2*pi) .- pi
    return new_engs
end

##################################################################################################
##################################################################################################
#=
This function will take b_n and 
make a tridiagonal matrix. 
Fine tuned for free case.. so don't use 
for interacting...use ED.make_M(bs)
=#
function makeM(L,bs)
    tmp = zeros(2*L,2*L)
    for i in 1:length(bs)
        tmp[i,i+1] = bs[i]
        tmp[i+1,i] = bs[i]
    end
    # if Krylov subspace is < 2*L-1
    for i in length(bs)+2:2*L
        tmp[i,i] = 1
    end
    return tmp
end


##################################################################################################
##################################################################################################
#=
Function to generate A_inf. Not used in this file...
L = length of system
T = period of drive 
jx,jz,g = Hamiltonian parameters
periods = number of periods for A_inf to evaluate 
=#
function get_Ainf(L,T,jx,jz,g,periods)
    site = 1  
    sigx = ED.sigX(L,site) # sigma_1^x
    # sigy = im.*ED.sigZ(L,site)*ED.sigX(L,site)
    # sigz = ED.sigZ(L,site)
    
    # constructing U, eigen not necessary, just use exp... not critical
    u1s,v1s = eigen(ED.H_Delta(L,0.0,0.0,g,1.0))
    U1 = v1s*Diagonal(exp.(-im.*u1s.*T/3))*v1s'
    u2s,v2s = eigen(ED.H_Delta(L,jx,0.0,0.0,1.0))
    U2 = v2s*Diagonal(exp.(-im.*u2s.*T/3))*v2s'
    u3s,v3s = eigen(ED.H_Delta(L,0.0,jz,0.0,1.0))
    U3 = v3s*Diagonal(exp.(-im.*u3s.*T/3))*v3s'
    U = U3*U2*U1

    ainfs = zeros(periods)
    # messing around with operators
    # op1 = (sigx + exp(im*pi/6)*sigy)/sqrt(2) 
    # op2 = (sigx + exp(im*pi/6)*sigy)/sqrt(2)
    # op1, op2 if you want to generalize A_inf...
    op1 = sigx
    op2 = sigx
    for ti in 1:periods
        print("ti = $ti\r")
        ainfs[ti] = real(1/(2^L)*tr(op2'*op1))
        op2 = U'*op2*U
    end
    fig,ax = plt.subplots(1,1)
    ax.plot(1:periods,ainfs)
    ax.set_xscale("log")
    plt.show()
end

##################################################################################################
##################################################################################################
#=
Alternate function for generating A_inf
L=  system size
U = propagator (spin basis)
times = times to evaluate
=#
function gen_Ainf(L,U,times)
    ainf = zeros(length(times))
    sigx = ED.sigX(L,1)
    quasi_phases,quasi_vecs = eigen(U)
    sigxv = quasi_vecs'*sigx*quasi_vecs
    for (ti,time) in enumerate(times)
      phases = Diagonal(quasi_phases .^ time)
      ainf[ti] = 1/(2^L)*real(tr(phases'*sigxv*phases*sigxv))
    end
    return ainf
end

##################################################################################################
##################################################################################################
#= 
Returns M matrix from Fabian paper.
L = length of system
T = period of drive
g = transverse magnetic field
=#
function get_M(L,T,g)
    a = cos(T*g)*cos(T)
    b = sin(T*g)*sin(T)
    c = -sin(T*g)*cos(T) 
    d = -cos(T*g)*sin(T)
    ap = cos(T*g)
    cp = -sin(T*g)
    M = zeros(2*L,2*L)
    for i in 2:L
        M[i,i] = a
        M[i,i-1] = b
        M[i,i+L] = c
        M[i,i+L-1]=-d
    end
    for i in 1:L-1
        M[i+L,i+L] = a
        M[i+L,i+L+1] = b
        M[i+L,i] = -c
        M[i+L,i+1] =  d
    end
    M[1,1] = ap
    M[1,L+1] = cp
    M[2*L,L] = -cp
    M[2*L,2*L] = ap

    return M
end

##################################################################################################
##################################################################################################
#=
This will construct the W matrix and basis "vectors" (operators)
in ED, thus valid for interactions.
L = length of system
N = number of iterations
U = propagator (spin basis)
op = seed operator (usually sigma_1^x)
threshold = minimum value of b_n, stopping condition (only need for j_z very close to zero)
=#
function arnoldi_ED(L,N,U,op,threshold=1e-10)
    vecs = [] # store basis vectors here
    push!(vecs,op)
    W = zeros(Complex{Float64},N,N) # define empty W matrix
    # Standard Arnoldi in ED, nothing too fancy.
    for n in 1:N
        tmp = op_conj(U,vecs[end]) # apply U
        for j in 1:length(vecs)
            overlap = op_dot(L,vecs[j],tmp) # find overlaps with previous basis vectors
            W[j,n] = overlap # store overlap values
            tmp = tmp - overlap.*vecs[j] # Gram-Schmidt orthogonalize w.r.t previous basis vectors
        end
        tmp_h = real(sqrt(op_dot(L,tmp,tmp))) # find norm of resulting state
        if abs(tmp_h) > threshold && n <N # if not at end
            W[n+1,n] = tmp_h # store this norm as sub-diagonal entry (simple to show by hand)
            push!(vecs, tmp ./ tmp_h) # finally normalize basis vector and store
        else
            # if we are at the end, lop off what we calculated and break out
            # this could probably be written better..
            W = W[1:n,1:n] 
            break
        end
    end
    return W,vecs
end

##################################################################################################
##################################################################################################
#=
This produces the approximate edge modes of W. 
Just solve for the LEFT eigenvector of W,
the LEFT null vector of (W -/+ 1)
and start with psi_1 = 1. The upper Hessenberg 
form of W will immediately lead to this result. 
One way this is an approximation is that we plug 
in for the eigenvalue, rather than solving for it.

W = W matrix for solving for edge modes
lambda = +/- 1, which type of edge mode you want to solve for

=#
function gen_wavefunction2(W,lambda)
    L = size(W)[1]
    t = zeros(L); t[1] = 1
    for k in 1:L-1
        tmp = (W[k,k] - lambda)*t[k]
        for j in 1:k-1
            tmp += W[k-j,k]*t[k-j]
        end
        t[k+1] = -tmp/W[k+1,k]
    end
    t = t ./ sqrt(t'*t)
    return t
end

##################################################################################################
##################################################################################################
#= 
This will perform Arnoldi on a Matrix in the free 
Majorana basis so no ED is required here. 
Main use is for checking analytic results.
Same Arnoldi algorith as before.

M = free single particle propagator, typically M matrix from Fabian paper
L = length of system
=#
function gen_wmat(M,L)
    vec = zeros(2*L)
    vec[1] = 1
    vecs = [vec]
    N = zeros(2*L,2*L)
    N[:,1] = vecs[end]
    for n in 2:2L
        tmp = M*vecs[end]
        for k in 1:length(vecs)
            tmp = tmp - (tmp'*vecs[k]).*vecs[k]
        end
        tmp = tmp./sqrt(tmp'*tmp)
        push!(vecs,tmp)
        N[:,n] = tmp
    end
    M2 = N'*M*N
    return M2,N
end



##################################################################################################
##################################################################################################
#= 
This function will produce b_n, given 
the starting "operator" (vector, usually (1,0,0,...))
and given an matrix M to tridiagonalize.
Similar to above, but this is Lanczos not Arnoldi
This does not require ED resources. 

M = matrix to tridiagonalize
op = starting vector
threshold = ending condition
=#
function lanczos_steps2(M,op,threshold = 1e-10)
    N = size(M)[1]
    Ops = [op]
    tmp = M*Ops[end]
    bs = [sqrt(tmp'*tmp)]
    tmp = tmp ./ bs[end]
    push!(Ops,tmp)
    for i in 3:N
        tmp1 = M*Ops[end] - bs[end].*Ops[end-1] 
        for k in 1:length(Ops)
            tmp1 -= Ops[k] .*(Ops[k]'*tmp1)
        end
        tmp2 = sqrt(tmp1'*tmp1)
        if abs(tmp2)>threshold
            push!(Ops,tmp1 ./ tmp2)
            push!(bs,tmp2)
        else 
            break
        end
    end
    return Ops,bs
end

##################################################################################################
##################################################################################################
#=
Similar to lanczos_steps2, just written slightly different...
L = system size 
H_F = floquet Hamiltonian
op = seed operator
N = number of states to calculate
threshold = end condition
=#
function lanczos_steps(L,H_F,op,N,threshold = 1e-10)
    ops = [op .* (1.0 + 0.0*im)]
    tmp = op_H(H_F,ops[end])
    bs = [abs(sqrt(op_dot(L,tmp,tmp)))]
    tmp = tmp ./ bs[end]
    push!(ops,tmp)
    for i in 3:N
        tmp1 = op_H(H_F,ops[end]) - bs[end] .* ops[end-1]
        for k in 1:length(ops)
            tmp1 -= ops[k] .* (op_dot(L,ops[k],tmp1))
        end
        tmp2 = abs(sqrt(op_dot(L,tmp1,tmp1)))
        if tmp2 > threshold
            push!(ops,tmp1 ./ tmp2)
            push!(bs,tmp2)
        else 
            break
        end
    end
    return ops,bs
end


##################################################################################################
##################################################################################################
#= 
This function takes the M matrix from the 
Fabian paper and moves it into the basis 
where everything is on the main diagonal
(just a reordering of the basis).
M = operator to change
=#
function change_basis(M)
    L = div(size(M)[1],2)
    ret = zeros(2*L,2*L)
    for i in 1:L
        ret[:,2*i-1] = M[:,i]
        ret[:,2*i] = M[:,i+L]
    end
    Mcopy = copy(ret)
    for j in 1:L
        ret[2*j-1,:] = Mcopy[j,:]
        ret[2*j,:] = Mcopy[j+L,:]
    end
    return ret
end

##################################################################################################
##################################################################################################
#= 
This produces the W matrix for the free case 
with the results determined from sympy/by hand.
L = system size
T,g = drive parameters
=#
function free_W(L,T,g)
    c1 = cos(T*g)
    c2 = cos(T)
    s1 = sin(T*g)
    s2 = sin(T)
    # will return the matrix ret
    ret = zeros(2*L,2*L)
    # define row1 and row2, will use to fill in ret
    row1 = zeros(2*L+1)
    row2 = zeros(2*L+1)
    row1[1] = s1
    row1[2] = c1*c2
    row1[3] = -s2*c1^2
    row2[1] = s2
    row2[2] = c1*c2
    row2[3] = -s1*c2^2
    for k in 4:2:2L
        row1[k] = row1[k-2]*(s1*s2)
        row1[k+1] = row1[k-1]*(s1*s2)
        
        row2[k] = row2[k-2]*(s1*s2)
        row2[k+1] = row2[k-1]*(s1*s2)
    end
    # now filling in ret
    ret[1,:] = row2[2:end]
    for k in 2:2:2*L
        ret[k,k-1:end] = row1[1:end-(k-1)]
        if k < 2*L
            ret[k+1,k:end] = row2[1:end-(k)]
        end
    end
    ret[1,:] = ret[1,:]/c2
    ret[:,end] = ret[:,end]/c2
    return ret
end

##################################################################################################
##################################################################################################
#= 
This will construct the W matrix in the free limit, 
majorana basis. It compares W constructed by the Arnoldi 
method defined above to the analytical answer found through 
sympy/by hand. It also constructs a simple approximation to 
the edge mode. The edge mode form is compared to the exact 
(for both zero and pi). Also the edge modes are time evolved 
and we see how it compares to A_inf(n) = W_{1,1}^n.
=#


function free_W_edge_modes()
    g = .3
    Tps = [2.0,3.5,5.5,8.25] # T = 3/2*Tp
    L = 10
    fig,ax = plt.subplots(3,length(Tps))
    for (ti,Tp) in enumerate(Tps)
        M = get_M(L,Tp,g) # M from Fabian paper
        display(M) 
        R = change_basis(M) # reorder basis vectors for M
        display(R)
        M2,N = gen_wmat(R,L) # find W via Arnoldi
        display(M2)
        display(free_W(L,Tp,g)) # analytical answer for W
        println("difference between computer and analytic")
        println(maximum(abs.(abs.(M2)-abs.(free_W(L,Tp,g)))))

        # spectrum of M2 
        u,v = eigen(M2)
        # potential edge states
        v1 = v[:,end]
        v2 = v[:,end-1]
        v3 = v[:,1]
        v4 = v[:,2]
        # approximate edge states
        t3 = gen_wavefunction2(M2,1)
        t4 = gen_wavefunction2(M2,-1)
        # plotting
        ax[1,ti].plot(imag.(log.(u)),".")
        ax[2,ti].plot(abs.(v1).^2,"-^",lw = .7,ms = 2,label="zero 1")
        ax[2,ti].plot(abs.(v2).^2,"-+",lw = .7,label="zero 2")
        ax[2,ti].plot(abs.(v3).^2,"--",lw = .7,label="pi 1")
        ax[2,ti].plot(abs.(v4).^2,":",lw = .7, label="pi 2")
        ax[2,ti].plot(t3.^2,"-o",lw = .7,ms = 2,label="+1")
        ax[2,ti].plot(t4.^2,"-s",lw = .7,ms = 2,label="-1")
        ax[2,ti].set_yscale("log")
        ax[2,ti].legend()
        # setting up A_inf calculation
        ainf1 = [1.0]
        ainf2 = [t3[1]]
        ainf3 = [t4[1]]
        vec = zeros(2*L); vec[1] = 1
        tend = 100
        # time evolving
        for ti in 1:tend
            vec = M2*vec
            t3 = M2*t3
            t4 = M2*t4
            push!(ainf1,vec[1])
            push!(ainf2,t3[1])
            push!(ainf3,t4[1])
        end
        # more plotting
        ax[3,ti].plot(0:tend,ainf1,"-o",alpha = .3,label = "\$ W_{1,1}^n \$")
        ax[3,ti].plot(0:tend,ainf2,".",label = "\$\\langle 1 | W^n|\\psi_{+} \\rangle\$")
        ax[3,ti].plot(0:tend,ainf3,".",label = "\$\\langle 1 | W^n|\\psi_{-} \\rangle\$")
        ax[3,ti].set_xlabel("time")
        ax[3,ti].legend()
    end
    fig.show()
end

##################################################################################################
##################################################################################################
#= 
Extension of previous method, checks overlap of 
<\psi|W|\psi> where \psi is our approximate edge 
state. Calculation is performed for many values of 
T.
=#
function free_W_edge_modes2()
    g = .3
    Tps = range(1,10,length = 1500)
    L = 10 # system size

    # storing results here..
    data1 = zeros(2*L,length(Tps)) # quasienergy
    data2 = zeros(2*L,length(Tps)) # not used....
    data3 = zeros(2*L,length(Tps)) # abs overlap zero mode
    data4 = zeros(2*L,length(Tps)) # abs overlap pi mode
    data5 = zeros(2*L,length(Tps)) # real overlap zero mode
    data6 = zeros(2*L,length(Tps)) # real overlap pi mode
    # iterate through Ts
    for (ti,Tp) in enumerate(Tps)
        M = get_M(L,Tp,g)
        M2,N = gen_wmat(M,L)
        u,v = eigen(M)
        data1[:,ti] = imag.(log.(u))

        # approximate solutions
        psi1 = gen_wavefunction2(M2,1) # possible zero mode
        psi2 = gen_wavefunction2(M2,-1) # possible pi mode 

        # cut off approximate solutions at different lengths (k)
        # and see how they compare
        for k in 1:2*L
            psi_tmp = zeros(2*L)
            psi_tmp[1:k] = psi1[1:k]
            psi_tmp = psi_tmp ./ sqrt(psi_tmp'*psi_tmp) # normalizing
            tmp = psi_tmp'*M2*psi_tmp # finding overlaps
            data3[k,ti] = abs(tmp) 
            data5[k,ti] = real(tmp)
            # now repeat for pi mode
            psi_tmp = zeros(2*L)
            psi_tmp[1:k] = psi2[1:k]
            psi_tmp = psi_tmp ./ sqrt(psi_tmp'*psi_tmp)
            tmp = psi_tmp'*M2*psi_tmp
            data4[k,ti] = abs(tmp)
            data6[k,ti] = real(tmp)
        end
    end
    
    # plotting
    fig,ax = plt.subplots(5,sharex = true)
    fig.subplots_adjust(hspace = 0.0)
    ax[1].plot(Tps,data1',".",color = "black",alpha = .4,ms = 1.0)
    for k in 1:2:2*L
        ax[2].plot(Tps,data3[k,:],".",label = "$(k)",ms = 2.0)
        ax[3].plot(Tps,data5[k,:],".",label = "$(k)",ms = 2.0)
        ax[4].plot(Tps,data4[k,:],".",label = "$(k)",ms = 2.0)
        ax[5].plot(Tps,data6[k,:],".",label = "$(k)",ms = 2.0)
    end
    ax[2].legend(ncol = 2)
    ax[3].legend(ncol = 2)
    ax[4].legend(ncol = 2)
    ax[5].legend(ncol = 2)
    ax[5].set_xlabel("T")
    ax[1].set_ylabel("quasienergy of W")
    ax[2].set_ylabel("\$|\\langle\\psi_+|W|\\psi_+\\rangle|^2\$")
    ax[3].set_ylabel("\$re(\\langle\\psi_+|W|\\psi_+\\rangle)\$")
    ax[4].set_ylabel("\$|\\langle\\psi_-|W|\\psi_-\\rangle|^2\$")
    ax[5].set_ylabel("\$re(\\langle\\psi_-|W|\\psi_-\\rangle)\$")
    fig.show()
end
##################################################################################################
##################################################################################################
#=
This produces exact Ainf along side time evolution of approximate 
edge modes.
=#
function interacting_W_edge_modes()
    Ls = [4,6] # different lengths
    # Hamiltonian parameters
    jx = 1.0
    g = .3
    jz = 0.05
    # different periods
    Tps = [2.0,3.5,5.5,8.25]
    Ts = 3/2 .* Tps
    # length of W matrix/ number of basis vectors
    N = 20
    Ns = [6,12,18] # different cut off lengths for approximate edge states (must be <= N)
    t_len = 5000 # length of ainf, in periods
    # store results here
    data_ainf1 = zeros(t_len,length(Ts),length(Ls)) # exact ainf
    data_ainf2 = zeros(t_len,length(Ts),length(Ls),length(Ns)) # time evolution for approximate zero mode
    data_ainf3 = zeros(t_len,length(Ts),length(Ls),length(Ns)) # time evoltuion for approximate pi mode
    data_decay1 = zeros(length(Ts),length(Ls),length(Ns)) # approx decay constants, see notes, zero mode
    data_decay2 = zeros(length(Ts),length(Ls),length(Ns)) # approx decay constants, pi mode
    data_decay1b = zeros(length(Ts),length(Ls),length(Ns)) # second decay constant, zero mode
    data_decay2b = zeros(length(Ts),length(Ls),length(Ns)) # second decay constant, pi mode
    for (Li,L) in enumerate(Ls) # step through system sizes
        sigx = ED.sigX(L,1) # sigma_1^x
        for (Ti,T) in enumerate(Ts) # step through different periods
            println(L,T)
            # propagator
            U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,  g,1.0))
            U2 = exp(-im * T/3 .* ED.H_Delta(L, jx,0.0,0.0,1.0))
            U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0, jz,0.0,1.0))
            U = U3*U2*U1
            # calculate ainf
            u,v = eigen(U)
            for ti in 1:t_len
                U_tmp = v*Diagonal(u .^ ti)*v'
                data_ainf1[ti,Ti,Li] = real(tr(sigx*U_tmp*sigx*(U_tmp'))/(2^L))
            end
            # step through different cut offs for approximate edge state of W
            for (Ni,N) in enumerate(Ns)
                W,qs = arnoldi_ED(L,N,U,sigx,1e-15) # find W 
                W = real(W) # not really needed..
                psi_p = gen_wavefunction2(W, 1) # zero approx 
                psi_m = gen_wavefunction2(W,-1) # pi approx
                data_decay1[Ti,Li,Ni] = -log(abs(psi_p'*W*psi_p)) # decay constant zero mode
                data_decay2[Ti,Li,Ni] = -log(abs(psi_m'*W*psi_m)) # decay constant pi mode
                data_decay1b[Ti,Li,Ni] = abs(psi_p[end])^2 # second decay constant zero  mode
                data_decay2b[Ti,Li,Ni] = abs(psi_m[end])^2 # second decay constant pi mode
                # ainfs for approximate edge modes
                for ti in 1:t_len
                    psi_p = W*psi_p
                    psi_m = W*psi_m
                    data_ainf2[ti,Ti,Li,Ni] = real(psi_p[1])
                    data_ainf3[ti,Ti,Li,Ni] = real(psi_m[1])
                end
            end
        end
    end
    # plotting
    fig,ax = plt.subplots(length(Ls),length(Ts),
                          sharey = true,sharex = true)
    fig.subplots_adjust(wspace = 0.0,hspace = 0.0)
    for (Li,L) in enumerate(Ls)
        for (Ti,T) in enumerate(Ts)
            ax[Li,Ti].plot(1:t_len,data_ainf1[:,Ti,Li],".",label = "exact")
            for (Ni,N) in enumerate(Ns)
                tmp_plt = data_ainf2[1,Ti,Li,Ni]*data_ainf2[:,Ti,Li,Ni] +
                    data_ainf3[1,Ti,Li,Ni]*data_ainf3[:,Ti,Li,Ni]
                tmp_plt2 = abs(tmp_plt[1]) .* exp.(-data_decay1[Ti,Li,Ni].*(1:t_len))
                tmp_plt3 = abs(tmp_plt[1]) .* exp.(-data_decay2[Ti,Li,Ni].*(1:t_len))
                
                tmp_plt2b = abs(tmp_plt[1]) .* exp.(-data_decay1b[Ti,Li,Ni].*(1:t_len))
                tmp_plt3b = abs(tmp_plt[1]) .* exp.(-data_decay2b[Ti,Li,Ni].*(1:t_len))
                
                ax[Li,Ti].plot(1:t_len,tmp_plt,
                               ".",
                               ms = 2.0,
                               label = "approx, N = $(N)")
                if Ni == 3
                    ax[Li,Ti].plot(1:t_len,tmp_plt2,
                                   "-|",
                                   lw = 2.0,
                                   label = "approx decay ASZM, N = $(N)")
                    ax[Li,Ti].plot(1:t_len,tmp_plt3,
                                   "-^",
                                   lw = 2.0,
                                   label = "approx decay ASPM, N = $(N)")

                    ax[Li,Ti].plot(1:t_len,tmp_plt2b,
                                   "-",
                                   lw = 2.0,
                                   label = "approx2 decay ASZM, N = $(N)")
                    ax[Li,Ti].plot(1:t_len,tmp_plt3b,
                                   "-",
                                   lw = 2.0,
                                   label = "approx2 decay ASPM, N = $(N)")
                                  
                end
                ax[Li,Ti].set_xlabel("\$t\$")
                ax[Li,Ti].set_xscale("log")
            end
        end
        ax[Li,1].set_ylabel("\$A_{\\infty}(t) \$")
    end
    ax[1,1].legend()
    fig.show()
end

##################################################################################################
##################################################################################################
#=
All free and in majorana basis.
Displays many plots of W matrix, log(W), M matrix, log(M),
simplifications of log(W), and the resulting reexponentiation 
of these simplified log(W), and their error with W.
=#
function free_W_overview()
    g = .3 # mag field
    Tps = [2.0,3.5,5.5,8.25] # T = 3/2*Tp # periods 
    L = 6 # system size
    # for plotting
    fig,ax = plt.subplots(7,length(Tps),
                          sharex = true,
                          sharey = true)
    fig.subplots_adjust(hspace = 0.0)
    # step through different periods
    for (ti,Tp) in enumerate(Tps)
        println("T = $(Tp) ##################################################")
        M = get_M(L,Tp,g) # get M from Fabian paper
        # get spectrum and display it in terminal
        u,v = eigen(M)
        display(u)
        display(abs.(v.^2))
        println("M")
        display(round.(M,sigdigits = 2))
        # plot M
        im1 = ax[1,ti].imshow(M)
        plt.colorbar(im1,ax = ax[1,ti])
        # get W from M, called M2.....
        M2,N = gen_wmat(M,L)
        # display spectrum in terminal
        u,v = eigen(M2)
        display(u)
        display(abs.(v.^2))
        println("W")
        display(round.(M2,sigdigits = 2))
        # plot W
        im2 = ax[2,ti].imshow(M2)
        plt.colorbar(im2,ax = ax[2,ti])
        println("Hm = log(M)")
        # log of M
        H1 = real.(log(M))
        display(round.(H1,sigdigits = 2))
        im3 = ax[3,ti].imshow(H1)
        # plot log of M
        plt.colorbar(im3,ax = ax[3,ti])
        # log of W (M2...)
        H2 = real.(log(M2))
        println("Hw = log(W)")
        display(round.(H2,sigdigits = 2))
        # plot log of W (M2...)
        im4 = ax[4,ti].imshow(H2)
        plt.colorbar(im4,ax = ax[4,ti])

        # simplify log of W, cutting off long range matrix elements
        # except for corners, where pi modes like to have large values
        H2b = deepcopy(H2)
        bands = 1
        for k in 1:2*L
            for j in 1:2*L
                if abs(k-j)>bands && !((k >= 2*L-2 && j <= 2) || (k <= 2 && j >= 2*L-2))
                    H2b[k,j] = 0
                end
            end
        end
        println("Hw'")
        display(round.(H2b,sigdigits = 2))
        # plot simplified log of W
        im5 = ax[5,ti].imshow(H2b)
        plt.colorbar(im5,ax = ax[5,ti])
        # reexponentiate this simplified matrix
        Wb = exp(H2b)
        println("W'")
        display(round.(Wb,sigdigits = 2))
        # plot it
        im6 = ax[6,ti].imshow(Wb)
        plt.colorbar(im6,ax = ax[6,ti])
        # now we compare it with the exact W
        println("|W - W'|")
        display(round.(abs.(Wb - M2),sigdigits = 2))
        im7 =ax[7,ti].imshow(abs.(Wb - M2))
        plt.colorbar(im7,ax = ax[7,ti])
        ax[end,ti].set_xlabel("T = $(Tp)")        
    end
    ax[1,1].set_ylabel("M")
    ax[2,1].set_ylabel("W")
    ax[3,1].set_ylabel("iHm = log(M)")
    ax[4,1].set_ylabel("iHw = log(W)")
    ax[5,1].set_ylabel("iHw'")
    ax[6,1].set_ylabel("W'")
    ax[7,1].set_ylabel("|W - W'|")
end

##################################################################################################
##################################################################################################
#=
Same as free_W_overview, but for interacting, 
We step through different J_z values and compare 
W with log(W) and see how the matrix elements change.
=#
function interacting_W_overview()
    L = 6 # system size
    jx = 1.0
    jzs = [0.0,.001,.01,.1]
    g = .3
    N = 4*L+1 # number of basis vectors to calculate
    Tps = [2.0,3.5,5.5,8.25]; Ts = 3/2 .*Tps
    sigx = ED.sigX(L,1)
    fig,ax = plt.subplots(2*length(Tps),length(Tps),
                          sharex = true,
                          sharey = true)
    fig.subplots_adjust(hspace = 0.0)
    # step through different Ts
    for (ti,T) in enumerate(Ts)
        println("T = $(T) ##################################################")
        for (jzi,jz) in enumerate(jzs) # run through different Jz
            # construct propagator.
            U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,  g,1.0))
            U2 = exp(-im * T/3 .* ED.H_Delta(L, jx,0.0,0.0,1.0))
            U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0, jz,0.0,1.0))
            U = U3*U2*U1
            u,v = eigen(U)
            # produce W, ED
            W,qs = arnoldi_ED(L,N,U,sigx,1e-15)
            W = real(W)
            println("W")
            println("here ----------------------------")
            println("W2")
            # checking trick to make W unitary by messing with last column
            # and matrix element N+1,N.... not critical
            W2 = W[1:N-1,1:N-1]
            println(det(W2))
            display(round.(W2,sigdigits = 2))
            W3 = copy(W2)
            W3[:,end] = W3[:,end] ./sqrt(1 - W[N,N-1]^2)
            display(round.(W3,sigdigits = 2))
            println(det(W3))
            W= copy(W3)
            println("here ----------------------------")
            # plotting matrix elements of W
            im1 = ax[2*jzi-1,ti].imshow(W)
            plt.colorbar(im1,ax = ax[2*jzi-1,ti])
            Hw = real.(log(W))
            # plotting log of W
            im2 = ax[2*jzi,ti].imshow(Hw)
            plt.colorbar(im2,ax = ax[2*jzi,ti])
            if ti == 1
                ax[2*jzi-1,1].set_ylabel("W($(jz))",fontsize = 5)
                ax[2*jzi,1].set_ylabel("Hw($(jz))",fontsize = 5)
            end
        end
        ax[end,ti].set_xlabel("T' = $(2/3*T)")
    end
    plt.show()
end

##################################################################################################
##################################################################################################
#=
This will produce plots for the b_n 
resulting from Lanczos on the Floquet Hamiltonian 
in the free case in the Majorana basis. 
For the propagator we use the "M" matrix from the 
Fabian paper.
=#
function floquet_ham_free_lanczos()
    L = 20
    g = .3
    jx = 1.0
#    Ts = [.5,1.5,3.6,5.,8.25]
    Ts = [2.0,3.5,5.5,8.25]

    data_engs = zeros(length(Ts),2*L)
    data_bs = zeros(length(Ts),2*L-1)

    for (Ti,T) in enumerate(Ts)
        tmpM = get_M(L,T,g) # Fabian matrix
        M = -im*log(tmpM) # antisymmetric, pure imaginary
        ops = zeros(Complex{Float64},size(M)[1]); ops[1] = 1.0
        ops,bs = lanczos_steps2(M,ops,1e-16)
        # build Krylov Hamiltonian and diagonalize, get quasienergy
        HK = makeM(L,bs) # plop bns into tri-diagonal matrix
        u,v = eigen(HK) # diagonalize
        data_engs[Ti,:] = u
        data_bs[Ti,:] = bs
    end
    # rest is plotting
    fig,ax = plt.subplots(2,length(Ts),sharey = "row")
    fig.subplots_adjust(wspace = 0.0)
    for (Ti,T) in enumerate(Ts)
        ax[1,Ti].plot(data_bs[Ti,:],"-o")
        ax[2,Ti].plot(data_engs[Ti,:],"-o")
        ax[2,Ti].axhline(y = -pi,color = "red",alpha = .5)
        ax[2,Ti].axhline(y =  pi,color = "red",alpha = .5)
        ax[2,Ti].axhline(y =  0 ,color = "red",alpha = .5)
        ax[1,Ti].set_xlabel("n")
        ax[2,Ti].set_xlabel("index")
    end
    ax[1,1].set_ylabel("\$ b_n \$")
    ax[2,1].set_ylabel("energy")
    fig.show()
end


##################################################################################################
##################################################################################################
#=
This will produce the b_n and H_K (spectrum) for the interacting 
drive. Will do the double Lanczos thing where we first do Lanczos 
on H_F, plug the b_n into H_K, alter the spectrum so that it is 
within the FBZ, then do Lanczos on the resulting new H_K',  and then 
plot those new b_n....

=#
function floquet_ham_interacting_lanczos()
    L = 60 # number of b_n
    # parameters of Hamiltonian
    g = .3
    jx = 1.0
    Ts = [.5,1.5,3.6,5.,8.25]

    # save quasienergy and b_ns
    data_engs = zeros(length(Ts),2,L)
    data_bs = zeros(length(Ts),2,L-1)
    # step through periods
    for (Ti,T) in enumerate(Ts)
        bs = Bs_5[Ti] # using the data from the beginning of this file
        HK = ED.make_M(bs) # put bns into tri-diagonal matrix
        u,v = eigen(HK) # diagonalize
        data_engs[Ti,1,:] = u
        data_bs[Ti,1,:] = bs
        # shifting energies of quasienergies in Krylov Hamiltonian
        # somewhat sketchy, as we terminated Lanczos algo early,
        # so these eigenvalues(vectors) are not the real values
        HK = v*Diagonal(shift_energies(u))*v'; u,v = eigen(HK)
        data_engs[Ti,2,:] = u
        # now we perform a second Lanczos on the new HK
        # and find the new b_ns and plot/save those.
        op = zeros(length(u)); op[1] = 1.0
        ops,bs_new = lanczos_steps2(HK,op)
        data_bs[Ti,2,:] = bs_new
    end
    # rest is plotting
    fig,ax = plt.subplots(2,length(Ts),sharey = "row")
    fig.subplots_adjust(wspace = 0.0)
    for (Ti,T) in enumerate(Ts)
        ax[1,Ti].plot(data_bs[Ti,1,:],"-o",ms = 4.)
        ax[1,Ti].plot(data_bs[Ti,2,:],"-x",ms = 4)
        ax[2,Ti].plot(data_engs[Ti,1,:],"-o",ms = 4)
        ax[2,Ti].plot(data_engs[Ti,2,:],"-x",ms = 4)
        ax[2,Ti].axhline(y = -pi,color = "red",alpha = .5)
        ax[2,Ti].axhline(y =  pi,color = "red",alpha = .5)
        ax[2,Ti].axhline(y =  0 ,color = "red",alpha = .5)
        ax[1,Ti].set_xlabel("n")
        ax[2,Ti].set_xlabel("index")
    end
    ax[1,1].set_ylabel("\$ b_n \$")
    ax[2,1].set_ylabel("energy")
    fig.show()
end

##################################################################################################
##################################################################################################
#=
Helper function for pseudospectrum calculation. 
See Dillon's paper.
=#
function pseudospectrumB(A,B,la,lb,eta)
    L = size(A)[1]
    b = zeros(Complex{Float64},2*L,2*L)
    Xa = eta.*(A - Diagonal([la for i in 1:L]))
    Xb = B - Diagonal([lb for i in 1:L])
    b[1:L,L+1:end] = Xa - im .* Xb
    b[L+1:end,1:L] = Xa + im .* Xb
    return b
end

##################################################################################################
##################################################################################################
#=
Interacting system (using Bs_5 from above). 
Perform double Lanczos proceedure and then 
uses the resulting b_n and calculates psuedospectrum
=#

function testing_pseudospectrum()
    # system parameters
    L = 8; jx = 1.0; g = .3; jz = 0.05;
    # different periods
#    Ts = 3/2 .*[.5,1.5,3.6,5.0,8.25]
      Ts = 3/2 .*[2.0,3.5,5.5,8.25]
    # parameter for pseudospectrum
    epsilon = .05
    # positions for pseudospectrum
    pos = range(0,1,length = 150)
    # nu is another parameter for pseudospectrum
    nus = 10 .^ range(-3,0,length = 5)
    gaps_pi = []
    gaps_0  = []
    Bs = []
    # step through periods
    for (Ti,T) in enumerate(Ts)
        println(Ti)
        bs1 = Bs_5[Ti] # get bns
        # make HK
        HK1 = ED.make_M(bs1)
        u1,v1 =  eigen(HK1)
        # shift energies
        u1b = shift_energies(u1)
        HK1b = v1*Diagonal(u1b)*(v1')
        # do a second Lanczos
        opB = zeros(Complex{Float64},size(HK1b)[1]); opB[1] = 1.0
        # get new b_ns
        opsB,bs1b = lanczos_steps2(HK1b,opB,1e-16)
        push!(Bs,bs1b)

        # set up pseudospectrum calculation
        M = ED.make_M(bs1b)
        L2 = size(M)[1]
        X = Diagonal(range(0,L2-1,length = L2))
        X = X ./ L2
        println(maximum(abs.(M*X - X*M).^2))
        
        data_pi = zeros(length(pos),length(nus))
        data_0 = zeros(length(pos),length(nus))
        # step through nu values and caluculate pseudospectrum
        for (nui,nu) in enumerate(nus), (xi,x) in enumerate(pos)
            b = pseudospectrumB(X,M,x,pi*1.0,nu)
            data_pi[xi,nui] = minimum(abs.(eigvals(b)))
            b = pseudospectrumB(X,M,x,0,nu)
            data_0[xi,nui] = minimum(abs.(eigvals(b)))
        end
        push!(gaps_0,data_0)
        push!(gaps_pi,data_pi)
    end
    # rest is plotting
    fig,ax = plt.subplots(3,length(Ts),sharey = "row",figsize = [10,8])
    fig.subplots_adjust(wspace = 0.0)
    for (Ti,T) in enumerate(Ts)
        data_0 = gaps_0[Ti]
        data_pi = gaps_pi[Ti]
        for (nui,nu) in enumerate(nus)
            ax[2,Ti].plot(pos,data_0[:,nui],"-o",ms = 3,label = "\$ \\nu = $(round(nu,digits = 4))\$")
            ax[3,Ti].plot(pos,data_pi[:,nui],"-o",ms = 3,label = "\$ \\nu = $(round(nu,digits = 4))\$")
        end
        ax[1,Ti].text(.5,.7,fontsize = "small",
            "\\begin{align*}L &= $(L) \\\\[-5pt] g &= $(g)\\\\[-5pt] J_z &= $(jz)\\\\[-5pt] T &=$(T) \\end{align*}",
                      transform = ax[1,Ti].transAxes, bbox = Dict("facecolor"=>"white"))
        ax[3,Ti].set_xlabel("\$\\lambda_x \$")
        ax[2,Ti].set_xlabel("\$\\lambda_x \$")
        ax[1,Ti].set_xlabel("\$n\$")
        ax[1,Ti].plot(Bs[Ti],"-o",ms = 5,color = "black",label = "\$ T = $(T) \$")
    end
    ax[2,1].set_yscale("log")
    ax[3,1].set_yscale("log")
    ax[2,1].set_ylabel("pseudospectrum gap at 0")
    ax[3,1].set_ylabel("pseudospectrum gap at \$\\pi\$")
    ax[1,1].set_ylabel("\$ b_n\$")
    ax[2,4].legend(loc = "lower right")
    ax[3,4].legend(loc = "lower right")
    fig.show()
end

##################################################################################################
##################################################################################################
#=
This produces the time evolution via Kyrlov methods
(KrylovKit) from the b_n produced by the Lanczos on the 
Floquet Hamiltonian. It produces time evolution for 
continuous time, but really only valid at stroboscopic times.
=#

function krylov_time_evolve()
    # system parameters
    L = 8
    jx = 1.0
    g = .3
    jz = 0.05
    Ts = 3/2 .*[.5,1.5,3.6,5.0,8.25]

    # times = 10 .^ range(-1,3,length = 40) # log scale time points
    times = range(1,1e1,length = 40) # linear scale time points
    signal = zeros(length(times),5)
    Bs = Bs_5 # using data from top of file (small system size)
    n2 = 2000 # plateau length, increase if times get longer
    # step through different periods (i)
    for i in 1:5
        bs = Bs[i] # get b_n data
        # add a plateau of length n2
        bs2 = Float64.(vcat(bs,[bs[end]  for i in 1:n2]))
        steps2 = length(bs2)
        # push!(Bs,bs2)
        # make a HK -> M
        M = spdiagm(1 => bs2,-1=>bs2)
        # start with wavefunction on site 1
        op = zeros(steps2+1); op[1] = 1.0
        # time evolve (KrylovKit)
        for (ti,time) in enumerate(times)
            println(ti)
            tmp4,info = exponentiate(M,im*time,op,
                                     tol = (1e-2)/times[end],
                                     maxiter = 1000)
            println(info)
            tmp5  = real.(op'*tmp4)
            # if tmp5 <.05 || info.converged == 0
            if info.converged == 0
                break
            else 
                signal[ti,i] = tmp5
            end
        end
    end
    # rest is plotting
    fig,ax = plt.subplots(1,1)
    markers = ["-o","-s","-^","-*","->"]
    for i in 1:5
        ax.plot(times,signal[:,i],markers[i],ms= 3.,lw = .5,color = "C$(i-1)",label = "\$T = $(Ts[i]) \$")
    end
    ax.set_xscale("log")
    ax.set_xlabel("\$ t \$")
    ax.set_ylabel("\$ A_\\infty(t) \$")
    ax.legend()
    plt.show()
end



##################################################################################################
##################################################################################################
#=
testing for pi modes in L = 4 system
=#
function small_mat(a,b,c)
    return sqrt((a^2 + b^2 + c^2)/2 + sqrt((a^2 + b^2 + c^2)^2 - 4*a^2*c^2)/2)
end


##################################################################################################
##################################################################################################
#= 
testing for 0 modes in L =4 system, just a change of a minus sign
=#
function small_mat2(a,b,c)
    return sqrt((a^2 + b^2 + c^2)/2 - sqrt((a^2 + b^2 + c^2)^2 - 4*a^2*c^2)/2)
end


##################################################################################################
##################################################################################################
#=
This code will look into modifications to the bns to see 
if we can find a simple model that produces pi modes while 
still looking like an SSH. 

The take home message is that we can of course produce 
an effective SSH model that hosts states that have +/- pi 
energy, but currently it is unclear how we "pin" the states 
to +/- energy. For the zero mode it is obvious (the usual 
edge mode calculation)...
=#

function modified_ssh_edge_mode()
    # system parameters
    L = 14
    g = .3
    jx = 1.0
    # solve for many period lengths
    Ts = range(.5,8.25,length = 200)
    # store results in these empty lists.
    engs= []
    Bs= []
    engs2 = []
    Bs2 = []
    engs2B = []
    Bs2B = []
    engs3 = []
    Bs3 = []
    Bs3B = []
    engs3B =[]

    cutoff = 24
    # we will save eigenvalue of first 3 b_n
    # (first 4 sites) in this array
    foursite_thing = []
    # step through period lengths
    for (Ti,T) in enumerate(Ts)
        # construct M matrix from Fabian paper
        tmpM =  get_M(L,T,g)
        # take long
        M = -im*log(tmpM) # antisymmetric, pure imaginary
        # perform Lanczos on M
        ops = zeros(Complex{Float64},size(M)[1]); ops[1] = 1.0
        ops,bs = lanczos_steps2(M,ops,1e-16)
        # get eigenvalue (pi?) of first 4  sites, save it
        # calling small_mat function
        push!(foursite_thing,abs(pi - small_mat(bs[1:3]...)))
        # save b_ns
        push!(Bs,real.(bs))
        HK = ED.make_M(bs)
        push!(engs,real.(eigvals(HK)))
        # save very small set of b_n
        bs2 = real.(bs)[1:end-cutoff]
        push!(Bs2,bs2)
        HK2 = ED.make_M(bs2)
        push!(engs2,real.(eigvals(HK2)))
        # now extend the very small subsystem above
        # by a few more sites..
        bs2B = real.(bs)[1:end-cutoff+2]
        push!(Bs2B,bs2B)
        HK2B = ED.make_M(bs2B)
        push!(engs2B,real.(eigvals(HK2B)))
        # now extend the very small (Bs2) set above by repeating
        # last two value to somewhat replicate
        # bulk of exact b_n 
        bs3 = copy(bs2)
        a = bs3[end-1]
        b = bs3[end]
        for i in 1:cutoff/2
            push!(bs3,a)
            push!(bs3,b)
        end
        push!(Bs3,bs3)
        HK3 = ED.make_M(bs3)
        push!(engs3,real.(eigvals(HK3)))
        # now do the same thing, extending
        # the small (Bs2B) subsystem
        # faking the bulk by repeating the
        # last two hoppings
        bs3B = copy(bs2B)
        a = bs3B[end-1]
        b = bs3B[end]
        for i in 1:cutoff/2
            push!(bs3B,a)
            push!(bs3B,b)
        end
        push!(Bs3B,bs3B)
        HK3 = ED.make_M(bs3B)
        push!(engs3B,real.(eigvals(HK3)))
    end
    # rest is plotting
    fig,ax = plt.subplots(5,2)
    brg = matplotlib.cm.get_cmap("brg",length(Ts))
    for i in 1:5
        ax[i,1].axhline(y = pi,color = "red",alpha = .3)
        ax[i,1].axhline(y = -pi,color = "red",alpha = .3)
    end
    for (Ti,T) in enumerate(Ts)
        xvals = [T for k in 1:length(engs[Ti])]
        xvals2 = [T for k in 1:length(engs2[Ti])]
        xvals3 = [T for k in 1:length(engs3[Ti])]
        xvals4 = [T for k in 1:length(engs2B[Ti])]
        xvals5 = [T for k in 1:length(engs3B[Ti])]
        ax[1,1].plot(xvals,engs[Ti],".",ms = 3.,color = brg(Ti*1.0/length(Ts)))
        ax[2,1].plot(xvals2,engs2[Ti],".",ms = 3.,color = brg(Ti*1.0/length(Ts)))
        ax[3,1].plot(xvals3,engs3[Ti],".",ms = 3.,color = brg(Ti*1.0/length(Ts)))
        ax[4,1].plot(xvals4,engs2B[Ti],".",ms = 3.,color = brg(Ti*1.0/length(Ts)))
        ax[5,1].plot(xvals5,engs3B[Ti],".",ms = 3.,color = brg(Ti*1.0/length(Ts)))
        if Ti%50 == 0
            println("Ti = $Ti")
            display(Bs[Ti])
            ax[1,2].plot(Bs[Ti],color = brg(Ti*1.0/length(Ts)),"-o",label = "\$T = $(T)\$")
            ax[2,2].plot(Bs2[Ti],color = brg(Ti*1.0/length(Ts)),"-o",label = "\$T = $(T)\$")
            ax[3,2].plot(Bs3[Ti],color = brg(Ti*1.0/length(Ts)),"-o",label = "\$T = $(T)\$")
            ax[4,2].plot(Bs2B[Ti],color = brg(Ti*1.0/length(Ts)),"-o",label = "\$T = $(T)\$")
            ax[5,2].plot(Bs3B[Ti],color = brg(Ti*1.0/length(Ts)),"-o",label = "\$T = $(T)\$")
            println(Bs3B[Ti])
        end
    end
    fig.show()
end


##################################################################################################
##################################################################################################
#= 
A_inf from floquet Hamiltonian, not really used...
kept as an example

=#
function ainf_example_HF()
    L = 6 # system size
    sigx = ED.sigX(L,1) # sigma_1^x
    # system parameters
    jx = 1.0
    g = .3
    jz = 0.05
    T = 3/2*8.25
    Ts = 3/2 .* [.5,1.5,3.6,5,8.25]
    # times to calculate
    times = range(1,1e2,length = 1000)
    data = []
    # step through different periods
    for (Ti,T) in enumerate(Ts)
        # construct propagator
        U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,g,1.0))
        U2 = exp(-im * T/3 .* ED.H_Delta(L,jx,0.0,0.0,1.0))
        U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0,jz,0.0,1.0))
        U = U3*U2*U1
        # get floquet Hamiltonian
        H_F = im/T*log(U)
        # construct A_inf
        ainf = []
        u,v = ED.diagonalize(H_F)
        sigxv = v'*sigx*v
        for (ti,time) in enumerate(times)
            print(ti, "\r")
            D = Diagonal(exp.(-im * time .* u))
            push!(ainf,1/(2^L)*tr(D'*sigxv*D*sigxv))
        end
        push!(data,ainf)
    end
    # rest is plotting
    fig,ax = plt.subplots(1,1)
    markers = ["-o","-s","-^","-*","->"]
    for (Ti,T) in enumerate(Ts)
        ax.plot(times,data[Ti],markers[Ti],ms= 3.,lw = .5,color = "C$(Ti-1)",label = "\$T = $(T) \$")
    end
    ax.set_xscale("log")
    ax.set_xlabel("\$ t \$")
    ax.set_ylabel("\$ A_\\infty(t) \$")
    ax.legend()
    plt.show()

end

##################################################################################################
##################################################################################################
#=
This generates the majoranas in the spin basis, 
uses an ED function.
L = system size
=#
function gen_majoranas(L)
    majoranas = []
    for i in 1:L
        tmp = ED.majorana(L,i)
        push!(majoranas,tmp) #sigx majorana
        sigz = ED.sigZ(L,i)
        push!(majoranas,im .* tmp*sigz) #sigy majorana
    end
    return majoranas
end
##################################################################################################
##################################################################################################
#=
decompose operator into pauli weights
L = system size 
op = operator
=#
function pauli_decomp(L,op,threshold = 1e-8)
    tmp_vec = []
    tmp_tot = 0.0
    for k in 0:(2^(2*L)-1)
        tmp_string = ""
        tmp = I
        tmp2 = k
        for l in 1:L
            tmp3 = tmp2&3
            if tmp3 == 0 
                tmp_string = tmp_string*"-"
                tmp = tmp
            elseif tmp3 == 1
                tmp_string = tmp_string*"x"
                tmp = tmp*ED.sigX(L,l)
            elseif tmp3 == 2
                tmp_string = tmp_string*"y"
                tmp = tmp*ED.sigX(L,l)*ED.sigZ(L,l)*im
            else # tmp3 == 3
                tmp_string = tmp_string*"z"
                tmp = tmp*ED.sigZ(L,l)
            end
            tmp2 = tmp2 >> 2
        end
        overlap = op_dot(L,tmp,op)
        if abs(overlap)>threshold
            tmp_tot += abs(overlap)^2
            println("$(tmp_string), overlap = $(abs(overlap)^2) ")
        end
    end
    println("total = $(tmp_tot)")
    println("")
end

##################################################################################################
##################################################################################################
#=
Get overlap of operator with all possible majoranas
L = system size
op = operator
=#
function majorana_decomp(L,op)
    tmp_vec = []
    majoranas = gen_majoranas(L)
    for l in 1:length(majoranas)
        tmp = op_dot(L,majoranas[l],op)
        push!(tmp_vec,abs(tmp)^2)
    end
    return tmp_vec
end





##################################################################################################
##################################################################################################
#=
Free Floquet drive, playing around with b_n from Lanczos on H_F.
This compares the b_n that result from the Majorana basis and 
the M matrix in the Fabian paper to the the b_n that result from 
the spin basis and turning off J_z. The latter requires performing a 
second Lanczos proceedure in order to have the b_n agree with the former. 
This second Lanczos proceedure hinges on an eigenspectrum calculation on 
HK, which only makes sense if the full HK is obtained. This is the case 
in the free limit. In the interacting limit, we stop the calculation early
and thus it is unclear what the spectrum of an incomplete HK represents.. 
Again, we are in the free case here, so everything is totally fine..
=#
function free_lanczos()
    # system parameters
    L = 4
    g = .3
    jx = 1.0
    Ts = [.5,1.5,3.6,5.,8.25]
    sigx = ED.sigX(L,1)
    # this generates the majoranas in the spin basis
    majoranas = gen_majoranas(L)
    # generate Hamiltonians for time evolution in binary drive
    mat1 = ED.H_Delta( L, 0.0, 0.0,   g, 1.0)
    mat2 = ED.H_Delta( L,  jx, 0.0, 0.0, 1.0)
    # store solutions here..
    engs1,engs2,engs1b = [],[],[]
    Bs1,Bs2,Bs1b = [],[],[]
    ainfs0, ainfs1,ainfs2,ainfs1b = [],[],[],[]
    # times = range(1e-3,1e2,length = 200)
    times = 0:100
    # step through periods
    for (Ti,T) in enumerate(Ts)
        # generate propagator
        U1 = exp(-im * T/2 .* mat1); U2 = exp(-im * T/2 .* mat2);
        U = U2*U1
        # find HF
        H_F = im*log(U)
        # decompose H_F into its pauli representation..
        # this is just for the terminal
        pauli_decomp(L,H_F)
        # go through all majoranas
        for p in 1:length(majoranas)
            # conjugate each one
            tmp_op = U'*majoranas[p]*U
            # find their overlap with each other
            tmp_vec = majorana_decomp(L,tmp_op)
            # print it out to the terminal..
            println("p = $p, vec = $(tmp_vec), total = $(sum(tmp_vec))")
        end
        # perform Lanczos on H_F
        # note that this method is slightly different from
        # lanczos_steps2.... core algorithm is the same,
        ops_0,bs1 = lanczos_steps(L,H_F,sigx,4*L)
        push!(Bs1,bs1)
        # get M matrix
        tmpM = get_M(L,T,g) # Fabian matrix #free_mat(L,g,T)
        M = -im*log(tmpM) # antisymmetric, pure imaginary
        # perform lanczos on that
        ops = zeros(Complex{Float64},size(M)[1]); ops[1] = 1.0 
        ops,bs2 = lanczos_steps2(M,ops,1e-16)
        push!(Bs2,bs2)
        # get eigenvalues of both to compare
        HK1 = ED.make_M(bs1); HK2 = makeM(L,bs2)
        u1,v1 =  eigen(HK1); u2,v2 =  eigen(HK2)
        # shift energies of HK1 since the many particle
        # spectrum has produced energies that are outside
        # the FBZ (after performing Lanczos).
        u1b = shift_energies(u1)
        HK1b = v1*Diagonal(u1b)*(v1')
        # perform Lanczos again after the shift.
        opB = zeros(Complex{Float64},size(HK1b)[1]); opB[1] = 1.0 
        opsB,bs1b = lanczos_steps2(HK1b,opB,1e-16)
        # printing out majorana decomposition, colors refer to plots
        println("decomposition of blue")
        for m in 1:length(ops_0)
            tmp_vec = majorana_decomp(L,ops_0[m])
            println("m = $m, vec = $(tmp_vec), total = $(sum(tmp_vec))")
        end
        # printing out majorana decomposition, colors refer to plots
        println("decomposition of green")
        for m in 1:length(opsB)
            tmp_op = zeros(Complex{Float64},2^L,2^L)
            for k in 1:length(opsB[m])
                tmp_op += opsB[m][k]*ops_0[k]
            end
            tmp_vec = majorana_decomp(L,tmp_op)
            println("m = $m, vec = $(tmp_vec), total = $(sum(tmp_vec))")
        end

        push!(Bs1b,bs1b)
        push!(engs1,real.(u1))
        push!(engs2,real.(u2))
        push!(engs1b,real.(u1b))

        # set up ainf calculation
        tmp_ainfs0 = zeros(length(times))
        tmp_ainfs1 = zeros(length(times))
        tmp_ainfs2 = zeros(length(times))
        tmp_ainfs1b = zeros(length(times))
        quasi_phases,quasi_vecs = eigen(U)
        sigxv = quasi_vecs'*sigx*quasi_vecs
        # step through time and fill in results
        for (ti,time) in enumerate(times)
            tmp_ainfs1[ti] = real((v1*Diagonal(exp.(im*time.*u1 ))*(v1'))[1,1])
            tmp_ainfs2[ti] = real((v2*Diagonal(exp.(im*time.*u2 ))*(v2'))[1,1])
            tmp_ainfs1b[ti] = real((v1*Diagonal(exp.(im*time.*u1b ))*(v1'))[1,1])
            phases = Diagonal(quasi_phases .^ time)
            tmp_ainfs0[ti] = 1/(2^L)*real(tr(phases'*sigxv*phases*sigxv))
        end
        push!(ainfs1,tmp_ainfs1); push!(ainfs2,tmp_ainfs2); push!(ainfs1b,tmp_ainfs1b)
        push!(ainfs0,tmp_ainfs0)
    end
    # rest is plotting
    fig,ax = plt.subplots(3,length(Ts),sharey = "row")
    fig.subplots_adjust(wspace = 0)
    for (Ti,T) in enumerate(Ts)
        ax[1,Ti].plot(Bs1[Ti],"-o",label = "spin")
      ax[1,Ti].plot(Bs2[Ti],"-x",ms = 10,label = "free")
      ax[1,Ti].plot(Bs1b[Ti],"--o",label = "spin, modified")
      ax[2,Ti].plot(engs1[Ti],"-o",label = "spin")
      ax[2,Ti].plot(engs2[Ti],"-x",label = "free")
      ax[2,Ti].plot(engs1b[Ti],"--o",label = "spin, modified")
      ax[2,Ti].axhline(y = pi,alpha = .2,color = "red")
      ax[2,Ti].axhline(y = -pi,alpha = .2,color = "red")
      ax[3,Ti].plot(times ,ainfs1[Ti],"-o",label = "spin")
      ax[3,Ti].plot(times ,ainfs2[Ti],"-x",label = "free")
      ax[3,Ti].plot(times ,ainfs1b[Ti],"-s",label = "spin modified")
      ax[3,Ti].plot(times , ainfs0[Ti],"-^",label = "exact")
        ax[3,Ti].set_xlabel("\$t/T \$")
        ax[2,Ti].set_xlabel("index")
        ax[1,Ti].set_xlabel("\$n \$")

        ax[1,Ti].text(.4,.75,
                     fontsize = "small",
                     "\\begin{align*}g &= $(g) \\\\[-5pt] T &= $(T)\\end{align*}",
                     transform = ax[1,Ti].transAxes,
                     bbox = Dict("facecolor"=>"white"))
        ax[2,Ti].text(.4,.05,
                     fontsize = "small",
                     "\\begin{align*}g &= $(g) \\\\[-5pt] T &= $(T)\\end{align*}",
                     transform = ax[2,Ti].transAxes,
                     bbox = Dict("facecolor"=>"white"))
    end
    ax[1,1].legend()
    ax[3,1].legend()
    ax[3,1].set_ylabel("\$A_\\infty\$")
    ax[2,1].set_ylabel("spectrum of \$ H_K \$")
    ax[1,1].set_ylabel("\$b_n \$")
    fig.show()
end

##################################################################################################
##################################################################################################
#=
Produces b_n from Lanczos on Floquet Hamiltonian 
for interacting case.
=#
function gen_bn_lanczos_HF()
    Ls = [6,8] # system sizes
    # system parameters
    g = .3
    jz = 0.05
    jx = 1.0
    T = 3/2 * 8.25
    N = 20 # number of states to calculate
    Bs = []
    # step through system sizes
    for (Li,L) in enumerate(Ls)
        # sigma_1^x and propagator
        sigx = ED.sigX(L,1)
        U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,g,1.0))
        U2 = exp(-im * T/3 .* ED.H_Delta(L,jx,0.0,0.0,1.0))
        U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0,jz,0.0,1.0))
        U = U3*U2*U1
        # floquet Hamiltonian
        H_F = im*log(U)
        # perform lanczos on H_F and sigx
        # note it is slightly different from
        # lanczos_steps2
        ops,bs = lanczos_steps(L,H_F,sigx,N)
        println(bs)
        push!(Bs,bs)
    end
    println("T = $(2/3*T)")
    # plotting
    fig,ax = plt.subplots(1,1)
    for (Li,L) in enumerate(Ls)
        ax.plot(abs.(Bs[Li]),"-o",label = "\$ L = $L\$")
    end
    ax.legend()
    fig.show()
end

##################################################################################################
##################################################################################################
#=
Takes b_ns and N and make a platue of length N.
Used for lanczos_interacting() below.
=#
function make_res(bs,N)
    bs_tmp = [bs[end] for i in 1:N]
    return vcat(bs,bs_tmp)
end

##################################################################################################
##################################################################################################
#=
Compares b_n from Lanczos before and after energy shift 
in interacting limit. The energy shift involves eigenspectrum 
calculation, so it is not clear what the "energies" or eigenvectors 
mean if we stop our calculation early since the Krylov subspace 
for sigma_1^x is too large.

=#
function lanczos_interacting()
    # system parameters
    L = 8; jx = 1.0; g = .3; jz = 0.05
    Ts = 3/2 .*[.5,1.5,3.6,5.0,8.25]
    # number of states to calculate
    Ns = [0,100,200]
    # save results in these empty arrays
    old_engs_Ts = []; new_engs_Ts = []
    old_Bs_Ts = []; new_Bs_Ts = []
    exact_ainf_Ts = []
    periods =0:100 # number of periods to time evolve
    ainfs_old_Ts = []; ainfs_new_Ts = []
    # step through periods
    for (Ti,T) in enumerate(Ts)
        
        old_engs = []; new_engs = []
        old_Bs = []; new_Bs =[]
        ainfs_old = []; ainfs_new = []
        # steps through number of b_n to calculate
        for (Ni,N) in enumerate(Ns)
            # make flat plateau for time evolve
            bs_old = make_res(Bs_5[Ti],N)
            push!(old_Bs,bs_old)
            # make M and M_tmp (shifted energies)
            M_old = ED.make_M(bs_old); e_old,v_old = eigen(M_old); push!(old_engs,e_old)
            M_tmp = v_old*Diagonal(shift_energies(e_old))*v_old'
            # display(maximum(abs.(shift_energies(e_old) - e_old)))
            # perform second lanczos
            op = zeros(length(e_old)); op[1] = 1.0
            ops,bs_new = lanczos_steps2(M_tmp,op)
            # update M_new with new bs after second lanczos
            push!(new_Bs,bs_new)
            M_new = ED.make_M(bs_new); e_new,v_new = eigen(M_new); push!(new_engs,e_new)
            
            steps =length(bs_old)
            # set up ainf calculation
            tmp_ainfs_old = zeros(length(periods))
            tmp_ainfs_new = zeros(length(periods))
            # step through time and perform ainf
            for (ti,time) in enumerate(periods)
                tmp_ainfs_old[ti] = real((v_old*Diagonal(exp.(im*time.*e_old ))*(v_old'))[1,1])
                tmp_ainfs_new[ti] = real((v_new*Diagonal(exp.(im*time.*e_new ))*(v_new'))[1,1])
            end
            push!(ainfs_old,tmp_ainfs_old)
            push!(ainfs_new,tmp_ainfs_new)
        end
        # save results
        push!(old_Bs_Ts,old_Bs)
        push!(new_Bs_Ts,new_Bs)
        push!(old_engs_Ts,old_engs)
        push!(new_engs_Ts,new_engs)
        push!(ainfs_old_Ts,ainfs_old)
        push!(ainfs_new_Ts,ainfs_new)
        U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,g,1.0))
        U2 = exp(-im * T/3 .* ED.H_Delta(L,jx,0.0,0.0,1.0))
        U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0,jz,0.0,1.0))
        U = U3*U2*U1
        # get exact ainf
        push!(exact_ainf_Ts, gen_Ainf(L,U,periods))
    end
    # rest is plotting
    fig,ax = plt.subplots(5,length(Ts))
    fig.subplots_adjust(hspace = .25,wspace = .24)
    for (Ti,T) in enumerate(Ts)
        println("Ti = $Ti")
        tmp_p,tmp_m,tmp_time = make_ainf_plus_minus(exact_ainf_Ts[Ti],periods)
        ax[4,Ti].plot(tmp_time,abs.(tmp_p),"k",alpha = .3,lw = 4,label = "exact")
        ax[5,Ti].plot(tmp_time,abs.(tmp_m),"k",alpha = .3,lw = 4,label = "exact")
        for (Ni,N) in enumerate(Ns)
            println(length(old_Bs_Ts[Ti]))
            ax[1,Ti].plot(old_Bs_Ts[Ti][Ni],"-o",ms = 3,label = "\$old: $N\$")
            ax[2,Ti].plot(old_Bs_Ts[Ti][Ni][1:20],"-o",ms = 3,label = "\$old: $N\$")
            ax[1,Ti].plot(new_Bs_Ts[Ti][Ni],"-o",ms = 3,label = "\$new: $N\$")
            ax[2,Ti].plot(new_Bs_Ts[Ti][Ni][1:20],"-o",ms = 3,label = "\$new: $N\$")
            xvals = range(0,1,length = length(old_engs_Ts[Ti][Ni]))
            ax[3,Ti].plot(xvals,old_engs_Ts[Ti][Ni],"-o",ms = 3.,label = "\$old: $N\$")
            ax[3,Ti].plot(xvals,new_engs_Ts[Ti][Ni],"-o",ms = 3.,label = "\$new: $N\$")

            tmp_p,tmp_m,tmp_time = make_ainf_plus_minus(ainfs_old_Ts[Ti][Ni],periods)
            ax[4,Ti].plot(tmp_time,abs.(tmp_p),"-o",ms = 1)
            ax[5,Ti].plot(tmp_time,abs.(tmp_m),"-o",ms = 1)
            tmp_p,tmp_m,tmp_time = make_ainf_plus_minus(ainfs_new_Ts[Ti][Ni],periods)
            ax[4,Ti].plot(tmp_time,abs.(tmp_p),"-o",ms = 1)
            ax[5,Ti].plot(tmp_time,abs.(tmp_m),"-o",ms = 1)
        end
        ax[1,Ti].set_xlabel("\$n\$",labelpad =0.0)
        ax[2,Ti].set_xlabel("\$n\$",labelpad = 0.0)
        ax[3,Ti].set_xlabel("\$i\$",labelpad = 0.0)
        ax[4,Ti].set_xlabel("\$t\$",labelpad = 0.0)
        ax[5,Ti].set_xlabel("\$t\$",labelpad = 0.0)
        ax[1,Ti].set_ylabel("\$b_n\$",labelpad = 0.0)
        ax[2,Ti].set_ylabel("\$b_n\$",labelpad = 0.0)
        ax[3,Ti].set_ylabel("spectrum",labelpad = 0.0)
        ax[4,Ti].set_ylabel("\$|A_\\infty^+|\$",labelpad = 0.0)
        ax[5,Ti].set_ylabel("\$|A_\\infty^-|\$",labelpad = 0.0)
    end
    ax[3,1].legend(loc = "lower right",framealpha = 1.0,handlelength =  1.)
    ax[5,1].legend(loc = "upper right",framealpha = 1.0,handlelength =  1.)
    fig.show()
end

##################################################################################################
##################################################################################################
#=
Take ainf and make ainf+/- as described in paper with Fabian

=#
function make_ainf_plus_minus(ainf,times)
    ainf_plus = (ainf[2:end] + ainf[1:end-1]) ./ 2
    ainf_minus= (ainf[2:end] - ainf[1:end-1]) ./ 2
    ret_times = (times[2:end] + times[1:end-1]) ./ 2
    return ainf_plus,ainf_minus,ret_times
end
##########################################################################################
#########################################################################################
#=
Computes W and takes log of W= W0+W1, for some appropriately chosen W0.
Also note that T/3 of interacting drive = T/2 of the free drive so that they agree with
each other when Jz=0. Thus T_int = (3/2) T_free
=#
function dothing()
    g =.3
#    T = .1 #=zero mode=#
#    T = pi/(g)*(3/2) + .4 #=Pi-mode=#

#= Parameters of Free drive in Fig 1: [2.0,3.5,5.5,8.25], g=0.3
This is respectively [0, 0-Pi,Trivial,Pi]
=#

    T =  (3/2)*8.25
    println("T=", T)
    L = 10
    jx = 1.0
    jz = .01
    println("jz=", jz)
    N = 10 #= For free case N <=2*L =#

    Wfree= free_W(L,(2/3)*T,g)
    u,v = eigen(Wfree)
    println("spectrum of Wfree")
    display(u)

    #= propagator. In paper written as exp(-i T/2..) exp(-i T/2..) exp(-iT/2..)
    in order to avoid introducing T'.
    =#
    
    U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,  g,1.0))
    U2 = exp(-im * T/3 .* ED.H_Delta(L, jx,0.0,0.0,1.0))
    U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0, jz,0.0,1.0))
    U = U3*U2*U1
   
    sigx = ED.sigX(L,1) # sigma_1^x
    W,qs = arnoldi_ED(L,N,U,sigx,1e-15) # find W
    println("|W|")
    display(round.(abs.(W),sigdigits = 2))
    println("|W'*W|")
    display(round.(abs.(W'*W),sigdigits = 2))
    println("|W*W'|")
    display(round.(abs.(W*W'),sigdigits = 2))

    u,v = eigen(W)
    println("spectrum of W")
    display(u)

  #  println("I*log W")
  #  display(round.(im*log(W),sigdigits = 2))

    retM = log(W) 
   
   println("main diagonal terms of |log(W)|")
    tmp0M = []
    for i in 1:(size(retM)[1])
        println(round.(abs.(retM[i,i]),sigdigits = 2))
        push!(tmp0M,retM[i,i])
    end

    println("first diagonal terms of |log(W)|")
    tmp1M = []
    for i in 1:(size(retM)[1]-1)
        println(round.(abs.(retM[i,i+1]),sigdigits = 2))
        push!(tmp1M,retM[i,i+1])
    end
    
    tmp2M = []
    println("second diagonal terms of |log(W)|")
    for i in 1:(size(retM)[1]-2)
        println(round.(abs.(retM[i,i+2]),sigdigits = 2))
        push!(tmp2M,retM[i,i+2])
    end


    println("|W_0|")
    W0 = copy(W)
    for i in 1:size(W0)[1]
        for j in 1:size(W0)[2]
            if abs(i-j) > 2 #Change how many diagonal terms in W you would like to keep
                W0[i,j] = 0.0 + 0.0*im
            end
        end
    end
    display(round.(abs.(W0),sigdigits = 2))

    println("|W0'*W0|")
    display(round.(abs.(W0'*W0),sigdigits = 2))
    println("|W0*W0'|")
    display(round.(abs.(W0*W0'),sigdigits = 2))

    u,v = eigen(W0)
    println("spectrum of W0")
    display(u)
   


   V1 = W - W0 #Here again you can manipulate things based on how you would like to split up W.
    #  V1=0
    println("|dW| = W - W_0")
    display(round.(abs.(V1),sigdigits = 2))
    println("|log W_0|")
    display(round.(abs.(log(W0)),sigdigits = 2))
#    println("I*log W_0")
#    display(round.(im*log(W0),sigdigits = 2))
    println("|W0^(-1)|")
    display(round.(abs.(W0^(-1)),sigdigits = 2))
    println("|log(W0) + W0^(-1)*dW|")
    ret = log(W0) + (W0^(-1))*V1  #Here make changes to how you would like to split the code
   
    display(round.(abs.(ret),sigdigits = 2))

    u,v = eigen(im*ret)
    println("spectrum of ret")
    display(u)


    println("main diagonal terms of |log(W0) + W0^(-1)*dW|")
    tmp0 = []
    for i in 1:(size(ret)[1])
        println(round.(abs.(ret[i,i]),sigdigits = 2))
        push!(tmp0,ret[i,i])
    end

    println("first diagonal terms of |log(W0) + W0^(-1)*dW|")
    tmp1 = []
    for i in 1:(size(ret)[1]-1)
        println(round.(abs.(ret[i,i+1]),sigdigits = 2))
        push!(tmp1,ret[i,i+1])
    end
    tmp2 = []
    println("second diagonal terms of |log(W0) + W0^(-1)*dW|")
    for i in 1:(size(ret)[1]-2)
        println(round.(abs.(ret[i,i+2]),sigdigits = 2))
        push!(tmp2,ret[i,i+2])
    end

    tmp3 = []
    println("third diagonal terms of |log(W0) + W0^(-1)*dW|")
    for i in 1:(size(ret)[1]-3)
        println(round.(abs.(ret[i,i+3]),sigdigits = 2))
        push!(tmp3,ret[i,i+3])
    end

    fig,ax = plt.subplots(1,1)
    ax.plot(abs.(tmp0),"-o",label = "main diagonal log(W0) +W0inv*dW ")
    ax.plot(abs.(tmp1),"-o",label = "n.n hopping log(W0) + W0inv *dW")
    ax.plot(abs.(tmp2),"-o",label = "n.n.n hopping log(W0) + W0inv *dW")
    ax.plot(abs.(tmp3),"-o",label = "n.n.n hopping log(W0) + W0inv *dW")
    plt.legend()
    plt.show()
   

    fig,ax = plt.subplots(1,1)
    ax.plot(abs.(tmp0M),"-o",label = "main diagonal log(W)")
    ax.plot(abs.(tmp1M),"-o",label = "n.n hopping log(W)")
    ax.plot(abs.(tmp2M),"-o",label = "n.n.n hopping log(W)")
    plt.legend()
    plt.show()

end



#=
Computes Ainf for different Jzs and Ls
=#

function simple_ainfFL()

    Ls = [4,6] # different lengths
    # Hamiltonian parameters
    jx = 1.0
    g = .3
    jz = 0.05
    T = 3/2*8.25

    t_len = 5000 # length of ainf, in periods
    # store results here
    data_ainf1 = zeros(t_len,length(Ls)) # exact ainf
    	       	 for (Li,L) in enumerate(Ls) # step through system sizes
               	 sigx = ED.sigX(L,1) # sigma_1^x
               	 println("L=",L," and ","T=",T)
               	 # propagator
               	   U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,  g,1.0))
              	   U2 = exp(-im * T/3 .* ED.H_Delta(L, jx,0.0,0.0,1.0))
               	   U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0, jz,0.0,1.0))
               	   U = U3*U2*U1
               	  # calculate ainf
               	  u,v = eigen(U)
               	  for ti in 1:t_len
                      U_tmp = v*Diagonal(u .^ ti)*v'
                      data_ainf1[ti,Li] = real(tr(sigx*U_tmp*sigx*(U_tmp'))/(2^L))
                   end
                end	

    		fig,ax = plt.subplots(1,1)
		ax.plot(1:t_len,data_ainf1[:,1],"*")
		ax.plot(1:t_len,data_ainf1[:,2],"+")
    		#for (Li,L) in enumerate(Ls)
     		#ax.plot(1:t_len,data_ainf1[:,Li])
     		#end
     		ax.set_xscale("log")
    		plt.show()
    		end

#= Another (faster) way of getting Ainf=#

function simple_ainfFL2()

    Ls = [10,12] # different lengths
    # Hamiltonian parameters
    jx = 1.0
    g = .3
    jz = 0.05
    T = 3/2*8.25

    t_len = 5000 # length of ainf, in periods
    # store results here
    data_ainf1 = zeros(t_len,length(Ls)) # exact ainf
    	       	 for (Li,L) in enumerate(Ls) # step through system sizes
               	 sigx = ED.sigX(L,1) # sigma_1^x
               	 println("L=",L," and ","T=",T)
                  # propagator
               	   U1 = exp(-im * T/3 .* ED.H_Delta(L,0.0,0.0,  g,1.0))
              	   U2 = exp(-im * T/3 .* ED.H_Delta(L, jx,0.0,0.0,1.0))
               	   U3 = exp(-im * T/3 .* ED.H_Delta(L,0.0, jz,0.0,1.0))
               	   U = U3*U2*U1
               	  # calculate ainf
		  quasi_phases,quasi_vecs = eigen(U)
		  sigxv = quasi_vecs'*sigx*quasi_vecs
    		  for ti in 1:t_len
      		  phases = Diagonal(quasi_phases .^ ti)
      		  data_ainf1[ti,Li] = 1/(2^L)*real(tr(phases'*sigxv*phases*sigxv))
    		  end
    		  end
               	  	

    		fig,ax = plt.subplots(1,1)
		ax.plot(1:t_len,data_ainf1[:,1],"*")
		ax.plot(1:t_len,data_ainf1[:,2],"+")
    		#for (Li,L) in enumerate(Ls)
     		#ax.plot(1:t_len,data_ainf1[:,Li])
     		#end
     		ax.set_xscale("log")
    		plt.show()
    		end
    


##################################################################################################
##################################################################################################
##################################################################################################
# run code down here.
# @timev free_W_edge_modes() # makes W and finds approximate edge modes for free case, quick to run
# @timev free_W_edge_modes2() # checks approximate edge modes for free case, quick to run
# @timev interacting_W_edge_modes() # finds approximate W edge modes, time evolves, quick to run
# @timev free_W_overview() # plots matrix elements of free majorana basis W and related matrices, quick to run
# @timev interacting_W_overview() # similar to above, but for interacting
# @timev floquet_ham_interacting_lanczos() # b_ns and shifting energies, interacting, taking data from top.
# @timev floquet_ham_free_lanczos() # b_ns and Krylov Hamiltonian spectrum, quick to run
# @timev testing_pseudospectrum() # trying to calculate clifford spectrum gap, not dumb
# @timev krylov_time_evolve() # krylov time evolve, left in as an example..
# @timev modified_ssh_edge_mode() # taking real b_n, replacing bulk seeing what happens, not dumb
# @timev ainf_example_HF() # ainf example using HF, so micromotion will be off, not really needed..
# @timev free_lanczos() # studying difference between b_n of majorana and spin basis Lanczos on HF
# @timev gen_bn_lanczos_HF()
# @timev lanczos_interacting()
# @timev simple_ainfFL()
 @timev simple_ainfFL2()
# @timev dothing()

end
