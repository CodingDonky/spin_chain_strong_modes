
#= Basis operators will be stored as combinations of x's and z's and
their coefficients and powers of i.

For example: z_1 x_2 y_3 = z_1 x_2 (i x_3 z_3).  

This information can be represented with binary vectors for the powers
of x's and z' so the above would be: (1 is on the right, and 3 is on
the left) 

x_powers = 1 1 0 
z_powers = 1 0 1 

This can be combined into a single integer:

110101 (= 53). 

Thus we only need to store binary integers of length 2*L to represent
all the possible Pauli-matrices. This code will assume that the |O_n)
basis operators will have support on essentially all possible
Pauli-strings, so we store the coefficients in an array of length
2^(2*L), where the index +1 represents the Pauli-matrix string itself.

****For integrable and free cases this is a very bad idea! Use
libNewComm3.jl instead! There we use a dictionary and store only the
basis operators that are needed. So if only 2L are needed, as in the
free case, then that is vastly superior. If 2^(2L) are needed, then
this code should be faster.****

** This also means we are essentially capped to L <= 14. It should 
mean we can go further in b_n space. If we want to work with 
semi-infinite systems, then use libNewComm3.jl**

Now if one is given the powers of x's and z's, then one does not 
need to know the nuber of i's as that we restrict our attention to 
float64 coefficients (Hermitian operators) and thus whereever x and 
z coincide adds a power to i. However, it is important to keep track 
of the i's for the commutation, for example:

[x1,yy] = 2i zy = 2 z_1 x_2 z_2
[zz,1x] = 2i zy = 2 z_1 x_2 z_2

The above two commutations yield the same result. If we did not keep
track of the phase factors we get:

[x_1, x_1 z_1 x_2 z_2] = z_1 x_2 z_2 - x_1 z_1 x_2 z_2 x_1 = 2 z_1 x_2 z_2 

[z_1 z_2, x_2] = z_1 z_2 x_2 - x_2 z_1 z_2 = -2 z_1 x_2 z_2

There is a sign error in the second line. We need to keep track of the
phases to avoid this.

So the operators will be represented as Arrays of length 2^(2L):

1 for the coefficients (Float64) 

2 for the powers of i (each are Bool)
one for (-1) factors and one i factors

The Hamiltonian is an operator that will not grow and thus only 
have support on ~L pauli operators. Thus it makes no sense to 
store it like the basis operators, instead we will have a 4th 
array that stores the operator strings (Integers) themselves, 
instead of using them as indices. 

=#




function make_H(L,jx,jy,jz,g)
    H_strings ::Array{Int64} = []
    H_coeffs  ::Array{Float64} = []
    H_deltas ::Array{Bool} = []
    H_epsilons ::Array{Bool} = []
    for i in 0:(L-2)
        z_pos  = Int64(1<<i)
        z_bond = Int64(3<<i)
        x_bond = Int64(3<<(i+L))
        y_bond = Int64(x_bond+z_bond)

        if jx != 0
            push!(H_strings,x_bond)
            push!(H_coeffs,jx)
            push!(H_deltas, 0)
            push!(H_epsilons,0)
        end 
        if jy != 0
            push!(H_strings,y_bond)
            push!(H_coeffs,jy)
            push!(H_deltas,0)
            push!(H_epsilons,1)
        end 
        if jz != 0 
            push!(H_strings,z_bond)
            push!(H_coeffs,jz)
            push!(H_deltas,0)
            push!(H_epsilons,0)
        end 
        if g != 0
            push!(H_strings,z_pos)
            push!(H_coeffs,g)
            push!(H_deltas,0)
            push!(H_epsilons,0)
        end 
    end 
    if g != 0
        push!(H_strings,1<<(L-1))
        push!(H_coeffs,g)
        push!(H_deltas,0)
        push!(H_epsilons,0)
    end 
    return H_strings,H_coeffs,H_deltas,H_epsilons
end


function hammingWeight(L::Int64,x::Int64)
    # https://en.wikipedia.org/wiki/Hamming_weight
    tmp = 0
    for i in 0:L-1
        tmp += (x>>i) %2
    end
    return tmp
    
    # m1  ::UInt64 = 0x5555555555555555 #010101...
    # m2  ::UInt64 = 0x3333333333333333 #00110011...
    # m4  ::UInt64 = 0x0f0f0f0f0f0f0f0f #0000111100..
    # m8  ::UInt64 = 0x00ff00ff00ff00ff #0000000011...
    # m16 ::UInt64 = 0x0000ffff0000ffff #000000000000000011..
    # m32 ::UInt64 = 0x00000000ffffffff

    # x = (x&m1) + (x>>1)&m1;
    # x = (x&m2) + (x>>2)&m2;
    # x = (x&m4) + (x>>4)&m4;
    # x = (x&m8) + (x>>8)&m8;
    # x = (x&m16)+ (x>>16)&m16
    # x = (x&m32)+ (x>>32)&m32
    # return x
end


function do_comm(L,Hs,Hc,Hd,He,op_c,op_d,op_e)
    new_op_c :: Array{Float64} = zeros(Float64,2^(2*L))
    new_op_d :: Array{Bool} = zeros(Bool,2^(2*L))
    new_op_e :: Array{Bool} = zeros(Bool,2^(2*L))

    ND = 2^(L) -1
    for i in 1:length(Hs)
        hstring::Int64 = Hs[i]
        hcoeff::Float64 = Hc[i]; hcoeff = 2*hcoeff # probably not neccessary
        hdelta::Bool = Hd[i]
        hepsilon::Bool = He[i]
        
        hstringX = hstring>>>L
        hstringZ = hstring & ND
  
        for j in 0:2^(2*L)-1
            if op_c[j+1] != 0.0
                ostringX::Int64 = (j)>>>L
                ostringZ::Int64 = (j) & ND
                dot1::Bool = hammingWeight(L,hstringX & ostringZ)%2
                dot2::Bool = hammingWeight(L,hstringZ & ostringX)%2
                
                if dot1 != dot2 # if they are equal, they commute
                    
                    new_coeff::Float64 = hcoeff*op_c[j+1]
                    new_string::Int64 = xor(hstring,j)
                    new_delta::Bool = xor(hdelta,op_d[j+1])
                    new_epsilon::Bool = xor(hepsilon,op_e[j+1],
                                            hdelta & op_d[j+1],
                                            dot1)
                    
                    new_op_c[new_string+1] = new_coeff + new_op_c[new_string+1]*(-1)^xor(new_epsilon,new_op_e[new_string+1])
                    new_op_d[new_string+1] = new_delta
                    new_op_e[new_string+1] = new_epsilon
                end
            end
        end
    end 
    return new_op_c, new_op_d, new_op_e
end

function do_b(op_c)
    return sqrt(op_c'*op_c)
end

function sub_ops!(L,coeff,
                  op1_c,op1_d,op1_e,
                  op2_c,op2_d,op2_e)
    # op1 - coeff*op2, modifies op1
    for j in 0:2^(2*L)-1
        # we assume op1_d == op2_d, which is typically true,
        # unless either op1_c == 0 or op2_c == 0, in that case
        # we want the delta corresponding to the non-zero part
        if op1_c[j+1] != 0
            op1_c[j+1] =op1_c[j+1] -  coeff*op2_c[j+1]*(-1)^xor(op1_e[j+1],op2_e[j+1])
        else
            op1_c[j+1] = op2_c[j+1]*(-coeff)
            op1_d[j+1] = op2_d[j+1]
            op1_e[j+1] = op2_e[j+1]
        end
    end
end


function update_step(L,H_s,H_c,H_d,H_e,
                     op1_c,op1_d,op1_e,
                     op2_c,op2_d,op2_e,b)
    tmp_c,tmp_d,tmp_e = do_comm(L,H_s,H_c,H_d,H_e,op1_c,op1_d,op1_e)

    sub_ops!(L,b,
             tmp_c,tmp_d,tmp_e,
             op2_c,op2_d,op2_e)

    return tmp_c,tmp_d,tmp_e
end

function overlap(L,
                 op1_c,op1_d,op1_e,
                 op2_c,op2_d,op2_e)
    tmp_sum ::Float64 = 0.0
    for k in 1:2^(2*L)
        tmp_sum += op1_c[k]*op2_c[k]*(-1)^xor(op1_e[k],op2_e[k])
    end
    return tmp_sum
end

function print_helper(L,xstring,zstring)
    xtmp = ""; ztmp = ""
    for k in 0:(L-1)
        xval = (xstring>>k)%2
        if xval == 0
            xtmp = xtmp*"_"
        else # xval == 1
            xtmp = xtmp*"x"
        end
        xval = (zstring>>k)%2
        if xval == 0
            ztmp = ztmp*"_"
        else
            ztmp = ztmp*"z"
        end
    end
    return xtmp, ztmp
end
        

function print_op(L,op)
    for (ki,k) in enumerate(op)
        if k != 0
            ostringX = (ki-1)>>>L
            ostringZ = (ki-1) & (2^(L)-1)
            xstr,zstr = print_helper(L,ostringX,ostringZ)
            println("x = $(xstr), z = $(zstr), coeff = $(k)")
        end
    end
end

function get_bs(L,j,jz,g,gamma,threshold,steps)
    jx = j*(1 + gamma)/2
    jy = j*(1 - gamma)/2
    Hs,Hc,Hd,He = make_H(L,jx,jy,jz,g)
    DH = 2^(2*L)

    O_0_c::Array{Float64} =zeros(Float64,DH)
    O_0_d::Array{Bool} = zeros(Bool,DH)
    O_0_e::Array{Bool} = zeros(Bool,DH)
    O_0_c[Int64(1<<L)+1] = 1

    
    Os_c = [O_0_c]
    Os_d = [O_0_d]
    Os_e = [O_0_e]
    
    
    O_1_c,O_1_d,O_1_e = do_comm(L,
                                Hs,Hc,Hd,He,
                                O_0_c,O_0_d,O_0_e) # A_1
    bs = [do_b(O_1_c)] # b_1
    if bs[end]<1e-14
        println("breaking out...")
        return bs
    end
    O_1_c = O_1_c ./ bs[end]
    push!(Os_c,O_1_c)
    push!(Os_d,O_1_d)
    push!(Os_e,O_1_e)

    for n in 2:steps
        println("n = $n")
        O_n_c,O_n_d,O_n_e = update_step(L,
                                        Hs,Hc,Hd,He,
                                        Os_c[end],Os_d[end],Os_e[end],
                                        Os_c[end-1],Os_d[end-1],Os_e[end-1],
                                        bs[end])
        push!(bs,do_b(O_n_c))
        if bs[end]<1e-14
            println("breaking out...")
            break
        end
        O_n_c = O_n_c ./bs[end]
        push!(Os_c,O_n_c)
        push!(Os_d,O_n_d)
        push!(Os_e,O_n_e)
        deleteat!(Os_c,1)
        deleteat!(Os_d,1)
        deleteat!(Os_e,1)
        
        # # reorthogonalize
        # op1_c = Os_c[end]
        # op1_d = Os_d[end]
        # op1_e = Os_e[end]
        # for m in 1:length(Os_c)-1
        #     op2_c = Os_c[m]
        #     op2_d = Os_d[m]
        #     op2_e = Os_e[m]
        #     tmp = overlap(L,
        #                   op1_c,op1_d,op1_e,
        #                   op2_c,op2_d,op2_e)
        #     sub_ops!(L,tmp,
        #              op1_c,op1_d,op1_e,
        #              op2_c,op2_d,op2_e)
        # end
        # tmpb = do_b(op1_c)
        # op1_c = op1_c ./ tmpb
        # Os_c[end] = op1_c
        # Os_d[end] = op1_d
        # Os_e[end] = op1_e
    end
    return bs
end



