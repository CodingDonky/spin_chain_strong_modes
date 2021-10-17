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

mutable struct my_op
    L :: Int64
    coeffs :: Array{Float64,1}
    deltas :: Array{Bool,1}
    epsilons :: Array{Bool,1}
end

function make_Op(L)
    L :: Int64 = L
    coeffs::Array{Float64,1} = zeros(Float64,2^(2*L))
    deltas::Array{Bool,1} = zeros(Bool,2^(2*L))
    epsilons::Array{Bool,1} = zeros(Bool,2^(2*L))
    return my_op(L,coeffs,deltas,epsilons)
end

struct my_H
    L :: Int64
    strings :: Array{Int64}
    coeffs :: Array{Float64}
    deltas :: Array{Bool}
    epsilons :: Array{Bool}
end

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
    return my_H(L,H_strings,H_coeffs,H_deltas,H_epsilons) 
end


function hammingWeightParity(L::Int64,x::Int64)
    # https://en.wikipedia.org/wiki/Hamming_weight
    # we will use the simple version since
    # L can change. testing  shows this is faster for us
    tmp::Int64 = 0
    for i in 0:L-1
        tmp += (x>>i) %2
    end
    return tmp%2
end


function do_comm!(L,H,op,newOp)
    # newOp is the workspace operator
    ND = 2^(L) -1
    for i in 1:length(H.strings)
        hstring::Int64 = H.strings[i]
        hcoeff::Float64 = H.coeffs[i]; hcoeff = 2*hcoeff # probably not neccessary
        hdelta::Bool = H.deltas[i]
        hepsilon::Bool = H.epsilons[i]
        
        hstringX = hstring>>>L
        hstringZ = hstring & ND
  
        for j in 0:2^(2*L)-1
            if op.coeffs[j+1] != 0.0
                ostringX::Int64 = (j)>>>L
                ostringZ::Int64 = (j) & ND
                dot1::Bool = hammingWeightParity(L,hstringX & ostringZ)
                dot2::Bool = hammingWeightParity(L,hstringZ & ostringX)
                if dot1 != dot2 # if they are equal, they commute
                    new_coeff::Float64 = hcoeff*op.coeffs[j+1]
                    new_string::Int64 = xor(hstring,j)
                    new_delta::Bool = xor(hdelta,op.deltas[j+1])
                    new_epsilon::Bool = xor(hepsilon,op.epsilons[j+1],
                                            hdelta & op.deltas[j+1],
                                            dot1)
                    newOp.coeffs[new_string+1] = new_coeff +
                        newOp.coeffs[new_string+1]*(-1)^xor(new_epsilon,newOp.epsilons[new_string+1])
                    newOp.deltas[new_string+1] = new_delta
                    newOp.epsilons[new_string+1] = new_epsilon
                end
            end
        end
    end 
end

function do_b(op)
    return sqrt(op.coeffs'*op.coeffs)
end

function sub_ops!(L,coeff,op1,op2)
    # op1 - coeff*op2, modifies op1
    for j in 0:2^(2*L)-1
        # we assume op1_d == op2_d, which is typically true,
        # unless either op1_c == 0 or op2_c == 0, in that case
        # we want the delta corresponding to the non-zero part
        if op1.coeffs[j+1] != 0
            op1.coeffs[j+1] =op1.coeffs[j+1] -  coeff*op2.coeffs[j+1]*(-1)^xor(op1.epsilons[j+1],op2.epsilons[j+1])
        else
            op1.coeffs[j+1] = op2.coeffs[j+1]*(-coeff)
            op1.deltas[j+1] = op2.deltas[j+1]
            op1.epsilons[j+1] = op2.epsilons[j+1]
        end
    end
end


function update_step!(L,H,op1,op2,op3,b)
    # op1, op2, op3 : current, past, new (older)
    # L|op1) -> op3
    do_comm!(L,H,op1,op3)
    # L|op1) - b|op2) -> op3
    sub_ops!(L,b,op3,op2)
end

function overlap(L,op1,op2)
    tmp_sum ::Float64 = 0.0
    for k in 1:2^(2*L)
        tmp_sum += op1.coeffs[k]*op2.coeffs[k]*(-1)^xor(op1.epsilons[k],op2.epsilons[k])
    end
    return tmp_sum
end

function normOp!(Op,b)
    Op.coeffs  = Op.coeffs ./ b
end

function clearOp!(Op)
    # clear out operator, save memory
    Op.coeffs .= 0.0
    Op.deltas .= false
    Op.epsilons .= false
end

function reorthogonalize!(L,Os)
    for m in 1:length(Os)-1
        tmp = overlap(L,Os[end],Os[m])
        sub_ops!(L,tmp,Os[end],Os[m])
        tmpb = do_b(Os[end])
        normOp!(Os[end],tmpb)
    end
end

function print_helper(L,xstring,zstring)
    xtmp = ""; ztmp = ""
    for k in 0:(L-1)
        xval = (xstring>>k)%2
        if xval == 0
            xtmp = xtmp*"_"
        else
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
        

function print_op()
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
    H = make_H(L,jx,jy,jz,g)
    OpA = make_Op(L)
    OpA.coeffs[Int64(1<<L)+1] = 1
    OpB = make_Op(L)
    do_comm!(L,H,OpA,OpB)
    bs = [do_b(OpB)]
    normOp!(OpB,bs[end])

    # reorthogonalize ~ n^2
    Os = [OpA,OpB]
    for n in 2:steps
        println("n = $n")
        push!(Os,make_Op(L))
        update_step!(L,H,Os[end-1],Os[end-2],Os[end],bs[end])
        push!(bs,do_b(Os[end]))
        normOp!(Os[end],bs[end])
        reorthogonalize!(L,Os)
    end

    # no reorthog.. yolo
    # OpC = make_Op(L) # empty for now
    # Os = [OpB,OpA,OpC]
    # for n in 2:steps
    #     println("n = $n")
    #     OpOld = Os[2]
    #     OpCur = Os[1]
    #     OpNew = Os[3]
    #     update_step!(L,H,OpCur,OpOld,OpNew,bs[end])
    #     push!(bs,do_b(OpNew))
    #     normOp!(OpNew,bs[end])
    #     clearOp!(OpOld)
    #     Os = [OpNew,OpCur,OpOld]
    # end
    
    return bs
end

