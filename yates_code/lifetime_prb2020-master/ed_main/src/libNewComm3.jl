#=
This is where commutator code goes.
removing the "myOp" struct, just now a dictionary
=#


# function Base.hash(x::Int64)
#     return UInt64(x)
# end


struct pauli_coeffs
    #=
    basically a tuple..
    Make sure it is not mutable, this helps 
    the compiler a great deal
    =#
    coeff::Float64
    delta::Bool
    epsilon::Bool
end 

function print_help(string::Int64)
    #=
    Helper function for printing the operator.
    =#
    xs = string>>>32
    zs = string & (2^(32)-1)
    ret = ""
    for i in 1:32
        tmpX = Bool(xs & 1)
        tmpZ = Bool(zs & 1)
        if tmpZ
            if tmpX
                ret = ret*"Y"
            else
                ret = ret*"Z"
            end 
        else
            if tmpX
                ret = ret*"X"
            else
                ret = ret*"_"
            end
        end 
        xs = xs>>>1
        zs = zs>>>1
    end
    return ret
end 

function print_op(strings::Dict{Int64,pauli_coeffs},num::Int64)
    #=
    Prints the operator.
    =#
    n = min(length(strings),num)
    for (i,key) in enumerate(keys(strings))
        coeffs = strings[key]
        pretty_string = print_help(key)
        @printf("% 2.1e, %d, %d, %s %d \n",
                coeffs.coeff,
                coeffs.delta,
                coeffs.epsilon,
                pretty_string,
                key)
        if i == num
            break
        end
    end 
    println("")
end

function make_H(L,jx,jy,jz,g)
    #=
    Make the Hamiltonian
    =#
    H::Dict{Int64,pauli_coeffs} = Dict() 
    for i in 0:(L-2)
        z_pos  = Int64(1<<i)
        z_bond = Int64(3<<i)
        x_bond = Int64(3<<(i+32))
        y_bond = Int64(x_bond+z_bond)
        if jx != 0
            H[x_bond] =  pauli_coeffs(jx,0,0)
        end 
        if jy != 0
            H[y_bond] = pauli_coeffs(jy,0,1)
        end 
        if jz != 0 
            H[z_bond] = pauli_coeffs(jz,0,0)
        end 
        if g != 0
            H[z_pos ] = pauli_coeffs(g ,0,0)
        end 
    end 
    if g != 0
        H[Int64(1<<(L-1))] = pauli_coeffs(g,0,0) 
    end 
    return H
end



function hammingWeight(x::Int64)
    # finding the number of 1's in a binary representation
    # the obvious way of doing it (iterating through the string)
    # might be too slow for me..
    # wikipedia has a good article on this..
    # https://en.wikipedia.org/wiki/Hamming_weight
    # tmp = 0 
    # for i in 0:32
    #     tmp += (x>>i) %2
    # end
    # return tmp

    m1  ::UInt64 = 0x5555555555555555 #010101...
    m2  ::UInt64 = 0x3333333333333333 #00110011...
    m4  ::UInt64 = 0x0f0f0f0f0f0f0f0f #0000111100..
    m8  ::UInt64 = 0x00ff00ff00ff00ff #0000000011...
    m16 ::UInt64 = 0x0000ffff0000ffff #000000000000000011..
    m32 ::UInt64 = 0x00000000ffffffff

    x = (x&m1) + (x>>1)&m1;
    x = (x&m2) + (x>>2)&m2;
    x = (x&m4) + (x>>4)&m4;
    x = (x&m8) + (x>>8)&m8;
    x = (x&m16)+ (x>>16)&m16
    x = (x&m32)+ (x>>32)&m32
    return x
end

function do_b(strings)
    tmp = 0
    for (string,coeffs) in strings
        tmp += abs(coeffs.coeff)^2
    end 
    return sqrt(tmp)
end

function do_mult(oldop,scalar)
    newOp ::Dict{Int64,pauli_coeffs} = deepcopy(oldop)
    for (string,coeffs) in newOp
        newOp[string] = pauli_coeffs(coeffs.coeff*scalar,coeffs.delta,coeffs.epsilon)
    end 
    return newOp
end 


function do_mult!(oldop,scalar)
    for (string,coeffs) in oldop
        oldop[string] = pauli_coeffs(coeffs.coeff*scalar,coeffs.delta,coeffs.epsilon)
    end 
end 

function norm_op!(op,b,threshold)
    if b> threshold
        for (string,coeffs) in op
            op[string] = pauli_coeffs(coeffs.coeff/b,coeffs.delta,coeffs.epsilon)
        end 
    else
        println("breaking out... b is too small")
        # this will likely throw an error... need to handle it better
        return 
    end
end



function clean_op!(op,threshold)
    for (string, coeffs) in op
        if abs(coeffs.coeff)<threshold
            delete!(op,string)
        end
    end
end

function do_subtract(op1,op2)
    newOp = deepcopy(op1)
    for (string2,coeffs2) in op2
        if haskey(newOp,string2)
            tmp = newOp[string2]
            newcoeff = tmp.coeff - coeffs2.coeff*(-1)^(xor(tmp.epsilon,coeffs2.epsilon))#, tmp.delta & coeffs2.delta))
            newOp[string2] = pauli_coeffs(newcoeff,tmp.delta, tmp.epsilon)
        else
            newOp[string2] = pauli_coeffs(coeffs2.coeff,coeffs2.delta, Bool(xor(1,coeffs2.epsilon)))
        end 
    end 
    return newOp
end 


function do_subtract!(op1,op2)
    for (string2,coeffs2) in op2
        if haskey(op1,string2)
            tmp = op1[string2]
            newcoeff = tmp.coeff - coeffs2.coeff*(-1)^(xor(tmp.epsilon,coeffs2.epsilon))
            op1[string2] = pauli_coeffs(newcoeff,tmp.delta, tmp.epsilon)
        else
            op1[string2] = pauli_coeffs(coeffs2.coeff,coeffs2.delta, Bool(xor(1,coeffs2.epsilon)))
        end 
    end 
end 

function do_comm(H,op)
    newOp ::Dict{Int64,pauli_coeffs} = Dict()
    # sizehint!(newOp.strings,2^(2*(L-1)))
    for (hstring,hcoeffs) in H
        hstringX = hstring>>>32
        hstringZ = hstring & (2^(32)-1)

        tmp_coeffs = hcoeffs.coeff
        tmp_delta = hcoeffs.delta
        tmp_epsilon = hcoeffs.epsilon
        
        for (ostring,ocoeffs) in op
            ostringX = ostring>>>32
            ostringZ = ostring & (2^(32)-1)
            dot1 = Bool(hammingWeight(hstringX & ostringZ)%2)
            dot2 = Bool(hammingWeight(hstringZ & ostringX)%2)
            if dot1 != dot2 # if they are equal, then they commute, don't add anything
                newcoeff = 2*tmp_coeffs*ocoeffs.coeff # 2*hcoeffs.coeff*ocoeffs.coeff
                newstring = xor(hstring,ostring)
                newdelta = xor(tmp_delta,ocoeffs.delta) # xor(hcoeffs.delta, ocoeffs.delta)
                newepsilon = xor(tmp_epsilon,ocoeffs.epsilon,
                                 tmp_delta & ocoeffs.delta,dot1)
                                 # xor(hcoeffs.epsilon,ocoeffs.epsilon, 
                                # hcoeffs.delta & ocoeffs.delta, dot1)
                if haskey(newOp,newstring)
                    tmp = newOp[newstring]
                    newnewcoeff = tmp.coeff + newcoeff*(-1)^xor(newepsilon,tmp.epsilon)
                    newOp[newstring] = pauli_coeffs(newnewcoeff, tmp.delta,tmp.epsilon)
                else
                    newOp[newstring] = pauli_coeffs(newcoeff, newdelta,newepsilon)
                end 
            end
        end 
    end 
    return newOp
end


function get_bs(L,j,jz,g,gamma,threshold,steps;start = "x")
    jx = j*(1 + gamma)/2
    jy = j*(1 - gamma)/2
    H =  make_H(L,jx,jy,jz,g)
    Psi_0 ::Dict{Int64,pauli_coeffs}  = Dict()
    if start == "y"
        Psi_0[Int64(1<<32 + 1)] = pauli_coeffs(1,1,0) # sig_y start
    elseif start == "z"
        Psi_0[Int64(1)] = pauli_coeffs(1,0,0) # sig_z start
    else
        Psi_0[Int64(1<<32)] = pauli_coeffs(1,0,0) # sig_x start, DEFAULT
    end
    # Psi_0 ::Dict{Int64,pauli_coeffs} = Dict(1<<32 => pauli_coeffs(1,0,0))
    bs = []
    result = do_comm(H,Psi_0) # b_1 |O_1)
    push!(bs,do_b(result)) # b_1
    norm_op!(result,bs[end],threshold) # |O_1)

    Gamma ::Dict{Int64,pauli_coeffs} = result # not copied over
    Xi ::Dict{Int64,pauli_coeffs} = Dict()
    if start == "y"
        Xi[Int64(1<<32 + 1)] = pauli_coeffs(bs[end],1,0) # xi = b_1 |O_0)
    elseif start == "z"
        Xi[Int64(1)] = pauli_coeffs(bs[end],0,0) # xi = b_1 |O_0)
    else
        Xi[Int64(1<<32)] = pauli_coeffs(bs[end],0,0) # xi = b_1 |O_0)
    end

    result = do_comm(H,Gamma)
    do_subtract!(result,Xi) # b_2 |O_2)
    push!(bs,do_b(result)) # b_2
    norm_op!(result,bs[end],threshold) # |O_2)
    do_mult!(Psi_0,bs[end]/bs[end-1])
    do_subtract!(result,Psi_0) # updated Psi_0
    Psi_0 = result

    for n in 3:2:steps
        result = do_comm(H,Psi_0)
        push!(bs,do_b(result)) #  b_{odd}
        norm_op!(result,bs[end],threshold) # |O_{odd})
        do_mult!(Gamma,bs[end]/bs[end-1])
        do_subtract!(result, Gamma)
        Gamma = result
        do_mult!(Xi,-bs[end]/bs[end-1])

        # tmp =  abs(Xi[1<<32].coeff)

        # clean_op!(Gamma,threshold*bs[end]/tmp)
        println("n = $n, length = $(length(Gamma))")

        result = do_comm(H,Gamma)
        do_subtract!(result,Xi)
        push!(bs,do_b(result))# b_{even}
        norm_op!(result,bs[end],threshold)
        do_mult!(Psi_0,bs[end]/bs[end-1])
        do_subtract!(result,Psi_0)
        Psi_0 = result

        # tmp =  abs(Xi[1<<32].coeff)
        # clean_op!(Psi_0,threshold*bs[end]/tmp)
        println("n = $(n+1), length = $(length(Psi_0))")
    end
    return bs
end
