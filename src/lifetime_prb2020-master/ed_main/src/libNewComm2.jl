#=
This is where commutator code goes.
This is an update to libNewComm.. will see if its faster.
=#

mutable struct myOp
    #=
    The default dictionary may not be the best option..
    We definitely want to be smart about this..
    =#
    L::Int64
    coeffs::Dict{Int64,Float64}
    deltas::Dict{Int64,Bool}
    epsilons::Dict{Int64,Bool}
end 

function print_help(string::Int64)
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

function print_op(op::myOp,num::Int64)
    n = min(length(op.coeffs),num)
    for (i,key) in enumerate(keys(op.coeffs))
        coeff = op.coeffs[key]
        delta = op.deltas[key]
        epsilon= op.epsilons[key]
        pretty_string = print_help(key)
        @printf("% 2.1e, %d, %d, %s %d \n",
                coeff,
                delta,
                epsilon,
                pretty_string,
                key)
        if i == num
            break
        end
    end 
    println("")
end

function make_H(L,jx,jy,jz,g)
    H = myOp(L,Dict(),Dict(),Dict())
    for i in 0:(L-2)
        z_pos  = Int64(1<<i)
        z_bond = Int64(3<<i)
        x_bond = Int64(3<<(i+32))
        y_bond = Int64(x_bond+z_bond)
        if jx != 0
            H.coeffs[x_bond] =  jx
            H.deltas[x_bond] =  0
            H.epsilons[x_bond] =  0
        end 
        if jy != 0
            H.coeffs[y_bond] =  jy
            H.deltas[y_bond] =  0
            H.epsilons[y_bond] =  1
        end 
        if jz != 0 
            H.coeffs[z_bond] =  jz
            H.deltas[z_bond] =  0
            H.epsilons[z_bond] =  0
        end 
        if g != 0
            H.coeffs[z_bond] =  g
            H.deltas[z_bond] =  0
            H.epsilons[z_bond] =  0
        end 
    end 
    if g != 0
        H.coeffs[Int64(1<<(L-1))] = g
        H.deltas[Int64(1<<(L-1))] = 0
        H.epsilons[Int64(1<<(L-1))] = 0
    end 
    return H
end


function hammingWeight(x::Int64)
    # finding the number of 1's in a binary representation
    # the obvious way of doing it (iterating through the string)
    # might be too slow for me..
    # wikipedia has a good article on this..
    # https://en.wikipedia.org/wiki/Hamming_weight
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

function do_b(op)
    return sqrt(sum(abs.(values(op.coeffs)).^2))
end

function do_mult(op,scalar)
    # newOp = myOp(op.L,deepcopy(op.coeffs),deepcopy(op.deltas),deepcopy(op.epsilons))
    newOp = myOp(op.L,copy(op.coeffs),copy(op.deltas),copy(op.epsilons))
    for (key,value) in newOp.coeffs
        newOp.coeffs[key] = value*scalar
    end 
    return newOp
end 

function do_subtract(op1,op2)
    # newOp = myOp(op1.L,deepcopy(op1.coeffs),deepcopy(op1.deltas),deepcopy(op1.epsilons))
    newOp = myOp(op1.L,copy(op1.coeffs),copy(op1.deltas),copy(op1.epsilons))
    for (key2, value2) in op2.coeffs
        if haskey(newOp.coeffs,key2)
            tmp = newOp.coeffs[key2] - value2*(-1)^(xor(newOp.epsilons[key2],op2.epsilons[key2]))
            newOp.coeffs[key2] = tmp
            # if abs(tmp) >threshold
            #     newOp.coeffs[key2] = tmp
            # else
            #     delete!(newOp.coeffs,key2)
            #     delete!(newOp.deltas,key2)
            #     delete!(newOp.epsilons,key2)
            # end
        else
            newOp.coeffs[key2] = value2
            newOp.deltas[key2] = op2.deltas[key2]
            newOp.epsilons[key2] = Bool(xor(1,op2.epsilons[key2]))
        end 
    end 
    return newOp
end 


function do_comm(H,op)
    L = H.L
    newOp = myOp(L,Dict(),Dict(),Dict()) 
    for (hstring,h_coeffs) in H.coeffs
        hstringX = hstring>>>32
        hstringZ = hstring & (2^(32)-1)
        for (ostring,o_coeffs) in op.coeffs
            ostringX = ostring>>>32
            ostringZ = ostring & (2^(32)-1)
            dot1 = Bool(hammingWeight(hstringX & ostringZ)%2)
            dot2 = Bool(hammingWeight(hstringZ & ostringX)%2)
            if dot1 != dot2 # if they are equal, then they commute, don't add anything
                new_coeff = 2*h_coeffs*o_coeffs
                new_string = xor(hstring,ostring)
                new_delta = xor(H.deltas[hstring],op.deltas[ostring])
                new_epsilon = xor(H.epsilons[hstring],op.epsilons[ostring],
                    H.deltas[hstring] & op.deltas[ostring],
                    dot1)
                if haskey(newOp.coeffs,new_string)
                    tmp = ( newOp.coeffs[new_string] + new_coeff*(-1)^xor(new_epsilon,newOp.epsilons[new_string]) )
                    newOp.coeffs[new_string] = tmp
                    # if abs(tmp)>threshold
                    #     newOp.coeffs[new_string] = tmp
                    # else
                    #     delete!(newOp.coeffs,new_string)
                    #     delete!(newOp.deltas,new_string)
                    #     delete!(newOp.epsilons,new_string)
                    # end
                else
                    newOp.coeffs[new_string] = new_coeff
                    newOp.deltas[new_string] = new_delta
                    newOp.epsilons[new_string] = new_epsilon
                end 
            end
        end 
    end 
    return newOp
end

function clean_op!(op,threshold)
    for (string, coeff) in op.coeffs
        if abs(coeff)<threshold
            delete!(op.coeffs,string)
            delete!(op.deltas,string)
            delete!(op.epsilons,string)
        end
    end
end

function get_bs(L,j,jz,g,gamma,threshold,steps)
    jx = j*(1 + gamma)/2
    jy = j*(1 - gamma)/2
    H =  make_H(L,jx,jy,jz,g)
    Ops = []
    bs = []
    sigX = myOp(L,Dict(),Dict(),Dict())
    sigX.coeffs[(1<<32)] = 1
    sigX.deltas[(1<<32)] = 0
    sigX.epsilons[(1<<32)] = 0
    push!(Ops,sigX)
    op2 = do_comm(H,Ops[1])
    b = do_b(op2)
    if abs(b) > threshold
        op2 = do_mult(op2,1/b)
        push!(Ops,op2)
        push!(bs,b)
    else
        return
    end 

    for n in 1:steps
        tmpA = do_comm(H,Ops[end])
        tmpB = do_mult(Ops[end-1],bs[end])
        tmpC = do_subtract(tmpA,tmpB)
        tmp_b = do_b(tmpC)
        if abs(tmp_b)>threshold
            tmpD = do_mult(tmpC,1/tmp_b)
            deleteat!(Ops,1)
            clean_op!(Ops[end],threshold)
            push!(Ops,tmpD)
            push!(bs,tmp_b)
        else
            println("\n breaking out...\n")
            break
        end
        println("n = $n, length = $(length(Ops[end].coeffs))")
    end
    return bs,Ops
end


