#=
This is where commutator code goes.
=#

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

mutable struct myOp
    #=
    The default dictionary may not be the best option..
    We definitely want to be smart about this..
    The keys in the strings dictionary are the pauli strings
    themselves, I store their information as the two halfs of
    Int64. The values of the string dictionary are the coefficients,
    signs, and phase factors, contained in a puali_coeffs instance.
    =#
    L::Int64
    strings::Dict{Int64,pauli_coeffs}
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

function print_op(op::myOp,num::Int64)
    #=
    Prints the operator.
    =#
    n = min(length(op.strings),num)
    for (i,key) in enumerate(keys(op.strings))
        coeffs = op.strings[key]
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
    H = myOp(L,Dict())
    for i in 0:(L-2)
        z_pos  = Int64(1<<i)
        z_bond = Int64(3<<i)
        x_bond = Int64(3<<(i+32))
        y_bond = Int64(x_bond+z_bond)
        if jx != 0
            H.strings[x_bond] =  pauli_coeffs(jx,0,0)
        end 
        if jy != 0
            H.strings[y_bond] = pauli_coeffs(jy,0,1)
        end 
        if jz != 0 
            H.strings[z_bond] = pauli_coeffs(jz,0,0)
        end 
        if g != 0
            H.strings[z_pos ] = pauli_coeffs(g ,0,0)
        end 
    end 
    if g != 0
        H.strings[Int64(1<<(L-1))] = pauli_coeffs(g,0,0) 
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
    tmp = 0
    for (string,coeffs) in op.strings
        tmp += abs(coeffs.coeff)^2
    end 
    return sqrt(tmp)
end

function do_mult(op,scalar)
    newOp = myOp(op.L,copy(op.strings))
    for (string,coeffs) in newOp.strings
        # coeffs.coeff = coeffs.coeff*scalar
        newOp.strings[string] = pauli_coeffs(coeffs.coeff*scalar,coeffs.delta,coeffs.epsilon)
    end 
    return newOp
end 

function clean_op!(op,threshold)
    for (string, coeffs) in op.strings
        if abs(coeffs.coeff)<threshold
            delete!(op.strings,string)
        end
    end
end

function do_subtract(op1,op2)
    newOp = myOp(op1.L,deepcopy(op1.strings))
    for (string2,coeffs2) in op2.strings
        if haskey(newOp.strings,string2)
            tmp = newOp.strings[string2]
            # if tmp.delta != coeffs2.delta
            #     println("somethings wrong")
            # end
            newcoeff = tmp.coeff - coeffs2.coeff*(-1)^(xor(tmp.epsilon,coeffs2.epsilon))#, tmp.delta & coeffs2.delta))
            # if abs(newcoeff)> 1e-12
            #     newOp.strings[string2].coeff = newcoeff
            # else 
            #     delete!(newOp.strings,string2)
            #     # println("too small")
            # end

            # newOp.strings[string2].coeff = newcoeff
            newOp.strings[string2] = pauli_coeffs(newcoeff,tmp.delta, tmp.epsilon)
        else
            newOp.strings[string2] = pauli_coeffs(coeffs2.coeff,coeffs2.delta, Bool(xor(1,coeffs2.epsilon)))
        end 
    end 
    return newOp
end 


function do_comm(H,op)
    L = H.L
    newOp = myOp(L,Dict()) 
    # sizehint!(newOp.strings,2^(2*(L-1)))
    for (hstring,hcoeffs) in H.strings
        hstringX = hstring>>>32
        hstringZ = hstring & (2^(32)-1)
        for (ostring,ocoeffs) in op.strings
            ostringX = ostring>>>32
            ostringZ = ostring & (2^(32)-1)
            dot1 = Bool(hammingWeight(hstringX & ostringZ)%2)
            dot2 = Bool(hammingWeight(hstringZ & ostringX)%2)
            if dot1 != dot2 # if they are equal, then they commute, don't add anything
                newcoeff = 2*hcoeffs.coeff*ocoeffs.coeff
                newstring = xor(hstring,ostring)
                newdelta = xor(hcoeffs.delta, ocoeffs.delta)
                newepsilon = xor(hcoeffs.epsilon,ocoeffs.epsilon, 
                                hcoeffs.delta & ocoeffs.delta, dot1)

                if haskey(newOp.strings,newstring)
                    tmp = newOp.strings[newstring]
                    # if newdelta != tmp.delta
                    #     println("we have a problem.... nd = $newdelta, td = $(tmp.delta)")
                    #     # the following line will trigger an error...
                    #     tmp.coeff += newcoeff*(-1)^(newepsilon + tmp.epsilon)*im^(newdelta + tmp.delta)
                    #     # its not suppose to make sense, I want it to crash.. lazy
                    #     tmp.coeff += 2 + 3.0*im
                    # else    
                    #     newnewcoeff = tmp.coeff + newcoeff*(-1)^xor(newepsilon,tmp.epsilon)
                    #     if abs(newnewcoeff)>1e-15
                    #         newOp.strings[newstring].coeff = newnewcoeff
                    #     else
                    #         delete!(newOp.strings,newstring)
                    #         # println("too small 2")
                    #     end
                    # end 
                    newnewcoeff = tmp.coeff + newcoeff*(-1)^xor(newepsilon,tmp.epsilon)

                    # newOp.strings[newstring].coeff = newnewcoeff
                    newOp.strings[newstring] = pauli_coeffs(newnewcoeff, tmp.delta,tmp.epsilon)
                else
                    newOp.strings[newstring] = pauli_coeffs(newcoeff, newdelta,newepsilon)
                end 
            end
        end 
    end 
    return newOp
end


function get_bs(L,j,jz,g,gamma,threshold,steps)
    jx = j*(1 + gamma)/2
    jy = j*(1 - gamma)/2
    H =  make_H(L,jx,jy,jz,g)
    Ops = []
    bs = []
    sigX = myOp(L,Dict())
    sigX.strings[(1<<32)] =  pauli_coeffs(1,0,0)
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
            clean_op!(tmpD,threshold)
            push!(Ops,tmpD)
            push!(bs,tmp_b)
        else
            println("\n breaking out...\n")
            break
        end
        println("n = $n, length = $(length(Ops[end].strings))")
    end
    return bs,Ops
end



function get_bs_cutoff(L,j,jz,g,gamma,threshold,steps,cutoff)
    jx = j*(1 + gamma)/2
    jy = j*(1 - gamma)/2
    H =  make_H(L,jx,jy,jz,g)
    Ops = []
    bs = []
    sigX = myOp(L,Dict())
    sigX.strings[(1<<32)] =  pauli_coeffs(1,0,0)
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
            clean_op!(tmpD,threshold)
            push!(Ops,tmpD)
            push!(bs,tmp_b)
        else
            println("\n breaking out...\n")
            break
        end
        println("n = $n, length = $(length(Ops[end].strings))")
        if tmp_b > cutoff
            break
        end
    end
    return bs,Ops
end
