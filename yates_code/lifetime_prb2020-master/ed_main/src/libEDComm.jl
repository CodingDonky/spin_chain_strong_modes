#=
This is where the ED commutator routines will go for now.
=#


comm_op_dot(L,op1,op2) =  1/(2^L)*tr(op1'*op2)
comm_get_b(L,op) = sqrt(comm_op_dot(L,op,op))
comm_get_comm(H,op) = H*op - op*H
comm_get_A(H,bn1,op1,op2) = comm_get_comm(H,op1) - bn1.*op2
comm_set_zero(val,threshold) = abs(val) < threshold ?  0.0 : val

function comm_check_Ops(L,ops)
    threshold = 1e-5
    @printf("n = %3d: ", length(ops))
    for i in 1:length(ops)
        tmp1 = abs(comm_op_dot(L,ops[i],ops[end]))^2
        @printf("%2.1e ",tmp1)
    end 
end 

function comm_get_bs(L,j,jz,g,gamma,threshold,steps;start = "x")
    #=
    This will calculate the b_n's in the usual ED 
    setup. This is full of round-off errors and numerical stability
    issues. Also, \sig_1^x breaks inversion and is 
    off-diagonal in parity.. At the moment it is not worth 
    optimizing this since we have another way of calculating these
    quantities.
    =#
    H = H_Delta(L,j,jz,g,gamma)
    if start == "y"
        Os = [im .*(sigZ(L,1)*sigX(L,1))]
    elseif start == "z"
        Os = [sigZ(L,1)]
    else 
        Os = [sigX(L,1)]
    end
    bs = []
    tmp = comm_get_comm(H,Os[1])
    # tmp = comm_set_zero.(tmp,threshold)
    println("")
    println("")
    println("")
    println(size(tmp))
    println("")
    println("")
    println("")
    tmpB = comm_get_b(L,tmp)
    if abs(tmpB) > threshold
        push!(Os,tmp ./tmpB)
        push!(bs,tmpB)
    else
        return
    end 

    for n in 1:steps-1
        tmpA = comm_get_comm(H,Os[end]) - bs[end]*Os[end-1]
        # tmpA = comm_set_zero.(tmpA,threshold)
        tmpB = comm_get_b(L,tmpA)
        if abs(tmpB)>threshold
            push!(Os,tmpA ./ tmpB)
            push!(bs,tmpB)
        else
            println("\n breaking out...\n")
            break
        end
        
        # reorthogonalize
        tmp_new = zeros(Complex{Float64},2^L,2^L)
        for m in length(Os)-1:-1:1
            tmp_new += Os[m] .*comm_op_dot(L,Os[m],Os[end])
        end 
        Os[end] -= tmp_new
        norm_tmp = sqrt.(comm_op_dot(L,Os[end],Os[end]))
        Os[end] = Os[end] ./ norm_tmp

        # comm_check_Ops(L,Os)
        @printf("b_%d:  %2.1e \r",n+1,abs(tmpB))
    end
    return bs, Os
end 


function comm_get_bs_cutoff(L,j,jz,g,gamma,threshold,steps,cutoff)
    #=
    This will calculate the b_n's in the usual ED 
    setup. This is full of round-off errors and numerical stability
    issues. Also, \sig_1^x breaks inversion and is 
    off-diagonal in parity.. At the moment it is not worth 
    optimizing this since we have another way of calculating these
    quantities.
    =#
    H = H_Delta(L,j,jz,g,gamma)
    Os = [sigX(L,1)]
    # Os = [im .*sigZ(L,1)*sigX(L,1)]
    bs = []
    tmp = comm_get_comm(H,Os[1])
    tmp = comm_set_zero.(tmp,threshold)
    tmpB = comm_get_b(L,tmp)
    if abs(tmpB) > threshold
        push!(Os,tmp ./tmpB)
        push!(bs,tmpB)
    else
        return bs, Os
    end 
    for n in 1:steps-1
        tmpA = comm_get_comm(H,Os[end]) - bs[end]*Os[end-1]
        tmpA = comm_set_zero.(tmpA,threshold)
        tmpB = comm_get_b(L,tmpA)
        if abs(tmpB)>threshold
            push!(Os,tmpA ./ tmpB)
            push!(bs,tmpB)
        else
            println("\n breaking out...\n")
            break
        end
        # reorthogonalize
        tmp_new = zeros(Complex{Float64},2^L,2^L)
        for m in length(Os)-1:-1:1
            tmp_new += Os[m] .*comm_op_dot(L,Os[m],Os[end])
        end 
        Os[end] -= tmp_new
        norm_tmp = sqrt.(comm_op_dot(L,Os[end],Os[end]))
        Os[end] = Os[end] ./ norm_tmp
        # comm_check_Ops(L,Os)
        @printf("b_%d:  %2.1e \r",n+1,tmpB)
        if bs[end] > cutoff
            break
        end
    end
    return bs, Os
end 

