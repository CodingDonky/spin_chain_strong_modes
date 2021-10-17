# This library will hold the ED functions.
# This code was based off of the helpful notes by Prof. Sandvik
# https://physics.bu.edu/~sandvik
# https://physics.bu.edu/~sandvik/perimeter/


# Flip the bits of integer a at bit sites i and j
function bitflip(a::Int64,i::Int64,j::Int64)::Int64
    return xor(a, 1<<i + 1<<j)
end

# Flip the bits of integer a at bit site i
function bitflip(a::Int64,i::Int64)::Int64
    return xor(a, 1<<i)
end

# get sigma_z at site i, does no "flipping"
function bitsign(a::Int64,i::Int64)::Int64
    return (1-2*((a>>i)&1))
end

# sum the first L (least) bits of an Int64
function bitsum(x::Int64,L::Int64)::Int64
    num ::Int64 = 0
    for i in 0:(L-1)
        num += (x>>i) & 1
    end
    return num
end

# Reverse the L rightmost bits of x
function bitreverse(x::Int64,L::Int64)::Int64
    ans ::Int64= 0
    for i in 0:(L-1)
        ans += ((x>>i) & 1) << (L-1-i)
    end
    return ans
end

# cycle the bits to the "right" by one.
function bitcycle(x::Int64,L::Int64)::Int64
    return (x<<1)&(2^L-1) + (x>>(L-1))
end


#=
Locate (the index of) integer x in an ordered list.. .
Julia probably has a better version of this..
I tested it, findfirst() sees to actually be slower.
=#
function locateState(x::Int64, list::Array{Int64,1})::Int64
    bmin = 1
    bmax = length(list)
    while bmin <=  bmax
        b = bmin+ floor(Int, (bmax - bmin)/2)
        if list[b]<x
            bmin = b +1
        elseif list[b]>x
            bmax = b -1
        else
            return b
        end
    end
    return -1
end


#=
This will find the representative states for a given
sub-block in the spin parity, mirror symmetry basis.

We specify mirror symmetry because in cases like (L even)
p = 0, m = 1, it will not be the same set of states as
p = 0, m = 0, as the latter have states like 0110,
1001, 0000 (for L = 4), that vanish in the former.
=#
function findBasis(L::Int64,p::Int64,m::Int64)
    ans :: Array{Int64,1} = []
    if m == 0
        for a in 0:(2^L-1)
            pt ::Int64 = bitsum(a,L)%2
            bt ::Int64 = bitreverse(a,L)
            if p == pt && a<=bt # note the <=
                push!(ans,a)
            end
        end
    else
        for a in 0:(2^L-1)
            pt ::Int64 = bitsum(a,L)%2
            bt ::Int64 = bitreverse(a,L)
            if p == pt && a<bt # note the <
                push!(ans,a)
            end
        end
    end
    return ans
end


#=
Find periodicity of state.. just for inversion this
is simple, its either 1 or 2
(0110) -> R = 1, one mirror flip to return to itself
(0011) -> R = 2, two mirror flips to return to itself
=#
function findPeriod(x::Int64,L::Int64)::Int64
    y = bitreverse(x,L)
    if x != y
        return 2
    else
        return 1
    end
end


#=
Reverse bitstring->y, if y<x then
y is the rep. If x<y then x is the rep.
For sign flips returning 0,1 depending on if
x or y is returned
=#
function findRep(x::Int64,L::Int64)
    y = bitreverse(x,L)
    if y < x
        return y, 1
    else
        return x, 0
    end
end

