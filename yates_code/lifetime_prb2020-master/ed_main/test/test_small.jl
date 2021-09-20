push!(LOAD_PATH,pwd()*"/../src/")
using ED
using Test
using LinearAlgebra

L = 4
gs = range(-4,4,length = 31)
jx = 1.0
jzs = range(-2,2,length = 21)
tol = 1e-12

function spectrum_check(u1,u2,v1,v2)
    # check spectrum is the same
    @test maximum(abs.(u1 - u2)) < tol
    @test sum(abs.(u1 - u2)) < tol
    # check eigenvectors are the same
    # @test maximum(abs.([abs(v1[:,i]'*v2[:,i])^2-1 for i in 1:16])) < tol
    # @test sum(abs.([abs(v1[:,i]'*v2[:,i])^2-1 for i in 1:16])) < tol
end


for (gi,g) in enumerate(gs), (jzi,jz) in enumerate(jzs)
    # check obc
    ham1 = ED.H(L,jx,jz,g)
    ham2 = ED.smallH(jx,g,jz)
    u1,v1 = ED.diagonalize(ham1)
    u2,v2 = ED.diagonalize(ham2)
    spectrum_check(u1,u2,v1,v2)

    # now check periodic bc
    ham3 = ED.HP(L,jx,jz,g)
    ham4 = ED.smallH(jx,g,jz,1.0)
    u3,v3 = ED.diagonalize(ham3)
    u4,v4 = ED.diagonalize(ham4)
    spectrum_check(u3,u4,v3,v4)

end
