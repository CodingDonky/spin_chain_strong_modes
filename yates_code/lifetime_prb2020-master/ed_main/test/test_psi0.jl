push!(LOAD_PATH,pwd()*"/../src/")
using Test
using ED
using LinearAlgebra

Ls = [6,8]
gs = range(-.4,.4,length = 41)
jx = 1.0
tol = 1e-13

for (Li,L) in enumerate(Ls), (gi,g) in enumerate(gs)

    ham = ED.H(L,jx,0.0,g)
    psi0 = ED.Psi0(L,g)
    # @printf("g^L = %.2e, max of comm = %.2e ",g^L,maximum(abs.(ham*psi0 - psi0*ham)))
    tmp = psi0*psi0
    tmp2 = maximum(Diagonal(tmp))
    @test abs(tmp2 - 1.0) < tol
end
