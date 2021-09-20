push!(LOAD_PATH,pwd()*"/../src/")
using Test
using Random
using ED
using LinearAlgebra

Ls = [4,6]
tol = 1e-13
trials = 100
for (Li,L) in enumerate(Ls)
    for site in 1:L
        sigX = ED.sigX(L,site)
        for n in 1:trials
            vec = randn(Complex{Float64},2^L)
            @test maximum(abs.(sigX*vec - ED.sigX_V(L,site,vec)))<tol
        end
    end
end
