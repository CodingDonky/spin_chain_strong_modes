using LinearAlgebra
using Dates
using PyPlot
using FileIO

using KrylovKit
using SparseArrays
using JLD
using Roots

using LaTeXStrings
using Plots

# Machine name so I can hardcode stuff for varios desktops / laptops if need be
machine = gethostname()

# Directory is print working directory
root_dir = pwd()
# This cuts off the directory to the githubs home directory
# at the outermost folder called "spin_chain_strong_modes"
root_dir = root_dir[1:findfirst("spin_chain_strong_modes",root_dir)[end]]

# Define other relevant folders for this repository
local_path_to_lifetimeprb2020 = root_dir*"/src/lifetime_prb2020-master/"
local_path_to_ED_lib = local_path_to_lifetimeprb2020*"ed_main/src"

push!(LOAD_PATH, local_path_to_ED_lib)
Base.load_path()
using ED


### Global Variables ###

# Pauli spin matrices
# id = [[1 0] ; [0 1]].+0im;
# σx = [[0 1] ; [1 0]].+0im;
# σy = [[0 -im] ; [im 0]];
# σz = [[1 0] ; [0 -1]].+0im;

# Pauli spin matrices
id = [1 0 ; 0 1].+0im;
σx = [0 1 ; 1 0].+0im;
σy = [0 -im ; im 0];
σz = [1 0 ; 0 -1].+0im;
# Pauli raising / lowering matrices
σp = σx+im.*σy # σ+
σm = σx-im.*σy # σ-

### Global Functions ###

function print_symbols()
    println("id, σx, σy, σz, ⊗(x,y)")
end

function commute(A, B)
    return A*B - B*A
end

# Tensor product shorthand
⊗(x,y) = kron(x,y);