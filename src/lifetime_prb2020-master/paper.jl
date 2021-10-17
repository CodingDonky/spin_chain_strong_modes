#=
This file will generate the plots for the paper. 
Currently, it does not requrie the ED library
(load it anyway).
Will add some basic ED examples later on.
REMEMBER TO UNCOMMENT RELEVANT PARTS OF THE CODE BELOW
REMEMBER THERE IS A DIFFERENCE BETWEEN RUNNING THE CODE AS
command line: julia paper.jl
VS
GOING TO JULIA REPL
julia> include("paper.jl").
SOME THINGS DO NOT RUN IN THE FORMER. BETTER TO USE Julia repl WHEN INTERACTING WITH PYPLOT.
HOWEVER IF ALL THAT THE CODE IS DOING IS SAVING PLOTS RATHER THAN SHOWING IT THEN RUNNING
IT IN command line is fine.
This file makes the figures exactly as they appear in the 
paper. This file uses Julia's PyPlot package to leverage
Python's matplotlib plotting library. Further, we use 
matplotlib's ability to generate LaTeX for the axes and 
legends. 

If there are errors when running this code it can be 3 things
1.) Julia error, this should be short error message
2.) Matplotlib error, this will be a longer message.
For example we pass the variable 
"marker_size = 10" into the PyPlot (Julia) functions. 
They will accept the code and translate it into Python and then 
matplotlib's plotting library will throw an error because it 
does not know what "marker_size" is (instead want "markersize").
3.) LaTeX error, this should spit out a huge error message that 
should look immediately familiar to anyone who works with LaTeX often.
For example, if one wants the y label to be "$|\psi|^2$", you need to 
escape the backslash and the $ symbols need backslashes, like this 
"\$|\\psi|^2\$".


This is where PyPlot is hosted
https://github.com/JuliaPy/PyPlot.jl

This link will help with making sure Python knows where to find 
LaTeX on your system.
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html

=#

# we add /ed_main/src/ to the LOAD_PATH variable so
# julia knows about our code there.
push!(LOAD_PATH,pwd()*"/ed_main/src/")

# this entire file is a module named "paper"
module paper

# loading packages
using LinearAlgebra
using PyPlot
using ED
# using Roots
# # check with get_bs.jl and gen_bs.jl...
# push!(LOAD_PATH,pwd()*"/../ED/src/")
using KrylovKit
using SparseArrays
# This bug might have been fixed, but as of 2019-2020
# one should load JLD last among the packages
using JLD


# Now we will "include" a bunch of files, this function effectively
# copy-and-pastes the contents of the files into this module.
# Completely equivalent to writing all the code in a single file,
# but more readable this way. 


# util funcs------------------------------------------------------
# utility functions are imported into the current namespace
include("paper_util.jl")

# in paper--------------------------------------------------------
# fig1 function is imported here
include("paper_fig1.jl")
# fig2 function is imported here
include("paper_fig2.jl")
# fig3 function is imported here
include("paper_fig3.jl")
# fig4_5 function is imported here
include("paper_fig4_5.jl")
# supp1 function is imported here
include("paper_supp1.jl")
# supp2_run and supp2 functions are imported here
include("paper_supp2.jl")

# extras-----------------------------------------------------------
# supp3 and supp4 functions are imported here
include("paper_supp3_4.jl")
# simple_bn function is imported here
include("paper_simple_bn.jl")
# simple_ainf function is imported here 
include("paper_simple_ainf.jl")
# simple_kyrlov function is imported here
include("paper_simple_krylov.jl")

####################################################################
# running the above functions, uncomment what is to be ran #########
####################################################################
# in paper----------------------------------------------------------
# @timev fig1()
# @timev fig2()
# @timev fig3()
# @timev fig4_5()
# @timev supp1()
# generates krylov time data, run if data is not present
# @timev supp2_run() 
# @timev supp2()

# extras, does not save any figures---------------------------------
# @timev supp3()
# @timev supp4()

# test functions----------------------------------------------------
 @timev simple_bn()
 @timev simple_ainf()
 @timev simple_krylov()

end
