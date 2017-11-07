using JuMP
using DataFrames
using Requests

include("Input.jl")

#include("Modeladmmm.jl")

include("Modeladmmdeep.jl")


#for i=1:33
#	println("parent for bus $i = ", buses[i].ancestor)
    #println("tail for line $i = ", lines[i].tail)
   
 # println("")
#end

#=for i=1:33
	println("bus $i = ", buses[i].R*Sbase)	
end
=#
Dist_model(buses, lines)
println("END")