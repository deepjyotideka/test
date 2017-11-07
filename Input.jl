using JuMP
using DataFrames
using Requests
using Ipopt

const Sbase = 100
Vbase = 250
Zbase = (Vbase^2)/Sbase

type Bus
   nodeID::Int
   Pd::Float64
   Qd::Float64
   gp::Float64
   smax::Float64
   R::Float64
   X::Float64
   children::Vector{Int}         
   ancestor::Vector{Int}  
   fqglob::Float64  
   uqglob::Float64       
   function Bus(nodeID, Pd, Qd, R, X, gp, smax, fqglob, uqglob)  
      b = new(nodeID, Pd, Qd)
      b.gp = gp
      b.smax = smax
      b.R = R
      b.X = X
      b.children = Int[]
      b.ancestor = Int[]
      b.fqglob = fqglob
      b.uqglob = uqglob
      return b
   end
end

##################################################################
type Line
   arcID::Int
   head::Int # the "from" node
   tail::Int # the "to" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   uqglob::Float64
   function Line(arcID, head, tail, r, x, uqglob)
      line = new(arcID, head, tail, r, x, uqglob)
      return line
   end
end
#########################################
# Bus/Node Data
busmat = readtable("/NodeData.csv") # Path for the file


  buses = Bus[]
  busIDmap = Dict()
for i in 1:size(busmat,1)

	  nodeID = i
      busIDmap[busmat[i,1]] = i
     
      Pd = busmat[i,2] /Sbase
      Qd = busmat[i,3] /Sbase
      R = busmat[i,6] /Sbase
      X = busmat[i,7] /Sbase

      fqglob = 0

    uqglob = 0
      
      gp = busmat[i,4] 
      smax = busmat[i,5]


      b = Bus(nodeID, Pd, Qd, R, X, gp, smax, fqglob, uqglob)
      push!(buses, b)

end

numbuses = length(buses)

   ######################################################################################################
   ## branch data
   branchmat = readtable("/LineData.csv") # Path for the file

   lines = Line[]
   for i in 1:size(branchmat,1)
      fbus = busIDmap[branchmat[i,2]]
      tbus = busIDmap[branchmat[i,1]]
      abus = busIDmap[branchmat[i,1]]
      x = branchmat[i,4]/Sbase
      r = branchmat[i,3] /Sbase
            uqglob = 0
          #  if i ==1
           #         push!(buses[fbus].ancestor, fbus)#ancestor
            #      end
      push!(buses[fbus].children, abus)#children
      push!(buses[tbus].ancestor, fbus)#ancestor
      l = Line(i, tbus, fbus, r, x, uqglob)  
      push!(lines,l)
   end
#########################################