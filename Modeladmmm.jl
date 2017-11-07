using JuMP
using JuMPChance
using DataFrames
using Requests
using Gurobi
using Ipopt
using Gadfly
using Vega

function Dist_model( buses, lines)

    const e = 0.05
    const vo = 1
    const p = 1/vo^2
    const numiter = 150
    const pf = 0.9 #powe factor for solar power
    const phi = acos(pf)
    const epsilon = 0.05

#node location of PV resources
    const PV1 = 2
    const PV2 = 3
    const PV3 = 6
    const PV4 = 18
    const PV5 = 21
    const PV6 = 25
    const PV7 = 32

    numbuses = length(buses)
    numlines = length(lines)
    numPV = 7

  Pactive = readtable("/Users/alihassan/Desktop/lddfdata.csv") # Active power flows as auxilary constants

    ob2 = Float64[]
    qj_ob = Float64[]
    loss_ob1 = Float64[]
    loss_ob = Float64[]
    v_ob = Float64[]

#For Dual variables
l1= Array{Float64,1}(numbuses)
l2= Array{Float64,1}(numbuses)
l3= Array{Float64,1}(numbuses)
l4= Array{Float64,1}(numbuses)

    for i=1:numbuses
      l1[i] = 0
      l2[i] = 0
    end
    for i=1:numbuses
      l3[i] = 0
      l4[i] = 0
    end

for iter=1:numiter

m = ChanceModel(solver=IpoptSolver())

@indepnormal(m, gp[t=1:numbuses],mean=0.5,var=0.05) # Randon variable for active power
@variable(m, gq[t=1:numbuses] ) #variable for reactive power generation
@variable(m, fqo[t=1:numbuses] ) # variable for reactive power flow
@variable(m, fqn[t=1:numbuses] ) # variable for reactive power flow
@variable(m, vqo[t=1:numbuses] >= 0) #variable for voltage 
@variable(m, vqn[t=1:numbuses] >= 0) #variable for voltage 

@objective(m, Min, sum( p/2*( fqn[t] - sum(buses[z].fqglob for z in buses[t].children) )^2 + p/2*( fqo[t] - buses[t].fqglob )^2  
                        + l1[buses[t].nodeID]*( fqn[t] - sum(buses[z].fqglob for z in buses[t].children) ) + l2[buses[t].nodeID]*( fqo[t] - buses[t].fqglob ) for t=1:numbuses ) + sum(
                        buses[t].R * (fqo[t])^2/vo + 
                        + p/2*( vqn[t] - sum(buses[z].uqglob for z in buses[t].children) )^2  +  p/2*( vqo[t] - buses[t].uqglob )^2
                        + l3[t]*( vqn[t] - sum(buses[z].uqglob for z in buses[t].children) ) +  l4[t]*( vqo[t] - buses[t].uqglob ) for t=1:numbuses) ) #Objective function  

for i in 1:numbuses

  if  i==PV1||i==PV2||i==PV3||i==PV4||i==PV5||i==PV6||i==PV7 
@constraint(m, gp[i] <= buses[i].gp, with_probability=1-epsilon)
@constraint(m, gp[i] >= 0, with_probability=1-epsilon)
end
@constraint(m, vqn[i]  <= vo * (e^2+2*e) )
@constraint(m, vqn[i]  >= vo * (e^2-2*e) )

  if  i==PV1||i==PV2||i==PV3||i==PV4||i==PV5||i==PV6||i==PV7
@constraint(m, fqn[buses[i].nodeID] - fqo[buses[i].nodeID]  + buses[buses[i].nodeID].Qd - gp[i]*(1/tan(phi)) <=0 ,with_probability=1-epsilon)
@constraint(m, -fqn[buses[i].nodeID] + fqo[buses[i].nodeID]  - buses[buses[i].nodeID].Qd - gp[i]*(1/tan(phi)) <=0 ,with_probability=1-epsilon)
  else
@constraint(m, fqn[buses[i].nodeID] - fqo[buses[i].nodeID]  + buses[buses[i].nodeID].Qd  <=0 )
@constraint(m, -fqn[buses[i].nodeID] + fqo[buses[i].nodeID]  - buses[buses[i].nodeID].Qd  <=0 )
  end

#@constraint(m, fqn[buses[i].nodeID] - fqo[buses[i].nodeID]  + buses[buses[i].nodeID].Qd - sqrt((buses[i].smax)^2 - (buses[i].gp)^2) <=0 )
#@constraint(m, -fqn[buses[i].nodeID] + fqo[buses[i].nodeID]  - buses[buses[i].nodeID].Qd - sqrt((buses[i].smax)^2 - (buses[i].gp)^2) <=0 )
@constraint(m, vqo[i] == vqn[i] + 2*(buses[i].R*Pactive[i,2] + buses[i].X*fqo[i]) )
end


solve(m,method=:Reformulate)

#Averaging Step
for i=1:numbuses
  if  i==22||i==25||i==33||i==18 #Nodes with no children
    buses[i].fqglob = ( getvalue(fqo[i])  )
    buses[i].uqglob = ( (getvalue(vqo[i])) )
    elseif i==2||i==3||i==6
    buses[i].fqglob = 1/3*( getvalue(fqo[i]) + sum(getvalue(fqn[z]) for z in buses[i].children)  )
    buses[i].uqglob = 1/3*( (getvalue(vqo[i])) + sum(getvalue(vqn[z]) for z in buses[i].children)  )
    else
    buses[i].fqglob = 1/2*( getvalue(fqo[i]) + sum(getvalue(fqn[z]) for z in buses[i].children)  )
    buses[i].uqglob = 1/2*( (getvalue(vqo[i])) + sum(getvalue(vqn[z]) for z in buses[i].children)  )
  end
end


#Lagrange multipliers update step
for i=1:numbuses
  if  i==22||i==25||i==33||i==18 #Nodes with no children
       l1[i] = l1[i] 
       l3[i] = l3[i] 
  else
       l1[i] = l1[i] + p*( sum(getvalue(fqn[z]) for z in buses[i].children ) - sum(buses[z].fqglob for z in buses[i].children) )
       l3[i] = l3[i] + p*( sum(getvalue(vqn[z]) for z in buses[i].children) - sum(buses[z].uqglob for z in buses[i].children) )
  end
       l2[i] = l2[i] + p*( getvalue(fqo[i]) - buses[i].fqglob )
       l4[i] = l4[i] + p*( getvalue(vqo[i]) - (buses[i].uqglob))  
end

  push!(qj_ob, sum( getvalue(fqn[x]) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x in 1:numbuses ) *Sbase )
 # push!(loss_ob, sum(  lines[e].r* abs( getvalue(fqo[lines[e].head]) ) for e in 1:numlines)*Sbase ) 
  #push!(loss_ob, sum(  buses[e].R* ( getvalue(fqo[e]) )^2 for e in 1:numbuses)*Sbase ) 
    push!(loss_ob, sum(  buses[e].R* ( buses[e].fqglob )^2 for e in 1:numbuses)*Sbase*Sbase) 

println("Obj = ",getobjectivevalue(m))
for i=1:numbuses
println("qj $i = ",( getvalue(fqn[i]) - getvalue(fqo[i])  +  buses[buses[i].nodeID].Qd )*Sbase)
println("v $i  = ", vo+buses[i].uqglob  )
println("")
end
println("-------------")
println("qj totsl = ", sum( getvalue(fqn[i]) - getvalue(fqo[i])  +  buses[buses[i].nodeID].Qd for i=1:numbuses )*Sbase)
println("-------------")
for i=2:numbuses
println("From ", buses[buses[i].nodeID].ancestor[1], " to ", buses[i].nodeID) 
println("fqgolbal $(i-1) = ", buses[i].fqglob *Sbase) 
println("fqo $(i-1) = ", getvalue(fqo[i]))
println("fqn $(i-1) = ", getvalue(fqn[i]))
println("")
end

if iter==numiter
 for e=1:numbuses
  push!(ob2, ( getvalue(fqn[e]) - getvalue(fqo[e])  +  buses[buses[e].nodeID].Qd ) *Sbase )
  push!(v_ob, vo+(buses[e].uqglob)  )
end
end

end
display(plot(layer(x=1:numbuses, y=v_ob, Geom.line), Guide.XLabel("bsu"), Guide.YLabel("V"),Guide.Title("V $(epsilon)")))
display(plot(layer(x=1:numiter, y=loss_ob, Geom.line), Guide.XLabel("Iterations"), Guide.YLabel("losses (MW)"),Guide.Title("Total losses $(epsilon)")))
display(plot(layer(x=1:numiter, y=qj_ob, Geom.line), Guide.XLabel("Iterations."), Guide.YLabel("qj (MVar)"),Guide.Title("Chance constraints with epsilon $(epsilon)")))

 for e=1:numbuses
println("q at $e = ", ob2[e])
println("")
end



#println("finish")
end 