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
    const numiter = 250 #number of iterations
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
l1= Array{Float64,1}(numbuses) # for flow from parent to node at parent
l2= Array{Float64,1}(numbuses) # for flow from parent to the node at node
l3= Array{Float64,1}(numbuses) # voltage at self
l4= Array{Float64,1}(numbuses) # voltage for parent

    for i=1:numbuses
      l1[i] = 0
      l2[i] = 0
      l3[i] = 0
      l4[i] = 0
    end

for iter=1:numiter

m = ChanceModel(solver=IpoptSolver())

@indepnormal(m, gp[t=1:numbuses],mean=0.5,var=0.05) # Randon variable for active power
@variable(m, gq[t=1:numbuses] ) #variable for reactive power generation

############# Local Variables ################
@variable(m, fqo[t=1:numbuses] ) # flow from parent to the node at node
@variable(m, fqn[t=1:numbuses] ) # flow to ndoe at parent
@variable(m, vqo[t=1:numbuses] >= 0) # voltage from parent (tail)
@variable(m, vqn[t=1:numbuses] >= 0) #variable to children (head)
##############################################

############ Objective Function (Augmented Lagrangian)##################
@objective(m, Min, sum(  buses[t].R * (fqo[t])^2/vo
                        +p/2*sum((fqn[z] - buses[z].fqglob)^2 for z in buses[t].children) + p/2*( fqo[t] - buses[t].fqglob )^2
                        + p/2*( vqn[t] - buses[t].uqglob )^2  +  p/2*( vqo[t] - sum(buses[z].uqglob for z in buses[t].ancestor))^2
                        + sum(l1[buses[z].nodeID]*(fqn[z] - buses[z].fqglob) for z in buses[t].children) + l2[buses[t].nodeID]*(fqo[t] - buses[t].fqglob )
                        + l3[t]*( vqn[t] - buses[t].uqglob) +  l4[t]*( vqo[t] - sum(buses[z].uqglob for z in buses[t].ancestor)) for t=2:numbuses)
                        + p/2*sum((fqn[z] - buses[z].fqglob)^2 for z in buses[1].children) + p/2*(fqo[1] - buses[1].fqglob)^2
                        + p/2*(vqn[1] - buses[1].uqglob )^2  + sum(l1[buses[z].nodeID]*(fqn[z] - buses[z].fqglob) for z in buses[1].children)
                        + l2[buses[1].nodeID]*(fqo[1] - buses[1].fqglob) + l3[1]*(vqn[1] - buses[1].uqglob) ) #Objective function

#Here buses[t].fqglob and buses[z].uqglob represents the global variables at the node for flow into and voltage at the node
############ Minimization Step ################
for i in 1:numbuses

  if  i==PV1||i==PV2||i==PV3||i==PV4||i==PV5||i==PV6||i==PV7  #nodes with PV
@constraint(m, gp[i] <= buses[i].gp, with_probability=1-epsilon)
@constraint(m, gp[i] >= 0, with_probability=1-epsilon)
end
@constraint(m, vqn[i]  <= vo * (e^2+2*e) )
@constraint(m, vqn[i]  >= vo * (e^2-2*e) )

  if  i==PV1||i==PV2||i==PV3||i==PV4||i==PV5||i==PV6||i==PV7

@constraint(m, sum(fqn[z] for z in buses[i].children) - fqo[buses[i].nodeID]  + buses[buses[i].nodeID].Qd - gp[i]*(1/tan(phi)) <=0 ,with_probability=1-epsilon)
@constraint(m, -sum(fqn[z] for z in buses[i].children) + fqo[buses[i].nodeID]  - buses[buses[i].nodeID].Qd - gp[i]*(1/tan(phi)) <=0 ,with_probability=1-epsilon)
  else
@constraint(m, sum(fqn[z] for z in buses[i].children) - fqo[buses[i].nodeID]  + buses[buses[i].nodeID].Qd  <=0 )
@constraint(m, -sum(fqn[z] for z in buses[i].children) + fqo[buses[i].nodeID]  - buses[buses[i].nodeID].Qd  <=0 )
  end

  if i!=1
@constraint(m, vqo[i] == vqn[i] + 2*(buses[i].R*Pactive[i,2] + buses[i].X*fqo[i]) )
  end

end

###########################################

solve(m,method=:Reformulate)

############Update Step ###################
for i=1:numbuses
  if  i==22||i==25||i==33||i==18 #Nodes with no children
    buses[i].fqglob =  1/2( getvalue(fqo[i])+getvalue(fqn[i]) )
    buses[i].uqglob =  (getvalue(vqn[i]))

    elseif i==1
    buses[i].fqglob =  (getvalue(fqo[i]) )
    buses[i].uqglob =  1/2*( (getvalue(vqn[i])) + sum(getvalue(vqo[z]) for z in buses[i].children)  )

    elseif i==2||i==3||i==6 #Nodes with 2 children (A bit of hardcoding for now) and parent => that's why divide by 3
    buses[i].fqglob = 1/2( getvalue(fqo[i])+getvalue(fqn[i]) )
    buses[i].uqglob = 1/3*( (getvalue(vqn[i])) + sum(getvalue(vqo[z]) for z in buses[i].children)  )

    else #Nodes with a children and parent
    buses[i].fqglob = 1/2( getvalue(fqo[i])+getvalue(fqn[i]) )
    buses[i].uqglob = 1/2*( (getvalue(vqn[i])) + sum(getvalue(vqo[z]) for z in buses[i].children)  )
  end
end
############################################

##################Lagrange multipliers update step #################
for i=1:numbuses
  if  i==1
       l1[i] = l1[i]
       l2[i] = l2[i] + p*( getvalue(fqo[i]) - buses[i].fqglob )
       l3[i] = l3[i] + p*( getvalue(vqn[i]) - buses[i].uqglob )
       l4[i] = l4[i]
  else

       l1[i] = l1[i] + p*(getvalue(fqn[i]) - buses[i].fqglob)
       l2[i] = l2[i] + p*( getvalue(fqo[i]) - buses[i].fqglob )
       l3[i] = l3[i] + p*( getvalue(vqn[i]) - buses[i].uqglob )
       l4[i] = l4[i] + p*( getvalue(vqo[i]) - sum(buses[z].uqglob for z in buses[i].ancestor))
  end

end
######################################


  push!(qj_ob, sum( sum(getvalue(fqn[z]) for z in buses[buses[x].nodeID].children) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x in 1:17 )
    + sum( (getvalue(fqn[x]) ) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x =18 ) #node with no children
    + sum( sum(getvalue(fqn[z]) for z in buses[buses[x].nodeID].children) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x in 19:21 )
    + sum( (getvalue(fqn[x]) ) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x =22 ) #node with no children
    + sum( sum(getvalue(fqn[z]) for z in buses[buses[x].nodeID].children) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x in 23:24 )
    + sum( (getvalue(fqn[x]) ) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x =25 )  #node with no children
    + sum( sum(getvalue(fqn[z]) for z in buses[buses[x].nodeID].children) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x in 26:32)
    + sum( (getvalue(fqn[x]) ) - getvalue(fqo[x])  +  buses[buses[x].nodeID].Qd for x =33 )  #node with no children
   *Sbase ) # for plotting the reactive power injections
  push!(loss_ob, sum(  buses[e].R* ( buses[e].fqglob )^2 for e in 1:numbuses)*Sbase*Sbase) # for plotting the losses

println("Obj = ",getobjectivevalue(m))

println("-------------")
println("qj total = ", sum( getvalue(fqn[i]) - getvalue(fqo[i])  +  buses[buses[i].nodeID].Qd for i=1:numbuses )*Sbase)
println("-------------")


if iter==numiter #for steady state value (last iteration)
 for e=1:numbuses
#  push!(ob2, ( getvalue(fqn[e]) - getvalue(fqo[e])  +  buses[buses[e].nodeID].Qd ) *Sbase )
  push!(v_ob, vo+(buses[e].uqglob)  ) # for plotting voltage
end
end

end
######### For test plotting ########
display(plot(layer(x=1:numbuses, y=v_ob, Geom.line), Guide.XLabel("bsu"), Guide.YLabel("V"),Guide.Title("V $(epsilon)")))
display(plot(layer(x=1:numiter, y=loss_ob, Geom.line), Guide.XLabel("Iterations"), Guide.YLabel("losses (MW)"),Guide.Title("Total losses $(epsilon)")))
display(plot(layer(x=1:numiter, y=qj_ob, Geom.line), Guide.XLabel("Iterations."), Guide.YLabel("qj (MVar)"),Guide.Title("Chance constraints with epsilon $(epsilon)")))

#=for e=1:numbuses
println("q at $e = ", ob2[e])
println("")
end=#

println("finish")
end
