module DCOPF
include("ScenarioReduction.jl")
include("NetworkBuild.jl")
using JuMP
using Gurobi
using .ScenarioReduction
using .NetworkBuild

export OPF

function OPF(Network, LoadH, NScenarios)
    BusLength = parse(Int, Network["BusLen"])
    Line_conc = Array(Vector{Int64}, BusLength)
    Line_x = Array(Vector{Float64}, BusLength)
    gen = Array(Vector{Int64}, BusLength)
    Loaddata = zeros(BusLength,NScenarios^length(Network["load"]),24)
    Exp_DAPrice = zeros(BusLength,NScenarios^length(Network["load"]),24)
    for k in 1:NScenarios
        [Loaddata[Network["load"]["1"]["bus"],NScen+(k-1)*NScenarios,t] = LoadH[1][1][NScen,t] for NScen in 1:NScenarios, t in 1:24]
        [Loaddata[Network["load"]["2"]["bus"],k+(NScen-1)*NScenarios,t] = LoadH[2][1][NScen,t] for NScen in 1:NScenarios, t in 1:24]
    end
    genInd = []
    for j in 1:BusLength
        Line_conc[j] = []
        Line_x[j] = []
        gen[j] = []
    end
    for i in 1:length(Network["branch"])
        push!(Line_conc[Network["branch"]["$i"]["fbus"]], Network["branch"]["$i"]["tbus"]) # adjacent node
        push!(Line_conc[Network["branch"]["$i"]["tbus"]], Network["branch"]["$i"]["fbus"])

        push!(Line_x[Network["branch"]["$i"]["fbus"]], Network["branch"]["$i"]["x"]) # edge reactance of connected node
        push!(Line_x[Network["branch"]["$i"]["tbus"]], Network["branch"]["$i"]["x"])
    end
    [push!(gen[Network["gen"]["$i"]["bus"]], Network["gen"]["$i"]["gen_ind"]) for i in 1:length(Network["gen"])] # each bus gen index
    [push!(genInd, length(gen[i])) for i in 1:BusLength]
    genMaxInd = maximum(genInd)

    for NScen in 1: NScenarios^length(Network["load"])
      for t in 1:24
          DAOPF = Model(solver = GurobiSolver(MIPGap = 1e-4, TimeLimit = 800, OutputFlag=0))
          @defVar(DAOPF, pgen[1:BusLength, 1:genMaxInd]>=0)
          @defVar(DAOPF, -400 <= delta[1:BusLength]<=400)
          @setObjective(DAOPF, Min, sum(Network["gen"]["$i"]["a"]*pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]^2+
                        Network["gen"]["$i"]["b"]*pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]] for i in 1:length(Network["gen"])))
          @addConstraint(DAOPF, Check_CSTR[i in 1:BusLength, j in 1:genMaxInd], pgen[i,j]*floor((genInd[i]-1)/100)>=0)
          @addConstraint(DAOPF, PB_CSTR[i in 1:BusLength], sum(pgen[i,j] for j in 1:genInd[i]) - sum(1/Line_x[i][j]*(delta[i] - delta[Line_conc[i][j]])
          for j in 1:length(Line_conc[i])) - Loaddata[i,NScen,t]==0)
              @addConstraint(DAOPF, delta[1]==0)
              @addConstraint(DAOPF, GenLimMin_CSTR[i in 1:length(Network["gen"])], pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]
                            >=Network["gen"]["$i"]["pmin"])
              @addConstraint(DAOPF, GenLimMax_CSTR[i in 1:length(Network["gen"])], pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]
                            <=Network["gen"]["$i"]["pmax"])
              @addConstraint(DAOPF, LineFlowF_CSTR[i in 1:length(Network["branch"])], 1/Network["branch"]["$i"]["x"]*
                            (delta[Network["branch"]["$i"]["fbus"]] - delta[Network["branch"]["$i"]["tbus"]]) <= Network["branch"]["$i"]["rateA"])
              @addConstraint(DAOPF, LineFlowB_CSTR[i in 1:length(Network["branch"])], 1/Network["branch"]["$i"]["x"]*
                            (delta[Network["branch"]["$i"]["fbus"]] - delta[Network["branch"]["$i"]["tbus"]]) >= -Network["branch"]["$i"]["rateA"])
              solve(DAOPF)
              [Exp_DAPrice[i,NScen,t] = getdual(PB_CSTR[i]) for i in 1:BusLength]
          end
      end
      return (Exp_DAPrice)
    end
end
