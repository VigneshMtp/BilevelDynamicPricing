#Top Layer
module DayAhead
include("ScenarioReduction.jl")
include("NetworkBuild.jl")
include("DCOPF.jl")

using JuMP
using Gurobi
using DataFrames
using Distributions

export LoadH, NScenarios, E_DAQuantity, E_DAPrice, Exp_DAprice, RTprice, BusLength, TLN_Scenarios, Nloadbus, loadbus
NScenarios = 5

LoadH = []
for i in 1:2
    if i ==1
        LoadData = convert(DataFrame, readtable("PJM_DAY.csv", header = true))
    else
        LoadData = convert(DataFrame, readtable("PJM_AE.csv", header = true))
    end
    HourlyLoad = (LoadData[1:100, 3:26])
    HourlyLoadSel,index,Prob = ScenarioReduction.ScenarioSelection(convert(Array, HourlyLoad), NScenarios)
    push!(LoadH, tuple(HourlyLoadSel,index,Prob))
end

Network = NetworkBuild.ReadMatData("NetworkData.m")
BusLength = parse(Int, Network["BusLen"])
E_DAQuantity = zeros(24,BusLength)
E_DAPrice = zeros(24,BusLength)
TLN_Scenarios = NScenarios^length(Network["load"])
Scen_Prob = zeros(TLN_Scenarios)
Loaddata = zeros(BusLength,TLN_Scenarios,24)
for k in 1:NScenarios
    [Loaddata[Network["load"]["1"]["bus"],NScen+(k-1)*NScenarios,t] = LoadH[1][1][NScen,t] for NScen in 1:NScenarios, t in 1:24]
    [Loaddata[Network["load"]["2"]["bus"],k+(NScen-1)*NScenarios,t] = LoadH[2][1][NScen,t] for NScen in 1:NScenarios, t in 1:24]
    [Scen_Prob[NScen+(k-1)*NScenarios] = LoadH[1][3][NScen]*LoadH[2][3][k] for NScen in 1:NScenarios]
end

Exp_DAprice = DCOPF.OPF(Network, LoadH, NScenarios)

function RTPrice(Exp_DAprice,spike)
    RTprice = zeros(BusLength,TLN_Scenarios,24)
    for i in 1:BusLength
        for NScen in 1: TLN_Scenarios
            for t in 1: 24
                M = ceil(rand()-(1-spike[t]))
                RTprice[i,NScen,t] = Exp_DAprice[i,NScen,t]*(1+(1-M)*rand(Normal(0,0.1),1)[1]+M*rand(Cauchy(0.25,0.05),1))[1]
            end
        end
    end
    return RTprice
end
SpikesDA = [0	0	0	0	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.25	0.25	0.25	0.25	0.25	0.25	0.1	0.1	0.1	0	0	0]
Exp_RTprice = RTPrice(Exp_DAprice,SpikesDA)
Line_conc = Array(Vector{Int64}, BusLength)
Line_x = Array(Vector{Float64}, BusLength)
gen = Array(Vector{Int64}, BusLength)

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

loadbus = Vector{Int64}(length(Network["load"]));
Nloadbus = Vector{Int64}(BusLength);
[loadbus[i] = Network["load"]["$i"]["bus"] for i in 1:length(Network["load"])]
[findfirst(loadbus,i)!= 0 ? Nloadbus[i] = 0 : Nloadbus[i] = i for i in 1:BusLength]
filter!(x->x!=0,Nloadbus)


for t in 1:24
  DAOPF = Model(solver = GurobiSolver(MIPGap = 1e-4, TimeLimit = 800, OutputFlag=0))
  @defVar(DAOPF, pgen[1:BusLength, 1:genMaxInd]>=0)
  @defVar(DAOPF, pgen_P[1:BusLength, 1:genMaxInd, 1:TLN_Scenarios]>=0)
  @defVar(DAOPF, pgen_N[1:BusLength, 1:genMaxInd, 1:TLN_Scenarios]>=0)
  @defVar(DAOPF, E_DA[1:BusLength]>=0)
  @defVar(DAOPF, E_DA_P[1:BusLength, 1:TLN_Scenarios]>=0)
  @defVar(DAOPF, E_DA_N[1:BusLength, 1:TLN_Scenarios]>=0)
  @defVar(DAOPF, -400 <= delta[1:BusLength]<=400)
  @defVar(DAOPF, -400 <= delta_w[1:BusLength, 1:TLN_Scenarios]<=400)
  @setObjective(DAOPF, Min, sum(Network["gen"]["$i"]["a"]*pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]^2+
                Network["gen"]["$i"]["b"]*pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]] for i in 1:length(Network["gen"]))
                +sum(Scen_Prob[NScen]*sum(Exp_RTprice[Network["gen"]["$i"]["bus"],NScen,t]*
                (pgen_P[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen]
                -pgen_N[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen]) for i in 1:length(Network["gen"])) for NScen in 1:TLN_Scenarios))

  # First Stage Constraints
  @addConstraint(DAOPF, Check_CSTR[i in 1:BusLength, j in 1:genMaxInd], pgen[i,j]*floor((genInd[i]-1)/100)>=0)
  @addConstraint(DAOPF, PB_CSTR[i in 1:BusLength], sum(pgen[i,j] for j in 1:genInd[i])-
                sum(1/Line_x[i][j]*(delta[i] - delta[Line_conc[i][j]]) for j in 1:length(Line_conc[i])) - E_DA[i]==0)
  @addConstraint(DAOPF, delta[1]==0)
  @addConstraint(DAOPF, GenLimMin_CSTR[i in 1:length(Network["gen"])], pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]
                >=Network["gen"]["$i"]["pmin"])
  @addConstraint(DAOPF, GenLimMax_CSTR[i in 1:length(Network["gen"])], pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]
                <=Network["gen"]["$i"]["pmax"])
  @addConstraint(DAOPF, LineFlowF_CSTR[i in 1:length(Network["branch"])], 1/Network["branch"]["$i"]["x"]*
                (delta[Network["branch"]["$i"]["fbus"]] - delta[Network["branch"]["$i"]["tbus"]]) <= Network["branch"]["$i"]["rateA"])
  @addConstraint(DAOPF, LineFlowB_CSTR[i in 1:length(Network["branch"])], 1/Network["branch"]["$i"]["x"]*
                (delta[Network["branch"]["$i"]["fbus"]] - delta[Network["branch"]["$i"]["tbus"]]) >= -Network["branch"]["$i"]["rateA"])

  # Second Stage Constraints
  @addConstraint(DAOPF, Check_CSTRII_P[i in 1:BusLength, j in 1:genMaxInd, NScen in 1:TLN_Scenarios],
                pgen_P[i,j,NScen]*floor((genInd[i]-1)/100)>=0)
  @addConstraint(DAOPF, Check_CSTRII_N[i in 1:BusLength, j in 1:genMaxInd, NScen in 1:TLN_Scenarios],
                pgen_N[i,j,NScen]*floor((genInd[i]-1)/100)>=0)

  @addConstraint(DAOPF, PB_CSTRII[i in 1:BusLength, NScen in 1:TLN_Scenarios], sum(pgen_P[i,j,NScen] - pgen_N[i,j,NScen]  for j in 1:genInd[i])
                -sum(1/Line_x[i][j]*(delta_w[i,NScen] - delta_w[Line_conc[i][j],NScen] -delta[i]+delta[Line_conc[i][j]]) for j in 1:length(Line_conc[i]))
                -(E_DA_P[i, NScen]-E_DA_N[i, NScen]) ==0)

  @addConstraint(DAOPF, GenLimMin_CSTRII[i in 1:length(Network["gen"]), NScen in 1:TLN_Scenarios],
                pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]+pgen_P[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen]
                -pgen_N[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen] >=Network["gen"]["$i"]["pmin"])

  @addConstraint(DAOPF, GenLimMax_CSTRII[i in 1:length(Network["gen"]), NScen in 1:TLN_Scenarios],
                pgen[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"]]+pgen_P[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen]
                -pgen_N[Network["gen"]["$i"]["bus"],Network["gen"]["$i"]["gen_ind"],NScen] <=Network["gen"]["$i"]["pmax"])
  @addConstraint(DAOPF, LineFlowF_CSTRII[i in 1:length(Network["branch"]), NScen in 1:TLN_Scenarios], 1/Network["branch"]["$i"]["x"]*
                (delta_w[Network["branch"]["$i"]["fbus"], NScen] - delta_w[Network["branch"]["$i"]["tbus"], NScen]) <= Network["branch"]["$i"]["rateA"])
  @addConstraint(DAOPF, LineFlowB_CSTRII[i in 1:length(Network["branch"]), NScen in 1:TLN_Scenarios], 1/Network["branch"]["$i"]["x"]*
                (delta_w[Network["branch"]["$i"]["fbus"], NScen] - delta_w[Network["branch"]["$i"]["tbus"], NScen]) >= -Network["branch"]["$i"]["rateA"])
  @addConstraint(DAOPF, DeltaConst[NScen in 1:TLN_Scenarios], delta_w[1,NScen]==0)
  @addConstraint(DAOPF, [i in 1:length(Network["load"]), NScen in 1:TLN_Scenarios],E_DA[Network["load"]["$i"]["bus"]]+
                E_DA_P[Network["load"]["$i"]["bus"],NScen]-E_DA_N[Network["load"]["$i"]["bus"],NScen] == Loaddata[Network["load"]["$i"]["bus"],NScen,t])

  @addConstraint(DAOPF, LoadZero[i in 1:length(Nloadbus), NScen in 1:TLN_Scenarios], E_DA_N[Nloadbus[i],NScen] == 0)
  @addConstraint(DAOPF, E_DA_Min[i in 1:length(Network["load"])], E_DA[Network["load"]["$i"]["bus"]] >= minimum(Loaddata[Network["load"]["$i"]["bus"],:,t]))
  @addConstraint(DAOPF, E_DA_Max[i in 1:length(Network["load"])], E_DA[Network["load"]["$i"]["bus"]] <= maximum(Loaddata[Network["load"]["$i"]["bus"],:,t]))
  solve(DAOPF)
  [E_DAQuantity[t,i] = getvalue(E_DA[i]) for i in 1:BusLength]
    [E_DAPrice[t,i] = getdual(PB_CSTR[i]) for i in 1:BusLength]
end
end
