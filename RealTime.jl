#Bottom Layer
include("DayAhead.jl")
include("ElectricVehicle.jl")
include("InputDataFile.jl")
using JuMP
using Gurobi
using Distributions
using .DayAhead
using .EVmod

#Declarations
N_Bus = DayAhead.BusLength
loadbus = DayAhead.loadbus
N_loadbus = length(DayAhead.loadbus)
DA_Price = DayAhead.E_DAPrice
DA_Quantity = DayAhead.E_DAQuantity

OutputFile  = Array{Any}(N_loadbus)
for i in 1:N_loadbus
    OutputFile[i] = open("ModelOutput$i.csv", "w")
    print(OutputFile[i],"Time\tDA_Quantity\tDA_Price\tDynamicPrice\tTotalConsumption\tPLConsumption\t")
    print(OutputFile[i],"ResidentConsumer\tRealTimePrice\tRealTimePriceUB\tRealTimePriceLB\tSOC\tGamma\tSpikes")
    print(OutputFile[i],"\n")
    close(OutputFile[i])
end

PLMax = zeros(N_loadbus,24)
[PLMax[i,t] = min(300,sum(EVmod.CEV[i][1,k].Time[t] for k in 1:length(EVmod.CEV[i][1,:]))) for i in 1:N_loadbus, t in 1:24]

Spikes = [0	0	0	0	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.25	0.25	0.25	0.25	0.25	0.25	0.1	0.1	0.1	0	0	0];
Gamma = [0.5 0.3 0]

# Creating RT-price distributions " 25*25 = 625 price for each bus and at each time
N_Replicaions = 25
M = 10000
E = zeros(24,N_Bus); P_cnsr = zeros(24,N_loadbus); DynPrice = zeros(N_Bus,24)
E,P_cnsr,DynPrice = InputDataFile.readData(E,P_cnsr,DynPrice,loadbus)

RT_priceDist = zeros(N_Bus, N_Replicaions*DayAhead.TLN_Scenarios, 24)
LB_RTprice = zeros(N_Bus, 24); UB_RTprice = zeros(N_Bus, 24)
for i in 1 : N_Replicaions
    RT_priceDist[:,(i-1)*DayAhead.TLN_Scenarios+1:i*DayAhead.TLN_Scenarios,:] = DayAhead.RTPrice(DayAhead.Exp_DAprice,Spikes)
end

[LB_RTprice[i,t] = quantile(RT_priceDist[i,:,t],0.1) for i in 1:N_Bus, t in 1:24]
[UB_RTprice[i,t] = quantile(RT_priceDist[i,:,t],0.9) for i in 1:N_Bus, t in 1:24]

CRT_price = zeros(24,N_Bus)
[CRT_price[t,i] = DA_Price[t,i]*(1+(1-ceil(rand()-(1-Spikes[t])))*rand(Normal(0,0.1),1)[1]+ceil(rand()-(1-Spikes[t]))*rand(Cauchy(0.25,0.05),1))[1] for t in 1:24, i in 1:N_Bus]

for GammaInd in 1: 3
  for adhoc in 1:2    #1 - optimal #2 - adhoc
    DynPrice = zeros(N_Bus,24)
    prevSOC_cpl=zeros(N_loadbus,EVmod.N_cpl)
    # KKT Problem
    for tau in 1:23
      SOprob = Model(solver = GurobiSolver(MIPGap = 1e-4, TimeLimit = 800, OutputFlag=0))
      @defVar(SOprob, zeta[1:N_Bus]>=0)
      @defVar(SOprob, NL_L[1:N_loadbus]>=0)
      @defVar(SOprob, DP[1:N_Bus, 1:24]>=0)
      @defVar(SOprob, Z[1:N_Bus]>=0)
      @defVar(SOprob, eta[1:N_Bus, 1:24]>=0)
      @defVar(SOprob, yrb[1:N_Bus, 1:24]>=0)
      @defVar(SOprob, P_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, Epl[i in 1:N_loadbus, t in 1:24]>=0)
      @defVar(SOprob, SOC_cpl[i in 1:N_loadbus, t in 0:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, rho_pl[i in 1:N_loadbus, 1:24]>=0)
      @defVar(SOprob, Lmu_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, Umu_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, Lnu_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, Unu_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, alpha_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl])
      @defVar(SOprob, gamma_cpl[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl]>=0)
      @defVar(SOprob, y1[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl],Bin)
      @defVar(SOprob, y2[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl],Bin)
      @defVar(SOprob, y3[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl],Bin)
      @defVar(SOprob, y4[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl],Bin)
      @defVar(SOprob, y5[i in 1:N_loadbus, t in 1:24, c in 1:EVmod.N_cpl],Bin)

      if adhoc ==1
        @setObjective(SOprob, Min, sum((DA_Price[tau,loadbus[i]])*DA_Quantity[tau,loadbus[i]] +
        CRT_price[tau,loadbus[i]]*((Epl[i,tau])-DA_Quantity[tau,loadbus[i]])+sum(DA_Price[t,loadbus[i]]*DA_Quantity[t,loadbus[i]]+
        0.5*(UB_RTprice[loadbus[i],t]+LB_RTprice[loadbus[i],t])*((Epl[i,t])-DA_Quantity[t,loadbus[i]])+ eta[loadbus[i],t] for t in tau+1:24)
        + round(Gamma[GammaInd]*(24-tau))*Z[loadbus[i]] for i in 1:N_loadbus))
      else
        @setObjective(SOprob, Min, 0)
      end
      #Upper Level Constraints
      @addConstraint(SOprob, DP_CSTR[i in 1:N_loadbus], (DA_Price[tau,loadbus[i]])*DA_Quantity[tau,loadbus[i]] -
      NL_L[i]  +
      CRT_price[tau,loadbus[i]]*((Epl[i,tau])-DA_Quantity[tau,loadbus[i]]) <=0)

      if adhoc == 1
        @addConstraint(SOprob, NonL_Lin[i in 1:N_loadbus], NL_L[i] >= (sum(rho_pl[i,t]*P_cnsr[t,i] + sum( - 11.5*PLMax[i,t]/1000*Umu_cpl[i,t,c] -
        sum(EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].Time[t] for k in 1:length(EVmod.CEV[i][c,:]))/1000*Unu_cpl[i,t,c] +
        (sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[h]/1000 for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t)-
        sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[h]/1000 for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t))*gamma_cpl[i,t,c] for c in 1:EVmod.N_cpl)
        for t in tau:24)+
        sum((sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:])) -
        sum(EVmod.CEV[i][c,k].beta*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:])))*alpha_cpl[i,t,c]
        for c in 1:EVmod.N_cpl,  t in tau+1:24)+
        sum((prevSOC_cpl[i,c]+sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:])) -
        sum(EVmod.CEV[i][c,k].beta*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:])))*alpha_cpl[i,t,c]
        for c in 1:EVmod.N_cpl,  t in tau)) -sum(DA_Price[t,loadbus[i]]*Epl[i,t] for t in tau+1:24))
      end
    @addConstraint(SOprob, Robust_CSTR[i in 1:N_loadbus, t in tau+1:24], Z[loadbus[i]]+eta[loadbus[i],t]
                   >=0.5*(UB_RTprice[loadbus[i],t]-LB_RTprice[loadbus[i],t])*yrb[loadbus[i],t])
    @addConstraint(SOprob,Robust_CSTR1[i in 1:N_loadbus, t in tau+1:24],yrb[loadbus[i],t]>=((Epl[i,t])-DA_Quantity[t,loadbus[i]]))
    @addConstraint(SOprob,Robust_CSTR1[i in 1:N_loadbus, t in tau+1:24],yrb[loadbus[i],t]>=-((Epl[i,t])-DA_Quantity[t,loadbus[i]]))
    #Lower Level Constraints
    @addConstraint(SOprob, PL[i in 1:N_loadbus, t in tau:24], Epl[i,t] == P_cnsr[t,i]+sum(P_cpl[i,t,c] for c in 1:EVmod.N_cpl))
    @addConstraint(SOprob, SOC_prev[i in 1:N_loadbus, c in 1:EVmod.N_cpl],SOC_cpl[i,tau-1,c] == prevSOC_cpl[i,c])
    @addConstraint(SOprob, SOC_CPL_t[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], SOC_cpl[i,t,c] == SOC_cpl[i,t-1,c]+
                   sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:]))-
                   sum(EVmod.CEV[i][c,k].beta*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[t]/1000 for k in 1:length(EVmod.CEV[i][c,:]))+
                   0.95*P_cpl[i,t,c])
    # KKT formulation for follower problem
    @addConstraint(SOprob, P_pl_DP[i in 1:N_loadbus], rho_pl[i,tau]<=DP[loadbus[i],tau])
    @addConstraint(SOprob, P_pl_EA[i in 1:N_loadbus, t in tau+1:24], rho_pl[i,t]<=DA_Price[t,loadbus[i]])
    @addConstraint(SOprob, Pcpl[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], - rho_pl[i,t] + Lmu_cpl[i,t,c]-Umu_cpl[i,t,c] - 0.95*alpha_cpl[i,t,c] <=0)
    @addConstraint(SOprob, SOCcpl[i in 1:N_loadbus, t in tau:23, c in 1:EVmod.N_cpl], +Lnu_cpl[i,t,c]-Unu_cpl[i,t,c]+alpha_cpl[i,t,c]-alpha_cpl[i,t+1,c]+gamma_cpl[i,t,c]<=0)
    @addConstraint(SOprob, SOCcpl2[i in 1:N_loadbus, c in 1:EVmod.N_cpl], +Lnu_cpl[i,24,c] -Unu_cpl[i,24,c]+alpha_cpl[i,24,c]+gamma_cpl[i,24,c]<=0)
    # Complementary Conditions
    @addConstraint(SOprob, CCSTR_y1[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], Lmu_cpl[i,t,c] <= M*(y1[i,t,c]))
    @addConstraint(SOprob, CCSTR_y12[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], P_cpl[i,t,c] <= M*(1-y1[i,t,c]))
    @addConstraint(SOprob, CCSTR_y2[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], Umu_cpl[i,t,c] <= M*(y2[i,t,c]))
    @addConstraint(SOprob, CCSTR_y22[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], 11.5*PLMax[i,t]/1000 - P_cpl[i,t,c]<= M*(1-y2[i,t,c]))
    @addConstraint(SOprob, CCSTR_y22_1[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], 11.5*PLMax[i,t]/1000 - P_cpl[i,t,c]>= 0)
    @addConstraint(SOprob, CCSTR_y3[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], Lnu_cpl[i,t,c] <= M*(y3[i,t,c]))
    @addConstraint(SOprob, CCSTR_y32[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], SOC_cpl[i,t,c] <= M*(1-y3[i,t,c]))
    @addConstraint(SOprob, CCSTR_y4[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], Unu_cpl[i,t,c] <= M*(y4[i,t,c]))
    @addConstraint(SOprob, CCSTR_y42[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], sum(EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].Time[t] for k in 1:length(EVmod.CEV[i][c,:]))/1000-SOC_cpl[i,t,c] <= M*(1-y4[i,t,c]))
    @addConstraint(SOprob, CCSTR_y42_2[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], sum(EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].Time[t] for k in 1:length(EVmod.CEV[i][c,:]))/1000-SOC_cpl[i,t,c] >=0)
    @addConstraint(SOprob, CCSTR_y5[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], gamma_cpl[i,t,c] <= M*(y5[i,t,c]))
    @addConstraint(SOprob, CCSTR_y52[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], SOC_cpl[i,t,c] -
                   (sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[h] for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t)-
                   sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[h] for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t))/1000<=M*(1-y5[i,t,c]))
    @addConstraint(SOprob, CCSTR_y52_2[i in 1:N_loadbus, t in tau:24, c in 1:EVmod.N_cpl], SOC_cpl[i,t,c] -
                   (sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolArrtime[h] for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t)-
                   sum(EVmod.CEV[i][c,k].alpha*EVmod.CEV[i][c,k].BC*EVmod.CEV[i][c,k].BoolDeptime[h] for k in 1:length(EVmod.CEV[i][c,:]), h in 1:t))/1000>=0)
    solve(SOprob)
    prevSOC_cpl = getvalue(SOC_cpl[:,tau,:])
    [DynPrice[loadbus[i],tau] = getvalue(NL_L[i])/(getvalue(Epl[i,tau])) for i in 1:N_loadbus]
    for i in 1:N_loadbus
      OutputFile[i] = open("ModelOutput$i.csv", "a")
      print(OutputFile[i],tau,"\t",round(DA_Quantity[tau,loadbus[i]],2),"\t",round(DA_Price[tau,loadbus[i]],2),"\t",round(DynPrice[loadbus[i],tau],2),"\t")
      print(OutputFile[i],round(getvalue(Epl[i,tau]),2),"\t",round(sum(getvalue(P_cpl[i,tau,:])),2),"\t",round(P_cnsr[tau,i],2),"\t")
      print(OutputFile[i],round(CRT_price[tau,loadbus[i]],2),"\t",round(UB_RTprice[loadbus[i],tau],2),"\t",round(LB_RTprice[loadbus[i],tau],2),"\t")
      print(OutputFile[i],sum(getvalue(SOC_cpl[i,tau,:])),"\t",Gamma[GammaInd],"\t",Spikes[tau],"\t")
      print(OutputFile[i],"\n")
      close(OutputFile[i])
    end
  end
end   #Adhoc
for i in 1:N_loadbus
    OutputFile[i] = open("ModelOutput$i.csv", "a")
    print(OutputFile[i],"24\t",round(DA_Quantity[24,loadbus[i]],2),"\t",round(DA_Price[24,loadbus[i]],2),"\t\t\t")
    print(OutputFile[i],"\t\t\t\t",round(P_cnsr[24,i],2),"\t")
    print(OutputFile[i],round(CRT_price[24,loadbus[i]],2),"\t",round(UB_RTprice[loadbus[i],24],2),"\t",round(LB_RTprice[loadbus[i],24],2),"\t")
    print(OutputFile[i],"\t\t",Gamma[GammaInd],"\t",Spikes[24],"\t")
    print(OutputFile[i],"\n")
    close(OutputFile[i])
end
end
