module ScenarioReduction
using DataFrames
export HourlyLoad, HourlyLoadSel, Prob
export ScenarioSelection

function ScenarioSelection(HourlyLoad1::Matrix, NScenarios::Int64)
    LoadActualLen = size(HourlyLoad1)[1]
    LoadLen = size(HourlyLoad1)[1] ; indsel = zeros(Int,NScenarios);
    Prob =fill(1/LoadActualLen,NScenarios) ; Cp = Array{Float64}(NScenarios)
    SelectedSubset = zeros(NScenarios,24)
    C_Old = fill(100000, LoadLen, LoadLen);
    for l in 1 : NScenarios
        C = Array{Float64}(LoadLen,LoadLen)
        Z = Array{Float64}(LoadLen)
        for k in 1 : LoadLen, u in 1 : LoadLen
            if l == 1
               k<=u ? C[k,u] = norm(HourlyLoad1[k,:]-HourlyLoad1[u,:], 2) : C[k,u] = C[u,k]
            else
               k<=u ? C[k,u] = min(C_Old[k,u],C_Old[k,indsel[l-1]]) : C[k,u] = min(C_Old[u,k],C_Old[k,indsel[l-1]])
            end
        end
        Z = sum(C,1)
        C_Old = C
        SelectedSubset[l,:]=convert(Array{Float64},HourlyLoad1[indmin(Z),:])
        HourlyLoad1 = HourlyLoad1[1:size(HourlyLoad1,1).!= indmin(Z),: ]
        LoadLen = size(HourlyLoad1)[1]
        indsel[l] = Int(indmin(Z))
    end
    for k in 1: LoadLen
        for u in 1:size(SelectedSubset)[1]
            Cp[u] = norm(HourlyLoad1[k,:]-SelectedSubset[u,:], 2)
        end
        Prob[indmin(Cp)] = Prob[indmin(Cp)] +1/LoadActualLen
    end
    return(SelectedSubset, indsel, Prob)
end
end
