module EVmod
export CEV, SEV, N_cpl, N_spl

type EV
    Type::String
    Model:: String
    BC::Float64
    alpha::Float64
    Time::Array{Float64,1}
    beta::Float64
    BoolArrtime::Array{Float64,1}
    BoolDeptime::Array{Float64,1}
end

CEV = []

function EVmapping(EVData,NPL,PLtype)
structEV = Matrix{EV}(NPL,length(EVData))
for i in 1: NPL
    for k in 1:length(EVData)
        EVData1 = split(EVData[k])
        structEV[i,k] = EV(EVData1[1], EVData1[2], parse(Float64,EVData1[3]), parse(Float64,EVData1[4])/100,
        map(x->parse(Float64,x),EVData1[5:28]),
            parse(Float64,EVData1[29])/100,zeros(24), zeros(24))
        if PLtype == 1
           structEV[i,k].BoolArrtime[findfirst(structEV[i,k].Time)] = 1
           structEV[i,k].BoolDeptime[findlast(structEV[i,k].Time)] = 1
        else
           structEV[i,k].BoolDeptime[findfirst(structEV[i,k].Time, 0)-1] = 1
           structEV[i,k].BoolArrtime[findlast(structEV[i,k].Time, 0)+1] = 1
        end
    end
end
return structEV
end

N_cpl = 100 ; N_spl = 20

EVData = readlines("EV_DataL1.txt") ; push!(CEV,EVmapping(EVData,N_cpl,1))
EVData = readlines("EV_DataL2.txt") ; push!(CEV,EVmapping(EVData,N_cpl,1))

end
