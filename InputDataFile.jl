module InputDataFile
function readData(E,P_cnsr,DynPrice,loadbus)
Inputfile = readlines("ResidentialLoad.txt")
  for i in 1:length(Inputfile)
      Data = split(Inputfile[i])
      ind = Data[1][end]
    	if Data[1] == "Load$ind"
         [E[t,loadbus[parse(Int64,ind)]] = parse(Data[t+2]) for t in 1:24]
         elseif Data[1] == "P_cnsr$ind"
         [P_cnsr[t,parse(Int64,ind)] = parse(Data[t+2]) for t in 1:24]
         elseif Data[1] == "DynPrice$ind"
         [DynPrice[loadbus[parse(Int64,ind)],t] = parse(Data[t+2]) for t in 1:24]
      end
  end
  ret = tuple(E,P_cnsr,DynPrice)
  return ret
end
end
