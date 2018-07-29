module NetworkBuild
export Network
export ReadMatData
#############################
# Reading Network data through *.csv file format
#############################

gen = ["bus", "gen_ind", "a", "b", "c", "pmin", "pmax"]
branch = ["fbus", "tbus", "r", "x", "b", "rateA"]
load = ["bus"]
Network = Dict()

function ReadMatData(file::String)
    if endswith(file, ".m")
       NetworkData = parse_data(file)
    end
    return NetworkData
end

function Str_to_Var(s::String)
   s = Symbol(s)
   @eval (($s))
 end

function parse_data(datafile :: String)
  datalines = readlines(datafile)
  Index = 1
  while Index <= length(datalines)
    if length(strip(datalines[Index]))<=0 || strip(datalines[Index])[1] == '%'
      Index = Index +1
      continue
    end
    if contains(datalines[Index], "function mpc")
       Network["Name"] = strip(split(datalines[Index], "=")[2], ['\r','\n', ';', '"',' '])
    elseif contains(datalines[Index], "mpc.baseMVA")
       Network["Base MVA"]= strip(split(datalines[Index], "=")[2], ['\r','\n', ';', '"', ' '])
    elseif contains(datalines[Index], "mpc.BusLen")
       Network["BusLen"]= strip(split(datalines[Index], "=")[2], ['\r','\n', ';', '"',' '])
    elseif contains(datalines[Index], "mpc.") # && contains(datalines[Index], "[")
            if contains(datalines[Index], "load") || contains(datalines[Index], "gen") || contains(datalines[Index], "branch")
                  SParameterType = strip(replace(split(datalines[Index])[1], "mpc.", ""))
                  ParameterType = Str_to_Var(SParameterType)
                  DataValue, Index = parse_MatrixData(datalines, Index)
                  MakeDict(DataValue,ParameterType,SParameterType)
            end
      end
    Index=Index + 1
  end
    return Network
end

function parse_MatrixData(datalines::Array{String,1}, Index::Int)
    MatrixData = Array[];
    while Index <= length(datalines)
      if length(strip(datalines[Index]))<=0 || contains(datalines[Index],"%") || contains(datalines[Index],"[")
        Index = Index + 1
      continue
      end
      contains(datalines[Index],"]") ? break : ""
      push!(MatrixData, split(replace(datalines[Index], ";","")))
      Index = Index + 1
    end
    return MatrixData, Index
end

function MakeDict(ParameterValue, ParameterType, SParameterType)
    Parameters=Dict();
    for irow in 1:length(ParameterValue)
        TypeDict = Dict()
        for icol in 1:length(ParameterType)
            icol ==1 || icol ==2 ? TypeDict[ParameterType[icol]]= parse(Int, ParameterValue[irow][icol]):
            TypeDict[ParameterType[icol]]= parse(Float64, ParameterValue[irow][icol])
        end
        Parameters["$irow"]=TypeDict
    end
    Network[SParameterType]=Parameters
end
end
