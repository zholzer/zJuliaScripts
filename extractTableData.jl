using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
inFile = "Data/Runs/beastRun_2024-05-27/ncbi_fluAH3N2_2024-05-07_Thailand_table.tsv" 

# output from header -> download date, country
# output sequence length, # samples, earliest date, latest date, # unique months
















noDateCount = 0 # initialize counts
noPlaceCount = 0
input = inFileLocation*inFile*".fasta"
df = DataFrame(ID = String[], Year = Int64[], Month = Int64[]) # initialize dataframe
FastaReader(input) do fr # parse fasta file
    for (header, seq) in fr
        splitHeader = split(header, "_") # split header
        splitDate = split(splitHeader[end], "-") # last input is date
        if splitDate[1] == "" # if no date skip
            global noDateCount += 1
            continue
        else # otherwise push to dataframe 
            if length(splitDate) >= 2
                yearAndMonth = [splitDate[1], splitDate[2]]
            elseif length(splitDate) == 1
                yearAndMonth = [splitDate[1], noMonth] # where to put samples missing months?
            end

            year = tryparse(Int64, yearAndMonth[1])
            month = tryparse(Int64, yearAndMonth[2])
            push!(df, [splitHeader[1], year, month])
        end
    end
end


minYear = minimum(df.Year) # get minimum and maximum year in dataset. can change for figures
maxYear = maximum(df.Year)