using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
# input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd
inFile = "Data/Intermediate/Trim/ncbi_fluB_2024-05-24_allDates_USA_trim.fasta"

splitFile = split(inFile, "/")
splitFile = split(splitFile[end], ".")
fileName = splitFile[1]
outFile = "Data/Oversampled/"*fileName*"-os.fasta"

n = 1 # partition by 1, 2, 3, 4, 6, and 12 month groups -> 1=each month, 12=each year
minSamplesNeeded = 10
# fixed number for all sampling
noMonth = "01"
noDay = "01"
oversampleBool = true
undersampleBool = true

df = DataFrame(ID = String[], Year = Int64[], Month = Int64[], Day = Int64[], Seq = String[]) # initialize dataframe
FastaReader(inFile) do fr # parse fasta file
    for (header, seq) in fr
        splitHeader = split(header, "_") # split header
        splitDate = split(splitHeader[end], "-") # last input is date
        if splitDate[1] == "" # if no date skip
            global noDateCount += 1
            continue
        else # otherwise push to dataframe 
            if length(splitDate) == 3
                adjustedDate = [splitDate[1], splitDate[2], splitDate[2]]
            elseif length(splitDate) == 2
                adjustedDate = [splitDate[1], splitDate[2], noDay]
            elseif length(splitDate) == 1
                adjustedDate = [splitDate[1], noMonth, noDay] # where to put samples missing months?
            end

            year = tryparse(Int64, adjustedDate[1])
            month = tryparse(Int64, adjustedDate[2])
            day = tryparse(Int64, adjustedDate[3])
            push!(df, [splitHeader[1], year, month, day, seq])
        end
    end
end

dfYearMonth = Array{Any}(undef,0)
dfYear = groupby(df, :Year)
conditions = Vector{BitVector}(undef, length(1:n:12))
for frame in dfYear
    let iter = 1
        for i in 1:n:12
            conditions[iter] = (frame.Month .>= i) .& (frame.Month .< i+n)
            iter = iter + 1
        end
    end
    for i in 1:length(conditions)
        if size(frame[conditions[i], :], 1) != 0
            push!(dfYearMonth, frame[conditions[i], :])
        end
    end
end

if (oversampleBool == false) && (undersampleBool == false)
    display("Neither oversampling or undersampling was selected.")
    exit()
else
    write(outFile, "")
end

for frame in dfYearMonth
    numRow = size(frame, 1)

    if oversampleBool == true
        let copyNum = 1
            while size(frame, 1) < minSamplesNeeded 
                let iter = rand(1:numRow)
                    editedID = frame[iter, "ID"] * "/copy" * string(copyNum) # needs to be random
                    push!(frame, [editedID, frame[iter, "Year"], frame[iter, "Month"], frame[iter, "Day"], frame[iter, "Seq"]])
                    #duplicate row of frame[i] and add indicator
                    copyNum = copyNum + 1
                end
            end
        end
    end

    if undersampleBool == true
        let copyNum = 0
            while size(frame, 1) > minSamplesNeeded 
                let iter = rand(1:numRow - copyNum)
                    delete!(frame, iter)
                    #duplicate row of frame[i] and add indicator
                    copyNum = copyNum + 1
                end
            end
        end
    end

    for row in eachrow(frame)
        header = row["ID"]*"_"*string(row["Year"])*"-"*string(row["Month"])*"-"*string(row["Day"])
        seq = row["Seq"]
        writefasta(outFile, [(header, seq)], "a")
    end
end

display("Done!")
