using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
# input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd
inFile = "ncbi_canineParvovirus_2024-04-15_China_VP2-clean.fasta"
outFileLocation = ""
n = 1 # partition by 1, 2, 3, 4, 6, and 12 month groups -> 1=each month, 12=each year
minSamplesNeeded = 10
# fixed number for all sampling
monthIfNoMonth = "01"
dayIfNoDay = "01"
oversample = true
undersample = true

outFileName = get_outfile_name_from_infile(inFile, outFileLocation)

fastaDF = DataFrame(ID = String[], Year = Integer[], Month = Integer[], Day = Integer[], Seq = String[]) # initialize dataframe
read_fasta_file(inFile::string, fastaDF::DataFrame, [monthIfNoMonth, dayIfNoDay])

function get_outfile_name_from_infile(inFile::String, outFile::String, missingDateOptions::AbstractVector)
    splitFile = split(inFile, "/")
    splitFile = split(splitFile[end], ".")
    fileName = splitFile[1]
    outFileName = outFileLocation*fileName*"-os.fasta"
    return outFileName
end

function read_fasta_file(inFile::string, fastaDF::DataFrame)
    FastaReader(inFile) do fr # parse fasta file
        for (header, seq) in fr

            splitHeader = split(header, "_") # split header
            splitDate = split(splitHeader[end], "-") # last input is date

            if splitDate[1] == "" # if no date skip
                global noDateCount += 1
                continue
            else # otherwise push to dataframe 
                push_sequence_to_df(splitDate)
            end

        end
    end
end

function push_sequence_to_df(splitDate, missingDateOptions)
    if length(splitDate) == 3
        adjustedDate = [splitDate[1], splitDate[2], splitDate[2]]
    elseif length(splitDate) == 2
        adjustedDate = [splitDate[1], splitDate[2], missingDateOptions[2]]
    elseif length(splitDate) == 1
        adjustedDate = [splitDate[1], missingDateOptions[1], missingDateOptions[2]] # where to put samples missing months?
    end

    year = tryparse(Integer, adjustedDate[1])
    month = tryparse(Integer, adjustedDate[2])
    day = tryparse(Integer, adjustedDate[3])
    push!(fastaDF, [splitHeader[1], year, month, day, seq])
end

dfYearMonth = Array{Any}(undef,0)
dfYear = groupby(fastaDF, :Year)
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

if (oversample == false) && (undersample == false)
    display("Neither oversampling or undersampling was selected.")
    exit()
else
    write(outFile, "")
end

for frame in dfYearMonth
    numRow = size(frame, 1)

    if oversample == true
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

    if undersample == true
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
