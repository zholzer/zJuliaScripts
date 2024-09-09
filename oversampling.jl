using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

function check_sampling_selection_and_clear_file(oversampleOption::Bool, undersampleOption::Bool, outFile::String)
    if (oversampleOption == false) && (undersampleOption == false)
        display("Neither oversampling or undersampling was selected.")
        exit()
    else
        write(outFile, "")
    end
end

function get_outfile_name_from_infile(inFile::String, outFile::String)
    splitFile = split(inFile, "/")
    splitFile = split(splitFile[end], ".")
    fileName = splitFile[1]
    outFileName = outFileLocation*fileName*"-os.fasta"
    return outFileName
end

function adjust_dates_for_missing_values(splitDate::Vector, missingDateOptions::Vector)
    if length(splitDate) == 3
        adjustedDate = [splitDate[1], splitDate[2], splitDate[2]]
    elseif length(splitDate) == 2
        adjustedDate = [splitDate[1], splitDate[2], missingDateOptions[2]]
    elseif length(splitDate) == 1
        adjustedDate = [splitDate[1], missingDateOptions[1], missingDateOptions[2]] # where to put samples missing months?
    end
    return adjustedDate
end

function read_fasta_file(inFile::String, fastaDF::DataFrame, missingDateOptions::Vector{String})
    FastaReader(inFile) do fr # parse fasta file
        for (header, seq) in fr

            splitHeader = split(header, "_") # split header
            ID = splitHeader[1]

            splitDate = split(splitHeader[end], "-") # last input is date

            if splitDate[1] == "" # if no date skip
                global noDateCount += 1
                continue
            else # otherwise push to dataframe 
                adjustedDate = adjust_dates_for_missing_values(splitDate, missingDateOptions)
                push!(fastaDF, [ID, adjustedDate[1], adjustedDate[2], adjustedDate[3], seq])
            end

        end
    end
end

function set_conditions(frame::SubDataFrame, conditions::Vector, monthsPerGroup::Integer)
    frameMonth = tryparse.(Int, frame.Month)
    let iter = 1
        for i in 1:monthsPerGroup:12
            conditions[iter] = (frameMonth .>= i) .& (frameMonth .< i+monthsPerGroup)
            iter = iter + 1
        end
    end
end

function push_to_dfsByGroupSize(frame::SubDataFrame, dfsByGroupSize::Array, conditions::Vector)
    for i in 1:length(conditions)
        if size(frame[conditions[i], :], 1) != 0
            push!(dfsByGroupSize, frame[conditions[i], :])
        end
    end
end

function split_fastaDF_by_year_and_month(fastaDF::DataFrame, dfsByGroupSize::Array, monthsPerGroup::Integer)
    fastaGroupedByYearDF = groupby(fastaDF, :Year)
    conditions = Vector{BitVector}(undef, length(1:monthsPerGroup:12))
    for frame in fastaGroupedByYearDF
        set_conditions(frame, conditions, monthsPerGroup)
        push_to_dfsByGroupSize(frame, dfsByGroupSize, conditions)
    end
end

function oversample_undersample(minSamplesNeeded::Integer)
    for frame in dfsByGroupSize
        numRow = size(frame, 1)

        if oversampleOption == true
            oversample(frame, minSamplesNeeded, numRow)
        end
        if undersampleOption == true
            undersample(frame, minSamplesNeeded, numRow)
        end
        
        write_to_os_file(frame::DataFrame, outFile::String)
        
    end
end

function oversample(frame::DataFrame, minSamplesNeeded::Integer, numRow::Integer)
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

function undersample(frame::DataFrame, minSamplesNeeded::Integer, numRow::Integer)
    let copyNum = 0
        while size(frame, 1) > minSamplesNeeded 
            let iter = rand(1:numRow - copyNum)
                delete!(frame, iter)
                copyNum = copyNum + 1
            end
        end
    end
end

function write_to_os_file(frame::DataFrame, outFile::String)
    for row in eachrow(frame)
        header = row["ID"]*"_"*string(row["Year"])*"-"*string(row["Month"])*"-"*string(row["Day"])
        seq = row["Seq"]
        writefasta(outFile, [(header, seq)], "a")
    end
end

# main

# script inputs
# input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd
inFile = "ncbi_canineParvovirus_2024-04-15_China_VP2-clean.fasta"
outFileLocation = ""
monthsPerGroup = 1 # partition by 1, 2, 3, 4, 6, and 12 month groups -> 1=each month, 12=each year
minSamplesNeeded = 10 # fixed number for all sampling
monthIfNoMonth = "01"
dayIfNoDay = "01"
oversampleOption = true
undersampleOption = true

fastaDF = DataFrame(ID = String[], Year = String[], Month = String[], Day = String[], Seq = String[]) # initialize dataframe
dfsByGroupSize = Array{DataFrame}(undef,0)

outFile = get_outfile_name_from_infile(inFile, outFileLocation)
check_sampling_selection_and_clear_file(oversampleOption, undersampleOption, outFile)

read_fasta_file(inFile, fastaDF, [monthIfNoMonth, dayIfNoDay])
split_fastaDF_by_year_and_month(fastaDF, dfsByGroupSize, monthsPerGroup)
oversample_undersample(minSamplesNeeded)

display("Done!")
