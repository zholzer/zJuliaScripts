using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
tableFile = "../Data/Runs/beastRun_2024-08-16_fluBVictoria_Japan/gisaid_fluBVictoria-HA_2024-08-16_Japan_final.tsv"
#tableFile = "Exploratory/Results/gisaid_dengueV1_2024_04_05_China_final_zVBSKY-table.tsv" 
# input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd
inFile = "../Data/Final/gisaid_fluBVictoria-HA_2024-08-16_Japan_final.fasta"

splitFile = split(tableFile, "/")
splitFile = split(splitFile[end], ".")
fileName = splitFile[1]
outFile = "Figures/"*fileName*".png"

n = 1 # partition by 1, 2, 3, 4, 6, and 12 month groups
noMonth = "0" # if 0, missing months are excluded from exploration. If 1-12, will set missing dates in that month (1 or 6 reccomended)


savePlot = false

noDateCount = 0 # initialize counts
noPlaceCount = 0

df = DataFrame(ID = String[], Year = Int64[], Month = Int64[]) # initialize dataframe
FastaReader(inFile) do fr # parse fasta file
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
numMonths = (maxYear - minYear + 1) * 12/n # number of months depending on n
yearMonthHeaders = Array{String}(undef, Int(numMonths)) # create the yyyy-mm headers
let iter = 1
    for y in minYear:maxYear
        for m in 1:n:12
            yearMonthHeaders[iter] = string(y, "-", m)
            iter += 1
        end
    end
end

countMatrix = zeros(Int, 1, Int(numMonths)) # get counts for each country at each year/month pair
for sample in eachrow(df)
    let iter = 1
        for y in minYear:maxYear
            for m in 1:n:12
                if ((sample.Year == y) && (sample.Month in m:m+n-1))
                    countMatrix[1, iter] += 1
                end
                iter += 1
            end
        end
    end
end

# Define the function to categorize values
categories = ["<2", "2-5", "6-10", "11-100", "101-1000", "1001+"]
function categorize(value)
    if value < 2
        return "<2"
    elseif value < 11
        return "2-5"
    elseif value < 101
        return "6-10"
    elseif value < 1001
        return "11-100"
    elseif value < 1001
        return "101-1000"
    else
        return "1001+"
    end
end
qualitativeMatrix = [categorize(count) for count in countMatrix] # categorize the count matrix

#dfCounts = DataFrame(countMatrix, yearMonthHeaders)
dfQuantCounts = DataFrame(countMatrix, yearMonthHeaders)
# for j in 1:Int(numMonths) # Add qualitative columns to dfCounts
#     dfCounts[yearMonthHeaders[j], 1] = qualitativeMatrix[:, 1]
#     dfQuantCounts[yearMonthHeaders[j], 1] = countMatrix[:, 1]
# end

skylineDf = DataFrame(CSV.File(tableFile, delim = '\t', header = 2))
row_data = Array(dfQuantCounts)'
yearMonths = names(dfQuantCounts)
# Create the bar chart

decimalYears = zeros(Float64, length(yearMonths))
# Split the date string into year and month
let iter = 1
    for yearMonth in yearMonths
        parts = split(yearMonth, "-")
        year = parse(Int, parts[1])
        month = parse(Int, parts[2])

        # Convert to decimal year
        decimal_year = year + (month - 1) / 12
        decimalYears[iter] = decimal_year
        iter = iter + 1;
    end
end

x = skylineDf[:, "time"]
y = skylineDf[:, "mean"]
lb = skylineDf[:, "lower"]
ub = skylineDf[:, "upper"]

p1 = plot(decimalYears, row_data, markersize=1.5, marker=:circle, markercolor=:purple,
label = "Relative Sample Size", xticks=(1:length(x), fill("", length(x))), 
ylabel = "Sample Size", legend = false, linecolor = "magenta", 
title = "Plot of Sample Size versus Population Size over Time")

p2 = plot(x, y, marker = :circle, markersize = 1.5, 
xlabel = "Year\nfile="*fileName, ylabel = "Population Size", markercolor = "blue", linecolor = "black",
label = "Population Size", legend = false)
plot!(x, ub, linestyle = :dash, linecolor = "darkslategray")
plot!(x, lb, fillrange = ub, fillalpha = 0.25, c = 1, linestyle = :dash, linecolor = "darkslategray")
display(plot(p1, p2, layout = (2, 1), framestyle=:box))
if savePlot
    savefig(outFile)
end
# 2009 - H1N1 swine flu breakout