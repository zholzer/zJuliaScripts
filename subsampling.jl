using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
inFile = "~/../Exploratory/ncbi_dengueV1-E_global_2024-06-19.fasta" # input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd 
n = 12 # partition by 1, 2, 3, 4, 6, and 12 month groups
noMonth = "0" # if 0, missing months are excluded from exploration. If 1-12, will set missing months in that month (1 or 6 reccomended)

lowTime = (2016, 01) # lower bound of dates wanted. set to 0,0 for all
highTime = (3000, 01) # upper bound of dates wanted set to 3000, 30 or other large numbers for all
chosenCountry = "India" # country to sample from

heatmapBool = true

out = false # if true outputs a subsample
outFile = "~/../Data/Sources/ncbi_dengueV1-E_India_2024-06-19.fasta" # name of outputfile

noDateCount = 0 # initialize counts
noPlaceCount = 0
df = DataFrame(ID = String[], Country = String[], Region = String[], Year = Int64[], Month = Int64[]) # initialize dataframe
FastaReader(inFile) do fr # parse fasta file
    for (header, seq) in fr
        splitHeader = split(header, "_") # split header
        splitPlace = split(splitHeader[2], ":") # second input is location as country:region
        splitDate = split(splitHeader[end], "-") # last input is date
        if splitDate[1] == "" # if no date skip
            global noDateCount += 1
            continue
        elseif splitPlace[1] == "" # if no location skip
            global noPlaceCount += 1
            continue
        else # otherwise push to dataframe
            if length(splitPlace) >= 2
                place = [splitPlace[1], splitPlace[2]]
            elseif length(splitPlace) >= 1 # case where no region
                place = [splitPlace[1], "NA"]
            end
    
            if length(splitDate) >= 2
                yearAndMonth = [splitDate[1], splitDate[2]]
            elseif length(splitDate) == 1
                yearAndMonth = [splitDate[1], noMonth] # where to put samples missing months?
            end

            year = tryparse(Int64, yearAndMonth[1])
            month = tryparse(Int64, yearAndMonth[2])
            push!(df, [splitHeader[1], place[1], place[2], year, month])
        end
    end
end

print("There are ", noDateCount, " samples with no date and ", noPlaceCount, " samples with no place.\n") # print number of nodate/noplace samples

minYear = lowTime[1] #minimum(df.Year) # get minimum and maximum year in dataset. can change to lowTime[1] for figures
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

dfGrouped = groupby(df, :Country) # group dataframe by country
numCountries = length(dfGrouped)
countMatrix = zeros(Int, numCountries, Int(numMonths)) # get counts for each country at each year/month pair
for i in 1:numCountries
    for sample in eachrow(dfGrouped[i])
        let iter = 1
            for y in minYear:maxYear
                for m in 1:n:12
                    if ((sample.Year == y) && (sample.Month in m:m+n-1))
                        countMatrix[i,iter] += 1
                    end
                    iter += 1
                end
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

countryNames = [group[1, :Country] for group in dfGrouped] # Create a DataFrame with countries and qualitative counts
dfCounts = DataFrame(Country = countryNames)
dfQuantCounts = dfCounts
for j in 1:Int(numMonths) # Add qualitative columns to dfCounts
    dfCounts[!, yearMonthHeaders[j]] = qualitativeMatrix[:, j]
    dfQuantCounts[!, yearMonthHeaders[j]] = countMatrix[:, j]
end

if heatmapBool
    # code for heatmap figure
    selected_xticks = 1:Int(2*12/n):length(yearMonthHeaders) # Select specific labels for the x-axis
    selected_xtick_labels = yearMonthHeaders[selected_xticks]

    catLength = length(categories) # Display the heatmap with the adjusted xticks
    colors = cgrad(:dense, catLength, categorical = true)
    p1 = heatmap(qualitativeMatrix, c = colors,
        xticks = (selected_xticks, selected_xtick_labels),  # Adjusted xticks
        yticks = (1:size(qualitativeMatrix, 1), countryNames),
        size = (1200, 1200), grid = true, cbar=false)
    M, N = size(qualitativeMatrix)
    vline!(selected_xticks, c=:black, legend = false, linewidth=.5) # makes a grid
    hline!(0.5:(M+0.5), c=:black, legend = false, linewidth=.5)

    # https://discourse.julialang.org/t/set-custom-colorbar-tick-labels/69806/4
    clrticks = categories
    yt = range(0,1,catLength+1)[1:catLength] .+ 0.5/catLength
    l = @layout [ a{.85w} [b{0.1h}; c] ]
    p2 = plot([NaN], lims=(0,1), framestyle=:none, legend=false);
    xx = range(0,1,100)
    zz = zero(xx)' .+ xx
    p3 = heatmap(xx, xx, zz, ticks=false, ratio=1, legend=false, fc=colors, lims=(0,1), # create figure with categorical color legend
                framestyle=:box, right_margin=20mm);
    [annotate!(1.25, yi, text(ti, 10, "Computer Modern")) for (yi,ti) in zip(yt,clrticks)]
    display(plot(p1, p2, p3, layout=l, margins=0mm))
end

dfCountry = dfGrouped[(Country=chosenCountry,)]
print("For ", chosenCountry, " the minimum year is ", minimum(dfCountry[:,"Year"]), " and the maximum is ", maximum(dfCountry[:,"Year"]), ".\n")
function outputSubsample(inFile, outFile, chosenCountry, lowTime, highTime) # function to output subsample at certain years for specific country
    countryIDs = dfCountry[:, 1]
    noDateCountLocal = 0
    FastaReader(inFile) do fr
        for (header, seq) in fr
            splitHeader = split(header, "_")
            splitDate = split(splitHeader[end], "-")

            if length(splitDate) >= 2
                yearAndMonth = [splitDate[1], splitDate[2]]
            elseif length(splitDate) == 1
                yearAndMonth = [splitDate[1], "0"] # where to put samples missing months?
            end

            if splitHeader[1] in countryIDs
                if splitDate[1] == ""
                    global noDateCountLocal += 1
                    continue
                else
                    year = tryparse(Int64, yearAndMonth[1])
                    month = tryparse(Int64, yearAndMonth[2])
                end
                if (year > lowTime[1]) & (year < highTime[1])
                    writefasta(outFile, [(header, seq)], "a")
                elseif year == lowTime[1]
                    if (month >= lowTime[2]) || (month == 0)
                        writefasta(outFile, [(header, seq)], "a")
                    end
                elseif year == highTime[1]
                    if (month <= highTime[2]) || (month == 0)
                        writefasta(outFile, [(header, seq)], "a")
                    end
                end
            end
        end
    end
    print("In ", chosenCountry, " there are ", noDateCountLocal, " samples with no date.\n")
end

if out # if wanted, output subsample
    outputSubsample(inFile, outFile, chosenCountry, lowTime, highTime)
end

# 2009 - H1N1 swine flu breakout