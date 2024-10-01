using FastaIO
using DataFrames
using Statistics
using Measures, Plots; gr()
using CSV

# script inputs
tableFile = "~/../Data/Runs/beastRun_2024-08-23_dengueV2-E_Japan-os/gisaid_dengueV2-E_2024_08_16_Japan_final-os_run0.tsv"
# input a fasta file with ID, Country:Region, [other stuff], date yyyy-mm-dd
inFile = "~/../Data/Oversampled/gisaid_dengueV2-E_2024_08_16_Japan_final-os.fasta"
n = 1 # partition by 1, 2, 3, 4, 6, and 12 month groups
noMonth = "0" # if 0, missing months are excluded from exploration. If 1-12, will set missing dates in that month (1 or 6 reccomended)
savePlot = false
plotTogether = false
noDateCount = 0 # initialize counts

function read_fasta_year_month(inFile::String, fastaDF::DataFrame)
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
                push!(fastaDF, [year, month])
            end
        end
    end
end

yearMonthDF = DataFrame(Year = Int64[], Month = Int64[]) # initialize dataframe
read_fasta_year_month(inFile, yearMonthDF)

minYear = minimum(yearMonthDF.Year) # get minimum and maximum year in dataset. can change for figures
maxYear = maximum(yearMonthDF.Year)
numMonths = Int((maxYear - minYear + 1) * 12/n) # number of months depending on n

function create_plot_headers(minYear::Int, maxYear::Int, n::Int)
    let iter = 1
        for y in minYear:maxYear
            for m in 1:n:12
                yearMonthHeaders[iter] = string(y, "-", m)
                iter += 1
            end
        end
    end
end

function get_count_of_each_year_month(fastaDF::DataFrame, minYear::Int, maxYear::Int, n::Int, numMonths::Int)
    for sample in eachrow(fastaDF)
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
end

yearMonthHeaders = Array{String}(undef, Int(numMonths)) # create the yyyy-mm headers
create_plot_headers(minYear, maxYear, n)

countMatrix = zeros(Int, 1, Int(numMonths)) # get counts for each country at each year/month pair
get_count_of_each_year_month(yearMonthDF, minYear, maxYear, n, numMonths)
dfQuantCounts = DataFrame(countMatrix, yearMonthHeaders)

skylineDf = DataFrame(CSV.File(tableFile, delim = '\t', header = 2))
row_data = Array(dfQuantCounts)'
yearMonths = names(dfQuantCounts)
# Create the bar chart

function convert_date_to_decimal_format(yearMonths::Vector{String})
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
end

decimalYears = zeros(Float64, length(yearMonths))
convert_date_to_decimal_format(yearMonths)

function get_main_filename(tableFile::String)
    splitFile = split(tableFile, "/")
    splitFile = split(splitFile[end], ".")
    fileName = splitFile[1]
    return fileName
end

y = skylineDf[:, "mean"]

if plotTogether
    x = skylineDf[:, "time"]
    fileName = get_main_filename(tableFile)

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
end

display(histogram2d(row_data,y,bins=100,label=""))
# https://discourse.julialang.org/t/plots-reduce-size-of-scatter-plots-etc/44445/3

if savePlot
    outFile = "Figures/"*fileName*".png"
    savefig(outFile)
end