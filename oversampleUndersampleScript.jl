using DataFrames
include("oversampleUndersampleSRC.jl")

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
