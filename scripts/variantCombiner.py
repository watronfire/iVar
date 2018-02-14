import argparse
import numpy as np
from csv import reader
from variant import Variant
from variant_manager import VariantManager
from collections import OrderedDict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def variantParser( variantFile ):
    returnDict = dict()
    firstLine = True
    for ro in reader( variantFile, delimiter="\t" ):
        if firstLine:
            firstLine = False
        else:
            # Calculate and assign the variant determinants.
            iposition = ro[1]
            iancestral = ro[2]
            icounts = ro[5]
            iTotCounts = ro[4]
            ifrequency = ro[6]
            isubstitution = ro[18]


            # Create a new dictionary entry with the variant information mapped to its bp position.
            returnDict[iposition] = Variant( iposition, iancestral, isubstitution, icounts, ifrequency, iTotCounts )
    return returnDict

# Generate a translation table dictionary
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
aminoAcids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codonTable = dict( zip( codons, aminoAcids ) )

# Parse parameters
parser = argparse.ArgumentParser()
parser.add_argument( "-r", "--reference", help="Path to virus map", type=str, required=True )
parser.add_argument( "-o", "--output", help="Path to output file", type=str, required=True )
parser.add_argument( "-f", "--frequency", help="Minimum variant frequency allowed", type=float, required=True )
parser.add_argument( "-p", "--pvalue", help="Maximum p-value for a variant to be accepted", type=float, required=False )
parser.add_argument( "-i", "--input", nargs="+", help="Path to input variant files", type=str, required=True )
parser.add_argument( "-t", "--test", help="Statistical test used to confirm variants. 1-Fisher's Exact Test, 2-Strict Frequency Limit, 3-Average Frequency Limit", type=int, required=False )
args = parser.parse_args()
pvLimit = 0.05

if args.test:
    # Assign the test. Furthermore, if fishers exact test is used then a maximum p-value is required.
    test = args.test
    if test == 1:
        pvLimit = args.pvalue
# The default program will accept variants based on the average frequency.
else:
    test = 3

reference = args.reference
pathToOutput = args.output
frequencyLimit = args.frequency
fileList = args.input
replicates = len( fileList )

inclusionList = [ "variance" ]


########################
##  Combine Variants  ##
#########################################################################################

# Generates variantDictionary with replicates if necessary.
variantDict = dict()
for i, file in enumerate( fileList ):
    with open( file, "r" ) as inputFile:
        if i == 0:
            variantDict = variantParser( inputFile )
        else:
            tempVariantDict = variantParser( inputFile )
            for pos in variantDict.keys() & tempVariantDict.keys():
                if variantDict[pos].substitution == tempVariantDict[pos].substitution:
                    variantDict[pos].addReplicate( tempVariantDict[pos].coverage, tempVariantDict[pos].frequency, tempVariantDict[pos].totalCounts )


########################
##  Confirm Variants  ##
#########################################################################################

for entry in list( variantDict.keys() ):
    if len( variantDict[entry].frequency ) != replicates:
        del variantDict[entry]

vm = VariantManager( variantDict, False )

histogramList, manhattanDict, highlights = vm.compareFreq( frequencyLimit, test, pvLimit)

variantDict = vm.getManagedVariants()


################
##  Graphing  ##
#########################################################################################

# Generates a histogram of p-values
plt.subplot( 2, 1, 1 )

# Calculates the lower end of the histogram so adequate resolution is seen.
lowestNum = min( histogramList )

n, bins, patches = plt.hist( histogramList, facecolor="red", alpha=0.75, bins=np.logspace(np.log10(lowestNum/10),np.log10(1.0), 50) )
plt.gca().set_xscale( "log" )
plt.gca().set_ylabel( "Frequency", weight = "bold" )
if test == 1:
    plt.gca().set_xlabel( "p-values", weight="bold" )
    plt.gca().set_title( "Histogram of p-values", weight="bold" )
else:
    plt.gca().set_xlabel( "Frequency of Variant", weight = "bold" )
    plt.gca().set_title( "Histogram of Variant Frequency", weight = "bold" )
plt.gca().axvline( x=0.05, color="black", linestyle="dashed", label="0.05" )
plt.gca().axvline( x=0.01, color="black", linestyle=":", label="0.01" )
plt.legend()
imageOutput = "/".join( pathToOutput.split("/")[:-1] ) + "/" + pathToOutput.split("/")[-1].split(".")[0]

# Generates a manhattan plot. Also colors values which would have been accepted by a flat threshold, instead of fisher
# exact test.
plt.subplot( 2, 1, 2 )
plt.plot( list( manhattanDict.keys() ), list( manhattanDict.values() ), "ro" )
plt.gca().set_yscale( "log" )
plt.gca().set_xlabel( "Genome Position (bp)", weight="bold" )
plt.gca().set_title( "Manhattan Plot", weight="bold" )
plt.gca().set_ylabel( "Frequency", weight = "bold" )
histLabel = "Above P-value Threshold" if test == 1 else "Above Frequency Threshold"
plt.plot( list( highlights.keys() ), list( highlights.values() ), "bo", label = histLabel )
plt.gca().axhline( y=frequencyLimit, color="black", linestyle="dashed", label=str( frequencyLimit ) )
plt.gcf().set_figheight( 12 )
plt.gcf().set_figwidth( 12 )

# Prevents duplicate entries in the legend.
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

# Saves the graphs
plt.savefig( imageOutput + ".png" )


##############
##  Output  ##
#########################################################################################

# Load virus map from file
virusMap = pd.read_csv( reference )

# Replace empty cells will empty strings.
virusMap = virusMap.replace( np.nan, "", regex=True )

# Introduce some new columns
virusMap["var_site"] = ""
virusMap["variant"] = virusMap["nucleotide"]
virusMap["change"] = ""
for i in range( 1, replicates + 1 ):
    virusMap["coverage_%s" % i] = ""
for i in range( 1, replicates + 1 ):
    virusMap["var_freq_%s" % i] = ""
if replicates > 1:
    virusMap["var_freq_ave"] = ""
virusMap["var_3nt_aa"] = ""
virusMap["var_aa"] = ""
virusMap["N_S"] = ""
virusMap["Sn"] = ""
for item in inclusionList:
    virusMap[item] = ""

for i in range( len( virusMap ) ):
    bp = str( virusMap.iloc[i,0] )
    variantFound = False
    if bp in variantDict:
        virusMap.iloc[i,7:] = variantDict[bp].getList( inclusionList )
        variantFound = True

    if "UTR" not in virusMap.iloc[i,1]:
        codon = ""
        codonPos = virusMap.iloc[i,2]

        if codonPos == 1 :
            codon = virusMap.iloc[i,8] + virusMap.iloc[i + 1, 8] + virusMap.iloc[i + 2, 8]
        elif codonPos == 2 :
            codon = virusMap.iloc[i - 1, 8] + virusMap.iloc[i, 8] + virusMap.iloc[i + 1, 8]
        else :
            codon = virusMap.iloc[i - 2, 8] + virusMap.iloc[i - 1, 8] + virusMap.iloc[i, 8]

        virusMap.iloc[i, virusMap.columns.get_loc( "var_3nt_aa" )] = codon
        virusMap.iloc[i, virusMap.columns.get_loc( "var_aa" )] = codonTable[codon]
        if variantFound:
            virusMap.iloc[i, virusMap.columns.get_loc( "N_S" )] = "S" if virusMap.iloc[i, 6] == virusMap.iloc[i, virusMap.columns.get_loc( "var_aa" )] else "N"


##################
##  Statistics  ##
#########################################################################################

# Need to output the statistics
virusMap["Sn"] = pd.to_numeric( virusMap["Sn"] )
if replicates == 1:
    virusMap["var_freq_1"] = pd.to_numeric( virusMap["var_freq_1"] )
    diversityNT = virusMap["var_freq_1"].sum()
    try:
        distanceN = virusMap.loc[ virusMap["N_S"] == "N", "var_freq_1" ].iloc[0]
    except IndexError:
        distanceN = 0
    try:
        distanceS = virusMap.loc[ virusMap["N_S"] == "S", "var_freq_1" ].iloc[0]
    except IndexError:
        distanceS = 0
else:
    virusMap["var_freq_ave"] = pd.to_numeric( virusMap["var_freq_ave"] )
    diversityNT = virusMap["var_freq_ave"].sum()
    try:
        distanceN = virusMap.loc[ virusMap["N_S"] == "N", "var_freq_ave" ].iloc[0]
    except IndexError:
        distanceN = 0
    try:
        distanceS = virusMap.loc[ virusMap["N_S"] == "S", "var_freq_ave" ].iloc[0]
    except IndexError:
        distanceS = 0
distance = diversityNT
diversityNT = diversityNT / 10272

complexitySn = virusMap["Sn"].sum() / 10272

richness = virusMap["Sn"].count()

try:
    selectionPN = distanceN / (distanceN + distanceS)
except ZeroDivisionError:
    selectionPN = 0

# Construct a dataframe with stats
stats_data = { "Test" : [ "Complexity_Sn", "Diversity_nt", "Richness", "Distance", "Distance","Distance", "Selection_pN" ],
               "region" : [ "CDS", "CDS", "CDS", "CDS", "CDS_N", "CDS_S", "CDS" ],
               "result" : [ complexitySn, diversityNT, richness, distance, distanceN, distanceS, selectionPN ] }

statsDF = pd.DataFrame( stats_data, columns=[ "Test", "region", "result" ] )
statsDF.to_csv( imageOutput + ".stats.csv", index=False )

virusMap.to_csv( pathToOutput, index=False )