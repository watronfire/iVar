import argparse
import numpy as np
from csv import reader
from math import log
from collections import OrderedDict
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Data Structure for holding the variant information
class Variant:
    def __init__( self, position, ancestral, substitution, coverage, frequency, totalCounts ):
        self.position = position
        self.ancestral = ancestral
        self.substitution = substitution
        self.coverage = list()
        self.coverage.append( coverage )
        self.frequency = list()
        self.totalCounts = list()
        self.totalCounts.append( totalCounts )
        self.freqA = 0.0
        self.Sn = 0.0
        if "/" not in frequency:
            self.frequency.append( float( frequency[:-1] ) / 100.0 )
        else:
            self.frequency.append( float( frequency.split( "/" )[0][:-1] ) / 100.0 )


    def addReplicate( self, coverage, freq, totalCounts ):
        self.coverage.extend( coverage )
        self.frequency.extend( freq )
        self.totalCounts.extend( totalCounts )


    def getAverage( self ):
        return sum( self.frequency ) / len( self.frequency )

    def __str__(self):
        frequencyList = list()
        for k in self.frequency:
            frequencyList.append( str( k * 100 ) + "%" )

        self.freqA = sum( self.frequency ) / len( self.frequency )
        freqAStr = str( self.freqA * 100 ) + "%"

        if self.freqA != 1.0:
            self.Sn = -( ( ( 1 - self.freqA ) * log( 1 - self.freqA ) ) + ( self.freqA * log( self.freqA ) ) ) / log( 2 )
        else:
            self.Sn = 0.0

        attributes = [ str( self.position ), self.substitution, self.ancestral + " -> " + self.substitution,
                        ",".join(self.coverage), ",".join( frequencyList ) ]

        if len( self.frequency ) > 1:
            attributes.append( freqAStr )

        return ",".join( attributes )

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
            iTotCounts = int( ro[4] )
            ifrequency = ro[6]
            isubstitution = ro[18]


            # Create a new dictionary entry with the variant information mapped to its bp position.
            returnDict[iposition] = Variant( iposition, iancestral, isubstitution, icounts, ifrequency, iTotCounts )
    return returnDict

# Calculates whether the variants are greater
def fisherTest( variant, ep ):
    outcomes = list()
    pvs = list()
    for j, number in enumerate( variant.coverage ):
        number = int( number )
        remainder = variant.totalCounts[j]
        expectedPercent = ep

        # Values which can be changed are the fisher-exact test expected values, and the outcomes threshold.
        oddsratio, pvalue = fisher_exact( [ [number, remainder], [expectedPercent * 100, (1-expectedPercent) * 100] ], alternative="greater" )
        pvs.append( pvalue )
        outcomes.append( pvalue < 0.05 )

    return outcomes, pvs

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
parser.add_argument( "-i", "--input", nargs="+", help="Path to input variant files", type=str, required=True)
args = parser.parse_args()

reference = args.reference
pathToOutput = args.output
frequencyLimit = args.frequency
fileList = args.input
replicates = len( fileList )

# Iterate throught the virus map, and load into list. As well, saving the sequence.
lines = list()
with open( reference, "r" ) as virusMap:
    rd = reader( virusMap )
    seq = "n"
    fl = True
    for row in rd:
        lines.append( row )
        if not fl:
            seq += row[ 5 ]
        else:
            fl = False

# Generates variantDictionary with replicates if necessary.
variantDict = dict()
for i, file in enumerate( fileList ):
    with open( file, "r" ) as inputFile:
        if i == 0:
            variantDict = variantParser( inputFile )
        else:
            tempVariantDict = variantParser( inputFile )
            for pos in variantDict.keys():
                if pos in tempVariantDict.keys() and variantDict[pos].substitution == tempVariantDict[pos].substitution:
                    variantDict[pos].addReplicate( tempVariantDict[pos].coverage, tempVariantDict[pos].frequency, tempVariantDict[pos].totalCounts )

graphDict = dict()
pvalues = list()
acceptedValues = dict()

# Clean up dictionary...
#sortedKeys = sorted( variantDict.keys(), key=int )

for entry in list( variantDict.keys() ):
    oc, pv = fisherTest( variantDict[entry], frequencyLimit )

    if replicates > 1:

        # Because of the way this program assembles the variantDict, entries in which replicates weren't found need to
        # be removed.
        if len( variantDict[entry].coverage ) < replicates:
            del variantDict[entry]
        else:
            # Converts entry to an int so that it can be sorted. If left uncoverted, then matplotlib will sort the variants
            # alphabetically, instead of numerically. Graph isn't the best...
            graphDict[int(entry)] = pv
            pvalues.extend( pv )

            # Creates a subset of graphDict which includes values which would have been accepted by a threshold test
            # rather than the fisher exact test being used.
            average = sum( variantDict[entry].frequency ) / len( variantDict[entry].frequency )
            if average > frequencyLimit :
                acceptedValues[int(entry)] = graphDict[int(entry)]

            # Deletes entries from variant dict unless all replicates pass fisher exact test.
            if not all( oc ):
                del variantDict[entry]

# Generates a histogram of p-values
plt.subplot( 2, 1, 1 )

lowestNum = 0.01
for i in list( graphDict.values() ):
    for k in i:
        if k < lowestNum:
            lowestNum = k
n, bins, patches = plt.hist( pvalues, facecolor="red", alpha=0.75, bins=np.logspace(np.log10(lowestNum/10),np.log10(1.0), 50) )
plt.gca().set_xscale( "log" )
plt.gca().set_ylabel( "Frequency", weight="bold" )
plt.gca().set_xlabel( "p-values", weight="bold" )
plt.gca().set_title( "Histogram of p-values", weight="bold" )
plt.gca().axvline( x=0.05, color="black", linestyle="dashed", label="0.05" )
plt.gca().axvline( x=0.01, color="black", linestyle=":", label="0.01" )
plt.legend()
imageOutput = "/".join( pathToOutput.split("/")[:-1] ) + "/" + pathToOutput.split("/")[-1].split(".")[0]

# Generates a manhattan plot. Also colors values which would have been accepted by a flat threshold, instead of fisher
# exact test.
plt.subplot( 2, 1, 2 )
plt.plot( list(graphDict.keys()), list(graphDict.values()), "ro" )
plt.plot( list( acceptedValues.keys() ), list( acceptedValues.values() ), "bo", label="Above Frequency Threshold" )
plt.gca().set_yscale( "log" )
plt.gca().set_ylabel( "p-value", weight="bold" )
plt.gca().set_xlabel( "Genome Position (bp)", weight="bold" )
plt.gca().set_title( "Manhattan Plot", weight="bold" )
plt.gca().invert_yaxis()
plt.gca().axhline( y=0.01, color="black", linestyle=":", label="0.01" )
plt.gca().axhline( y=0.05, color="black", linestyle="dashed", label="0.05" )
plt.gcf().set_figheight( 12 )
plt.gcf().set_figwidth( 12 )

# Prevents duplicat entries in legend.
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

# Saves the graphs
plt.savefig( imageOutput + ".png" )

# Creates, or overrides, an output file.
with open( pathToOutput, "w" ) as outputFile :

    # List to hold lines before writing to file.
    writingBuffer = list()

    # Some statistical measures which will be calculated along the way.
    complexitySn = 0
    diversityNT = 0
    richness = 0
    distanceN = 0
    distanceS = 0

    for i, line in enumerate( lines ):

        # String which will hold the updated line in document.
        tempLine = ",".join( line )

        # If reading the first line then append additional column names.
        if i == 0:
            if replicates > 1:
                coverages = [ "coverage_%s" % i for i in range( 1, replicates + 1 ) ]
                var_freqs = [ "var_freq_%s" % i for i in range( 1, replicates + 1 ) ]
                tempLine += ",var_site,variant,Change," + ",".join( coverages ) + "," + ",".join( var_freqs ) + ",var_freq_ave,var_3nt_aa,var_aa,N_S,Sn"
            else:
                tempLine += ",var_site,variant,change,coverage,var_freq,var_3nt_aa,var_aa,N_S,Sn,"
            writingBuffer.append( tempLine )
            continue

        # If the position is found in the variant dictionary, then we add the information we have for the variant to the line.
        variantFound = False
        if line[0] in variantDict:           #and variantDict[line[0]].freqA > 0.03:
            tempLine += "," + str( variantDict[line[0]] ) + ","
            variantFound = True

        # Else add the necessary spaces and the variants sequence.
        else :
            if replicates > 1:
                tempLine += ",," + line[5] + ",,," + ( ",," * ( replicates ) )
            else:
                tempLine += ",," + line[5] + ",,,,"

        # Next translate variant sequence.
        codon = ""
        if "UTR" not in line[1] :

            pos = int( line[0] )
            codonPos = int( line[2] )

            # If a variant was found then we use its substitution.
            var = ""
            if variantFound :
                var = variantDict[line[0]].substitution
            else :
                var = seq[pos]

            # Create the Codon at each position.
            if codonPos == 1 :
                codon = var + seq[pos + 1] + seq[pos + 2]
            elif codonPos == 2 :
                codon = seq[pos - 1] + var + seq[pos + 1]
            else :
                codon = seq[pos - 2] + seq[pos - 1] + var

            # Determine whether mutation is synonymous or non-synonymous. Also going to calculate some statistics here.
            mutationType = ""
            if variantFound :

                # Calculate specified statistics.
                complexitySn += variantDict[line[0]].Sn
                diversityNT += variantDict[line[0]].freqA
                richness += 1

                if line[6] == codonTable[codon] :
                    mutationType = "S" + "," + str( variantDict[line[0]].Sn )
                    distanceS += variantDict[line[0]].freqA
                else :
                    mutationType = "N" + "," + str( variantDict[line[0]].Sn )
                    distanceN += variantDict[line[0]].freqA

                    # add the codon, and its translation to the line.
            tempLine += codon + "," + codonTable[codon] + "," + mutationType

        # write the line to file.
        # outputFile.write( tempLine  )
        writingBuffer.append( tempLine )

    distance = diversityNT
    complexitySn /= 10272
    diversityNT /= 10272
    try:
        selectionPN = distanceN / (distanceN + distanceS)
    except ZeroDivisionError:
        selectionPN = 0

    if replicates > 1:
        writingBuffer[0] += ",,Test,region,result"
    else:
        writingBuffer[0] += ",Test,region,result"
    writingBuffer[1] += ",,,,,Complexity_Sn,CDS," + str( complexitySn )
    writingBuffer[2] += ",,,,,Diversity_nt,CDS," + str( diversityNT )
    writingBuffer[3] += ",,,,,Richness,CDS," + str( richness )
    writingBuffer[4] += ",,,,,Distance,CDS," + str( distance )
    writingBuffer[5] += ",,,,,Distance,CDS_N," + str( distanceN )
    writingBuffer[6] += ",,,,,Distance,CDS_S," + str( distanceS )
    writingBuffer[7] += ",,,,,Selection_pN,CDS," + str( selectionPN )

    for entry in writingBuffer :
        outputFile.write( entry + "\n" )
