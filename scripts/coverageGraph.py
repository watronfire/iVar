import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pyplot as plt
import sys
from csv import reader

inputFile = sys.argv[1]
outputFile = sys.argv[2]
inputName = inputFile.split("/")[-1]

coverageDict = dict()

with open( inputFile, "r" ) as input:
	for line in reader( input, delimiter="\t" ):
		coverageDict[int( line[1] )] = line[2]

plt.plot( list( coverageDict.keys() ), list( coverageDict.values() ) )
plt.gca().set_yscale( "log" )
plt.gca().set_ylabel( "log Depth", weight="bold" )
plt.gca().set_xlabel( "Genome Position (bp)", weight="bold" )
plt.gca().set_title( "Coverage Map for " + inputName, weight="bold" )
plt.gcf().set_figheight( 6 )
plt.gcf().set_figwidth( 12 )


plt.savefig( outputFile )
