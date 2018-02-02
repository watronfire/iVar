from math import log
from statistics import stdev

class Variant( object ):
    def __init__( self, position, ancestral, substitution, coverage, frequency, totalCounts ):
        self.position = position
        self.ancestral = ancestral
        self.substitution = substitution

        self.coverage = list()
        self.coverage.append( coverage )

        self.frequency = list()
        self.frequency.append( float( frequency[:-1] ) / 100.0 )

        self.totalCounts = list()
        self.totalCounts.append( int( totalCounts ) )

        self.freqA = 0.0

        self.Sn = 0.0
        self.updateSn( self.frequency[0] )
        self.variance = 0.0

        # List which variant_manager will use to flag the variants.
        # Somehow the list needs to hold references to the attributes of a variant object.

    def updateSn( self, f ):
        if f != 1.0:
            self.Sn = -( ( ( 1 - f ) * log( 1 - f ) ) + ( f * log( f ) ) ) / log( 2 )
        else:
            self.Sn = 0.0

    def calculateVariance( self ):
        # If only two replicates, then variance is represented by the log distance of the two replicates.
        if len( self.frequency ) == 2:
            self.variance = abs( log( self.frequency[0] ) - log( self.frequency[1] ) )
        # Else we calculate the standard deviation of all frequencies.
        else:
            logFreq = [ log( i ) for i in self.frequency ]
            self.variance = stdev( logFreq )

    def addReplicate( self, coverage, freq, totalCounts ):
        self.coverage.extend( coverage )
        self.frequency.extend( freq )
        self.totalCounts.extend( totalCounts )

        # Update freqA
        self.freqA = sum( self.frequency ) / len( self.frequency )

        # Update Sn
        self.updateSn( self.freqA )

        self.calculateVariance()

    def getChange( self ):
        return self.ancestral + " -> " + self.substitution

    def getAverageFreq( self ):
        return sum( self.frequency ) / len( self.frequency )

    def getAverageCov( self ):
        return sum( self.coverage ) / len( self.coverage )

    def getList( self, inclusionList ):
        outputList = list()
        outputList.append( self.position )
        outputList.append( self.substitution )
        outputList.append( self.getChange() )
        outputList.extend( self.coverage )
        outputList.extend( self.frequency )
        if len( self.frequency ) > 1 :
            outputList.append( self.freqA )
        outputList.extend( ["", "", ""] )
        outputList.append( self.Sn )
        for item in inclusionList:
            outputList.append( eval( "self." + item ) )
        return outputList