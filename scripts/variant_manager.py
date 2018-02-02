from scipy.stats import fisher_exact

class VariantManager( object ):

    def __init__( self, variants, delete ):
        self.delete = delete
        self.variants = variants

    # Accepts variants based on whether or not all replicate frequencies are greater than some minimum.
    # returns a boolean list of length replicates which
    def freqLimit( self, variant, freqLim ) :
        outcomes = list()
        for number in variant.frequency :
            outcomes.append( number > freqLim )
        return outcomes, variant.frequency

    # Calculates whether the variants are greater than the expected percentage
    def fisherTest( self, variant, freqLim ) :
        outcomes = list()
        pvs = list()
        for j, number in enumerate( variant.coverage ) :
            number = int( number )
            remainder = variant.totalCounts[j]
            expectedPercent = freqLim

            # Values which can be changed are the fisher-exact test expected values, and the outcomes threshold.
            oddsratio, pvalue = fisher_exact( [[number, remainder], [expectedPercent * 100, (1 - expectedPercent) * 100]], alternative = "greater" )
            pvs.append( pvalue )
            outcomes.append( pvalue < pvLimit )

        return outcomes, pvs

    # Accepts variants based on whether or not the average frequency of all replicates is greater than some minimum.
    def averageFreqLimit( self, variant, freqLim ) :
        outcomes = list()
        outcomes.append( variant.getAverageFreq() > freqLim )
        return outcomes, variant.frequency

    statTest = { 1 : fisherTest, 2 : freqLimit, 3 : averageFreqLimit }

    def compareVariance( self, cutoff ):
        for entry in list( self.variants.keys() ):
            replicates = len( self.variants[entry] )
            if self.variants[entry].variance > cutoff and replicates > 1:
                if self.delete:
                    del self.variants[entry]

    def compareCoverage( self, cutoff, strict ):
        for entry in list( self.variants.keys() ):
            if strict:
                resultList = [x > cutoff for x in self.variants[entry].coverage ]
                if not all( resultList) and self.delete:
                    del self.variants[entry]
            else:
                if self.variants[entry].getAverageCov() < cutoff and self.delete:
                    del self.variants[entry]

    def compareFreq( self, frequencyLimit, test, pvLim ):

        # Data structures for graphing.
        # if test = 1, holds pvalues, else frequencies
        histogramList = list()
        # Holds frequencies mapped to genomic position
        manhattanDict = dict()
        # Holds frequencies which pass test.
        highlights = dict()

        for entry in list( self.variants.keys() ):
            if test != 1:
                oc, pv = self.statTest[test]( self, self.variants[entry], frequencyLimit )
            else:
                oc, pv = self.statTest[test]( self, self.variants[entry], frequencyLimit, pvLim )

            histogramList.extend( pv )

            manhattanDict[int( entry )] = self.variants[entry].frequency

            if all( oc ) :
                highlights[int( entry )] = self.variants[entry].frequency
            elif self.delete:
                del self.variants[entry]

        return histogramList, manhattanDict, highlights

    def getManagedVariants( self ):
        return self.variants
