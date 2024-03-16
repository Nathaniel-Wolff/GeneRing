#!/usr/bin/env python
import sys
from matplotlib import pyplot as plt
import scipy
from scipy import stats as stats
import argparse

import bisect
import os.path
import os
import numpy as np

import Trace
from Trace import Trace


__version__="01.00.00"
__author__ ="Nathaniel Wolff"

"""This module uses the same argparse params as TraceAnalyzer module to read in a list of trace files.
    The each raw trace file is converted to a generator list of tuples using a global parsing function.
    The Moleculify method uses the Chromatin module to returns a Molecule object. 
    This Molecule object can be queried to produce its bubbles. 
    """
class fileListParser:
    """
    """
    def __init__(self):
        self.params = {}
    def parseRawTrace(self, file):
        """Global function to parse a raw trace to a format compatible with either this module or TraceAnalyzer."""
        with open(file, "r") as file:
            openArray = np.array(np.array([line.rstrip().split("\t") for line in file.readlines()], dtype=float))
            return openArray
    def queryKmer(self, t, k, thisKmer):
        """ 1) Returns substring index (locations) of a character/kmer (thisKmer) of size k in text = t -> substring index
        """
        ''' Create index from all substrings of size 'length' '''

        t = t  #text to be searched...in the usage below, this is a directory
        k = k  #k-mer length (k)

        index = []
        # creates a substring index w/parameters t = text and k = substring size
        for offset in range(0, len(t) - k + 1):
            # enumerating kmers
            kmer = t[offset: offset + k]
            # saving them and their offsets as tuples
            index.append((kmer, offset))

        # sorting the substring index
        index.sort()

        #take sorted substring index, kmer and start binary search
        assert len(thisKmer) == k
        hits = []

        #Find first location of kmer in self.index (hint: use bisect.bisect_left function)

        # first needs to run through an alphabetized list of the substrings and find the first hit
        justSubStrings = [tuple[0] for tuple in index]
        firstHit = bisect.bisect_left(justSubStrings, thisKmer)

        # using the firstHit index (which gives the index in self.index) to find the first occurrence in the text
        firstOccurence = index[firstHit][1]
        hits.append(firstOccurence)

        # Iterate through self.index from first location of kmer to last adding matches to hits
        for tuple in index[firstHit:]:
            if tuple[0] == thisKmer:
                # add the occurrences to hits
                hits.append(tuple[1])
                # print(hits)
        return hits

    def parseFiles(self):
        """

        *Now with os.walk
        Based on Robert's TraceAnalyzer
        1) Receives a list of files name from argparse.
        2) Parses arguments and returns error message to the user if the user supplied no files.
        3) Saves a dict of the params to the instance of the class (contains a list of the raw trace file names under key "files")
        Returns a dict of all the files using vars
        """

        DEFAULT_LENGTH = 2246
        DEFAULT_THRESHOLD = 4
        DEFAULT_SMOOTH = 10

        parser = argparse.ArgumentParser(description='Analyze Chromatin Ring Trace Files.')
        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

        #Population only takes 1 arg, the outermost directoru.
        #This means the user must have all of the subdirectories containing raw traces in the same folder, whose name is population.
        parser.add_argument('populations',
                            type=str,
                            help='The outermost directory that contains the trace files you wished to analyze.')

        parser.add_argument('-t', '--threshold',
                            default=DEFAULT_THRESHOLD,
                            type=float,
                            help='The threshold distance required for two points on opposing strands to be considered as part of a linker.')

        parser.add_argument('-l', '--length',
                            default=DEFAULT_LENGTH,
                            type=int,
                            help='The expected length in basepairs for the DNA molecule.')

        parser.add_argument('-s', '--smooth',
                            default=DEFAULT_SMOOTH,
                            type=int,
                            help='The # of coordinates to include in a sliding window of the average distance between points on opposing strands.')

        parser.add_argument('-u', '--user',
                            default=None,
                            type=str,
                            help='The name of the person who completed the trace and is using this software: for record keeping.')

        parser.add_argument('-o', '--out_path',
                            default=None,
                            type=str,
                            help='The name of the folder you wish to save the output to.')

        parser.add_argument('-p', '--plotless',
                            action='store_true')

        parser.add_argument('-i', '--imgres',
                            default=scipy.NAN,
                            type=scipy.float16,
                            help='The image resolution of the raw image used to make the trace.')
        #now taes the names of the folders whose distributions are to be produced
        #f flag, for consistency
        parser.add_argument('-d', '--directory',
                            action='store_false',
                            default='true',
                            )

        args = parser.parse_args()
        currentDir = os.getcwd()
        #Consolidating all the text files from the innermost directory into one list
        goalDirectory = currentDir + "/"

        counter = 0
        x = 3

        #iterating over each population member and consolidating stuff
        for population in args.populations.split(","):
            files = []
            for (dirpath, dirnames, filenames) in os.walk(population, topdown=False):


                #goal here is to see if each dirpath contains "/" x amount of times (x is currently hardcoded; later it will be given by the user)

                #approach:
                #generate a substring index for each dirpath as list sorted in lexographic order
                # get the "/" hits from the substring index and check its length
                # if it == x, then read in the text files from that directory

                if len(self.queryKmer(t = dirpath, k = 1, thisKmer = "/")) == x:
                    filesAndDirs = [self.parseRawTrace(dirpath+"/"+file) for file in filenames if file.endswith(".txt")]
                    files += filesAndDirs

            parsedFiles = files
            self.params[population] = parsedFiles

            if len(parsedFiles) == 0:
                sys.stderr.write('No Trace Files Found. Check the file path.')
                sys.exit()

        return self.params
    def parseFilesTTest(self, windowStart, windowEnd):
        #did not use infoParse class function bc cleaning it might be more trouble than its worth
        """Takes the params dict and
        1) Opens a raw trace file and parses each one to a generator"""
        populations = list(self.params.keys())

        bubbleSizes = {}
        scaleFactors = []
        tracesPerPopulation = {}
        for population, parsedFiles in self.params.items():
            parsedFiles = parsedFiles
            population = population

            popBubbleSizes = []
            noTracesPerPopulation = 0

            for parsedFile in parsedFiles:
                #parsing trace file to list
                #Building Trace object -> Chromatin.Molecule object -> Bubble Sizes Generator -> Bubble Sizes
                # traces are about the double the length of PHO5, about 2246 BP
                thisTrace = Trace(parsedFile)
                scaleFactor = thisTrace.scale(thisTrace.length)
                scaleFactors.append(scaleFactor)
                noTracesPerPopulation += 1
                #getting set of bubble sizes and adding to the total bubble sizes
                popBubbleSizes += thisTrace.returnLocalBubblesOrLinkers((windowStart,windowEnd))

            #saving the number of traces in a population to the dictionary
            tracesPerPopulation[population] = noTracesPerPopulation

            #saving population's bubble sizes to the dictionary
            bubbleSizes[population] =  popBubbleSizes

        #MannWhitneyUTest
        #now is used
        mwuStat,pvalue = stats.mannwhitneyu(list(bubbleSizes.items())[1][1], list(bubbleSizes.items())[0][1])
        alpha = .05

        if pvalue < alpha:
            sys.stdout.write("null hypothesis rejected; significant difference between " + populations[0] + " and " + populations[1] + "."
                             + "\n" + "P value:" + str(pvalue))
        else:
            sys.stdout.write(
                "null hypothesis not rejected; no significant difference between " + populations[0] + " and " + populations[
                    1] + "over the window: " + str(windowStart) + "-" +  str(windowEnd) + "."
                + "\n" + "P value:" + str(pvalue))
        #print("samples per population", tracesPerPopulation)
        #calculating the average scale to be used in plotting
        # calculating the proper number of bins to use, that is, the max - min of the data divided by scale
        averageScale = int(np.average(scaleFactors))

        return bubbleSizes, averageScale

    def plotDistributions(self, distributionDict, averageScale):
        """Method for plotting all the distributions in the same figure."""

        #setting up figure
        figureHeight = 4
        figureWidth = 4
        figureSize = (figureHeight,figureWidth)
        ax = plt.figure(figsize = figureSize)

        # set up 2 panels of the same size
        panelWidth = 3
        panelHeight = 1

        firstLeft = .5/figureWidth
        secondLeft = .5/figureWidth

        firstBottom = .375 / figureHeight
        secondBottom = (.375+panelHeight+1)/figureHeight

        leftOffset = 0
        bottomOffset = .025


        firstMainPanel = plt.axes([firstLeft, firstBottom, panelWidth/figureWidth, panelHeight/figureHeight])
        #firstMainPanel.set_title(list(distributionDict.keys())[0])
        firstMainPanel.set_title("Activated")
        firstMainPanel.set_ylabel("Probability Density")
        firstMainPanel.set_xlabel("Bubble Size (BP)")

        secondMainPanel = plt.axes([secondLeft, secondBottom, panelWidth/figureWidth, panelHeight/figureHeight], sharex = firstMainPanel, sharey = firstMainPanel ) #label = list(distributionDict.keys())[1])
        secondMainPanel.set_title(str(list(distributionDict.keys())[1]))
        secondMainPanel.set_ylabel("Probability Density")
        secondMainPanel.set_xlabel("Bubble Size (BP)")

        panels = [firstMainPanel, secondMainPanel]

        #textAnchors = [(firstLeft-leftOffset, firstBottom-bottomOffset), (secondLeft-leftOffset,secondBottom-bottomOffset)]

        #plotting the distributions in the respective panels as normalized heatmaps
        for nameAndDistribution,panel in zip(list(distributionDict.items()), panels):
            population,bubbleSizes = nameAndDistribution[0], nameAndDistribution[1]
            #calculating the proper number of bins to use, that is, the max - min of the data divided by scale
            properBins = int((max(bubbleSizes) - min(bubbleSizes))/averageScale)
            oneHistogram = np.histogram(bubbleSizes)

            # trying to check what the distribution is in the range of small bubble sizes
            panel.hist(bubbleSizes, density = True, color = "green", bins = properBins, range = (30,500) )

            plt.xticks([0,50,100,150,200,250,300,350,400,450,500])
            panel.set_xlim(0,500)

        #plt.title("Bubble Size Distribution of PHO5 Active and Repressed Strains")
        plt.show()

def main():
    testFileName = "GB26_E1_3.txt"
    folder = "Activated"
    thisFileListParser = fileListParser()
    thisParams = thisFileListParser.parseFiles()
    #change the window later to command line; using the UCSC's genome browser's 0-468 promoter length for PHO5

    distributionsDict,scaleFactors = thisFileListParser.parseFilesTTest(0, 500)
    thisFileListParser.plotDistributions(distributionsDict,scaleFactors)


if __name__ == "__main__":
        main()
