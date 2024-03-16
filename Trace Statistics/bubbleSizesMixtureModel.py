#importing sklearn

from sklearn.model_selection import train_test_split as split
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot as plt

import numpy as np
from numpy import e as e
from multiBubbleDistributionStats import fileListParser as parser

def ashmansD(mus = [], stdevs = []):
    """Public function to calculate Ashamn's D, to determine how cleanly a
    distribution can be parsed into a bimodal one."""
    #idk if .5 works
    return np.abs(mus[1]-mus[0])/(np.sqrt(2*(  stdevs[1]**2  +  stdevs[0]**2))   )



#creating a set of subplots, one for the predicted distribution and one for the true one
#this is for both of 2 populations

figure, axes = plt.subplots(4)
figure.set_size_inches(4,4, forward = True)

plt.subplots_adjust(bottom=0.05,
                    top=.9,
                    wspace=0.4,
                    hspace=0.7)


#parse in data
thisParser = parser()

#get bubbleSizes dictionaries
thisParser.parseFiles()
theseBubbleSizes = thisParser.parseFilesTTest(0, windowEnd=460)

nameAndBubbles = theseBubbleSizes[0]

averageScale = theseBubbleSizes[1]

#preprocessing data from bubble size dictionary, natural log transform
logTransformDict = {}


for name,bubbles in nameAndBubbles.items():



    #saving dictionary value as a tuple of the log transform of the data and the proper amount of bins to use in log space
    logTransformDict[name] = bubbles

"""Important:
We are assuming that both distributions are bimodal...see paper/ask Hinrich for rationale. """
#fitting GM models for both populations using the same parameters
counter =  0
for population,sizes in logTransformDict.items():
    name = population
    sizes = sizes


    #reshaping array since there is only a single feature
    sizes = np.array(sizes).reshape(-1,1)
    # test/train split; 20/80
    train,test = split(sizes, test_size=.2, train_size = .8)

    # logTransform  = [np.log(bubble) for bubble in bubbles]
    # calculating the proper amount of bins to use from the split, using the test set
    properBins = int( (max(test) - min(test) ) / averageScale)
    theseBins = properBins

    #fitting to the GMM; we are assuming bimodal; that is, n_components = 2
    thisModel = GaussianMixture(n_components=2).fit(train)

    #plotting the predictions as PDFs again
    #prediction =  thisModel.predict(test)

    #getting the scores and converting to actual probabilities
    #taking the actual probabilities from the score_samples

    actualProbabilities = thisModel.predict_proba(test)


    #plotting actual probabilities with the density parameter off

    axes[counter].hist(actualProbabilities, bins=theseBins, density = True)
    print("these bins", theseBins)

    #axes[counter].plot([item[1] for item in actualProbabilities])


    #titling
    axes[counter].set_title(name)

    #print("first", np.histogram(actualProbabilities, bins =30))
    #trying to autoscale the axes
    #axes[counter].autoscale()

    # plotting the log transform of the test data as a histogram on one of the subplots
    
    axes[counter+1].hist(test, bins=theseBins, density = True)
    axes[counter+1].set_title(name + " Test Set")

    counter += 2

    #reporting the estimated paramaters, most notably, the modes, along with the population name
    thisModes = [item for item in thisModel.means_]


    #calculating Ashman's D for both of the averages to see if the distributions are truly bimodal
    #covariances matrix and taking the trace of the matrix to gets the variances
    covariances = thisModel.covariances_
    #print("modes", thisModes, "covariances", covariances)

    thisStdevs = [np.sqrt(np.trace(covariances[i]) / 2) for i in range(0, 2)]
    popAshman = ashmansD(thisModes,thisStdevs)
    print("Ashman",popAshman )


plt.show()