import argparse
import math
import numpy
import pandas



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="solFile", help="input the soluble sequence emissions here")
    parser.add_argument(dest="transMemFile", help="input transmembrane file here")
    parser.add_argument(dest="stateSeq", help="input the state sequence file here")
    parser.add_argument(dest = "analysisSeq", help="input the sequence file to be analyzed here")
    args = parser.parse_args()

    solubleFile = open(args.solFile, "r")
    transMembraneFile = open(args.transMemFile, "r")
    stateFrequenciesFile = open(args.stateSeq, "r")
    analysisFile = open(args.analysisSeq, "r")
    global solubleFrequencies
    solubleFrequencies = {}
    global transMembraneFrequencies
    transMembraneFrequencies = {}
    global stateFrequencies
    stateFrequencies = {}
    global solDict
    solDict = {}
    global transDict
    transDict = {}
    global stateDict
    stateDict = {}
    global solProb
    solProb = []
    global transProb
    transProb = []
    analysisSequence = analysisFile.read()
    analysisSequence = analysisSequence.rstrip('\r\n').replace("\n", '')
    global scoresMatrix
    scoresMatrix = numpy.zeros((2, len(analysisSequence)), dtype=float, order='C')
    global transProbFinal
    global solProbFinal
    global finalState

    global currentState
    currentState = "I"
    def getSolubleFrequencies():
        solubleCode = solubleFile.read()
        solubleCode = solubleCode.rstrip('\r\n').replace("\n",'')
        size = 1
        i = 0
        solNuc = ""
        while i<=(len(solubleCode)-1):
            solNuc = solubleCode[i:(i+size)]
            if solNuc in solubleFrequencies:
                solubleFrequencies[solNuc] +=1
            else:
                solubleFrequencies[solNuc] = 1
            i = i + 1
        solTotals = 0
        for i in solubleFrequencies:
            solTotals = solTotals + solubleFrequencies[i]
        for i in solubleFrequencies:
            temp = solubleFrequencies[i]
            probability = temp / solTotals
            solubleFrequencies[i] = probability
        for i in solubleFrequencies:
            solDict[i] = solubleFrequencies[i]
    getSolubleFrequencies()

    def getTransMembraneFrequencies():
        transMembraneCode = transMembraneFile.read()
        transMembraneCode = transMembraneCode.rstrip('\r\n').replace("\n",'')
        transMembraneFrequencies = {}
        size = 1
        i = 0
        transNuc = ""
        while i<=(len(transMembraneCode)-1):
            transNuc = transMembraneCode[i:(i+size)]
            if transNuc in transMembraneFrequencies:
                transMembraneFrequencies[transNuc] +=1
            else:
                transMembraneFrequencies[transNuc] = 1
            i = i + 1
        transTotals = 0
        for i in transMembraneFrequencies:
            transTotals = transTotals + transMembraneFrequencies[i]
        for i in transMembraneFrequencies:
            temp = transMembraneFrequencies[i]
            probability = temp / transTotals
            transMembraneFrequencies[i] = probability
        for i in transMembraneFrequencies:
            transDict[i] = transMembraneFrequencies[i]
    getTransMembraneFrequencies()

    def getStateFrequencies():
        stateCode = stateFrequenciesFile.read()
        for line in stateCode:
            temp = ""
            if(line[0:1] != '\n'):
                temp = "B" + line[0:1]
                if temp in stateFrequencies:
                    stateFrequencies[temp] +=1
                else:
                    stateFrequencies[temp] = 1
        for line in stateCode:
            temp = ""
            if(line[len(line)-1: len(line)] != '\n'):
                temp = line[len(line)-1: len(line)] + "E"
                if temp in stateFrequencies:
                    stateFrequencies[temp] +=1
                else:
                    stateFrequencies[temp] = 1
        stateCode = stateCode.rstrip('\r\n').replace("\n",'')

        size = 2
        i = 0
        stateNuc = ""
        while i<=(len(stateCode)-2):
            stateNuc = stateCode[i: (i+2)]
            if stateNuc in stateFrequencies:
                stateFrequencies[stateNuc] +=1
            else:
                stateFrequencies[stateNuc] = 1
            i = i + 1
        stateTotals = 0
        startTotals = 0
        endTotals = 0
        for i in stateFrequencies:
            if (i == 'BS' or i == 'BT'):
                startTotals = startTotals + stateFrequencies[i]
            else:
                if (i == 'SE' or i == 'TE'):
                    endTotals = endTotals + stateFrequencies[i]
                else:
                    stateTotals = stateTotals + stateFrequencies[i]
        for i in stateFrequencies:
            temp = stateFrequencies[i]
            if (i == 'BS' or i == 'BT'):
                startProb = temp / startTotals
                stateDict[i] = startProb
            else:
                if (i == 'SE' or i == 'TE'):
                    endProb = temp / endTotals
                    stateDict[i] = endProb
                else:
                    probability = temp / stateTotals
                    stateDict[i] = probability
    getStateFrequencies()
    def getStateSequence():
        startCode = analysisSequence[0:1]
        x = 0
        while x<len(analysisSequence):
            transProb.append(0)
            x = x +1
        transStartAcidProb = transDict[startCode]
        transProb[0] = math.log(transStartAcidProb,2)
        i = 1
        while i<len(analysisSequence):
            currentCode = analysisSequence[i:i+1]
            currentStateProb = math.log(transDict[currentCode],2)
            lastProb = math.log(transStartAcidProb,2)
            conditionalProb = (currentStateProb * lastProb)/lastProb
            transProb[i] = conditionalProb + transProb[i-1]
            lastProb = currentStateProb
            i +=1
        solStartAcidProb = math.log(solDict[startCode],2)
        n = 0
        while n<len(analysisSequence):
            solProb.append(0)
            n = n + 1
        solProb[0] = solStartAcidProb
        k = 0
        while k<len(analysisSequence):
            curCode = analysisSequence[k:k+1]
            curStateProb = math.log(solDict[curCode],2)
            if(k == 0):
                latestProb = solStartAcidProb
            if(k != 0):
                latestProb = solProb[k-1]
            condProb = (curStateProb * latestProb)/latestProb
            solProb[k] = condProb + solProb[k-1]
            latestProb = curStateProb
            k = k + 1
        a = 0
        while a<len(transProb):
            scoresMatrix[0,a] = transProb[a]
            a = a +1
        b = 0
        while b<len(solProb):
            scoresMatrix[1,b] = solProb[b]
            b = b + 1
        transProbFinal = transProb[len(transProb)-1]
        solProbFinal = solProb[len(solProb)-1]
    getStateSequence()
    def getFinalState():
        numberOfElements = len(solProb)
        z = numberOfElements -1
        finalState = numpy.zeros((1, numberOfElements), dtype= str, order='C')
        while z > -1:
            solTemp = scoresMatrix[0,z]
            transTemp = scoresMatrix[1,z]
            if(solTemp>transTemp):
                finalState[0,z] = 'S'
            if(transTemp>solTemp):
                finalState[0,z] = 'T'
            z = z-1
        print(finalState)
    getFinalState()
main()