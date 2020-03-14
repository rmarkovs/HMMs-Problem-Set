import argparse
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
    #print(solubleFile)
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
    global currentState
    currentState = "I"
    #Solublecode = file.read()

    #code = code.rstrip('\r\n').replace("\n", '')
    #code = code.capitalize()

    #nucleotides = {}
    #size = int(args.kmer)
    #i = 0
    #nuc = ""

    # while i <= (len(code) - 1 - size):
    #     nuc = code[i:(i + size)]
    #     if nuc in nucleotides:
    #         nucleotides[nuc] += 1
    #     else:
    #         nucleotides[nuc] = 1
    #     i = i + 1

    # outfile = open("results.txt", "w")
    # for key in nucleotides:
    #     outfile.write(key + "   " + str(nucleotides[key]))
    #     outfile.write("\n")
    def getSolubleFrequencies():
        solubleCode = solubleFile.read()
        solubleCode = solubleCode.rstrip('\r\n').replace("\n",'')
        #print(solubleCode)
        #solubleFrequencies = {}
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
            #print(solNuc)
            #print(i)
        #print(solubleFrequencies)
        #print(solubleFrequencies)
        return solubleFrequencies
        #print(solubleFrequencies)
    getSolubleFrequencies()

    def getTransMembraneFrequencies():
        transMembraneCode = transMembraneFile.read()
        transMembraneCode = transMembraneCode.rstrip('\r\n').replace("\n",'')
        #transMembraneFrequencies = {}
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
        #print(transMembraneFrequencies)
        return transMembraneFrequencies
    getTransMembraneFrequencies()

    def getStateFrequencies():
        #may need to go line by line
        #stateFrequencies = {}

        stateCode = stateFrequenciesFile.read()
        #added in to get start probabilities
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
        #stateCode = stateFrequenciesFile.rstrip('\r\n').replace("\n",'')
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
        #print(stateFrequencies)
        return stateFrequencies
    getStateFrequencies()
    def getProbabilities():
        #solDict = {}
        #print(getSolubleFrequencies())
        solDict = getSolubleFrequencies()
        #solDict = getSolubleFrequencies()
        #print(solDict)
        solTotals = 0
        for i in solDict:
            solTotals = solTotals + solDict[i]
        for i in solDict:
            temp = solDict[i]
            probability = temp/solTotals
            solDict[i] = probability
        #print(solDict)
        transDict = getTransMembraneFrequencies()
        #print(transDict)
        transTotals = 0
        for i in transDict:
            transTotals = transTotals + transDict[i]
        for i in transDict:
            temp = transDict[i]
            probability = temp/transTotals
            transDict[i] = probability
        #print(transDict)
        stateDict = getStateFrequencies()
        #print(stateDict)
        stateTotals = 0
        startTotals = 0
        endTotals = 0
        for i in stateDict:
            if(i == 'BS' or i == 'BT'):
                startTotals = startTotals + stateDict[i]
            else:
                if(i == 'SE' or i == 'TE'):
                    endTotals = endTotals + stateDict[i]
                else:
                    stateTotals = stateTotals + stateDict[i]
        for i in stateDict:
            temp = stateDict[i]
            if (i == 'BS' or i == 'BT'):
                startProb = temp/startTotals
                stateDict[i] = startProb
            else:
                if (i == 'SE' or i == 'TE'):
                    endProb = temp/endTotals
                    stateDict[i] = endProb
                else:
                    probability = temp/stateTotals
                    stateDict[i] = probability
        #print(stateDict)
    getProbabilities()
    def getStateSequence():
        analysisSequence = analysisFile.read()
        analysisSequence = analysisSequence.rstrip('\r\n').replace("\n",'')
        stateSequence = ""
        #startCode = ""
        startCode = analysisSequence[0:1]

        # startAcid = ""
        # startAcid = startCode
        probTransStart = 0
        transStartAcidProb = transDict[startCode]
        conditionalProbTransStart = (transStartAcidProb + stateDict["BT"] - (transStartAcidProb * stateDict["BT"]))/stateDict["BT"]
        probTransStart = stateDict["BT"] * conditionalProbTransStart
        solStartAcidProb = solDict[startCode]
        conditionalProbSolStart = (solStartAcidProb + stateDict["BS"] - (solStartAcidProb * stateDict["BS"]))/stateDict["BS"]
        probSolStart = stateDict["BS"] * conditionalProbSolStart



    getStateSequence()
    #def calculateStateProbability(code, )
    # def probabilityChain():
    #
    # probabilityChain()
main()