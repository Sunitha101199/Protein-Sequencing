"""
Protein Sequencing Project
Name: Sunitha
Roll Number: 2021501001
"""

from os import replace
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename,"r")
    dna = ""
    for i in f.read():
        for j in i.splitlines():
            dna+=j
    return dna


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    rna = []
    tempStr = str(dna[startIndex:]).replace("T","U")
    for i in range(0,len(tempStr),3):
        rna.append(tempStr[i:i+3])
        if tempStr[i:i+3]=="UAA" or tempStr[i:i+3]=="UAG" or tempStr[i:i+3]=="UGA":
            break
    return rna


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    codonToAmino = {}
    import json
    f = open(filename)
    data = json.load(f)
    for i in data:
        for j in data[i]:
            if "T" in j:
                j = str(j).replace("T","U")
            codonToAmino[j]=i
    return codonToAmino


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    aminoAcid = []
    for i in codons:
        if i=="AUG" and "Start" not in aminoAcid:
            aminoAcid.append("Start")
        elif i=="UAA" or i == "UAG" or i == "UGA":
            aminoAcid.append("Stop")
        else:
            aminoAcid.append(codonD[i])
    return aminoAcid


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna = readFile(dnaFilename)
    # bases, unused, synthesized = len(dna), 0, 0
    codonToAmino = makeCodonDictionary(codonFilename)
    proteins = []
    i=0
    while i!=len(dna):
        rna = []
        if dna[i:i+3]=="ATG":
            rna = dnaToRna(dna, i)
            proteins.append(generateProtein(rna, codonToAmino))
            i += 3*len(rna)
            # synthesized+=1
        else:
            # unused+=1
            i+=1
    # print(bases, unused, synthesized)
    return proteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common_Proteins = []
    for i in proteinList1:
        if i in proteinList2 and i not in common_Proteins:
                common_Proteins.append(i)
    return common_Proteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combine_Proteins = [j for i in proteinList for j in i]
    return combine_Proteins


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoAcidFrequency = {}
    for i in aaList:
        if i not in aminoAcidFrequency:
            aminoAcidFrequency[i]=0
        aminoAcidFrequency[i]+=1
    return aminoAcidFrequency


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    comPro1 = combineProteins(proteinList1)
    comPro2 = combineProteins(proteinList2)
    aminoAcidFrequency1 = aminoAcidDictionary(comPro1)
    aminoAcidFrequency2 = aminoAcidDictionary(comPro2)
    aminoAcidDifferences = []
    unique = []
    for i in comPro1+comPro2:
        if i not in unique and i!="Start" and i!="Stop":
            unique.append(i)
    print(unique)
    for i in unique:
        li = []
        if i in aminoAcidFrequency1 and i not in aminoAcidFrequency2:
            if abs(aminoAcidFrequency1[i]/len(comPro1))>cutoff:
                li.append(i)
                li.append(abs(aminoAcidFrequency1[i]/len(comPro1)))
                li.append(0)
        elif i in aminoAcidFrequency2 and i not in aminoAcidFrequency1:
            if abs(aminoAcidFrequency2[i]/len(comPro2))>cutoff:
                li.append(i)
                li.append(0)
                li.append(abs(aminoAcidFrequency2[i]/len(comPro2)))
        else:
            if abs((aminoAcidFrequency1[i]/len(comPro1)) - (aminoAcidFrequency2[i]/len(comPro2)))>cutoff:
                li.append(i)
                li.append(abs(aminoAcidFrequency1[i]/len(comPro1)))
                li.append(abs(aminoAcidFrequency2[i]/len(comPro2)))
        if len(li)==3:
            aminoAcidDifferences.append(li)
    return aminoAcidDifferences


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("These are the proteins in two genes with length more than two.")
    for i in commonalities:
        if len(i)>2:
            print("-".join(i[1:len(i)-1]))
    print("These are the amino acids that occurred at the most different rates.")
    for i in differences:
        print(i[0]+":",round(i[1],2),"%","in Sep1, ",round(i[2],2),"% ","in Sep2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
 
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
 
