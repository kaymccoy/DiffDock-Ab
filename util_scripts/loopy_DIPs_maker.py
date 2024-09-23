# imports

import sys
import os
import dill
import glob
sys.path.append('/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib')
import mstpython as mst

distCut = 10.0
loopyCut = 0.65
outName = "loopyDIPS-" + str(loopyCut)
outDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/" + outName
if not os.path.exists(outDir):
    os.mkdir(outDir)

outLines = []
outPath = outDir + "/splits.csv"
outStructPath = outDir + "/structures"
if not os.path.exists(outStructPath):
    os.mkdir(outStructPath)

baseDillPath = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP-old/datasets/DIPS/pairs_pruned/"
basePdbPath = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/DIPS_pruned_raw/"



# function to get the chains for the pair in the original training set

def getPairChains(dillFile):
    chains1 = []
    chains2 = []
    with open(dillFile, 'rb') as in_strm:
        atomPairPair = dill.load(in_strm)
    chainCol1 = atomPairPair[1].loc[:,"chain"]
    chainCol2 = atomPairPair[2].loc[:,"chain"]
    for i in chainCol1:
        if i not in chains1:
            chains1.append(i)
    for i in chainCol2:
        if i not in chains2:
            chains2.append(i)
    finalChains = [chains1,chains2]
    return(finalChains)


# function for dictionary of helix / beta / loop

def geomDicFunc(pdbPath):

    outDic = {}

    with open(pdbPath,"r") as pdbF:
        pdbLines = pdbF.readlines()

    for i in pdbLines:
        if i[0:5] == 'HELIX':
            chainStart = i[19]
            numStart = int(i[21:25])
            iCodeStart = i[25]
            chainEnd = i[31]
            numEnd = int(i[33:37])
            iCodeEnd = i[37]
            if chainStart != chainEnd:
                print("helix or sheet stretches across multiple chains? Skipping as unclear how to classify...")
                return 0
            else:
                addLine = [numStart,iCodeStart,numEnd,iCodeEnd]
                if chainStart in outDic:
                    outDic[chainStart].append(addLine)
                else:
                    outDic[chainStart] = [addLine]

        elif i[0:5] == 'SHEET': #***DIFF THAN HELIX, FIX!
            chainStart = i[21]
            numStart = int(i[22:26])
            iCodeStart = i[26]
            chainEnd = i[32]
            numEnd = int(i[33:37])
            iCodeEnd = i[37]
            if chainStart != chainEnd:
                print("helix or sheet stretches across multiple chains? Skipping as unclear how to classify...")
                return 0
            else:
                addLine = [numStart,iCodeStart,numEnd,iCodeEnd]
                if chainStart in outDic:
                    outDic[chainStart].append(addLine)
                else:
                    outDic[chainStart] = [addLine]

    return outDic

def isLoopy(rChain,rNum,rIcode,geomDic): # get whether or not a residue's loopy according to its chain / num, checked in the dict

     # check if chain is a key; if not, return loopy
    if rChain in geomDic:

        # if so, check if the number is within one of the ranges; if not, return loopy, if so, return not loopy, if right on the edge, check if the iCode is before / after appropriately to chose what to return!

        checkList = geomDic[rChain]
        for i in checkList:
            numStart = i[0]
            iCodeStart = i[1] # will be a space (' ') if no iCode
            numEnd = i[2]
            iCodeEnd = i[3]
            if (rNum >= numStart) and (rNum <= numEnd):
                if (rIcode >= iCodeStart) and (rIcode <= iCodeEnd):
                    return False
        return True

    else:
        return True

   

# open the csv used for testing / training...

with open("/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/DIPS/data_file.csv", "r") as df:
    dipsBase = df.readlines()

conResListB = []
conResListLarge = []

# iterate over each csv entry... and get the loopiness of contact residues

cl = mst.ContactList()

count = 0
for d in dipsBase:

    if count%1000 == 0:
        print(str(count) + " done...")
    count += 1

    loopNum = 0
    loopDenom = 0

    # find its dill path

    fileBase = d.split(",")[0]
    if fileBase == 'path':
        continue
    dillPath = baseDillPath + fileBase

    # use that to find the partner chains

    partnerChains = getPairChains(dillPath)

    # find its raw PDB path according to PDB ID and name of the dill file 

    pdbID = fileBase.split("/")[-1].split(".pdb")[0]
    pdbGlobStr = basePdbPath + pdbID + "*.pdb"

    pdbFiles = glob.glob(pdbGlobStr)
    pdbFile = ""

    foundFile = False
    for i in pdbFiles:
        chainNames = i.split("~")[-1].split(".pdb")[0]
        baseChainList = chainNames.split("_")
        chainList = []
        for ii in baseChainList:
            newList = ii.split("-")
            chainList.append(newList)
        if chainList == partnerChains:
            pdbFile = i
            foundFile = True
            break

    if not foundFile:
        print("failed at finding file...")
        print("testing fileBase:")
        print(fileBase)
        print("testing partnerChains:")
        print(partnerChains)
        print("testing pdbGlobStr:")
        print(pdbGlobStr)
        print("testing pdbFiles:")
        print(pdbFiles)
        continue

    # read in PDB file --> make a dictionary of whether the residue's helix, sheet, or loop (dict of dicts with chain as first key, residue number + iCode as second key)

    geomDic = geomDicFunc(pdbFile) 
    if geomDic == 0:
        continue

    # get A & B atoms from chain(s) listed in partnerChains

    AB = mst.Structure(pdbFile, "QUIET")
    A = mst.emptyStructure()
    B = mst.emptyStructure()

    for i in partnerChains[0]:
        aChain = AB.getChainByID(i)
        A.appendChain(aChain)

    for i in partnerChains[1]:
        bChain = AB.getChainByID(i)
        B.appendChain(bChain)

    aReses = A.getResidues()
    bReses = B.getResidues()

    #cfB = mst.ConFind("/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/testfiles/rotlib.bin", B)

    conResListA = []
    loopNumA = 0
    totNumA = 0
    conResListB = []
    loopNumB = 0
    totNumB = 0

    i = 0
    while i < len(aReses):

        conRes = aReses[i]
        i += 1

        try:
            conCA = conRes.findAtom("CA",True)
            conCoor = conCA.getCoor()
        except:
            continue

        ii = 0
        while ii < len(bReses):
            compRes = bReses[ii]
            ii += 1

            try:
                compConCA = compRes.findAtom("CA",True)
                compCoor = compConCA.getCoor()
            except:
                continue

            compDist = conCoor.distance(compCoor)
            if compDist <= distCut:

                if conRes not in conResListA:
                    conResListA.append(conRes)
                    conResChain = conRes.getChainID(True)
                    conResNum = conRes.num
                    conResIcode = conRes.iCode

                    totNumA += 1
                    if isLoopy(conResChain,conResNum,conResIcode,geomDic): 
                        loopNumA += 1

                if compRes not in conResListB:
                    conResListB.append(compRes)
                    bConResChain = compRes.getChainID(True)
                    bConResNum = compRes.num
                    bConResIcode = compRes.iCode

                    totNumB += 1
                    if isLoopy(bConResChain,bConResNum,bConResIcode,geomDic):
                        loopNumB += 1

    # make csv entry + another column that has loopy %

    if totNumA == 0:
        continue
    if totNumB == 0:
        continue

    loopyPerA = float(loopNumA)/totNumA    
    loopyPerB = float(loopNumB)/totNumB
    loopyPer = max(loopyPerA,loopyPerB)
    if loopyPer < loopyCut:
        continue
    else:
        
        chainStr = ""
        for chE in chainList:
            for ch in chE:
                chainStr += "_" + ch

        nameBase = pdbID + chainStr

        if loopyPerA >= loopyPerB:
            outPathA = outStructPath + "/" + nameBase + "_l_u.pdb"
            outPathB = outStructPath + "/" + nameBase + "_r_u.pdb"
        else:
            outPathA = outStructPath + "/" + nameBase + "_r_u.pdb"
            outPathB = outStructPath + "/" + nameBase + "_l_u.pdb"

        A.writePDB(outPathA,"QUIET")
        B.writePDB(outPathB,"QUIET")

        outLine = nameBase + ",test," + str(loopyPer) + "\n"
        outLines.append(outLine)

# at end, save csv file with loopy-ness % details!
  
with open(outPath,"w") as of:
    of.writelines(outLines)