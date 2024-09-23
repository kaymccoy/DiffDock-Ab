### REMINDER ###

#use "conda activate biopython"

#sometimes alignment segfaults; using a lazy list of things to skip b/c it's so infrequent



### IMPORTS ###

import argparse
import glob
import os
from Bio.PDB import PDBParser, PDBIO, CEAligner
io = PDBIO()
aligner = CEAligner()

import faulthandler
faulthandler.enable()



### FUNCTIONS ###

def load_structure(pdb_filename):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_filename)
    return structure

def write_structure(structure,out_name,io):
    io.set_structure(structure)
    io.save(out_name)

def extract_chains(structure, chain_ids):
    extracted_structure = structure.copy() 
    for model in extracted_structure:
        for chain_id in list(model):
            if chain_id.id not in chain_ids:
                model.detach_child(chain_id.id) 

    return extracted_structure

def alignTo(structMove, structStay,aligner):
    aligner.set_reference(structStay)
    aligner.align(structMove)
    return structMove



### ARGUMENTS ###

parser = argparse.ArgumentParser()
parser.add_argument("-b", dest="b", help="the base directory; should be an AADaM database dir that includes a summary of the db called lightDb.txt, has relaxed antigen structures in a dir called relaxedAntigens, and ImmuneBuilder predicted antibodies in a dir called immuneBuilderAbs") # 
parser.add_argument("-t", dest="t", help="use the true / bound structure for the antibody or antigen, if the antibody model / relaxed antigen failed for one but succeeded for the other (default = false, skip if either one wasn't succesfully made)")
parser.add_argument("-s", dest="s", help="only run on a subset of the structures, up to this number")
parser.add_argument("-o", dest="o", help="old lite database format")
arguments = parser.parse_args()

if arguments.s:
    maxNum = int(arguments.s)

outDir1 = arguments.b + "/alignedPairs"
if arguments.s:
    outDir1 += "_" + str(arguments.s)
if not os.path.exists(outDir1):
    os.mkdir(outDir1)
outDir2 = outDir1 + "/structures"
if not os.path.exists(outDir2):
    os.mkdir(outDir2)



### MAIN BODY ###
    
if arguments.o:
    segFaultSkipL = []
else:
    segFaultSkipL = ['7ssh','3q1s','6o2a']
    
csvFile = ['path,split\n']
extracted_filename = 'extracted_structure.pdb'

# read in summary file

abAgSum = []
abAgSumFile = arguments.b + "/lightDb.txt"
with open(abAgSumFile,"r") as sumFile:
    for line in sumFile:
        abAgSum.append(line.split(","))

header = True
pairsList = []
count = 0
agOnlyNotFound = 0
abOnlyNotFound = 0
for aSum in abAgSum:

    if header == True:
        header = False
        continue

    if arguments.s and count >= maxNum:
        break

    # grab the antigen chains
    agString = aSum[1].strip().replace(" | ","-")
    agChains = agString.split("-")

    # grab the ab chains

    if arguments.o:
        abString = aSum[2].strip().replace(" | ","-")
        abChains = abString.split("-")
    else:
        abChainH = aSum[2].strip()
        abChainL = aSum[3].strip()

        abChains = []
        if abChainH != "NA":
            abChains.append(abChainH)
        if abChainL != "NA":
            abChains.append(abChainL)
        abString = "-".join(abChains)

    # make out file names

    pdbID = aSum[0]
    if pdbID in segFaultSkipL:
        print("skipping because in segfault list...")
        continue
    csvName = pdbID + "_" + agString + "_" + abString
    agOutName = outDir2 + "/" + csvName + "_r_u.pdb"
    abOutName = outDir2 + "/" + csvName + "_l_u.pdb"

    if os.path.exists(agOutName) and os.path.exists(abOutName):
        print("skipping because out files already exist!")
        csvStr = csvName + ",test\n"
        csvFile.append(csvStr)
        continue

    # grab the true structure
    
    print("working on: " + pdbID,flush=True)
    if arguments.o:
        structPattern = arguments.b + "/structuresIMGT/" + pdbID + "*"
    else:
        structPattern = arguments.b + "/structures/" + pdbID + "*"
    structPath = glob.glob(structPattern)[0]
    trueStruct = load_structure(structPath)

    # get the true ag & ab structures

    agStruct = extract_chains(trueStruct,agChains)
    abStruct = extract_chains(trueStruct,abChains)
    # get the true ab structure

    # grab the ab model (if arguments.t and this isn't findable, use the true structure instead, so long as the ag was relaxed successfully)
        
    relaxedAgFound = False
    if arguments.o:
        relaxedAgPattern = arguments.b + "/afAntigens/" + pdbID + "*.pdb"
    else:
        relaxedAgPattern = arguments.b + "/relaxedAntigens/" + pdbID + "*" + agString + "_ag*.pdb" 
    relaxedAgList = glob.glob(relaxedAgPattern)
    if len(relaxedAgList) > 0:
        relaxedAgFound = True
        relaxedAg = load_structure(relaxedAgList[0])
    else:
        print("relaxed ag not found for pdb ID:" + pdbID + " with antigen chain(s): " + agString )
        agOnlyNotFound += 1

    # grab the relaxed ag (if arguments.t and this isn't findable, use the true structure instead, so long as the ab was modeled successfully)
    
    ibAbFound = False
    if arguments.o:
        ibAbPattern = arguments.b + "/ibAntibodies/" + pdbID + "*.pdb" 
    else:
        ibAbPattern = arguments.b + "/immuneBuilderAbs/" + pdbID + "*" + abString + "_.pdb" 
    ibAbList = glob.glob(ibAbPattern)
    if len(ibAbList) > 0:
        ibAbFound = True
        ibAb = load_structure(ibAbList[0])
    else:
        print("modeled ab not found for pdb ID:" + pdbID + " with antibody chain(s): " + abString )
        abOnlyNotFound += 1

    # continue / skip this PPI set, if necessary by t flag
        
    if (not relaxedAgFound) and (not ibAbFound):
        print("skipping because ag and ab not found...")
        agOnlyNotFound -= 1
        abOnlyNotFound -= 1
        continue
    elif (not arguments.t) and ((not relaxedAgFound) or (not ibAbFound)):
        continue
    
    count += 1 

    # align the ag struct to the original, if found

    agFail = False
    if relaxedAgFound:
        try:
            agToSave = alignTo(relaxedAg,agStruct,aligner)
            write_structure(agToSave,agOutName,io)
        except:
            agFail = True
    else:
        write_structure(agStruct,agOutName,io)

    # align the ab struct to the original, if found

    abFail = False
    if ibAbFound:
        try:
            abToSave = alignTo(ibAb,abStruct,aligner)
            write_structure(abToSave,abOutName,io)
        except:
            abFail = True
    else:
        write_structure(abStruct,abOutName,io)

    if abFail and agFail:
        continue
    elif abFail:
        if arguments.o:
            print("ab alignment failed...")
            quit()
        abOnlyNotFound += 1
        write_structure(abStruct,abOutName,io)
    elif agFail:
        if arguments.o:
            print("ag alignment failed...")
            quit()
        agOnlyNotFound += 1
        write_structure(agStruct,agOutName,io)

    # add a line to the csv file that summarizes them as test files
        
    csvStr = csvName + ",test\n"
    csvFile.append(csvStr)

print("agOnlyNotFound:")
print(agOnlyNotFound)
print("abOnlyNotFound:")
print(abOnlyNotFound)

# save the csv file too!

outFile = outDir1 + "/splits.csv"     
with open(outFile,'w') as of:
    of.writelines(csvFile)

### python3 /dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ab_ag_aligner.py -b /dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_earlyStructureDB_Jan9th2024dl_Mar16th2022endDate_9595seqIDs_5cutoff_NOX_strict_plus_Jan17th2023afterDate_8080seqIDs_5cutoff_NOX_strict -t 1 -s 250
    
### python3 /dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ab_ag_aligner.py -b /dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_earlyStructureDB_Jan9th2024dl_Mar16th2022endDate_9595seqIDs_5cutoff_NOX_strict_plus_Jan17th2023afterDate_8080seqIDs_5cutoff_NOX_strict -t 1 
    
### python3 /dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ab_ag_aligner.py -b /dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_earlyStructureDB_Jan9th2024dl_Mar16th2022endDate_9090seqIDs_5cutoff_NOX_notstrict_plus_Jan17th2023afterDate_8080seqIDs_5cutoff_NOX_strict -t 1 
    


### python3 /dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ab_ag_aligner.py -b /dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_SAbDab_01-17-22_all_c180_c290_c380_r3pt5_d2022-03-16_handParsed -t 1 -o 1