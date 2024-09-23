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

def get_chains(structure):
    chainList = []
    for chain in structure.get_chains():
        chainList.append(chain.get_id())
    return chainList

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
parser.add_argument("-b", dest="b", help="the base directory; should be the DIPS partners, before being aligned") # 
parser.add_argument("-s", dest="s", help="only run on a subset of the structures, up to this number")
arguments = parser.parse_args()

if arguments.s:
    maxNum = int(arguments.s)

outDir1 = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65_aligned"
if arguments.s:
    outDir1 += "_" + str(arguments.s)
if not os.path.exists(outDir1):
    os.mkdir(outDir1)
outDir2 = outDir1 + "/structures"
if not os.path.exists(outDir2):
    os.mkdir(outDir2)



### MAIN BODY ###

segFaultList = ['1dd1_A_C','1dd1_B_C','3e47_A_G','3e47_O_U']
csvFile = ['path,split\n']
extracted_filename = 'extracted_structure.pdb'

# read in summary file

dipsSum = []
dipsSumFile = arguments.b + "/splits.csv"
with open(dipsSumFile,"r") as sumFile:
    dipsSum = sumFile.readlines()

nameBaseList = []

header = True
pairsList = []
count = 0
for dEntry in dipsSum:

    nameBase = dEntry.split(",")[0]
    if nameBase in segFaultList:
        continue
    if nameBase in nameBaseList: #prevents dupes
        continue
    else:
        nameBaseList.append(nameBase)

    pdbID = nameBase.split("_")[0]

    if header == True:
        header = False
        continue

    if arguments.s and count >= maxNum:
        break

    # get the ag-like structure

    agPath = arguments.b + "/structures/" + nameBase + "_r_u.pdb"
    agList = glob.glob(agPath)
    if len(agList) > 0:
        relaxedAg = load_structure(agList[0])
    else:
        print("relaxed ag-like partner not found for pdb ID:" + pdbID)
        print("agPath was: " + agPath)
        continue

    # get the ab-like structure
        
    abPath = arguments.b + "/structures/" + nameBase + "_l_u.pdb"
    abList = glob.glob(abPath)
    if len(abList) > 0:
        ibAb = load_structure(abList[0])
    else:
        print("ab-like partner not found for pdb ID:" + pdbID)
        continue

    # continue / skip this PPI set, if necessary by t flag

    agChains = get_chains(relaxedAg)
    agString = ",".join(agChains)
    abChains = get_chains(ibAb)
    abString = ",".join(abChains)

    # make out file names

    agOutName = outDir2 + "/" + nameBase + "_r_u.pdb"
    abOutName = outDir2 + "/" + nameBase + "_l_u.pdb"

    if os.path.exists(agOutName) and os.path.exists(abOutName):
        print("skipping because out files already exist!")
        count += 1 
        csvFile.append(dEntry)
        continue

    # grab the true structure
    
    print("working on: " + nameBase,flush=True)
    structPattern = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/DIPS_pruned_raw/" + pdbID + "*" + abString + "*.pdb"
    structPath = glob.glob(structPattern)[0]
    trueStruct = load_structure(structPath)

    # get the true ag & ab structures

    agStruct = extract_chains(trueStruct,agChains)
    abStruct = extract_chains(trueStruct,abChains)
    # get the true ab structure

    # align the ag struct to the original, if found

    try:
        agToSave = alignTo(relaxedAg,agStruct,aligner)
        write_structure(agToSave,agOutName,io)
    except:
        continue

    # align the ab struct to the original, if found

    try:
        abToSave = alignTo(ibAb,abStruct,aligner)
        write_structure(abToSave,abOutName,io)
    except:
        print("failed at writing ab - removing ag struct " + agOutName)
        os.remove(agOutName)
        continue

    # add a line to the csv file that summarizes them as test files

    count += 1     
    csvFile.append(dEntry)

# save the csv file too!

outFile = outDir1 + "/splits.csv"     
with open(outFile,'w') as of:
    of.writelines(csvFile)

### python3 /dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock_proj/loopy_DIPs_aligner.py -b /dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65 