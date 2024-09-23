import glob


inStr = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock_proj/AADAM_vanilla_structs_fix_1/*_ligand-gt.pdb"
pairStr = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock_proj/AADAM_vanilla_structs_fix_1/*_receptor.pdb"
outDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock_proj/AADAM_vanilla_structs_fixed/"

iList = glob.glob(inStr)
pairList = glob.glob(pairStr)
for i in iList:
    nameBase = i.split("/")[-1].split("_ligand-gt.pdb")[0]
    for ii in pairList:
        nameBase2 = ii.split("/")[-1].split("_receptor.pdb")[0]
        if nameBase == nameBase2:
            with open(i,"r") as f1:
                outFile = f1.readlines()
            outFile.pop() #removes END line
            with open(ii,"r") as f2:
                toAppend = f2.readlines()
            outFile += toAppend

            with open(outDir + nameBase + ".pdb","w") as of:
                of.writelines(outFile)


