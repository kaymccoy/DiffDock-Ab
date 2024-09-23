import glob
import shutil

# summary file 

abAgSumFile = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_SAbDab_01-17-22_all_c180_c290_c380_r3pt5_d2022-03-16_handParsed/lightDb.txt'
pdbList = []
outFile = ["path,split\n"]

with open(abAgSumFile,"r") as sumFile:
    header = True
    for line in sumFile:
        if header:
            header = False
            continue
        pdbList.append(line.split(",")[0])

print(pdbList)

abbDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/abbOutput/"

afAntigens = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldAntigens/ag_trimmed/"

for p in pdbList:
    abFile = glob.glob(abbDir + p + "*.pdb")[0]
    agFile = glob.glob(afAntigens + p + "*.pdb")[0]

    outDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/ADAAM_test_set/structures/" 
    ab_lig_out = outDir + p + "_l_u.pdb"
    ag_rec_out = outDir + p + "_r_u.pdb"
    
    shutil.copyfile(abFile, ab_lig_out)
    shutil.copyfile(agFile, ag_rec_out)

    outFile.append(p + ",test\n")
    
with open("/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/ADAAM_test_set.csv", "w") as of:
    of.writelines(outFile)