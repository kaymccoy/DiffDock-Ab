csvListFile = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/AADaM_notstrict/splits_allLD_val.csv"
baseDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/AADaM_notstrict/structures/"
outCSV = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/AADaM_notstrict/splits_allLD_val_filtered.csv"

with open(csvListFile, "r") as f:
    csvList = f.readlines()

i = 1
while i < len(csvList):
    pdbBase = csvList[i].split(",")[0]
    pdbL = baseDir + pdbBase + "_l_b.pdb"
    pdbR = baseDir + pdbBase + "_r_b.pdb"

    caCount = 0
    with open(pdbL, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("ATOM"):
            if line[13:15] == "CA":
                caCount += 1

    with open(pdbR, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("ATOM"):
            if line[13:15] == "CA":
                caCount += 1

    if caCount > 1000:
        print("Removing " + pdbBase + " from list")
        csvList.pop(i)
    
    i+=1

with open(outCSV, "w") as f:
    for line in csvList:
        f.write(line)