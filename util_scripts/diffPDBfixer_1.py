import glob
import prody

prody.confProDy(verbosity='critical')

def fix_pdb(infile_name, outfile_name):
    prot = prody.parsePDB(infile_name)
    prody.writePDB(outfile_name, prot)

inDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock/visualization_top/AADAM_test_1/"
outDir = "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbDiffDock_proj/AADAM_vanilla_structs_fixed/"

inList = glob.glob(inDir + "*.pdb")

for i in inList:
    baseName = i.split("/")[-1]
    outName = outDir + baseName
    fix_pdb(i,outName)
