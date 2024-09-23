import glob
import os

list1 = glob.glob('/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65_aligned/structures/*_l_u.pdb')

list2 = glob.glob('/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65_aligned/structures/*_r_u.pdb')

for i in list1:
    iBase = i.split("/")[-1].split("_l_u.pdb")[0]
    rMatch = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65_aligned/structures/' + iBase + "_r_u.pdb"
    if not os.path.exists(rMatch):
        print('missing partner: ')
        print(rMatch)

for i in list2:
    iBase = i.split("/")[-1].split("_r_u.pdb")[0]
    rMatch = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/DiffDock-PP/datasets/loopyDIPS-0.65_aligned/structures/' + iBase + "_l_u.pdb"
    if not os.path.exists(rMatch):
        print('missing partner: ')
        print(rMatch)