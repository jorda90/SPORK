from SPORK_GTFEntry import GTFEntry
import argparse
import pickle

parser = argparse.ArgumentParser()

parser.add_argument('-gtf', '--gtf', nargs='+', help='gtf input file')
parser.add_argument('-rn', '--reference_name', help='reference name which will determine the name of output file')

args = parser.parse_args()


gtfs = []
abs_gtf_file_path = args.gtf

for gtf_f in abs_gtf_file_path:
#sys.stdout.write("Reading in GTF file "+abs_gtf_file_path+"\n")
    print gtf_f
    with open(gtf_f, "r") as gtf_file:
        for gtf_line in gtf_file.readlines():

            if gtf_line[0] == '#' or gtf_line[0] == '@':
                continue
            gtf = GTFEntry(gtf_line)
            if gtf.feature in ["exon"]:
                gtfs.append(gtf)

chrom_gtfs = {}
for gtf in gtfs:
    if gtf.chromosome not in chrom_gtfs:
        chrom_gtfs[gtf.chromosome] = [gtf]
    else:
        chrom_gtfs[gtf.chromosome].append(gtf)

# Pre sort the gtfs into a donor and acceptor oriented list by chromosome and strand

gtfs_info = {}
chrom_gtfs_don = {}
chrom_gtfs_acc = {}
chrom_don_libs = {}
chrom_acc_libs = {}
for chrom in chrom_gtfs:
    chrom_gtfs_don[chrom+"+"] = []
    chrom_gtfs_don[chrom+"-"] = []
    chrom_don_libs[chrom+"+"] = []
    chrom_don_libs[chrom+"-"] = []
    for don_gtf in sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.donor):
        chrom_gtfs_don[chrom+don_gtf.strand].append(don_gtf)
        chrom_don_libs[chrom+don_gtf.strand].append(don_gtf.donor)

    chrom_gtfs_acc[chrom+"+"] = []
    chrom_gtfs_acc[chrom+"-"] = []
    chrom_acc_libs[chrom+"+"] = []
    chrom_acc_libs[chrom+"-"] = []
    for acc_gtf in sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.acceptor):
        chrom_gtfs_acc[chrom+acc_gtf.strand].append(acc_gtf)
        chrom_acc_libs[chrom+acc_gtf.strand].append(acc_gtf.acceptor)

#with open("chrom_gtfs_don_keys.txt", "w") as f:
#    f.write('GTF Finder don chroms: '+str(chrom_gtfs_don.keys())+'\n')

#with open("chrom_gtfs_acc_keys.txt", "w") as f:
#    f.write('GTF Finder acc chroms: '+str(chrom_gtfs_acc.keys())+'\n')

chrom_acc_libs_duplicates = {}
chrom_gtfs_acc_duplicates = {}
for chrom in chrom_acc_libs:
    print chrom
    if len(chrom_acc_libs[chrom]) == 0:
        chrom_acc_libs_duplicates[chrom] = []
        chrom_gtfs_acc_duplicates[chrom] = []
        continue
    chrom_acc_libs_duplicates[chrom] = [chrom_acc_libs[chrom][0]]
    chrom_gtfs_acc_duplicates[chrom] = [chrom_gtfs_acc[chrom][0]]
    prev_acc = chrom_acc_libs[chrom][0]
    curr_index = 0
    for index in range(1, len(chrom_acc_libs[chrom])):
        if prev_acc == chrom_acc_libs[chrom][index]:
            if chrom_gtfs_acc[chrom][index].gene_name != chrom_gtfs_acc_duplicates[chrom][curr_index].gene_name:
                chrom_gtfs_acc_duplicates[chrom][curr_index].synonyms.add(chrom_gtfs_acc[chrom][index].gene_name)
            continue
        prev_acc = chrom_acc_libs[chrom][index]
        curr_index += 1
        chrom_acc_libs_duplicates[chrom].append(prev_acc)
        chrom_gtfs_acc_duplicates[chrom].append(chrom_gtfs_acc[chrom][index])

chrom_don_libs_duplicates = {}
chrom_gtfs_don_duplicates = {}
for chrom in chrom_don_libs:
    print chrom
    if len(chrom_don_libs[chrom]) == 0:
        chrom_don_libs_duplicates[chrom] = []
        chrom_gtfs_don_duplicates[chrom] = []
        continue
    chrom_don_libs_duplicates[chrom] = [chrom_don_libs[chrom][0]]
    chrom_gtfs_don_duplicates[chrom] = [chrom_gtfs_don[chrom][0]]
    prev_don = chrom_don_libs[chrom][0]
    curr_index = 0
    for index in range(1, len(chrom_don_libs[chrom])):
        if prev_don == chrom_don_libs[chrom][index]:
            if chrom_gtfs_don[chrom][index].gene_name != chrom_gtfs_don_duplicates[chrom][curr_index].gene_name:
                chrom_gtfs_don_duplicates[chrom][curr_index].synonyms.add(chrom_gtfs_don[chrom][index].gene_name)
            continue
        prev_don = chrom_don_libs[chrom][index]
        curr_index += 1
        chrom_don_libs_duplicates[chrom].append(prev_don)
        chrom_gtfs_don_duplicates[chrom].append(chrom_gtfs_don[chrom][index])

gtfs_info["chrom_gtfs_don"] = chrom_gtfs_don_duplicates
gtfs_info["chrom_gtfs_acc"] = chrom_gtfs_acc_duplicates
gtfs_info["chrom_don_libs"] = chrom_don_libs_duplicates
gtfs_info["chrom_acc_libs"] = chrom_acc_libs_duplicates


pickle.dump(gtfs_info, open(str(args.reference_name) + "_gtfs_info_duplicates_removed","wb"), -1)
