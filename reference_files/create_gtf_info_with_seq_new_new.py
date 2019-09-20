import copy
from Bio import SeqIO
from Bio import SeqFeature
from BCBio import GFF
import pickle
import pandas as pd
from SPORK_GTFEntry import GTFEntry
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-gtf', '--gtf_file', help='annotation file')
parser.add_argument('-fa', '--fasta', help='fasta file')

args = parser.parse_args()


f_handle = open(args.fasta)
f_dict = SeqIO.to_dict(SeqIO.parse(f_handle, "fasta"))
f_handle.close()

chrom_length = {}
for chrom in f_dict:
    chrom_length[chrom] = len(f_dict[chrom].seq)

limit_info = dict(
        gff_id = f_dict.keys(),
        gff_type = ["exon"])


gtfs = []
a_handle = open(args.gtf_file, "rU")
for rec in GFF.parse(a_handle, base_dict=f_dict, limit_info=limit_info):
    print(rec.id)
    for ftr in rec.features:
        if len(ftr.sub_features) > 0:
            for sftr in ftr.sub_features:
                chromosome = rec.id
                if "source" in sftr.qualifiers:
                    source = sftr.qualifiers["source"][0]
                else:
                    source = "-"
                feature = "exon"
                start = int(sftr.location.start)
                stop = int(sftr.location.end)
                score = "."
                strand = "+" if int(sftr.location.strand) == 1 else "-"
                frame = "."
                if "gene_name" in sftr.qualifiers:
                    gene_name = sftr.qualifiers["gene_name"][0]
                elif "gene_id" in sftr.qualifiers:
                    gene_name = sftr.qualifiers["gene_id"][0]
                #new_location = SeqFeature.FeatureLocation(start = start - 1000001, end = stop + 1000001, strand = sftr.location.strand)
                #sftr.location = new_location
                #seq = str(sftr.extract(rec.seq))
                seq = ""
                gtf = GTFEntry(chromosome, source, feature, start, stop, score, strand, frame, gene_name, seq)
                gtfs.append(gtf)
        else:
            chromosome = rec.id
            if "source" in ftr.qualifiers:
                source = ftr.qualifiers["source"][0]
            else:
                source = "-"
            feature = "exon"
            start = int(ftr.location.start)
            stop = int(ftr.location.end)
            score = "."
            strand = "+" if int(ftr.location.strand) == 1 else "-"
            frame = "."
            if "gene_name" in ftr.qualifiers:
                gene_name = ftr.qualifiers["gene_name"][0]
            elif "gene_id" in ftr.qualifiers:
                gene_name = ftr.qualifiers["gene_id"][0]
            #new_location = SeqFeature.FeatureLocation(start = start - 1000001, end = stop + 1000001, strand = ftr.location.strand)
            #ftr.location = new_location
            #seq = str(ftr.extract(rec.seq))
            seq = ""
            gtf = GTFEntry(chromosome, source, feature, start, stop, score, strand, frame, gene_name, seq)
            gtfs.append(gtf)
            
chrom_gtfs = {}
for gtf in gtfs:
    if gtf.chromosome not in chrom_gtfs:
        chrom_gtfs[gtf.chromosome] = [gtf]
    else:
        chrom_gtfs[gtf.chromosome].append(gtf)
        
chrom_seq_don = {}
chrom_seq_acc = {}
for chrom in chrom_gtfs:
    chrom_seq_don[chrom + "+"] = {}
    chrom_seq_don[chrom + "-"] = {}
    
    positive_exons = [gtf for gtf in chrom_gtfs[chrom] if gtf.strand=="+"]
    previous_don = 0
    sorted_don_gtf = sorted(positive_exons,key=lambda gtf: gtf.donor)
    for index, don_gtf in enumerate(sorted_don_gtf):
        
        current_index = index
        next_donor = sorted_don_gtf[current_index]
        
        while next_donor.donor == don_gtf.donor and current_index < len(sorted_don_gtf) - 1:
            current_index += 1
            next_donor = sorted_don_gtf[current_index]
        
        
        if don_gtf.donor in chrom_seq_don[chrom+don_gtf.strand]:
            pass
            #if don_gtf.gene_name in chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] and don_gtf.seq not in chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name]:
                #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name].append(don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]])
            #else:
                #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name] = [don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]]]
        else:
            chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] = {"previous_don" : min(1000000, (don_gtf.donor - previous_don)/2)}
            if next_donor.donor == don_gtf.donor:
                chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"] = min(1000000, chrom_length[chrom] - don_gtf.donor)
            else:
                chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"] = min(1000000, (next_donor.donor - don_gtf.donor)/2)
            
            #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] = {don_gtf.gene_name : [don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]]]}
        
        if don_gtf.donor != previous_don:
            previous_don = don_gtf.donor
    
    negative_exons = [gtf for gtf in chrom_gtfs[chrom] if gtf.strand=="-"]
    previous_don = 0
    sorted_don_gtf = sorted(negative_exons,key=lambda gtf: gtf.donor)
    for index, don_gtf in enumerate(sorted_don_gtf):
        
        current_index = index
        next_donor = sorted_don_gtf[current_index]
        
        while next_donor.donor == don_gtf.donor and current_index < len(sorted_don_gtf) - 1:
            current_index += 1
            next_donor = sorted_don_gtf[current_index]
        
        if don_gtf.donor in chrom_seq_don[chrom+don_gtf.strand]:
            pass
            #if don_gtf.gene_name in chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] and don_gtf.seq not in chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name]:
                #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name].append(don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]])
            #else:
                #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor][don_gtf.gene_name] = [don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]]]
        else:
            chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] = {"previous_don" : min(1000000, (don_gtf.donor - previous_don)/2)}
            if next_donor.donor == don_gtf.donor:
                chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"] = min(1000000, chrom_length[chrom] - don_gtf.donor)
            else:
                chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"] = min(1000000, (next_donor.donor - don_gtf.donor)/2)
            
            #chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor] = {don_gtf.gene_name : [don_gtf.seq[1000001 - chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["previous_don"] : -1000001 + chrom_seq_don[chrom+don_gtf.strand][don_gtf.donor]["next_don"]]]}
        
        if don_gtf.donor != previous_don:
            previous_don = don_gtf.donor
    
    chrom_seq_acc[chrom + "+"] = {}
    chrom_seq_acc[chrom + "-"] = {}
    

    previous_acc = 0
    sorted_acc_gtf = sorted(positive_exons,key=lambda gtf: gtf.acceptor)
    for index, acc_gtf in enumerate(sorted_acc_gtf):
                                                                                         
        current_index = index
        next_acceptor = sorted_acc_gtf[current_index]
        
        while next_acceptor.acceptor == acc_gtf.acceptor and current_index < len(sorted_acc_gtf) - 1:
            current_index += 1
            next_acceptor = sorted_acc_gtf[current_index]
                                                                                         
        if acc_gtf.acceptor in chrom_seq_acc[chrom+acc_gtf.strand]:
            pass
            #if acc_gtf.gene_name in chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] and acc_gtf.seq not in chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name]:
                #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name].append(acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]])
            #else:
                #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name] = [acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]]]
        else:
            chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] = {"previous_acc" : min(1000000, (acc_gtf.acceptor - previous_acc)/2)}
            if next_acceptor.acceptor == acc_gtf.acceptor:
                chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"] = min(1000000, chrom_length[chrom] - acc_gtf.acceptor)
            else:
                chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"] = min(1000000, (next_acceptor.acceptor - acc_gtf.acceptor)/2)
            
            #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] = {acc_gtf.gene_name : [acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]]]}
        
        if previous_acc != acc_gtf.acceptor:
            previous_acc = acc_gtf.acceptor
            
    previous_acc = 0
    sorted_acc_gtf = sorted(negative_exons,key=lambda gtf: gtf.acceptor)
    for index, acc_gtf in enumerate(sorted_acc_gtf):
                                                                                         
        current_index = index
        next_acceptor = sorted_acc_gtf[current_index]
        
        while next_acceptor.acceptor == acc_gtf.acceptor and current_index < len(sorted_acc_gtf) - 1:
            current_index += 1
            next_acceptor = sorted_acc_gtf[current_index]
                                                                                         
        if acc_gtf.acceptor in chrom_seq_acc[chrom+acc_gtf.strand]:
            pass
            #if acc_gtf.gene_name in chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] and acc_gtf.seq not in chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name]:
                #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name].append(acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]])
            #else:
                #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor][acc_gtf.gene_name] = [acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]]]
        else:
            chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] = {"previous_acc" : min(1000000, (acc_gtf.acceptor - previous_acc)/2)}
            if next_acceptor.acceptor == acc_gtf.acceptor:
                chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"] = min(1000000, chrom_length[chrom] - acc_gtf.acceptor)
            else:
                chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"] = min(1000000, (next_acceptor.acceptor - acc_gtf.acceptor)/2)
            
            #chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor] = {acc_gtf.gene_name : [acc_gtf.seq[1000001 - chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["previous_acc"] : -1000001 + chrom_seq_acc[chrom+acc_gtf.strand][acc_gtf.acceptor]["next_acc"]]]}
        
        if previous_acc != acc_gtf.acceptor:
            previous_acc = acc_gtf.acceptor
            
a_handle = open(args.gtf_file, "rU")
for rec in GFF.parse(a_handle, base_dict=f_dict, limit_info=limit_info):
    print(rec.id)
    for ftr in rec.features:
        if len(ftr.sub_features) > 0:
            for sftr in ftr.sub_features:
                chromosome = rec.id
                strand = "+" if int(sftr.location.strand) == 1 else "-"
                
                donor = int(sftr.location.end) if strand == "+" else int(sftr.location.start)
                acceptor = int(sftr.location.start) if strand == "+" else int(sftr.location.end)
                
                if "gene_name" in sftr.qualifiers:
                    gene_name = sftr.qualifiers["gene_name"][0]
                elif "gene_id" in sftr.qualifiers:
                    gene_name = sftr.qualifiers["gene_id"][0]
                #new_location = SeqFeature.FeatureLocation(start = start - 1000001, end = stop + 1000001, strand = sftr.location.strand)
                #sftr.location = new_location
                #seq = str(sftr.extract(rec.seq))
                #seq = ""
                #gtf = GTFEntry(chromosome, source, feature, start, stop, score, strand, frame, gene_name, seq)
                #gtfs.append(gtf)
                start = int(sftr.location.start)
                stop = int(sftr.location.end)
                new_location_don = SeqFeature.FeatureLocation(start = start - chrom_seq_don[chromosome + strand][donor]["previous_don"], end = stop + chrom_seq_don[chromosome + strand][donor]["next_don"], strand = sftr.location.strand)
                sftr.location = new_location_don
                if gene_name in chrom_seq_don[chromosome + strand][donor]:
                    chrom_seq_don[chromosome + strand][donor][gene_name].append(str(sftr.extract(rec.seq)))
                else:
                    chrom_seq_don[chromosome + strand][donor][gene_name] = [str(sftr.extract(rec.seq))]
                
                
                new_location_acc = SeqFeature.FeatureLocation(start = start - chrom_seq_acc[chromosome + strand][acceptor]["previous_acc"], end = stop + chrom_seq_acc[chromosome + strand][acceptor]["next_acc"], strand = sftr.location.strand)
                sftr.location = new_location_acc
                
                if gene_name in chrom_seq_acc[chromosome + strand][acceptor]:
                    chrom_seq_acc[chromosome + strand][acceptor][gene_name].append(str(sftr.extract(rec.seq)))
                else:
                    chrom_seq_acc[chromosome + strand][acceptor][gene_name] = [str(sftr.extract(rec.seq))]
                    
        else:
            chromosome = rec.id
            strand = "+" if int(ftr.location.strand) == 1 else "-"
            
            donor = int(ftr.location.end) if strand == "+" else int(ftr.location.start)
            acceptor = int(ftr.location.start) if strand == "+" else int(ftr.location.end)
            
            if "gene_name" in ftr.qualifiers:
                gene_name = ftr.qualifiers["gene_name"][0]
            elif "gene_id" in ftr.qualifiers:
                gene_name = ftr.qualifiers["gene_id"][0]
            #new_location = SeqFeature.FeatureLocation(start = start - 1000001, end = stop + 1000001, strand = ftr.location.strand)
            #ftr.location = new_location
            #seq = str(ftr.extract(rec.seq))
            #seq = ""
            #gtf = GTFEntry(chromosome, source, feature, start, stop, score, strand, frame, gene_name, seq)
            #gtfs.append(gtf)
            
            start = int(ftr.location.start)
            stop = int(ftr.location.end)
            
            new_location_don = SeqFeature.FeatureLocation(start = start - chrom_seq_don[chromosome + strand][donor]["previous_don"], end = stop + chrom_seq_don[chromosome + strand][donor]["next_don"], strand = ftr.location.strand)
            ftr.location = new_location_don
            if gene_name in chrom_seq_don[chromosome + strand][donor]:
                chrom_seq_don[chromosome + strand][donor][gene_name].append(str(ftr.extract(rec.seq)))
            else:
                chrom_seq_don[chromosome + strand][donor][gene_name] = [str(ftr.extract(rec.seq))]
            
            new_location_acc = SeqFeature.FeatureLocation(start = start - chrom_seq_acc[chromosome + strand][acceptor]["previous_acc"], end = stop + chrom_seq_acc[chromosome + strand][acceptor]["next_acc"], strand = ftr.location.strand)
            ftr.location = new_location_acc
                
            if gene_name in chrom_seq_acc[chromosome + strand][acceptor]:
                chrom_seq_acc[chromosome + strand][acceptor][gene_name].append(str(ftr.extract(rec.seq)))
            else:
                chrom_seq_acc[chromosome + strand][acceptor][gene_name] = [str(ftr.extract(rec.seq))]
            
        
gtfs_info_seq = {}
gtfs_info_seq["chrom_seq_don"] = chrom_seq_don
gtfs_info_seq["chrom_seq_acc"] = chrom_seq_acc

pickle.dump(gtfs_info_seq, open("GRCh38" + "_gtfs_info_seq_10_new_new_new","wb"), 0)

