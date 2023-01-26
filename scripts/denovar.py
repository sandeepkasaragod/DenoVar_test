import read_file
import tryptic_peptide
import os
import re
import sys
from itertools import islice
import pandas as pd
from os.path import join as join

#prot_db_cleave = "GCF_000001405.39_GRCh38.p13_protein_refseq109_formated_validated_added_missing.fasta" # DB used for cleaving

nr_pep_list = {} # store tyrptic_peptide and corresponding protein accession
pep_len_list_prot_acc = {} # store tyrptic_peptide length as key and its protein accession as value
pep_len_list_peptide = {} # store tyrptic_peptide length as key and its peptide as value
variant_pep_details = {} # storing the variant peptide with protein accession

dicts_dbSNP_variants = {} # storing dbsnp raw variants information
dict_seq = {} # NP accession and protein sequence
gene_to_np = {} # NP accession and gene symbol

dicts_dbSNP_variant_protein = {}

def distance(a, b):
    f = (lambda x, y: 0 if x == y else 1)
    return sum(map(f, a, b))
    #sum(map(lambda (x, y): 0 if x == y else 1, zip(a, b)))

def create_single_line_fasta(db_dir, infile, refseq):
	strs = ""
	if not os.path.isfile(join(db_dir, infile)):
		print ("creating a single line fasta file")
		write_file = open(join(db_dir, infile), 'w')
		rds = read_file.read_fasta(join(db_dir, refseq))
		for rows in rds:
			strs = strs + rows[1].rstrip()
		write_file.write(strs.rstrip())
		write_file.close()
		for i in open(join(db_dir, infile)):
			strs = strs + i.rstrip()
		return strs
	else:
		for i in open(join(db_dir, 'singleline.fasta')):
			strs = i.rstrip()
	return strs
					
def create_nr_pep_list(enzyme, clevage, min_len, max_len, refseq):
    if enzyme =="Trypsin":
        rds = read_file.read_fasta(refseq)
        for rows in rds:
            for iter_clevage in range(clevage + 1):
                a = tryptic_peptide.trypsin(rows[1].rstrip(), iter_clevage, 7, 35)
                for i in a:
                    if i not in nr_pep_list:
                            nr_pep_list[i] = [rows[0].split('|')[3]]
                    else:
                            nr_pep_list[i].append(rows[0].split('|')[3])
                            
    elif enzyme=="LysC":
        rds = read_file.read_fasta(refseq)
        for rows in rds:
            for iter_clevage in range(clevage + 1):
                a = tryptic_peptide.tryptic_peptide_lysc(rows[1].rstrip(), iter_clevage, 7, 35)
                for i in a:
                    if i not in nr_pep_list:
                            nr_pep_list[i] = [rows[0].split('|')[3]]
                    else:
                            nr_pep_list[i].append(rows[0].split('|')[3])
                            
    elif enzyme=="Chymotrypsin":
        rds = read_file.read_fasta(refseq)
        for rows in rds:
            for iter_clevage in range(clevage + 1):
                a = tryptic_peptide.tryptic_peptide_lysc(rows[1].rstrip(), iter_clevage, 7, 35)
                for i in a:
                    if i not in nr_pep_list:
                            nr_pep_list[i] = [rows[0].split('|')[3]]
                    else:
                            nr_pep_list[i].append(rows[0].split('|')[3])

enzyme_cleave_peps = {} # for storing the peps, when protein accession is already known
def create_cleavege_peptide(prot_seq, enzyme, clevage, min_len, max_len):
    enzyme_cleave_peps = {}
    if enzyme == "Trypsin":
        for iter_clevage in range(int(clevage) + 1):
            a = tryptic_peptide.trypsin(prot_seq, iter_clevage, min_len, max_len)
            for i in a:
                if i not in enzyme_cleave_peps:
                    enzyme_cleave_peps[i] = [1]
                else:
                    enzyme_cleave_peps[i].append(1)

    elif enzyme =="LysC":
        pass
    elif enzyme == "Chymptrypsin":
        pass
    
    return enzyme_cleave_peps
# Creating peptide_lenght = [peptides] and peptide_lenght = [its corresponding NP accessions]
def create_pep_len_list():
    for k, v in nr_pep_list.items():
        if len(k) not in pep_len_list_peptide:
                pep_len_list_peptide[len(k)] = [k]
                pep_len_list_prot_acc[len(k)] = [v]
        else:
                pep_len_list_peptide[len(k)].append(k)
                pep_len_list_prot_acc[len(k)].append(v)


#Get the protein accession number for variant peptides
def get_variant_pep_prot_info(infile):
    for variant_pep in open(infile):
        if variant_pep.rstrip()[-1] == "K" or variant_pep.rstrip()[-1] == "R":
            if len(variant_pep.rstrip()) in pep_len_list_peptide:
                for iter_peptide in pep_len_list_peptide[len(variant_pep.rstrip())]:
                    dists = distance(variant_pep.rstrip(), iter_peptide)
                    if dists == 1:
                        peptide_idx = pep_len_list_peptide[len(variant_pep.rstrip())].index(iter_peptide) # peptided index to trace backs its NP_ID
                        get_np_acc = pep_len_list_prot_acc[len(variant_pep.rstrip())][peptide_idx] # NP accession for variant peptide
                        variant_pep_details[variant_pep.rstrip()] = get_np_acc #Variant peptide = [NP accession lists]

def get_variant_pep_prot_info_update(variant_pep):
    if variant_pep.rstrip()[-1] == "K" or variant_pep.rstrip()[-1] == "R":
        if len(variant_pep.rstrip()) in pep_len_list_peptide:
            for iter_peptide in pep_len_list_peptide[len(variant_pep.rstrip())]:
                dists = distance(variant_pep.rstrip(), iter_peptide)
                if dists == 1:
                    peptide_idx = pep_len_list_peptide[len(variant_pep.rstrip())].index(iter_peptide) # peptided index to trace backs its NP_ID
                    get_np_acc = pep_len_list_prot_acc[len(variant_pep.rstrip())][peptide_idx] # NP accession for variant peptide
                    variant_pep_details[variant_pep.rstrip()] = get_np_acc #Variant peptide = [NP accession lists]
                    #print (variant_pep, get_np_acc, variant_pep_details)
    return (variant_pep_details)

def normal_seq(db_dir, infile):
        rds = read_file.read_fasta(join(db_dir, infile))
        for rows in rds:
            dict_seq[rows[0].split(' ')[0]] = rows[1].rstrip()

def read_dbSNPDB(infile):
    for i in open(infile):
        split_i = i.split('\t')
        if split_i[0].rstrip() + "_" + split_i[1].rstrip() not in dicts_dbSNP_variants:
            dicts_dbSNP_variants[split_i[0].rstrip() + "_" + split_i[1].rstrip()] = i.rstrip()

#read_dbSNPDB()
#print(len(dicts_dbSNP_variants))
#print ("Reading db finished")

def create_dbSNPDB(accession): #protein accession
    normal_sequence = dict_seq[accession]
    if accession in dicts_dbSNP_variants:
        try:
            for i in dicts_dbSNP_variants[accession]:
                split_i = i.split('\t')
                get_mut_col = split_i[2] #A143T
                mut_loc = re.findall('\d+', get_mut_col)
                ref_aa = get_mut_col[0]
                alt_aa = get_mut_col[-1]
                try:
                    if normal_sequence[int(mut_loc[0]) -1 ] == ref_aa:
                        if ref_aa != alt_aa:
                            mut_seq = normal_sequence[:int(mut_loc[0]) - 1] + alt_aa + normal_sequence[int(mut_loc[0]):]
                            dicts_dbSNP_variant_protein[mut_seq] = i.rstrip()
                except:
                    pass
        except:
            pass

def read_input_header(infile):
    with open(infile) as fls:
        for i in islice(fls, 0, 1):
            return i

def get_variant_aa_pos(variant_pep, normal_pep, acc): #For getting the genome locationl
    normal_seq = dict_seq[acc]
    variant_pos = [(x, variant_pep[x]) for x in range(len(variant_pep)) if variant_pep[x] != normal_pep[x]]
    #print (variant_pos[0][0], variant_pos[0][1])
    position = int(normal_seq.index(normal_pep)) + int(variant_pos[0][0])
    #print (normal_seq, variant_pos[0][1], position, normal_seq[position])
    return normal_seq[position], position + 1, variant_pos[0][1]  #print (position, normal_seq[position])

#normal_seq_single_line = ""
#for i in open('Single_line.fasta'):
#    normal_seq_single_line = i.rstrip()

def get_header(infile):
	with open(infile) as f:
		for i in islice(f, 0, 1):
			return i.rstrip()
cnt = 1
if __name__=="__main__":
    if len(sys.argv) == 13:
        input_file = sys.argv[1]
        ## db_list = sys.argv[2]
        enzyme = sys.argv[2]
        aa_min_len = sys.argv[3]
        aa_max_len = sys.argv[4]
        max_num_clevage = sys.argv[5]
        is_prt_acc = sys.argv[6]
        db_path = sys.argv[7]
        db_dir = sys.argv[8]
        refseq_db = sys.argv[9]
        temp_dir = sys.argv[10] # blast processed files
        threads = sys.argv[11]
        output_dir = sys.argv[12]
        ##is_prt_acc = sys.argv[2]
        normal_seq_single_line = create_single_line_fasta(db_dir, "singleline.fasta", refseq_db)
        normal_seq(db_dir, refseq_db)
        read_dbSNPDB(db_path)
        output_file = open(join(output_dir, input_file.split('.')[0] + '.txt'), 'w')
        if is_prt_acc == "Yes":
            hdr = get_header(input_file)
            df = pd.read_csv(input_file, sep='\t')
            try:
                v_pep = df.columns.get_loc('peptides')
                prot_acc = df.columns.get_loc('protein_accession')
            except:
                  print ("column missing, check for the column names peptides and protein_accession")
                  exit(1)
            output_file.write(hdr + '\t' + "Annotation_details" + '\n')
            with open(input_file) as f:
                for line in islice(f, 1, None):
                    split_line = line.split('\t')
                    #normal_aa = split_line[1].rstrip()
                    variant_peptide = split_line[v_pep]
                    combine = ""
                    if variant_peptide not in normal_seq_single_line:
                        cnt+=1
                        for iter_prt_acc in split_line[prot_acc].split(':'):
                            #Calling the normal protein
                            enzyme_cleave_peps = {} #making the dictionary empty for every protein sequence
                            if iter_prt_acc.rstrip() in dict_seq: # Check in Normal sequence dictinary
                                normal_protein_seq = dict_seq[iter_prt_acc.rstrip()]
                                enzyme_cleave_db = create_cleavege_peptide(normal_protein_seq, enzyme, max_num_clevage, aa_min_len, aa_max_len)
                                for iter_pep, iter_prot_acc in enzyme_cleave_db.items():
                                    if len(iter_pep) == len(variant_peptide) and distance(iter_pep, variant_peptide) == 1:
                                        variant_info = get_variant_aa_pos(variant_peptide, iter_pep, iter_prt_acc.rstrip()) #(N, 99, H)
                                        if (iter_prt_acc.rstrip() + "_" + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])) in dicts_dbSNP_variants:
                                            variant_infos = dicts_dbSNP_variants[iter_prt_acc.rstrip() + "_" + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])]
                                            split_var_info = variant_infos.split('\t')
                                            if len(combine) >1:
                                                combine = combine + ";" + iter_prt_acc.rstrip() + ":= " +'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])+ "||Clinvar_info=" + split_var_info[6] + "||COSMIC_id=" + split_var_info[9] + "|Mutation_id=" + split_var_info[10] + "|Mutation_status=" + split_var_info[11] + "|Primary_site=" + split_var_info[12] + "|Pubmed=" + split_var_info[13] + "||" + "ICGC_id=" + split_var_info[15] + "|verification_status="+ split_var_info[18] + "|Biological_verification_status="+ split_var_info[19].rstrip()
                                            else:
                                                combine = combine + iter_prt_acc.rstrip() + ":= "+ 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])+ "||Clinvar_info=" + split_var_info[6] + "||COSMIC_id=" + split_var_info[9] + "|Mutation_id=" + split_var_info[10] + "|Mutation_status=" + split_var_info[11] + "|Primary_site=" + split_var_info[12] + "|Pubmed=" + split_var_info[13] + "||" + "ICGC_id=" + split_var_info[15] + "|verification_status="+ split_var_info[18] + "|Biological_verification_status="+ split_var_info[19].rstrip()
                                        else:
                                            if len(combine) > 1:
                                                combine = combine + ";" + iter_prt_acc.rstrip() + ":= " + 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2]) + '||' + 'NOT_AVAILABLE'
                                            else:
                                                combine = combine + iter_prt_acc.rstrip() + ":= " + 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2]) +'||' + 'NOT_AVAILABLE'

                            else:
                                #print (line.rstrip() + '\t' + "-" + '\t' + "-") 
                                pass
                        if len(combine) >1:
                            output_file.write(line.rstrip() + '\t' +  combine + '\n')
                        else:
                            output_file.write(line.rstrip() + '\t' + 'Unknown peptide/Peptide without clevage site' + '\n')
                    else:
                        output_file.write(line.rstrip() + '\t' + "Normal peptide" + '\n')
                	#print (line.rstrip() + '\t' + combine)
            output_file.close()
        else:
            dicts = {}
            df = pd.read_csv(input_file, sep='\t')
            df['id'] = df.index
            df.to_csv(join(temp_dir, input_file.split('.')[0] + '_id_added.txt'), sep='\t', index=False)
            df = pd.read_csv(join(temp_dir,input_file.split('.')[0] + '_id_added.txt'), sep='\t')
            header_cols = df.columns
            v_pep = df.columns.get_loc('peptides')
            idss = df.columns.get_loc('id')
            
            write_file = open(join(temp_dir, 'for_blastp.fasta'), 'w')
            for ids, pep in zip(df['id'], df['peptides']):
                write_file.write('>' + str(ids) + '\n' + pep.strip() + '\n')
            # running blastp
            write_file.close()
            if not os.path.isfile(join(db_dir,'refseq.fasta.phr')):
                os.system('makeblastdb -in ' + join(db_dir, 'refseq.fasta') + ' -dbtype prot -out ' + join(db_dir, 'refseq.fasta'))
                
            os.system('blastp -db ' + join(db_dir, 'refseq.fasta') +\
             ' -query ' + join(temp_dir, 'for_blastp.fasta') + ' -out ' + join(temp_dir, input_file.split('.')[0] + '_prot_acc.txt') + \
            ' -outfmt "7 qseqid sseqid pident lenght mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -num_threads ' + threads)
            #sleep(2)
            
            for i in open(join(temp_dir, input_file.split('.')[0] + '_prot_acc.txt')):
                if '#' not in i.rstrip():
                    split_i = i.split('\t')
                    if split_i[0] not in dicts:
                        dicts[split_i[0]] = [split_i[1]]
                    else:
                        dicts[split_i[0]].append(split_i[1]) #List of protein accessions

            output_file.write('\t'.join(df.columns) + '\t' + "Annontation_details" + '\n')
            for line in open(join(temp_dir, input_file.split('.')[0] + '_id_added.txt')):
                split_line = line.split('\t')
                variant_peptide = split_line[v_pep]
                ids = split_line[idss].rstrip()
                combine=""
                if variant_peptide not in normal_seq_single_line:
                    cnt+=1
                    enzyme_cleave_peps = {}
                    if ids in dicts:
                        for iter_acc in dicts[ids]:
                            if iter_acc in dict_seq:
                                normal_protein_seq = dict_seq[iter_acc]
                                enzyme_cleave_db = create_cleavege_peptide(normal_protein_seq, enzyme, max_num_clevage, aa_min_len, aa_max_len)
                                for iter_pep, iter_prot_acc in enzyme_cleave_db.items():
                                    if len(iter_pep) == len(variant_peptide) and distance(iter_pep, variant_peptide) == 1:
                                        variant_info = get_variant_aa_pos(variant_peptide, iter_pep, iter_acc)
                                        if iter_acc + "_" + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2]) in dicts_dbSNP_variants:
                                            variant_infos = dicts_dbSNP_variants[iter_acc + "_" + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])]
                                            split_var_info = variant_infos.split('\t')
                                            if len(combine) >1:
                                                combine = combine + ";" + iter_acc + ":= " +'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])+ "||Clinvar_info=" + split_var_info[6] + "||COSMIC_id=" + split_var_info[9] + "|Mutation_id=" + split_var_info[10] + "|Mutation_status=" + split_var_info[11] + "|Primary_site=" + split_var_info[12] + "|Pubmed=" + split_var_info[13] + "||" + "ICGC_id=" + split_var_info[14] + "|verification_status="+ split_var_info[17] + "|Biological_verification_status="+ split_var_info[18].rstrip()
                                            else:
                                                combine = combine + iter_acc + ":= "+ 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2])+ "||Clinvar_info=" + split_var_info[6] + "||COSMIC_id=" + split_var_info[9] + "|Mutation_id=" + split_var_info[10] + "|Mutation_status=" + split_var_info[11] + "|Primary_site=" + split_var_info[12] + "|Pubmed=" + split_var_info[13] + "||" + "ICGC_id=" + split_var_info[14] + "|verification_status="+ split_var_info[17] + "|Biological_verification_status="+ split_var_info[18].rstrip()
                                        else:
                                            if len(combine) > 1:
                                                combine = combine + ";" + iter_acc + ":= " + 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2]) + '||' + 'NOT_AVAILABLE'
                                            else:
                                                combine = combine + iter_acc + ":= " + 'p.' + str(variant_info[0]) + str(variant_info[1]) + str(variant_info[2]) +'||' + 'NOT_AVAILABLE'
                            
                        if len(combine) > 1:
                            output_file.write(line.rstrip() + '\t' +  combine + '\n')
                        else:
                            output_file.write(line.rstrip() + '\t' + 'Unknown peptide/Peptide without clevage site' + '\n')
                else:
                    output_file.write(line.rstrip() + '\t' + "Normal peptide" + '\n')
            output_file.close()                               
