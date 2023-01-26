import pandas as pd
import os
import gzip
from itertools import islice
import read_file
import re
from os.path import join as join
from urllib import request
import sys

dicts_np_gene = {} # gene = NP_accession
dicts_np_seq = {}
dict_nm_np = {}

def download_refseq(db_path, infile, url):
	os.system('wget -c ' + join(url, infile) + ' -O ' + join(db_path, 'refseq.fasta.gz'))
	os.system('gunzip ' + join(db_path, 'refseq.fasta.gz'))

def download_gtf(db_path, infile, url):
	os.system('wget -c ' + join(url, infile) + ' -O ' + join(db_path, 'genomic.gtf.gz'))

def read_gtf(db_path):
	with gzip.open(join(db_path, 'genomic.gtf.gz')) as f:
		for i in islice(f, 0, None):
			df = i.decode('utf8')
			if not df.startswith('#'):
				df_split = df.split('\t')
				if df_split[2] == 'CDS' and 'protein_id' in df_split[8]:
					gene_id = re.search('gene_id\ (.*?)\;', df_split[8]).group(1)
					protein_id = re.search('protein_id\ (.*?)\;', df_split[8]).group(1)
					dicts_np_gene[protein_id.strip('"')] = gene_id.strip('"')
		
def read_fasta_file(db_path):
	rds = read_file.read_fasta(join(db_path, 'refseq.fasta'))
	for rows in rds:
		#gene = rows[0].split('|')[4].split('#')[0].strip()
		np_acc = rows[0].split(' ')[0]
		#dicts_np_gene[np_acc] = gene
		dicts_np_seq[np_acc] = rows[1].rstrip()
	
def download_icgc(infile, outdir, url): #readme file which has the cancer type details
	os.makedirs(outdir, exist_ok=True)
	for i in open(infile):
		if i.startswith(' - '):
			cancer_type = re.search('\ -\ (.*?)\ ', i.rstrip()).group(1)
			if not os.path.isfile(join(outdir, cancer_type + '.tsv.gz')):
				cancer_type = re.search('\ -\ (.*?)\ ', i.rstrip()).group(1)
				os.system('wget -c ' + url + cancer_type + '/' + 'simple_somatic_mutation.open.' + cancer_type + '.tsv.gz' + ' -O ' + join(outdir, cancer_type + '.tsv.gz'))
			else:
				print (cancer_type, " Already downloaded")

def read_enst_acc(infile):
	dicts_enst_np = {}
	with gzip.open(infile) as f:
		for j in islice(f, 0, None): 
			df = j.decode('utf8').split('\t')
			dicts_enst_np[df[0].split('.')[0]] = df[2].rstrip()
	return dicts_enst_np

def read_np_acc(infile):
	with gzip.open(infile) as f:
		for j in islice(f, 0, None):
			df = j.decode('utf8').split('\t')
			dict_nm_np[df[0].split('.')[0]] = df[1].rstrip()
	return dict_nm_np
		
def match_variants(protein_acc, var_pos):
        if protein_acc in dicts_np_seq:
                sequence = dicts_np_seq[protein_acc]
                var_location = re.search('\d+', var_pos).group(0)
                #print (var_location, len(sequence))
                if len(sequence) >= int(var_location):
                        #print (sequence[int(var_location) - 1], var_pos[0])
                        if sequence[int(var_location)-1] == var_pos[0]:
                                return 1
                        else:
                                return 0

def read_cosmic_header(infile):
	with gzip.open(infile) as f:
		for i in islice(f, 0, 1):
			df = i.decode('utf8').split('\t')
			enst = df.index('Accession Number')
			cosmic_id = df.index('LEGACY_MUTATION_ID')
			genomic_mut_id = df.index('GENOMIC_MUTATION_ID')
			aa_mut = df.index('Mutation AA')
			status = df.index('Mutation somatic status')
			primary_site = df.index('Primary site')
			chrs = df.index('Mutation genome position')
			pubmed = df.index('Pubmed_PMID')
	return enst, cosmic_id, genomic_mut_id, aa_mut, status, primary_site, chrs, pubmed

def match_variants(protein_acc, var_pos):
	if protein_acc in dicts_np_seq:
		sequence = dicts_np_seq[protein_acc]
		var_location = re.search('\d+', var_pos).group(0)
		#print (var_location, len(sequence))
		if len(sequence) >= int(var_location):
			#print (sequence[int(var_location) - 1], var_pos[0])
			if sequence[int(var_location)-1] == var_pos[0]:
				return 1
			else:
				return 0

def parse_cosmic(infile, enst_acc_file, download_dir):
	dicts_enst_np = read_enst_acc(enst_acc_file)
	write_file = open(join(download_dir, 'cosmic_db.txt'), 'w')
	header = read_cosmic_header(join(download_dir, infile))
	with gzip.open(join(download_dir, infile)) as f:
		for i in islice(f, 1, None):
			split_i = i.decode('ISO-8859-1').split('\t')
			enst = split_i[header[0]].split('.')[0]
			cosmic_id = split_i[header[1]]
			genomic_mut_id = split_i[header[2]]
			aa_mut = split_i[header[3]]
			status = split_i[header[4]]
			primary_site = split_i[header[5]]
			chrs = [split_i[header[6]] if len(split_i[header[6]]) > 1 else "-"]
			pubmed = [split_i[header[7]] if len(split_i[header[7]]) > 1 else "-"]
			if "?" not in aa_mut and '=' not in aa_mut and "*" not in aa_mut and len(aa_mut) > 1 and "_" not in aa_mut:
				if enst in dicts_enst_np:
					for np_acc in dicts_enst_np[enst].split(','):
						if np_acc in dicts_np_gene:
							matching_var_loc = match_variants(np_acc, aa_mut.split('.')[1])
							if matching_var_loc != 0:
								write_file.write(np_acc + '\t' + enst + '\t' +  dicts_np_gene[np_acc] + '\t' + \
										cosmic_id + '\t' + genomic_mut_id + '\t' +  aa_mut + '\t' + \
										 status + '\t' +  primary_site + '\t' +  pubmed[0] +  '\t' + chrs[0] + '\n')
	write_file.close()		

def download_clinvar(url, clinvar_file, download_dir):
	os.system("wget -c " + url + '/' + clinvar_file + ' -O ' + join(download_dir, clinvar_file))

def annovar_db_status(annovar_path, genome_build_version, threads):
	'''look up for database and download them accordingly'''
	if not os.path.isfile(join(annovar_path, 'humandb', genome_build_version + '_cytoBand.txt')):
		print ("downloading ANNOVAR cytoband database..............", '\n')
		os.system('perl ' + join(annovar_path, 'annotate_variation.pl') + ' -buildver ' + genome_build_version + ' -downdb cytoBand ' +  join(annovar_path, 'humandb'))
	if not os.path.isfile(join(annovar_path, 'humandb', genome_build_version + '_refGene.txt')):
		print ("downloading ANNOVAR refGene database...............", '\n')
		os.system(join(annovar_path, 'annotate_variation.pl') + ' -buildver ' +  genome_build_version + ' -downdb -webfrom annovar refGene ' +  join(annovar_path, 'humandb')) 
	if not os.path.isfile(join(annovar_path, 'humandb', genome_build_version + '_avsnp150.txt')):
		print ("downloading ANNOVAR avsnp150 database...............", '\n')
		os.system(join(annovar_path, 'annotate_variation.pl') + ' -buildver ' + genome_build_version + ' -downdb -webfrom annovar avsnp150 ' +  join(annovar_path, 'humandb'))
		#run annovar download 

def run_annovar(annovar_path, infile, download_dir, genome_build_version, threads):
	annovar_db_status(annovar_path, genome_build_version, threads)
	if not os.path.isfile(join(download_dir, 'clinvar_for_annotation.vcf')):
		print ("uncompressing the clinvar vcf file.........", '\n')
		os.system('bcftools view ' + join(download_dir, infile) + ' > ' + join(download_dir, 'clinvar_for_annotation.vcf'))
	if not os.path.isfile(join(download_dir, 'clinvar_annotatation.' + genome_build_version + '_multianno.txt')):
		print ("running annovar annotation..................",'\n')
		os.system('perl ' + join(annovar_path, 'table_annovar.pl') + ' ' + join(download_dir, 'clinvar_for_annotation.vcf') + ' '  + join(annovar_path, 'humandb') + ' -buildver ' + genome_build_version + ' -remove -protocol refGene,cytoBand,avsnp150 -operation g,r,f -nastring . -vcfinput -polish -outfile ' + join(download_dir, 'clinvar_annotatation') + ' -thread ' + threads)
	return join(download_dir, 'clinvar_annotatation.' + genome_build_version + '_multianno.txt')

def read_clinvar_header(infile):
	with open(infile) as f:
		for i in islice(f, 0, 1):
			df = i.strip().split('\t')
			gene_and_acc_and_mut_loc = df.index('AAChange.refGene') #0
			chrs = df.index('Chr') #1
			chr_loc = df.index('Start')#2
			rs_id = df.index('avsnp150')#3
			clinvar_info = df.index('Otherinfo11')#4
			location_category = df.index('Func.refGene')#5
			mutation_type = df.index('ExonicFunc.refGene')#6
			#print (gene_and_acc_and_mut_loc, chrs, chr_loc, rs_id, clinvar_info, location_category, mutation_type)
			return gene_and_acc_and_mut_loc, chrs, chr_loc, rs_id, clinvar_info, location_category, mutation_type

def parse_clinvar(infile, np_acc_file, annovar_path, download_dir, genome_build_version, threads):
	#download annovar and annotate variants; this will be the input for clinvar
	try:
		annovar_run = run_annovar(annovar_path, infile, download_dir, genome_build_version, threads)
		dict_nm_np = read_np_acc(np_acc_file)
		write_file = open(join(download_dir, 'clinvar_db.txt'), 'w')
		header = read_clinvar_header(annovar_run)
		with open(annovar_run) as f:
			for k in islice(f, 1, None):
				split_j = k.strip().split('\t')
				for j in split_j[header[0]].split(','):
					if split_j[header[5]] =="exonic" and split_j[header[6]] =="nonsynonymous SNV":
						gene = j.split(':')[0]
						nm_acc = j.split(':')[1]
						mut_loc = j.split(':')[-1].split('.')[1]
						chrs = split_j[header[1]]
						chr_loc = split_j[header[2]]
						clinvar_info = split_j[header[4]]
						location_category = split_j[header[5]]
						mutation_type = split_j[header[6]]
						rs_id = split_j[header[3]]
						
						if nm_acc in dict_nm_np:
							protein_acc = dict_nm_np[nm_acc]
							matching_var_loc = match_variants(protein_acc, mut_loc)
							if matching_var_loc != 0:
								write_file.write(gene + '\t' + nm_acc + '\t' + mut_loc + '\t' + protein_acc + '\t' + chrs + '\t' +  chr_loc + '\t' + rs_id + '\t' + clinvar_info + '\n')
		write_file.close()
	except:
		# on erorr: delete the unfinished file
		print (sys.exc_info())#
		if os.path.isfile(join(download_dir, 'clinvar_db.txt')):
			os.remove(join(download_dir, 'clinvar_db.txt'))

def read_icgc_header(infile):
	try:
		with gzip.open(infile) as f:
			for i in islice(f, 0, 1):
				df = i.decode('utf8').split('\t')
				transcript = df.index('transcript_affected')
				variant_type = df.index('consequence_type')
				mut_loc = df.index('aa_mutation')
				verification_status = df.index('verification_status')
				biological_validation_status = df.index('biological_validation_status')
				chrs = df.index('chromosome')
				chr_loc = df.index('chromosome_start')
				mu_id = df.index('icgc_mutation_id')
		return transcript, variant_type, mut_loc, verification_status, biological_validation_status, chrs, chr_loc, mu_id
	except:
		pass

def parse_icgc(indir, enst_acc_file, download_dir):
	dicts_enst_np = read_enst_acc(enst_acc_file)
	write_file = open(join(download_dir, 'icgc_db.txt'), 'w')
	for lst_file in os.listdir(indir):
		try:
			header = read_icgc_header(join(indir, lst_file))
			with gzip.open(join(indir, lst_file)) as f:
				for j in islice(f, 1, None):
					split_j = j.decode('ISO=8859-1').split('\t')
					if split_j[header[0]] in dicts_enst_np and split_j[header[1]] == 'missense_variant':
						for iter_np in dicts_enst_np[split_j[header[0]]].split(','):
							if iter_np in dicts_np_gene:
								matching_var_loc = match_variants(iter_np, split_j[header[2]])
								if matching_var_loc != 0:
									write_file.write(split_j[header[7]] + '\t' + split_j[header[0]] + '\t' + split_j[header[1]] + '\t' + \
											split_j[header[2]] + '\t' + split_j[header[3]] + '\t' + split_j[header[4]] + '\t' + \
											split_j[header[5]] + '\t' + split_j[header[6]] + '\t' + iter_np + '\t' + dicts_np_gene[iter_np] + '\n')
		except:
			print ('File corrupted or truncated ', lst_file, '\n', sys.exc_info())
	write_file.close()
								
def create_db(download_dir, db_dir, outfile, icgc_readme_file, icgc_raw_file_dir, icgc_url, cosmic_file, clinvar_file, clinvar_url, enst_acc_file, np_acc_file, annovar_path, refseq_url, refseq_db, gtf_file, gtf_url, genome_build_version, threads):
	if not os.path.isfile(join(db_dir, outfile)):
		if not os.path.isfile(join(db_dir, 'refseq.fasta')):
			download_refseq(db_dir, refseq_url, refseq_db)

		#reading fasta file
		read_fasta_file(db_dir)

		if not os.path.isfile(join(db_dir, 'genomic.gtf.gz')):
			download_gtf(db_dir, gtf_file, gtf_url)
	
		#reading gtf file
		read_gtf(db_dir)
	
		db_file = []
		for lst_db in os.listdir(download_dir):
			if '_db.txt' in lst_db:
				db_file.append(lst_db)
		print ("available databases: ", ', '.join(db_file), '\n')

		if 'cosmic_db.txt' not in db_file:
			#Check for the raw cosmic_file exist [ADD condition here]
			if os.path.isfile(join(download_dir, cosmic_file)):
				print ("creating cosmic database", '\n')
				parse_cosmic(cosmic_file, enst_acc_file, download_dir)
			else:
				print ('Please download the COSMIC database from the COSMIC repository and add the file name to config.yml', '\n')
				exit(1)
		if 'clinvar_db.txt' not in db_file:
			if not os.path.isfile(join(download_dir, clinvar_file)):
				print ('creating clinvar database', '\n')
				download_clinvar(clinvar_url, clinvar_file, download_dir)
			parse_clinvar(clinvar_file, np_acc_file, annovar_path, download_dir, genome_build_version, threads)
		if 'icgc_db.txt' not in db_file:
			if not os.path.isfile(join(icgc_readme_file, icgc_raw_file_dir)):
				download_icgc(icgc_readme_file, icgc_raw_file_dir, icgc_url)
			print ('creating icgc database', '\n')
			parse_icgc(icgc_raw_file_dir, enst_acc_file, download_dir)
		lst_db = [x for x in os.listdir(download_dir) if '_db.txt' in x]
		if len(lst_db) < 3:
			print ("All three database should be avilable for creating a master database", "\n", "Please have a look at download folders")
			exit(1)

		clinvar_info = {}
		print ('creating clinvar dictionary', '\n')
		for i in open(join(download_dir, "clinvar_db.txt")):
			split_i = i.split('\t')
			clinvar_info[split_i[3].strip() + "_" + split_i[2]] = split_i[0] + '\t' + split_i[1] + '\t' + split_i[4] + '\t' + split_i[6] + '\t' + split_i[7].rstrip()

		cosmic_info = {}
		print ('creating cosmic dictionary', '\n')
		for j in open(join(download_dir, "cosmic_db.txt")):
			split_j = j.split('\t')
			cosmic_info[split_j[0].strip() + "_" + split_j[5].split('.')[1]] = split_j[1] + '\t' + split_j[2] + '\t' + split_j[3] + '\t' + split_j[4] + '\t' + split_j[6] + '\t' + split_j[7] + \
                                                 			'\t' + split_j[8].rstrip()
		icgc_info = {}
		print ('creating icgc dictionary', '\n')
		for k in open(join(download_dir, "icgc_db.txt")):
			split_k = k.split('\t')
			icgc_info[split_k[8].strip() + "_" + split_k[3]] = split_k[0] + '\t' + split_k[1] + '\t' + split_k[2] + '\t' + split_k[4] + '\t' + split_k[5] + '\t' + split_k[6] + '\t' \
									 + split_k[8] + '\t' + split_k[9].rstrip()
	
		#creating uniq_list
		print ('creating master list of unique variants', '\n')
		master_list = {}
		for k, v in clinvar_info.items():
			master_list[k] = 1

		for k, v in cosmic_info.items():
			master_list[k] = 1

		for k, v in icgc_info.items():
			master_list[k] = 1
		print ("unique variants available in the database: ", len(master_list), '\n')
	
		write_file = open(outfile, "w")
		for key, val in master_list.items():
			split_key = key.split('_')
			if key in clinvar_info:
				clinvar = clinvar_info[key]
			else:
				clinvar = '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-'

			if key in cosmic_info:
				cosmic = cosmic_info[key]
			else:
				cosmic = '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-'

			if key in icgc_info:
				icgc = icgc_info[key]
			else:
				icgc = '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-'
			write_file.write(split_key[0] + "_" + split_key[1] + '\t' + split_key[2].strip() + '\t' \
					 + clinvar + '\t' + cosmic + '\t' + icgc + '\n')

		write_file.close()
		print ('database created sucessfully', '\n')
		#download_icgc('./../downloads/README.txt', './../downloads/icgc_rawfiles')
	
	
#create_db('./../downloads', './../database/db.txt', './../downloads/README.txt', './../downloads/icgc_rawfiles/', 'https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/' ,'./../downloads/CosmicGenomeScreensMutantExport.tsv.gz', 'clinvar_20220328.vcf.gz', 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/', './../database/ENST_and_NP_acc.txt', './../database/NM_NP_acc.txt', './../tools/annovar/', './../database/Refseq109.fasta', 'hg38', '10')
	#Listing unique_variants:
	
if __name__ == "__main__":
	if len(sys.argv) >= 16:
		download_dir = sys.argv[1]
		db_dir = sys.argv[2]
		database_path = sys.argv[3]
		icgc_readme = sys.argv[4]
		icgc_rawfile_dir = sys.argv[5]
		icgc_url = sys.argv[6]
		cosmic_file = sys.argv[7]
		clinvar_file = sys.argv[8]
		clinvar_url = sys.argv[9]
		enst_acc = sys.argv[10]
		nm_acc = sys.argv[11]
		annovar_path = sys.argv[12]
		refseq_url = sys.argv[13]
		refseq_path = sys.argv[14]
		gtf_path = sys.argv[15]
		gtf_url = sys.argv[16]
		genome_build = sys.argv[17]
		threads = sys.argv[18]
		create_db(download_dir, db_dir, database_path, icgc_readme, icgc_rawfile_dir, icgc_url, cosmic_file, clinvar_file, clinvar_url, enst_acc, nm_acc, annovar_path, refseq_url,refseq_path, gtf_path, gtf_url, genome_build, threads)
	else:
		print (len(sys.argv))


