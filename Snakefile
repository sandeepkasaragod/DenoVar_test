import os
from os.path import join as join

configfile: 'config.yml'

os.makedirs(config['download_path'], exist_ok=True)
os.makedirs(config['temp'], exist_ok=True)

rule all:
	input:
		config['db_path'],
		join(config['output'], config['input_file'])

rule download_database:
	params: 
		threads = workflow.cores,
		db_path = config['db_path'],
		genome_build = config['genome_build'],
		script_path = "scripts/create_db.py",
		clinvar_url = config['clinvar_url'],
		icgc_url = config['icgc_url'],
		clinvar_file = config['clinvar_file'],
		icgc_dir = join(config['download_path'], config['icgc_dir_path']),
		cosmic_file = config['cosmic_file'],
		refseq_file = config['refseq_file'],
		refseq_url = config['refseq_url'],
		gtf_file = config['gtf_file'],
		gtf_url = config['gtf_url']
	input:
		download_path = config['download_path'],
		database_dir = config['db_dir'],
		icgc_readme = join(config['download_path'], config['icgc_readme']),
		enst_acc = config['enst_acc'],
		nm_acc = config['nm_acc'],
		annovar_path = config['annovar_path'],
	output:
		result = config['db_path']
	shell:
		"python {params.script_path} {input.download_path} {input.database_dir} {params.db_path} {input.icgc_readme} {params.icgc_dir} {params.icgc_url} {params.cosmic_file} \
					 {params.clinvar_file} {params.clinvar_url} {input.enst_acc} {input.nm_acc} {input.annovar_path}	\
					{params.refseq_file} {params.refseq_url} {params.gtf_file} {params.gtf_url} {params.genome_build} {params.threads}"

rule denovar:
	params:
		enzyme = config['enzyme'],
		aa_min_len = config['aa_min_len'],
		aa_max_len = config['aa_max_len'],
		max_num_cleavage = config['max_cleavage'],
		prot_acc = config['protein_accession'],
		script_path = "scripts/denovar.py",
		refseq_file = "refseq.fasta",
		temp = config['temp'],
		threads = config['threads'],
		output_dir = config['output'],
	input:
		db_path = config['db_path'],
		database_dir = config['db_dir'],
		input_file = config['input_file'],
	output:
		txt = join(config['output'], config['input_file'])
	shell:
		"python {params.script_path} {input.input_file} {params.enzyme} {params.aa_min_len} {params.aa_max_len} {params.max_num_cleavage} {params.prot_acc} {input.db_path} {input.database_dir} {params.refseq_file} {params.temp} {params.threads} {params.output_dir}"
