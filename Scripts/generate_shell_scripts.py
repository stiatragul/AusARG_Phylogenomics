#!/usr/bin/env python

"""
Ian Brennan
created 12 July 2021

"""

import pandas as pd
import subprocess
import os
import argparse
import re


def get_args():
	parser = argparse.ArgumentParser(
		description= "Generate shell scripts for each \n"
						"stage in the workflow",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# which script
	parser.add_argument(
		'--script',
		type=str,
		default=None,
		help='Which python script? Options are: '
				'dedupe_clean_filter, trinity_filtered_assembly, '
				'match_contigs, make_PRG, quality_2_sample, quality_2_lineage')

	# sample
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with sample info.'
		)	

	# parallel
	parser.add_argument(
		'--parallel',
		type=int,
		default=1,
		help='If --analyze *yes*, how many instances should be run in parallel?'
			' Keep in mind total memory usage will be --parallel x --mem'
			' and total CPU usage will be --parallel x --CPU')

	# CPUs
	parser.add_argument(
		'--CPU',
		type=int,
		default=1,
		help='# of CPUs to use for each sample'
	   )

	# Project name
	parser.add_argument(
		'--project',
		type=str,
		default=None,
		help='Name this project. This string will be added at the beginning of output files')

	# Input directory
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help='Base directory for pipeline.'
		)

	# Analyze?
	parser.add_argument(
		'--analyze',
		type=str,
		default="No",
		help='if *yes*, will run analyses from the shell script. '
				'If not included (no --analyze flag), will only output shell file')

	# memory
	parser.add_argument(
		'--mem',
		type=int,
		default=2,
		help='RAM to use for mapping/filtering in Gb'
	   	)

	# minid value
	parser.add_argument(
		'--minid',
		type=str,
		default=0.25,
		help='Minid cutoff for read acceptance '
			'see: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/'
		)	

	# reference file
	parser.add_argument(
		'--ref',
		type=str,
		default=None,
		help='File of reference sequences to filter against.'
		)

	# evalue
	parser.add_argument(
		'--evalue',
		type=float,
		default=1e-10,
		help="Minimum evalue req'd for match, given in "
				 "Xe-X format."
		)

	# database
	parser.add_argument(
		'--db',
		type=str,
		default=None,
		help="Database to use to search."
		)	

	# matches to keep
	parser.add_argument(
		'--keep',
		type=str,
		default=None,
		help='Tags to keep; comma-delimited. Options are '
			 ' easy_recip_match, complicated_recip_match '
		)

	# skip sample
	parser.add_argument(
		'--skip',
		type=str,
		default=None,
		help='Sample(s) to skip when making PRGs'
		)	

	# skip sample
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Name the output directory'
		)						

	return parser.parse_args()

def make_shells(args):

	# read in the sample file
	d = pd.read_csv(args.file)

	# get the names of all samples to run
	lineages = d['lineage'].unique().tolist()
	#lineages = d['lineage'].tolist()
	samps = d['sample'].tolist()

	# make a file to store the output
	sh_out = '%s_%s.sh' % (args.project,args.script)
	sh_out = os.path.join(args.dir,sh_out)
	print("*** shell script written to: %s ***" % sh_out)
	o = open(sh_out, 'w')

	# make a loop to create the scripts
	for ix, samps in enumerate(samps):

		# dedupe_clean_filter.py
		if args.script == "dedupe_clean_filter":
			dcf_call = "python Scripts/dedupe_clean_filter_reads.py --dir %s --file %s --sample %s --mem %s --CPU %s --minid %s --ref %s" % (args.dir, args.file, samps, args.mem, args.CPU, args.minid, args.ref)
			o.writelines(dcf_call + '\n')

		# trinity_filtered_assembly.py
		if args.script == "trinity_filtered_assembly":
			tri_call = "python Scripts/trinity_filtered_assembly.py --dir %s --sample %s --mem %s --CPU %s" % (args.dir, samps, args.mem, args.CPU)
			o.writelines(tri_call + '\n')

		# match_contigs_to_probes.py
		if args.script == "match_contigs":
			mat_call = "python Scripts/match_contigs_to_probes.py --dir %s --sample %s --evalue %s --db %s" % (args.dir, samps, args.evalue, args.db)
			o.writelines(mat_call + '\n')

		# quality_2_assembly.py
		if args.script == "quality_2_sample":
			qua_call = "python Scripts/quality_2_assembly.py --dir %s --ind %s --file %s --outdir %s" % (args.dir, samps, args.file, args.outdir)
			o.writelines(qua_call + '\n')


	# make a loop for the lineages
	for ix, lineages in enumerate(lineages):

		# make_PRG.py
		if args.script == "make_PRG":
			prg_call = "python Scripts/make_PRG.py --dir %s --lineage %s --file %s --keep %s" % (args.dir, lineages, args.file, args.keep)	
			o.writelines(prg_call + '\n')

		# quality_2_assembly.py
		if args.script == "quality_2_lineage":
			qua_call = "python Scripts/quality_2_assembly.py --dir %s --ind %s --file %s --outdir %s" % (args.dir, lineages, args.file, args.outdir)
			o.writelines(qua_call + '\n')

	# close the shell file		
	o.close()

	# return the name of the file
	return sh_out

def run_shell(shell_name, args):
	#print("running: %s" % shell_name)
	subprocess.call("parallel -j %s --bar :::: %s" % (args.parallel, shell_name), shell=True)


def main():
	# get arguments
	args = get_args()
	# make the shells
	shell_name = make_shells(args)
	# run the shell?
	if args.analyze == 'yes':
		run_shell(shell_name, args)
#		# make the function here
	

if __name__ == "__main__":
	main()
