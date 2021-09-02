import glob
import os
import pandas as pd
import numpy as np
from snakemake.utils import validate, min_version

#### GLOBAL scope functions ####
def get_resource(rule,resource) -> int:
	'''
	Attempt to parse config.yaml to retrieve resources available for a given
	rule. It will revert to default if a key error is found. Returns an int.
	with the allocated resources available for said rule. Ex: "threads": 1
	'''

	try:
		return config['resources'][rule][resource]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]

def get_params(rule,param) -> int:
	'''
	Attempt to parse config.yaml to retrieve parameters available for a given
	rule. It will crash otherwise.
	''' 
	try:
		return config['parameters'][rule][param]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve parameter for {rule}/{param}: Exiting...')
		sys.exit(1)


#### GLOBAL PARAMETERS ####

min_version('6.2.1')

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
DATA = config['samples']
CHR_LIST = config['chr_list']


#### Load rules ####
include: 'rules/vep_annotation.smk'


rule all:
	input:
		expand(f"{OUTDIR}/final/{{CHR}}.txt", CHR=CHR_LIST)
	threads:
		get_resource('default', 'threads')
	resources:
		mem=get_resource('default', 'mem'),
		walltime=get_resource('default', 'walltime')
	shell:
		f"rm -r {OUTDIR}/split; rm -r {OUTDIR}/annotate_no_header"
