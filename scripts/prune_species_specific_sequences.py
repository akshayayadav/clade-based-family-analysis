#!/usr/bin/python

import re
import sys
import os
from ete3 import Tree, PhyloTree


#get leaves under any given node
def get_node_leaves(node):
	leaf_node_arr = node.get_leaves()
	leaf_arr=list()
	for leaf in leaf_node_arr:
		leaf_arr.append(leaf.name)
	
	return(leaf_arr)


#reading the profile file to get the species id from which the sequences will be pruned
def read_profile_file(profile_fileName):
	profile_file=open(profile_fileName, "r")
	prune_regex_str = ''
	for line in profile_file:
		line=line.rstrip()
		linearr=re.split(r'\s+',line)
		species_id=linearr[0]
		seq_count=int(linearr[1])
		if(seq_count!=0):
			prune_regex_str = prune_regex_str + "|" + species_id

	prune_regex_str = prune_regex_str[1:]
	return(prune_regex_str)


#checks for one or more outgroup sequences in the list of nodes
def get_sequences_for_pruning(leaf_arr, prune_re):
	prune_seq_arr=list(filter(prune_re.search, leaf_arr))
	return (prune_seq_arr)


#get prune tree using the sequences in prune_seq_arr
def prune_tree (prune_seq_arr, fam_tree):
	fam_tree.prune(prune_seq_arr, preserve_branch_length=True)



def get_prune_regex(profile_fileName):
	prune_regex_str = read_profile_file(profile_fileName)
	prune_re = re.compile(prune_regex_str)
	return (prune_re)


def process_family_tree(fam_tree_fileName, prune_re, out_dirName):
	fam_tree = PhyloTree(fam_tree_fileName, format=1)
	fam_tree_id = os.path.splitext(os.path.basename(fam_tree_fileName))[0]
	print fam_tree_id

	leaf_arr = get_node_leaves(fam_tree)
	prune_seq_arr = get_sequences_for_pruning(leaf_arr, prune_re)
	try:
		prune_tree (prune_seq_arr, fam_tree)
	except:
		return 0

	fam_tree.write(format=1, outfile=out_dirName+"/"+fam_tree_id)


def process_fam_tree_files_from_directory(fam_tree_dirName, prune_re, out_dirName):
	for fam_tree_fileName in os.listdir(fam_tree_dirName):
		process_family_tree(fam_tree_dirName+"/"+fam_tree_fileName, prune_re, out_dirName)



fam_tree_dirName = sys.argv[1]
profile_fileName = sys.argv[2]
out_dirName = sys.argv[3]

prune_re=get_prune_regex(profile_fileName)
process_fam_tree_files_from_directory(fam_tree_dirName, prune_re, out_dirName)

