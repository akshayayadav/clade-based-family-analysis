#!/usr/bin/env python

from __future__ import division
import re
import sys
import os
from ete3 import PhyloTree

#this function detects and selects only the largest and ingroup clade
def check_for_clade_overlap(ingroup_clade, ingroup_clade_arr):
	new_clade_detection_flag=1
	for clade in ingroup_clade_arr:
		overlap=list(set(ingroup_clade).intersection(clade))
		if(len(overlap)>0):
			new_clade_detection_flag=0
			break
	if(new_clade_detection_flag==1):
		ingroup_clade_arr.append(ingroup_clade)
		return 1
	else:
		return 0


#get leaves under any given node
def get_node_leaves(node):
	leaf_node_arr = node.get_leaves()
	leaf_arr=list()
	for leaf in leaf_node_arr:
		leaf_arr.append(leaf.name)
	
	return(leaf_arr)

#checks for one or more sequences from exclude_species_re, in the list of nodes
def check_for_outgroups(leaf_arr, outgrp_re):
	outgrp_matches=list(filter(outgrp_re.search, leaf_arr))
	if(len(outgrp_matches)==0):
		return 1
	else:
		return 0

def print_clade_sequences(node, fam_tree_fileName, clade_counter):
	leaf_node_arr = node.get_leaves()
	leaf_arr=list()
	for leaf in leaf_node_arr:
		print '{0} {1} {2}'.format(os.path.splitext(os.path.basename(fam_tree_fileName))[0] ,clade_counter,leaf.name)
	
#get nodes in the family trees that form monophyletic clades for species not in the  
def get_ingroup_monophyletic_clade_nodes(fam_tree, outgrp_re, fam_tree_fileName, species_dict):
	ingroup_clade_arr=list()
	node_dict={}
	clade_counter=1
	for node in fam_tree.traverse("levelorder"):
		leaf_arr = get_node_leaves(node)
		if(check_for_outgroups(leaf_arr, outgrp_re)):
			flag=check_for_clade_overlap(leaf_arr, ingroup_clade_arr)
			if(flag==1):
				#node_dict[node]=1
				clade_flag=check_ingroup_clade_composition(node, species_dict)
				if (clade_flag==1):
					print_clade_sequences(node, fam_tree_fileName, clade_counter)
					clade_counter+=1
				#print "********************************"
	
	return(node_dict)


def check_ingroup_clade_composition(ingrp_clade_node, species_dict):
	ingrp_clade_leaves_arr=get_node_leaves(ingrp_clade_node)
	ingrp_clade_species_dict=ingrp_clade_species_count(ingrp_clade_leaves_arr)
	no_of_represented_species=compare_species_compositions(ingrp_clade_species_dict, species_dict)
	if(no_of_represented_species):
		return(1)
	else:
		return(0)

#gets species counts for the ingroup clade		
def ingrp_clade_species_count(ingrp_clade_leaves_arr):
	ingrp_clade_species_dict={}
        for seq_id in ingrp_clade_leaves_arr:
                species_id = re.split(r'\.',seq_id)[0]
		if(ingrp_clade_species_dict.has_key(species_id)):
			ingrp_clade_species_dict[species_id]+=1
		else:
			ingrp_clade_species_dict[species_id]=1

	return(ingrp_clade_species_dict)

#comparing ingroup clade species composition with the expected species composition
def compare_species_compositions(ingrp_clade_species_dict, species_dict):
	species_representation_cutoff=1.0
	no_of_represented_species=0
	for cl_sp in ingrp_clade_species_dict:
		if(ingrp_clade_species_dict[cl_sp]>=species_dict[cl_sp]):
			no_of_represented_species+=1
	if((no_of_represented_species/len(species_dict))>=species_representation_cutoff):
		return(no_of_represented_species)
	else:
		return(0)



def read_profile_file(profile_fileName):
	profile_file=open(profile_fileName, "r")
	outgrp_regex_str = ''
	species_dict = {}
	for line in profile_file:
		line=line.rstrip()
		linearr=re.split(r'\s+',line)
		species_id=linearr[0]
		seq_count=int(linearr[1])
		if(seq_count==0):
			outgrp_regex_str = outgrp_regex_str + "|" + species_id
		else:
			species_dict[species_id] = seq_count

	
	outgrp_regex_str = outgrp_regex_str[1:]
	return([outgrp_regex_str,species_dict])



def process_fam_tree(fam_tree_fileName, exclude_regex_str, species_dict):
	exclude_species_re = re.compile(exclude_regex_str)
	fam_tree = PhyloTree(fam_tree_fileName, format=1)
	get_ingroup_monophyletic_clade_nodes(fam_tree, exclude_species_re, fam_tree_fileName, species_dict)
##########################################################################


#input family tree file. First argument
fam_tree_fileName=sys.argv[1]
profile_fileName=sys.argv[2]


#regular expression that excludes species: exclude_species_re
#dictionary/hash with expected species composition. Any given clade will be selected if it contains all the species in the dictionary with counts >= counts given in the dictionary: species_dict

exclude_regex_str, species_dict =read_profile_file(profile_fileName)
process_fam_tree(fam_tree_fileName, exclude_regex_str, species_dict)
