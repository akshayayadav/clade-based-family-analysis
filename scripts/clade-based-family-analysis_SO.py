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

#checks for one or more outgroup sequences in the list of nodes
def check_for_outgroups(leaf_arr, outgrp_re):
	outgrp_matches=list(filter(outgrp_re.search, leaf_arr))
	if(len(outgrp_matches)==0):
		return 1
	else:
		return 0
#get nodes in the family trees that form monophyletic ingroup clades 
def get_ingroup_monophyletic_clade_nodes(fam_tree, outgrp_re):
	ingroup_clade_arr=list()
	node_dict={}
	for node in fam_tree.traverse("levelorder"):
		leaf_arr = get_node_leaves(node)
		if(check_for_outgroups(leaf_arr, outgrp_re)):
			flag=check_for_clade_overlap(leaf_arr, ingroup_clade_arr)
			if(flag==1):
				node_dict[node]=1
				#print node
				#print "********************************"
	
	return(node_dict)
############################################################################
def get_SO_duplication_events(fam_tree, node_dict, species_dict, fam_tree_fileName):
	SO_events = fam_tree.get_descendant_evol_events()
	clade_counter=0
	for ev in SO_events:
		if (node_dict.has_key(ev.node)):
			ingrp_clade_node=ev.node
			clad_flag=check_ingroup_clade_composition(ingrp_clade_node, species_dict)
			if not (clad_flag):
				continue
			clade_counter+=1
			if(ev.etype == "D"):		
				print '{0} {1} {2} {3}'.format(os.path.splitext(os.path.basename(fam_tree_fileName))[0], clade_counter, ingrp_clade_node.support,"D")
			else:
				print '{0} {1} {2} {3}'.format(os.path.splitext(os.path.basename(fam_tree_fileName))[0], clade_counter, ingrp_clade_node.support, "S")

#checks if the ingroup monophyletic clade contains expected species counts					
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
	species_representation_cutoff=0.7
	no_of_represented_species=0
	for cl_sp in ingrp_clade_species_dict:
		if(ingrp_clade_species_dict[cl_sp]>=species_dict[cl_sp]):
			no_of_represented_species+=1
	if((no_of_represented_species/len(species_dict))>=species_representation_cutoff):
		return(no_of_represented_species)
	else:
		return(0)

#ignoring unrooted trees - unrooted are multifucation at the root
def detect_multifurcation(fam_tree):
	root = fam_tree.get_tree_root()
	outgroups = root.get_children()
	if (len(outgroups) != 2):
		return 0
	else:
		return 1

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

#execute the workflow for a given family tree file
def process_family_tree(fam_tree_fileName, profile_fileName):

	outgrp_regex_str, species_dict = read_profile_file(profile_fileName)
	outgrp_re = re.compile(outgrp_regex_str)

	#for lgf5 families
	#outgrp_re = re.compile(r'arath|prupe|solly|cucsa|vitvi')

	#for lgf5 families
	#species_dict={"glyma":2, "phavu":2, "vigan":2, "vigra":2, "vigun":2, "cajca":2, "medtr":2, "cicar":2, "tripr":2, "lotja":2, "lupan":2, "aradu":2, "araip":2, "arahy":2}




	#for orthofinder families
	#outgrp_re = re.compile(r'prupe\_pep\_prupe|cucsa\_pep\_cucsa')

	#for orthofinder families
	#species_dict={"glyma_pep_glyma":2, "phavu_pep_phavu":2, "vigan_pep_vigan":2, "vigra_pep_vigra":2, "vigun_pep_vigun":2, "cajca_pep_cajca":2, "medtr_pep_medtr":2, "cicar_pep_cicar":2, "lotja_pep_lotja":2, "lupan_pep_lupan":2, "arahy_pep_arahy":2}


	#for orthofinder cerca-bauto
	#outgrp_re = re.compile(r'nissc\_nissc')

	#for orthofinder cerca-bauto
	#species_dict={"cerca_cerca":2, "bauto_bauto":2 }
	#species_dict={"cerca_cerca":2, "bauto_bauto":1 }
	#species_dict={"cerca_cerca":1, "bauto_bauto":2 }

	#for orthofinder chafa-mimpu
	#outgrp_re = re.compile(r'nissc\_nissc')

	#for orthofinder chafa-mimpu
        #species_dict={"chafa_chafa":2, "mimpu_mimpu":2 }
        #species_dict={"chafa_chafa":2, "mimpu_mimpu":1 }
        #species_dict={"chafa_chafa":1, "mimpu_mimpu":2 }


	
	#for cerca-bauto
	#outgrp_re = re.compile(r'adu|aip|ath|cfa|car|csa|gma|lja|lan|mtr|mpu|nsc|pvu|ppe|sly|tpr|vun|vvi')

	#for cerca-bauto
	#species_dict={"cca":2, "bto":2 }
	#species_dict={"cca":2, "bto":1 }
	#species_dict={"cca":1, "bto":2 }



	#for chafa-mimpu
	#outgrp_re = re.compile(r'adu|aip|ath|bto|cca|car|csa|gma|lja|lan|mtr|nsc|pvu|ppe|sly|tpr|vun|vvi')
	
	#for chafa-mimpu
	#species_dict={"cfa":2, "mpu":2}
	#species_dict={"cfa":2, "mpu":1}
	#species_dict={"cfa":1, "mpu":2}

	

	#for legumes
	#outgrp_re = re.compile(r'ath|ppe|sly|csa|vvi')

	#for legumes
	#species_dict={"gma":2, "pvu":2, "vun":2, "mtr":2, "car":2, "tpr":2, "lja":2, "lan":2, "adu":2, "aip":2, "bto":2, "cca":2, "cfa":2, "mpu":2, "nsc":2}


	#for chafa-mimpu-papil
	#outgrp_re = re.compile(r'ath|bto|cca|ppe|sly|csa|vvi')

	#for chafa-mimpu-papil
	#species_dict={"gma":2, "pvu":2, "vun":2, "mtr":2, "car":2, "tpr":2, "lja":2, "lan":2, "adu":2, "aip":2, "cfa":2, "mpu":2, "nsc":2}


	#for papil
	#outgrp_re = re.compile(r'ath|bto|cca|cfa|mpu|ppe|sly|csa|vvi')

	#for papil
	#species_dict={"gma":2, "pvu":2, "vun":2, "mtr":2, "car":2, "tpr":2, "lja":2, "lan":2, "adu":2, "aip":2, "nsc":2}
	


	fam_tree = PhyloTree(fam_tree_fileName, format=1)
	if not (detect_multifurcation(fam_tree)):
		return 0
	node_dict = get_ingroup_monophyletic_clade_nodes(fam_tree, outgrp_re)
	get_SO_duplication_events(fam_tree, node_dict, species_dict, fam_tree_fileName)



#takes name of directory containing family tree files and sends for processing
def read_fam_tree_files_from_directory(fam_tree_dirName, profile_fileName):
	for fam_tree_fileName in os.listdir(fam_tree_dirName):
		process_family_tree(fam_tree_dirName+"/"+fam_tree_fileName, profile_fileName)

fam_tree_dirName = sys.argv[1]
profile_fileName = sys.argv[2]
read_fam_tree_files_from_directory(fam_tree_dirName, profile_fileName)

