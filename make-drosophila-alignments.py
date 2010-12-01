#! /usr/local/bin/python

import sys, os, math, string, random, pickle
import stats, biofile, newick
from optparse import OptionParser

def readORFs(fname):
	orfs = []
	for line in file(fname,'r').readlines():
		if line[0] != '#':
			orf = line.strip().split('.')[0]
			orfs.append(orf)
	return orfs

def readFlybaseMapping(fname):
	id_map = {}
	# Ensembl Gene ID	Ensembl Transcript ID	Ensembl Protein ID	FlyBase gene id	FlyBase transcript id	FlyBase protein id
	for line in file(fname,'r').readlines():
		if line[0] != '#':
			flds = line.strip().split('\t')
			fb_gene_id = flds[3]
			fb_trans_id = flds[4]
			id_map[fb_trans_id] = fb_gene_id
	return id_map

def readOneToOneOrthologs(fname, master_spec, tree_species):
	ortho_dict = {}
	# cluster_id	classification	dana	dere	dgri	dmel	dmoj	dper	dpse	dsec	dsim	dvir	dwil	dyak
	specs = ['dana','dere','dgri','dmel','dmoj','dper','dpse','dsec','dsim','dvir','dwil','dyak']
	master_id = specs.index(master_spec)
	ids = [specs.index(x) for x in tree_species]
	for line in file(fname,'r').readlines():
		if line[0] != '#':
			flds = line.strip().split('\t')
			orf_list = flds[2:14]
			(dana,dere,dgri,dmel,dmoj,dper,dpse,dsec,dsim,dvir,dwil,dyak) = orf_list
			classif = [x for x in flds[1]]
			spec_class = ''.join([classif[x] for x in ids])
			if spec_class == '1'*len(tree_species): # All 1:1 orthologs
				orfs = [orf_list[x] for x in ids]
				master_orf = orf_list[master_id]
				spec_orf_list = zip(tree_species, orfs)
				ortho_dict[master_orf] = spec_orf_list
	return ortho_dict

if __name__ == "__main__":
	# Goal: turn directory of files into alignments
	# (al_len, spec_orf_pairs, aligned_prots) = alignment_dict[orf]


	parser = OptionParser(usage="%prog [options] <conservation type> <master species> <map file linking species ID to FASTA ORFeome> " + \
						  "<base directory for FASTA ORFeome> <FASTA ID format> <alignment pickle file> <Newick tree for amino-acid conservation>" + \
						  "<Newick tree for codon conservation>")
	parser.add_option("-n", "--n", dest="num_to_align", type="int", default=None, help="how many alignments to produce")
	(options, args) = parser.parse_args()
									
	master_spec = args[0]
	tree_fname = args[1]
	in_dir = os.path.expanduser(args[2])
	ortholog_fname = os.path.expanduser(args[3])
	id_map_fname = os.path.expanduser(args[4])
	alignment_out_fname = os.path.expanduser(args[5])
	ortholog_out_fname = os.path.expanduser(args[6])

	# Get tree
	if os.path.isfile(tree_fname):
		tree = newick.tree.readTree(file(tree_fname,'r'))
	else:
		# Interpret tree as string
		tree = newick.tree.parseTree(tree_fname)
	tree_species = [x.name for x in tree.leaves]
	assert master_spec in tree_species
	print "# Species: (%s)" % ','.join(tree_species)

	orfs = []
	dirlist = os.listdir(in_dir)
	for x in dirlist:
		if x.startswith("FB"):
			orfs.append(x.strip().split(".")[0])

	print "# Found %d alignments" % len(orfs)
	id_map = readFlybaseMapping(id_map_fname)
	ortho_dict = readOneToOneOrthologs(ortholog_fname, master_spec, tree_species)
	print "# Found %d 1:1 ortholog sets" % len(ortho_dict.keys())

	n_to_align = min(len(orfs), options.num_to_align)
	alignment_dict = {}
	ortholog_dict = {}
	n_written = 0
	n_failed = 0
	n_duplicates = 0
	for trans_id in orfs[0:n_to_align]:
		fname = os.path.join(in_dir,'%s.fasta' % trans_id)
		orf_alignment_dict = biofile.readFASTADict(fname)
		try:
			# Alignments are FBtr transcript IDs
			# Orthologs are FBgn gene IDs
			# ID map turns FBtr into FBgn
			gene_id = id_map[trans_id]
			#print trans_id, gene_id
			spec_orf_list = ortho_dict[gene_id]
			#print trans_id, spec_orf_list
			spec_orf_dict = dict(spec_orf_list)
			del spec_orf_dict[master_spec]
			spec_orf_dict[master_spec] = trans_id
			new_spec_orf_list = spec_orf_dict.items()
			#print trans_id, gene_id, new_spec_orf_list
			#print trans_id, orf_alignment_dict.keys()
			#print trans_id, n_written, orf_alignment_dict[master_spec][0:10]
			# Now fetch alignment 
			master_seq = orf_alignment_dict[master_spec]
			alignment_spec_orf_list = [(spec,orf) for (spec,orf) in new_spec_orf_list if spec in orf_alignment_dict.keys()]
			protal = [orf_alignment_dict[spec] for (spec, sorf) in alignment_spec_orf_list]
			al_len = len(protal)
			try:
				alignment_dict[gene_id].append((len(protal[0]),(al_len, alignment_spec_orf_list, protal)))
				n_duplicates += 1
			except KeyError:
				alignment_dict[gene_id] = [(len(protal[0]),(al_len, alignment_spec_orf_list, protal))]
			#print "# Wrote %s" % ','.join(["%s-%s" % (s,o) for (s,o) in alignment_spec_orf_list])
		except KeyError, ke:
			#print "# Can't find ID %s" % ke
			n_failed += 1
			continue

	final_alignment_dict = {}
	for gene_id in alignment_dict.keys():
		# For alignments with multiple hits -- multiple transcripts -- pick longest alignment
		align_list = alignment_dict[gene_id]
		if len(align_list) > 1:
			# More than one alignment
			align_list.sort()
		(prot_len,(al_len, alignment_spec_orf_list, protal)) = align_list[0]
		trans_id = dict(alignment_spec_orf_list)[master_spec]
		final_alignment_dict[trans_id] = (al_len, alignment_spec_orf_list, protal)
		ortholog_dict[trans_id] = alignment_spec_orf_list
		n_written += 1
	pickle.dump(final_alignment_dict, file(alignment_out_fname,'w'))
	print "# Wrote %d alignments (%d failures, %d duplicates) to %s" % (n_written, n_failed, n_duplicates, alignment_out_fname)
	pickle.dump(ortholog_dict, file(ortholog_out_fname,'w'))
	print "# Wrote %d 1:1 orthologs to %s" % (len(ortholog_dict.values()), ortholog_out_fname)
