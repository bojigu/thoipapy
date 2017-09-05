#! /usr/bin/python

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import sys

# Parses a BLAST output file in XML format and returns a tab separated table
#ARGV[1] is the filepath to blast output in XML format
#ARGV[2] is the parser output filepath
#ARGV[3] is the output filepath to a list of cluster names (from the representative set of sequences) and hit names (from the blast)

blast_records = NCBIXML.parse(open("my_blast.xml"))
parsed_out = open("parseout", "w")
hit_dict = open("hit_dict", "w")

#Write the header to the parsed output
parsed_out.write("Cluster_name\t" + "Query_name\t" + "Query_length\t" + "Hit_name\t" + "Query_start\t" + "Query_end\t" + "Frame\t" + "Hit_start\t" + "Hit_end\t" + "ID\n")

#Iterate through blast records
for record in blast_records:
	#hit_counter = 0
	#Get cluster name of the representative sequence
	cluster_name = record.query.split()[0]
	query_name = record.query.split()[1]
	print(record.query)
	for desc in record.descriptions:
		if desc.num_alignments < 2:
			for aln in record.alignments:
				hit_name = aln.hit_def
				for hsp in aln.hsps:
					#if hit_counter == 0:
						#if len(hsp.query) > 200:
					query_len = hsp.align_length
					query_start = hsp.query_start
					query_end = hsp.query_end
					hit_start = hsp.sbjct_start
					hit_end = hsp.sbjct_end
					frame = hsp.frame[1]
					hit_id = float(hsp.identities) / float(query_len)
					parsed_out.write(cluster_name + "\t" + query_name + "\t" + str(query_len) + "\t" + hit_name + "\t" + str(
				query_start) + "\t" + str(query_end) + "\t" + str(frame) + "\t" + str(hit_start) + "\t" + str(
				hit_end) + "\t" + str(hit_id) + "\n")
					#hit_counter += 1
hit_dict.close()
parsed_out.close()