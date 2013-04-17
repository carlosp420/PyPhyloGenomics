from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import subprocess
import sys
import time

def batchAlignment(files, out_path):
	'''
	Reads in sequences from files, group homologous sequences and aligns them.

	files = a list of file names. The first file in the list should
	be the master file (e.g. Bombyx mori).
	out_path = the folder where the alignment files shall be store.
	'''
	
	for idx, f in enumerate(files): # Read each file into a dictionary.
		exec("handle%d = %s" % (idx+1, "SeqIO.index(f,'fasta')")) 

	count = 0 # To count the number of homologous groups.
	print "DO NOT FREAK OUT!"
	print "Command line shells will pop up now."
	time.sleep(3)

	print "Parsing the master file"
	for ID1 in handle1: # Parsing the master file.
		seqs_list = []
		seqs_list.append(SeqRecord(handle1[ID1].seq, id = ID1))
		print "Pooling gene %s" % ID1;

		idx = len(files)
		while idx > 1: # Parsing the other files.
			script = '''for ID2 in handle%d:
	if ID1 == ID2.split('-')[0]:
		seqs_list.append(SeqRecord(handle%d[ID2].seq, id = ID2))'''
			exec(script % (idx, idx))
			idx -= 1

		#TODO: Add a command for PSAs, for cases where the homologous group
													# has only 2 sequences.

		if len(seqs_list) > 2: # Now perfom the MSA.
			count += 1

			muscle_cline = MuscleCommandline(out=out_path+"\\"+ID1+
											 ".aln", fasta=True)
		
			child = subprocess.Popen(str(muscle_cline),
									 stdin = subprocess.PIPE,
									 stderr = subprocess.PIPE,
									 shell = (sys.platform!="win32")) # In case the user's machine is not Windows.

			SeqIO.write(seqs_list, child.stdin, "fasta")
			child.stdin.close()

			time.sleep(3) # Wait 3 seconds to not fill off the screen with CMDs.

	print "%d alignments have been saved in %s." % (count, out_path)


#### USE ####

# The first file in the list below should be Bombyx mori sequences:
# files = ['Bmori_exons.fas', 'Danaus_exons.fas','Heliconius_exons.fas','Manduca_exons.fas']
# out_path = C:\Users\...\NewFolder  # Debes crear un nuevo folder para q se guarden los alineamientos.
# batchAlignment(files, out_path)
