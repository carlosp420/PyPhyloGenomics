'''
Created on Nov 29, 2011

@author: chriswheat
'''

from Bio.Blast.Applications import NcbiblastxCommandline
import sys
import subprocess
from Bio import SeqIO
import os


def Make_and_Blast_DB(datafolder,databasefolder,db,query,outfile,blast_type,softwarefolder):
    if blast_type == "tblastn" or blast_type == "blastp":
        seq=SeqIO.parse(datafolder+QueryFile, "fasta")
        outFile=open(datafolder+QueryFile+".AA.fas", "w")
        for i in seq:
            outFile.write(str(">"+i.description+"\n"+i.seq.translate()+"\n"))
        outFile.close() 
        blast_query = datafolder+QueryFile+".AA.fas"
    else:
        blast_query = datafolder+query    
    blast_db = databasefolder+db
    print blast_db
    blast_out = outfile
#    softwarefolder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/localapps/"
    blast_exe = softwarefolder+blast_type 
    print blast_exe
    format_exe = softwarefolder+"makeblastdb"  
    
    'uncomment to build database'    
    if blast_type == "blastp":
        cline =format_exe+" -in "+blast_db+" -dbtype prot" # this is a straightforward way of generating the cmd call
    else:
        cline =format_exe+" -in "+blast_db+" -dbtype nucl" # this is a straightforward way of generating the cmd call
    return_code = subprocess.call(str(cline), shell=(sys.platform!="linux"))
    print return_code, "\n"
    print "database created \n"
    
    print "blasting: ", blast_query
    print "against db = ",blast_db
    cline1 = NcbiblastxCommandline(cmd=blast_exe, query=blast_query, db=blast_db, evalue=0.00001, outfmt=6, out=blast_out, num_threads=2)
    return_code1 = subprocess.call(str(cline1), shell=(sys.platform!="linux")) # this is calling the blast
    print return_code1
    print "blast finished"
    
    
' input parameters for running'
softwarefolder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/thirdparty/ncbi-blast-2.2.23/bin/"
databasefolder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/databases/"
datafolder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/LongExons/"
#DBFile = "Dplex_geneset_OGS1_pep.fasta"
DBFile = "Dplex_genome_v2.fasta"
'blastn tblastx tblastn'
BlastType = "tblastn"
#AA_blastfile = "BmoriSingleCopyOrthologs.txt.AA.fas"
QueryFile = "BmoriLongExonsInFrame.fas"
Blast_outfile=QueryFile+"_"+BlastType+"_"+DBFile+".txt"

def main():
    os.chdir(datafolder)
    Make_and_Blast_DB(datafolder,databasefolder,DBFile,QueryFile,Blast_outfile,BlastType,softwarefolder)
    
main()