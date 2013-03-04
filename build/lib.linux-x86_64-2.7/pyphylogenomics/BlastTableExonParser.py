'''
Created on Aug 4, 2011

@author: chriswheat
'''


from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
import os
import sys
from operator import itemgetter, attrgetter  # for sorting


def QueryExonParser(blast_in, out_file,THRESHOLD_EVALUE,THRESHOLD_IDENT,THRESHOLD_LEN):
    inFile=open(blast_in, "r") 
    query_subj_dict={}
    MyQIDList=[]
    MySIDList=[]       
    for line in inFile:
        if 'Query' in line:
            continue # this makes the loop go back to the start again, skipping the rest of the operations below     
        line=line.strip()
        fragment=line.split()
# alignment_length is the length of the actual alignment, so since this is BLASTx, it is in amino acids
# Query    Subject    %_id    alignment_length    mistmatches    gap_openings    q.start    q.end    s.start    s.end    e-value    bit_score
#    header="Query [0],Subject[1],%_id[2],alignment_length[3],mistmatches[4],gap_openings[5],q.start[6],q.end[7],s.start[8],s.end[9],e-value[10],bit_score[11]\n"
        if float(fragment[10])<= THRESHOLD_EVALUE and int(fragment[3])>= THRESHOLD_LEN and float(fragment[2])>= THRESHOLD_IDENT:
            if fragment[1].find("|") > 0 :
                dbID = fragment[1].split("|")[2]
            else:
                dbID = fragment[1]

            if fragment[0] not in MyQIDList:
                query_subj_dict.setdefault(fragment[0],fragment)
                MyQIDList.append(fragment[0]) # this makes it only take the top hit in the list, so it filters out redundancy in the query hits
#                resultFile.write(line+"\n") #these outputs are ranked upon evalue
            if fragment[0] in MyQIDList:
                aln_len = int(fragment[3]) #this is the amino acid alignment length
                if aln_len > query_subj_dict[fragment[0]][3]:
                    query_subj_dict.setdefault(fragment[0],fragment)
    return query_subj_dict

def ScaffoldParser(QEP_dict,MinExonicDistance):
    exon_blasts = QEP_dict.values()
    scaffolds = []
    exon_scaffolds = []
    exon_start_dict = {}
    exon_scaffold_dict = {}
    for exon in exon_blasts:
        scaffolds.append(exon[1])
        exon_start_dict.setdefault(exon[0],exon[8])
        exon_scaffold_dict.setdefault(exon[0],exon[1])
        exon_scaffolds.append(exon[0]+" "+exon[1])
    unique_scaffolds=set(exon_scaffold_dict.values())
    IDs = exon_scaffold_dict.keys()
    sorted_exon_blasts = sorted(exon_blasts, key=itemgetter(1,0)) #sorting by item1, then item0
    all_space_filtered_ids=[]
    for scaffold in unique_scaffolds:
        ids_on_scaffold = []
        for id in IDs:
            if exon_scaffold_dict[id] == scaffold:
                ids_on_scaffold.append(id)
        IDs_exon_location = []
        for gene in ids_on_scaffold:
            IDs_exon_location.append(gene+" "+exon_start_dict[gene])
        space_filtered_ids = SetChromoDistanceMin(IDs_exon_location,MinExonicDistance)
        all_space_filtered_ids.append(space_filtered_ids)
    print "candidate exons: ", len(exon_blasts)
    print "spread across scaffolds: ", len(unique_scaffolds) 
    all_filt_ids = clean_line_object(all_space_filtered_ids)
    return all_filt_ids

def SetChromoDistanceMin(IDs_exon_location,MinDist):
    entries = []
    for item in IDs_exon_location:
        entries.append(item.split())
    sorted_list = sorted(entries, key=itemgetter(1))
    loci = []
    loci.append(sorted_list[0][0])
    lastlocusstart = int(sorted_list[0][1])
    for item in sorted_list:
        diff = int(item[1]) - int(lastlocusstart)
        if diff > MinDist:
            loci.append(item[0])
            lastlocusstart = item[1]
    return loci
    

def getAlignmentLength(query_start, query_end):  #this is the amino acid alignment length
    if query_start < query_end:
        #add +1 because it doesn't start from zero    
        return query_end - query_start + 1
    else:
        return query_start - query_end + 1

def from_query_end(query_start, query_end, length):
    #Returns hit distance from query sequence end
    if query_start < query_end:
        return length - query_end
    else:
        return length - query_start

def from_query_start(query_start, query_end):
    #Returns hit distance from query sequence start
    if query_start < query_end:
        return query_start -1
    else:
        return query_end -1

def aln_subject_start(subject_start, subject_end):
    if subject_start < subject_end:
        return subject_start -1
    else:
        return subject_end -1

def aln_subject_end(subject_start, subject_end):
    if subject_start < subject_end:
        return subject_end
    else:
        return subject_start

def GetExonsInFrame(QEP_dict,ids,wholefastafiles):
    exon_inframes={}
    for id in ids:
#        print id
#        print wholefastafiles[id].seq
#        print QEP_dict[id]        
#        print QEP_dict[id][6]
#        print QEP_dict[id][7]
#        difference = int(QEP_dict[id][7]) - int(QEP_dict[id][6])
#        print difference
#        if difference > 0:
#            direction = "FORWARD"
#        elif difference < 0:
#            direction = "REVERSE"        
#        print direction    
        exon_start = int(QEP_dict[id][6])
        exon_stop = int(QEP_dict[id][7])
#        print isinstance(exon_start, int )
#        if direction == "FORWARD":

        exon_inframe_start=[]
        exon_inframe_stop=[]
#        print "exon_start", float(exon_start)/3
#        print "exon_stop", float(exon_stop)/3
        if exon_start < 4:
            exon_inframe_start = 3
        elif exon_start > 3:
#            print float(exon_start)/3 - exon_start/3
#            print float(exon_start+1)/3 - (exon_start+1)/3            
#            print float(exon_start+2)/3 - (exon_start+2)/3
            if (float(exon_start)/3 - exon_start/3 )== 0:
                exon_inframe_start = exon_start
            elif (float(exon_start+1)/3 - (exon_start+1)/3)== 0:
                exon_inframe_start = exon_start+1
            elif (float(exon_start+2)/3 - (exon_start+2)/3 )== 0:
                exon_inframe_start = exon_start+2
                
        if (float(exon_stop)/3 - exon_stop/3 )== 0:
            exon_inframe_stop = exon_stop
        elif (float(exon_stop-1)/3 - (exon_stop-1)/3 )== 0:
            exon_inframe_stop = exon_stop-1
        elif (float(exon_stop-2)/3 - (exon_stop-2)/3 )== 0:
            exon_inframe_stop = exon_stop-2    
        
#        print exon_inframe_start, " ", exon_inframe_stop
        exon_inframes.setdefault(id,wholefastafiles[id].seq[int(exon_inframe_start):int(exon_inframe_stop)])
#    for id in ids:
#        print id, exon_inframes[id]
    return exon_inframes
         
    
def FastaGeneIDExtract(inffile,ids):
    outFile=open("temp.fas", "w")
    records_all = SeqIO.index(inffile, 'fasta')
    heads=[]
    for id in ids:
        heads.append(id)
    count = 1
    for head in heads:
        outFile.write(str(">"+head+"\n"+records_all[head].seq+"\n"))
        count = count +1
    outFile.close()
    print "Filtered fasta sequences extracted: ", count
    records_all = SeqIO.index("temp.fas", 'fasta')
    os.system('rm temp.fas -f ')
    return records_all    
    
def clean_line_object(dirtyline):
    dirtyline = str(dirtyline)
    newline = dirtyline.replace("[","")
    newline = newline.replace("]","")
    newline = newline.replace("'","")    
    newline = newline.replace(" ","")
    cleanline = newline.split(",")
    return cleanline

def Dict2Fasta(dict,outfilename):
    outFile=open(outfilename, "w")
    headers = dict.keys()
    for head in headers:
        outFile.write(str(">"+head+"\n"+dict[head]+"\n"))
    outFile.close()
    print "finished writing extracted exons into fasta file: ", outfilename

'running the program '
'filter blast results based on these threshold values'
'LEN is assumed to be the ungapped aligment length at the AA level since using blastX'
' input parameters for running'

Folder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/LongExons"
databasefolder = "/Users/chriswheat/Documents/workspace/WheatBioInfo/databases/"
BlastFile = "BmoriSingleCopyOrthologs.txt.DNA.fas_blastn_silkgenome.fas.txt"
outfile=BlastFile+".ExonData.txt"
THRESHOLD_EVALUE = 1e-6
THRESHOLD_IDENT = 98
THRESHOLD_LEN = 280 
MinExonicDistance = 810000
inffile = databasefolder+"silkcds.fas"

def main():
    os.chdir(Folder)
    QEP_dict = QueryExonParser(BlastFile, outfile,THRESHOLD_EVALUE,THRESHOLD_IDENT,THRESHOLD_LEN)
    ids = ScaffoldParser(QEP_dict,MinExonicDistance)
#    for id in ids:
#        print id
    print "Minimum Inter-Exonic distance: ", MinExonicDistance
    wholefastafiles = FastaGeneIDExtract(inffile,ids)
#    print wholefastafiles
    exon_inframes = GetExonsInFrame(QEP_dict,ids,wholefastafiles)
#    for id in ids:
#        print id, exon_inframes[id]
    Dict2Fasta(exon_inframes,"BmoriLongExonsInFrame.fas")
    
main()


