"""
------------------------------------------------------------------------------
execute BLAST commands to create BLAST database
------------------------------------------------------------------------------
"""

# It should work normally in Windows. See if it works as well in Unix!

import subprocess # Python recommends to use now subprocess.Popen() instead of os.popen()

def makeblastdb():
    command = 'makeblastdb -in db.fas -dbtype nucl -parse_seqids -input_type fasta' 
    p = subprocess.Popen(command, shell=True) 

    command = 'blastdb_aliastool -dblist "db.fas" -dbtype nucl -out db.fas -title "db" '
    p = subprocess.Popen(command, shell=True)

    command = 'makembindex -input db.fas -iformat fasta -output db'
    p = subprocess.Popen(command, shell=True)
