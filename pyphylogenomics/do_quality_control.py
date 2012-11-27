import sys;
import subprocess;
import shlex;



if len(sys.argv) < 2:
	print "Error. Enter FASTQ file to process";
	sys.exit(0);

fastqFile = sys.argv[1].strip();

#
print "\nConverting quality data from phred to solexa";
cmd = ["convert_phred_2_solexa.py", fastqFile, "tmp"];
subprocess.call(cmd);
print subprocess.check_output("grep -c '^@W' tmp", shell=True).strip() + "\t: Original number of reads.\n" 


#
print "fastq_quality_filter -q 20 -p 70";
cmd = "fastq_quality_filter -q 20 -p 70 -i tmp -o filter1.fastq";
cmd = shlex.split(cmd);
process = subprocess.Popen(cmd);
process.wait();
print subprocess.check_output("grep -c '^@W' filter1.fastq", shell=True).strip() + "\t: filter1 number of reads.\n" 


#
print "fastq_quality_trimmer -t 35 -l 40";
cmd = "fastq_quality_trimmer -t 35 -l 40 -i filter1.fastq -o filter2.fastq";
cmd = shlex.split(cmd);
process = subprocess.Popen(cmd);
process.wait();
print subprocess.check_output("grep -c '^@W' filter2.fastq", shell=True).strip() + "\t: filter2 number of reads.\n" 


# 
print "Removing indexes: fastx_trimmer -f 9"
cmd = "fastx_trimmer -f 9 -i filter2.fastq -o filter3.fastq";
cmd = shlex.split(cmd);
process = subprocess.Popen(cmd);
process.wait();

print "\nFile to process: filter3.fastq\nend;";
