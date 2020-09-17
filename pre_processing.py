import os 
import sys
import subprocess

def fastqc_run(filename, directory_in, directory_out):
	os.system('/programs/FastQC-0.11.5/fastqc -q -o ' + directory_out + ' ' + directory_in + filename)

def read_count_file(directory_rawfastqc,directory_trimfastqc,directory_al_bam):
	read_dict = {}
	for filename in os.listdir(directory_rawfastqc):
		if filename.split('_')[-1]=='fastqc.zip':
			inp_1 ='unzip -c ' + directory_rawfastqc + filename + ' ' +filename.split('.')[0]+'/fastqc_data.txt | sed -n 9p'
			out_1 = subprocess.check_output(inp_1,shell=True)
			filename_2 = 'trimmed_'+filename.split('Lexogen3_')[-1]
			inp_2 ='unzip -c ' + directory_trimfastqc + filename_2 + ' ' +filename_2.split('.')[0]+'/fastqc_data.txt | sed -n 9p'
			out_2 = subprocess.check_output(inp_2,shell=True)
			key = filename.split('Lexogen3_')[-1].split('_fastqc')[0]
			value_1 = out_1.split()[-1]
			value_2 = out_2.split()[-1]
			read_dict[key]=[value_1,value_2]
	for filename in os.listdir(directory_al_bam):
		if filename.split('R1')[-1] == 'Log.final.out':
			inp_3 = "grep 'Uniquely mapped reads number' " + directory_al_bam + filename
			inp_4 = "grep 'Uniquely mapped reads %' " + directory_al_bam + filename
			out_3 = subprocess.check_output(inp_3,shell=True)
			out_4 = subprocess.check_output(inp_4,shell=True)
			uniq_mapped_reads = out_3.split()[-1]
			uniq_percent_reads = out_4.split()[-1][0:4]
			key = filename.split("Log")[0].split('trimmed_')[-1]
			read_dict[key].append(uniq_mapped_reads)
			read_dict[key].append(uniq_percent_reads)
	out_file = open('mitoY_rnaseq_read_count.txt','w')
	out_file.write('Sample' + '\t' + 'Raw'+'\t'+ 'Trimmed'+'\t' + 'Mapped' + '\t' + 'Percent' +'\n')
	for key in read_dict:
		out_file.write(key+'\t'+read_dict[key][0]+'\t'+read_dict[key][1]+ '\t' + read_dict[key][2] + '\t' + read_dict[key][3]+'\n')
	out_file.close()

##############################################################################################################################
######################  RELEVANT DIRECTORES ##################################################################################
##############################################################################################################################

os.chdir('/workdir/mam737/mitoY_rnaseq_final/')

dir_1 = './raw_files/'
dir_2 = './raw_fastqc/'
dir_3 = './trimmed_files/'
dir_4 = './trimmed_fastqc/'
dir_5 = './aligned_bam/'
dir_6 = './htseq_counts/'

##############################################################################################################################
######################  RELEVANT DIRECTORES ##################################################################################
##############################################################################################################################



##############################################################################################################################
######################  QUALITY CONTROL ON RAW DATA ##########################################################################
##############################################################################################################################

num_files = next(os.walk(dir_1))[2]
if len(num_files) != 72:
	print('Missing Samples. Please Check Raw Data')

##Check that files have been unzipped##
for filename in os.listdir(dir_1):
	if filename.endswith('.gz'):
		print("Unzipping " + filename)
		os.system('gunzip ' + dir +filename)

##Check that FASTQC is done on all raw reads##
##If not, run FASTQC on all files##
if not os.path.exists(dir_2):
	os.makedirs(dir_2)
	for filename in os.listdir(dir_1):
		print(filename)
		fastqc_run(filename,dir_1,dir_2)

##Make MultiQC for Raw Reads#
if not os.path.exists('./raw_multiqc/'):
	#to run on our computers, need to set up the environment parameter
	# run these before running script
	#export PYTHONPATH=/programs/multiqc/lib/python3.6/site-packages
	#export PATH=/programs/multiqc/bin:$PATH

	#run multiqc
	os.system('multiqc -o raw_multiqc -n raw_multiqc_report.html ./raw_fastqc/8436_2804_54202_HMFLTBGX2_72xLexogen3_*_R1_fastqc.zip')

##############################################################################################################################
######################  QUALITY CONTROL ON RAW DATA ##########################################################################
##############################################################################################################################



##############################################################################################################################
######################  TRIMMING  ############################################################################################
##############################################################################################################################

 ##Trimming using Trimmomatic##
if not os.path.exists(dir_3):
	os.makedirs(dir_3)
	for filename in os.listdir(dir_1):
		os.system('java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE ' + dir_1+filename + " " + dir_3 + 'trimmed_' +filename.split('Lexogen3_')[-1] + ' ILLUMINACLIP:/workdir/mam737/ref_seqs/truseq.fa:2:30:10 LEADING:3 TRAILING:3 HEADCROP:10 SLIDINGWINDOW:4:20 MINLEN:20')#


###FastQC Trimmomatic ##
if not os.path.exists(dir_4):
	os.makedirs(dir_4)
	for filename in os.listdir(dir_3):
		print(filename)
		fastqc_run(filename,dir_3,dir_4)#

###MultiQC Trimmomatic##
if not os.path.exists('./trimmed_multiqc/'):
	os.system('multiqc -o trimmed_multiqc -n trimmed_multiqc_report.html ./trimmed_fastqc/*_R1_fastqc.zip')

##############################################################################################################################
######################  TRIMMING  ############################################################################################
##############################################################################################################################



##############################################################################################################################
######################  STAR ALIGNMENT #######################################################################################
##############################################################################################################################

#NOTE STAR INDEXES HAVE ALREADY BEEN GENERATED USING THE FOLLOWING COMMAND
#nohup /programs/STAR-2.6/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /workdir/mam737/ref_seqs/STARindex --genomeFastaFiles /workdir/mam737/ref_seqs/dmel_r6.18/dmel-all-chromosome-r6.18.fasta --sjdbGTFfile /workdir/mam737/ref_seqs/dmel_r6.18/updated_dmel_r6.18_gffread.gtf &#

## STAR on TRIMMOMATIC##
if not os.path.exists(dir_5):
	os.makedirs(dir_5)
	print('Running Star on Trimmomatic Files')
	for filename in os.listdir(dir_3):
		os.system('/programs/STAR-2.6/bin/Linux_x86_64/STAR --genomeDir /workdir/mam737/ref_seqs/STARindex/ --readFilesIn ' + dir_3 + filename + ' --runThreadN 12 --outFileNamePrefix ' + dir_5 + filename.split('.')[0] + ' --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts')

## ALIGNED MULTIQC TRIMMOMATIC ##
if not os.path.exists('./aligned_multiqc/'):
	os.system('multiqc -o aligned_multiqc -n aligned_multiqc_report.html ./trimmed_fastqc/*_R1_fastqc.zip ' + dir_5)

if not os.path.isfile('./mitoY_rnaseq_read_count.txt'):
	read_count_file(dir_2,dir_4,dir_5)

##############################################################################################################################
######################  READ COUNTING ########################################################################################
##############################################################################################################################

####### NOTE ######
#to run on our computers, need to set up the environment parameter
# run these before running script
#export PYTHONPATH=/programs/HTSeq-0.11.0/lib/python3.6/site-packages/
#export PATH=/programs/HTSeq-0.11.0/bin:$PATH
# you may get an error if you used the similar set of commands above to run multiqc
# if so, exit out, and relogin

if not os.path.exists(dir_6):
	os.makedirs(dir_6)
	for filename in os.listdir(dir_5):
		if filename.endswith('.bam'):
			os.system('python /programs/HTSeq-0.11.0/bin/htseq-count -f bam -i gene_id ' + dir_5 + filename + ' /workdir/mam737/ref_seqs/dmel_r6.18/updated_dmel_r6.18_gffread.gtf > ' + dir_6 +  filename.split("Aligned")[0] + '_output_basename.counts')


##############################################################################################################################
######################  READ COUNTING ########################################################################################
##############################################################################################################################
	