##########################################################################################################
#Title:
#seq2HLA - HLA typing from RNA-Seq sequence reads
#
#Release: 2.2
#
#Author:
#Sebastian Boegel
#TRON - Translational Oncology at the University Medical Center Mainz, 55131 Mainz, Germany
#University Medical Center of the Johannes Gutenberg-University Mainz, Mainz, Germany
#
#License:
#The MIT License (MIT)
#Copyright (c) 2012 Sebastian Boegel
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

#Contact:
#boegels@uni-mainz.de
#
#Synopsis:
#We developed an in-silico method "Seq2HLA", written in python and R, which takes standard RNA-Seq sequence reads in fastq format 
#as input, uses a bowtie index comprising all HLA alleles and outputs the most likely HLA class I and class II genotypes (in 4 digit resolution), 
#a p-value for each call, and the expression of each class 
#
#Usage: 
#python seq2HLA.py -1 <readfile1> -2 <readfile2> -r "<runname>" [-p <int>]* [-3 <int>]**
#*optional: number of parallel search threads for bowtie optional (Default:6)
#**optional: trim int bases from the low-quality end of each read
#readfile can be uncompressed or gzipped fastq file
#runname should contain path information, e.g. "folder/subfolder/..../run", in order to store all resulting into to folder and all filenames will have the suffix run-

#Output:
#The results are outputted to stdout and to textfiles. Most important are: 
#i) <prefix>-ClassI.HLAgenotype2digits => 2 digit result of Class I
#ii) <prefix>-ClassII.HLAgenotype2digits => 2 digit result of Class II
#iii) <prefix>-ClassI.HLAgenotype4digits => 4 digit result of Class I
#iv) <prefix>-ClassII.HLAgenotype4digits => 4 digit result of Class II
#v) <prefix>.ambiguity => reports typing ambuigities (more than one solution for an allele possible)
#vi) <prefix>-ClassI.expression => expression of Class I alleles
#vii) <prefix>-ClassII.expression => expression of Class II alleles
#
#Dependencies:
#0.) seq2HLA is a python script, developed with Python 2.6.8
#1.) bowtie must be reachable by the command "bowtie". seq2HLA was developed and tested with bowtie version 0.12.7 (64-bit). The call to bowtie is invoked with 6 CPUs. You can change that by paramter -p.
#2.) R must be installed, seq2HLA.py was developed and tested with R version 2.12.2 (2011-02-25)
#3.) Input must be paired-end reads in fastq-format
#4.) Index files must be located in the folder "references".
#5.) Packages: biopython (developed with V1.58), numpy (1.7.1)

#Version history:
#2.2: improved performance, automatic detection of read length (option -l no longer required), user can choose number of parralel search threads (-p), seq2HLA now works with automatic path recognition, so it can be invoked from every path (April 2014)
#2.1: supports gzipped fastq files as input
#2.0: 4-digit typing
#1.0: 2-digit typing
###########################################################################################################

from operator import itemgetter
import sys
import linecache
import ast
import os
from Bio import SeqIO
import numpy
import operator
from optparse import OptionParser
import fourdigits
import gzip
import subprocess

#These variables need to be global, as they are filled and used by different modules
readcount={}
readspergroup={}
allelesPerLocus={}
version="2.2."

def main(runName,readFile1,readFile2,fastaClassI,fastaClassII,bowtiebuildClassI,bowtiebuildClassII,trim3,threads):
	#print version number
	print "Now running seq2HLA version "+str(version)+"!"
	#determine if input is uncompressed or gzipped fastq files and determine number of lines (important for expression calculation
	f=gzip.open(readFile1, 'r')
	#if first line can be read without error, input is gzipped
	try:
		first_line = f.readline()
                gzipped=1
		#process_wc=subprocess.Popen(['bash','-c','zcat '+readFile1+' | sed \'2q;d\' - | wc -L'],stdout=subprocess.PIPE)		
		cmd = "zcat %s | sed '2q;d' | wc -L" % readFile1
		process_wc=subprocess.Popen( ['bash','-c',cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		readlength=process_wc.communicate()[0]
		print "Input is a gipped file ....."
        except Exception, e:
		cmd = "sed '2q;d' %s | wc -L" % readFile1
                process_wc=subprocess.Popen( ['bash','-c',cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                readlength=process_wc.communicate()[0]
                gzipped=0
		print "Input is uncompressed fastq file ...."
	#as shown in the publication HLA typing with RNA-Seq works best by allowing as less mismatches as necessary
	if int(readlength)<=50:
                mismatch=1
        else:
                if int(readlength)<=100:
                        mismatch=2
                else:
                        mismatch=3
	#concatenate mapping parameiters
	mapopt="-p "+str(threads)+" -a -v"+str(mismatch)
	print "The read length of your input fastq was determined to be "+str(int(readlength))+", so "+str(mismatch)+" mismatches will be allowed and "+str(threads)+" threads will be used by bowtie."
	
	#call HLA typing for Class I
	mainClassI(runName+"-ClassI",readFile1,readFile2,bowtiebuildClassI,fastaClassI,mapopt,trim3,gzipped)
	#call HLA typing for Class II
	mainClassII(runName+"-ClassII",readFile1,readFile2,bowtiebuildClassII,fastaClassII,mapopt,trim3,gzipped)	

#---------------Class I-------------------------------
def mainClassI(runName,readFile1,readFile2,bowtiebuild,hla1fasta,mapopt,trim3,gzip):
	#-------1st iteration-----------------------------------
	print "----------HLA class I------------"
	sam1=runName+"-iteration1.sam"
	iteration=1
	print "First iteration starts....\nMapping ......"
	mapping(sam1,runName,readFile1,readFile2,bowtiebuild,1,mapopt,trim3,gzip)
	medians=[]
	medians.append(0)
	medians.append(0)
	medians.append(0)
	medianflag=False
	#Calculation of first digital haplotype.....
	output1=runName+".digitalhaplotype1"	
	print "Calculation of first digital haplotype....."
	map=createRefDict(hla1fasta,"A","B","C")
	readMapping(map,sam1)
	fourDigitString1=predictHLA1(sam1,medians,output1,medianflag)
	print "1st iteration done.\nNow removing reads that mapped to the three top-scoring groups ......."
	try:
		removeReads(runName,createRemoveList(runName,map))
	except IOError:
                print "Nothing to remove\n"
	
	#------2nd iteration------------------------------------------
	print "Second iterations starts .....\n Mapping ......"
	medians=[]
	iteration=2
	sam2=runName+"-iteration2.sam"
	newReadFile1=runName+"-2nditeration_1.fq"
	newReadFile2=runName+"-2nditeration_2.fq"
	mapping(sam2,runName,newReadFile1,newReadFile2,bowtiebuild,2,mapopt,trim3,gzip)
	medianfile=runName+".digitalhaplotype1"
	medians.append(linecache.getline(medianfile, 2).split('\t', 3)[2])
	medians.append(linecache.getline(medianfile, 3).split('\t', 3)[2])
	medians.append(linecache.getline(medianfile, 4).split('\t', 3)[2])
	medianflag=True
	output2=runName+".digitalhaplotype2"
	finaloutput=runName+".HLAgenotype2digits"
	#Calculation of second digital haplototype
	print "Calculation of second digital haplotype....."
	map=createRefDict(hla1fasta,"A","B","C")
	readMapping(map,sam2)
	fourDigitString2=predictHLA1(sam2,medians,output2,medianflag)
	print "2nd iteration done."
	reportHLAgenotype(output1,output2,finaloutput,1)
	print "Calculation of locus-specific expression ..."
	try:
		expression("A","B","C",694,694,694,map,runName,readFile1)
	except IOError:
                print "A: 0 RPKM\nB: 0 RPKM\nC: 0 RPKM\n"
	
	#-----3rd iteration in case of at least one homozygous call-----------
	f1=fourDigitString1.split(",")
	f2=fourDigitString2.split(",")
	a1=f1[0]
	a1fourDigit=f1[1]
	if a1=="no":
		a1fourDigit_solutions=f1[2]
	else:
		a1fourDigit_solutions=int(f1[2])
	b1=f1[3]
	b1fourDigit=f1[4]
	if b1=="no":
		b1fourDigit_solutions=f1[5]
	else:
		b1fourDigit_solutions=int(f1[5])
	
	c1=f1[6]
	c1fourDigit=f1[7]
	if c1=="no":
		c1fourDigit_solutions=f1[8]
	else:
		c1fourDigit_solutions=int(f1[8])
	a2=f2[0]
	b2=f2[3]
	c2=f2[6]
	
	#try-catch block to prevent IO-Error in case of no expression
	try:
		if a2=="no" or b2=="no" or c2=="no":
			removeReads(runName,createRemoveListFourDigits(runName,map,a1fourDigit,b1fourDigit,c1fourDigit))
			sam3=runName+"-iteration3.sam"
			mapping(sam3,runName,newReadFile1,newReadFile2,bowtiebuild,2,mapopt,trim3,gzip)
			readMapping(map,sam3)
			output3=runName+".digitalhaplotype3"
			fourDigitString3=predictHLA1(sam3,medians,output3,medianflag)
			f3=fourDigitString3.split(",")
			a3=f3[0]
			b3=f3[3]
			c3=f3[6]
	
		#1st allele--------------------------------
		Aallele1=a1fourDigit
		if not a1fourDigit_solutions=="":
			if a1fourDigit_solutions>1:
				Aallele1+="'"
			final_a=Aallele1+"\t"+str(getMaxP(sam1+".4digitsA1.solutions"))+"\t"
		else:
			final_a="no\tNA\t"
	
		Ballele1=b1fourDigit
		if not b1fourDigit_solutions=="":
			if b1fourDigit_solutions>1:
				Ballele1+="'"
			final_b=Ballele1+"\t"+str(getMaxP(sam1+".4digitsB1.solutions"))+"\t"
		else:
			final_b="no\tNA\t"
		
		Callele1=c1fourDigit
		if not c1fourDigit_solutions=="":
			if c1fourDigit_solutions>1:
				Callele1+="'"
			final_c=Callele1+"\t"+str(getMaxP(sam1+".4digitsC1.solutions"))+"\t"
		else:
			final_c="no\tNA\t"
		#2nd allele--------------------------------	
		if a2=="no":
			if a3=="no" or a1fourDigit_solutions==1 or not a3.split(":")[0]==a1fourDigit.split(":"):
				if a1fourDigit_solutions=="":
					final_a+="no\tNA"
				elif a1fourDigit_solutions>1:
					final_a+=Aallele1[0:-1]+"\t"+str(getP("A",finaloutput))
				else:
					final_a+=Aallele1+"\t"+str(getP("A",finaloutput))
			else:
				Aallele2=f3[1]
				if int(f3[2])>1:
					Aallele2+="'"
				final_a+=Aallele2+"\t"+str(getMaxP(sam3+".4digitsA2.solutions"))
		else:
			Aallele2=f2[1]
			if int(f2[2])>1:
				Aallele2+="'"
			final_a+=Aallele2+"\t"+str(getMaxP(sam2+".4digitsA2.solutions"))
	
		if b2=="no":
			if b3=="no" or b1fourDigit_solutions==1 or not b3.split(":")[0]==b1fourDigit.split(":"):
				if b1fourDigit_solutions=="":
					final_b+="no\tNA"
				elif b1fourDigit_solutions>1:
					final_b+=Ballele1[0:-1]+"\t"+str(getP("B",finaloutput))
				else:
					final_b+=Ballele1+"\t"+str(getP("B",finaloutput))
			else:
				Ballele2=f3[4]
				if int(f3[5])>1:
					Ballele2+="'"
				final_b+=Ballele2+"\t"+str(getMaxP(sam3+".4digitsB2.solutions"))
		else:
			Ballele2=f2[4]
			if int(f2[5])>1:
				Ballele2+="'"
			final_b+=Ballele2+"\t"+str(getMaxP(sam2+".4digitsB2.solutions"))
	
		if c2=="no":
			if c3=="no" or c1fourDigit_solutions==1 or not c3.split(":")[0]==c1fourDigit.split(":"):
				if c1fourDigit_solutions=="":
					final_c+="no\tNA"
				elif c1fourDigit_solutions>1:
					final_c+=Callele1[0:-1]+"\t"+str(getP("C",finaloutput))
				else:
					final_c+=Callele1+"\t"+str(getP("C",finaloutput))
			else:
				Callele2=f3[7]
				if int(f3[8])>1:
					Callele2+="'"
				final_c+=Callele2+"\t"+str(getMaxP(sam3+".4digitsC2.solutions"))
		else:
			Callele2=f2[7]
			if int(f2[8])>1:
				Callele2+="'"
			final_c+=Callele2+"\t"+str(getMaxP(sam2+".4digitsC2.solutions"))
		
		finaloutput4digit=runName+".HLAgenotype4digits"
		#output the 4-digit type (stdout and to file)
		reportHLA4digitGenotype(final_a,final_b,final_c,finaloutput4digit,1)
	except IOError:
		#no expression"
                print "no expression"
		finaloutput4digit=runName+".HLAgenotype4digits"
		reportHLA4digitGenotype("no\tNA\tno\tNA","no\tNA\tno\tNA","no\tNA\tno\tNA",finaloutput4digit,1)
	
#--------------Class II----------------------------------------
def mainClassII(runName,readFile1,readFile2,bowtiebuild,hla2fasta,mapopt,trim3,gzip):
	#-------1st iteration----------------------------------------
	print "----------HLA class II------------" 
	sam1=runName+"-iteration1.sam"
	iteration=1
	print "ClassII: first iteration starts....\nMapping ......"
	mapping(sam1,runName,readFile1,readFile2,bowtiebuild,1,mapopt,trim3,gzip)
	medians=[]
	medians.append(0)
	medians.append(0)
	medians.append(0)
	medianflag=False
	output1=runName+".digitalhaplotype1"	
	print "ClassII: calculation of first digital haplotype....."
	map=createRefDict(hla2fasta,"DQA1","DQB1","DRB1")
	readMapping(map,sam1)
	fourDigitString1=predictHLA2(sam1,medians,output1,medianflag)
	print "1st iteration done.\nNow removing reads that mapped to the three top-scoring groups ......."
	try:
		removeReads(runName,createRemoveList(runName,map))
	except IOError:
		print "Nothing to remove\n"
	
	#------2nd iteration------------------------------------------
	print "Second iterations starts .....\n Mapping ......"
	medians=[]
	iteration=2
	sam2=runName+"-iteration2.sam"
	newReadFile1=runName+"-2nditeration_1.fq"
	newReadFile2=runName+"-2nditeration_2.fq"
	mapping(sam2,runName,newReadFile1,newReadFile2,bowtiebuild,2,mapopt,trim3,gzip)
	medianfile=runName+".digitalhaplotype1"
	medians.append(float(linecache.getline(medianfile, 2).split('\t', 3)[2])/2.0)
	medians.append(float(linecache.getline(medianfile, 3).split('\t', 3)[2])/2.0)
	medians.append(float(linecache.getline(medianfile, 4).split('\t', 3)[2])/2.0)
	medianflag=True
	output2=runName+".digitalhaplotype2"
	finaloutput=runName+".HLAgenotype"
	#Calculation of second digital haplototype
	print "ClassII: calculation of second digital haplotype....."
	map=createRefDict(hla2fasta,"DQA1","DQB1","DRB1")
	readMapping(map,sam2)
	fourDigitString2=predictHLA2(sam2,medians,output2,medianflag)
	print "2nd iteration done."
	reportHLAgenotype(output1,output2,finaloutput,2)
	print "Calculation of locus-specific expression ..."
	try:
		expression("DQA1","DQB1","DRB1",400,421,421,map,runName,readFile1)
	except IOError:
		print "DQB1: 0 RPKM\nDRB1: 0 RPKM\nDQA1: 0 RPKM\n"

	#-----3rd iteration in case of at least one homozygous call-----------
	f1=fourDigitString1.split(",")
        f2=fourDigitString2.split(",")
        a1=f1[0]
        a1fourDigit=f1[1]
        if a1=="no":
                a1fourDigit_solutions=f1[2]
        else:
                a1fourDigit_solutions=int(f1[2])
        b1=f1[3]
        b1fourDigit=f1[4]
        if b1=="no":
                b1fourDigit_solutions=f1[5]
        else:
                b1fourDigit_solutions=int(f1[5])

        c1=f1[6]
        c1fourDigit=f1[7]
        if c1=="no":
                c1fourDigit_solutions=f1[8]
        else:
                c1fourDigit_solutions=int(f1[8])
        a2=f2[0]
        b2=f2[3]
        c2=f2[6]
	
	#try-catch block to prevent IO-Error in case of no expression
	try:
		if a2=="no" or b2=="no" or c2=="no":
                	removeReads(runName,createRemoveListFourDigits(runName,map,a1fourDigit,b1fourDigit,c1fourDigit))
	                sam3=runName+"-iteration3.sam"
        	        mapping(sam3,runName,newReadFile1,newReadFile2,bowtiebuild,2,mapopt,trim3,gzip)
                	readMapping(map,sam3)
	                output3=runName+".digitalhaplotype3"
        	        fourDigitString3=predictHLA2(sam3,medians,output3,medianflag)
	                f3=fourDigitString3.split(",")
        	        a3=f3[0]
                	b3=f3[3]
	                c3=f3[6]	

		#1st allele-------------------------------
	        #HLA-A
		Aallele1=a1fourDigit
	        if not a1fourDigit_solutions=="":
        	        if a1fourDigit_solutions>1:
                	        Aallele1+="'"
	                final_a=Aallele1+"\t"+str(getMaxP(sam1+".4digitsDQA1.solutions"))+"\t"
        	else:
                	final_a="no\tNA\t"

	        #HLA-B
		Ballele1=b1fourDigit
		if not b1fourDigit_solutions=="":
        	        if b1fourDigit_solutions>1:
                	        Ballele1+="'"
	                final_b=Ballele1+"\t"+str(getMaxP(sam1+".4digitsDQB1.solutions"))+"\t"
        	else:
                	final_b="no\tNA\t"

	        #HLA-C
		Callele1=c1fourDigit
	        if not c1fourDigit_solutions=="":
        	        if c1fourDigit_solutions>1:
                	        Callele1+="'"
	                final_c=Callele1+"\t"+str(getMaxP(sam1+".4digitsDRB1.solutions"))+"\t"
        	else:
                	final_c="no\tNA\t"
	        #2nd allele
        	if a2=="no":
                	if a3=="no" or a1fourDigit_solutions==1 or not a3.split(":")[0]==a1fourDigit.split(":"):
                        	if a1fourDigit_solutions=="":
                                	final_a+="no\tNA"
	                        elif a1fourDigit_solutions>1:
        	                        final_a+=Aallele1[0:-1]+"\t"+str(getP("A",finaloutput))
                	        else:
                        	        final_a+=Aallele1+"\t"+str(getP("A",finaloutput))
	                else:
        	                Aallele2=f3[1]
                	        if int(f3[2])>1:
                        	        Aallele2+="'"
	                        final_a+=Aallele2+"\t"+str(getMaxP(sam3+".4digitsDQA2.solutions"))
		else:
                	Aallele2=f2[1]
	                if int(f2[2])>1:
        	                Aallele2+="'"
                	final_a+=Aallele2+"\t"+str(getMaxP(sam2+".4digitsDQA2.solutions"))

	        if b2=="no":
        	        if b3=="no" or b1fourDigit_solutions==1 or not b3.split(":")[0]==b1fourDigit.split(":"):
                	        if b1fourDigit_solutions=="":
                        	        final_b+="no\tNA"
	                        elif b1fourDigit_solutions>1:
        	                        final_b+=Ballele1[0:-1]+"\t"+str(getP("B",finaloutput))
                	        else:
                        	        final_b+=Ballele1+"\t"+str(getP("B",finaloutput))
	                else:
        	                Ballele2=f3[4]
                	        if int(f3[5])>1:
                        	        Ballele2+="'"
	                        final_b+=Ballele2+"\t"+str(getMaxP(sam3+".4digitsDQB2.solutions"))
        	else:
                	Ballele2=f2[4]
	                if int(f2[5])>1:
        	                Ballele2+="'"
                	final_b+=Ballele2+"\t"+str(getMaxP(sam2+".4digitsDQB2.solutions"))

	        if c2=="no":
        	        if c3=="no" or c1fourDigit_solutions==1 or not c3.split(":")[0]==c1fourDigit.split(":"):
                	        if c1fourDigit_solutions=="":
                        	        final_c+="no\tNA"
	                        elif c1fourDigit_solutions>1:
        	                        final_c+=Callele1[0:-1]+"\t"+str(getP("C",finaloutput))
                	        else:
                        	        final_c+=Callele1+"\t"+str(getP("C",finaloutput))
	                else:
        	                Callele2=f3[7]
                	        if int(f3[8])>1:
                        	        Callele2+="'"
	                        final_c+=Callele2+"\t"+str(getMaxP(sam3+".4digitsDRB2.solutions"))
        	else:
                	Callele2=f2[7]
	                if int(f2[8])>1:
        	                Callele2+="'"
                	final_c+=Callele2+"\t"+str(getMaxP(sam2+".4digitsDRB2.solutions"))

		finaloutput4digit=runName+".HLAgenotype4digits"
        	#output 4-digit type (to stdout and to file)
		reportHLA4digitGenotype(final_a,final_b,final_c,finaloutput4digit,2)
	except IOError:
		#no expression
		print "no expression"
		finaloutput4digit=runName+".HLAgenotype4digits"
		reportHLA4digitGenotype("no\tNA\tno\tNA","no\tNA\tno\tNA","no\tNA\tno\tNA",finaloutput4digit,2)

#performs the bowtie mapping for the 2 iterations using the given parameters
def mapping(sam,runName,readFile1,readFile2,bowtiebuild,iteration,mapopt,trim3,gzip):
	if iteration==1:
		if gzip==0:
			mapping_cmd="(bowtie -3 "+trim3+" -S "+mapopt+" --al "+runName+".aligned "+bowtiebuild+" -1 "+readFile1+" -2 "+readFile2+" |  awk -F \'\\t\' '$3 != \"*\"{ print $0 }' > "+sam+") 2> "+runName+".bowtielog"
		else:
			mapping_cmd="(bowtie -3 "+trim3+" -S "+mapopt+" --al "+runName+".aligned "+bowtiebuild+" -1 <(zcat "+readFile1+") -2 <(zcat "+readFile2+") |  awk -F \'\\t\' '$3 != \"*\"{ print $0 }' > "+sam+") 2> "+runName+".bowtielog"
	if iteration==2:
		mapping_cmd="bowtie -3 "+trim3+" -S "+mapopt+" "+bowtiebuild+" -1 "+readFile1+" -2 "+readFile2+" "+sam
	#execute bowtie
	process_mapping=subprocess.Popen(['bash','-c',mapping_cmd],stdout=subprocess.PIPE)
	out,err=process_mapping.communicate()
	
	if iteration==1:
		#print alignment stats
		printcommand="cat "+runName+".bowtielog"
		printcommand_proc=subprocess.Popen(['bash','-c',printcommand],stdout=subprocess.PIPE)
		out,err=printcommand_proc.communicate()
		print out

#create dictionary "map", that contains all IMGT/HLA-allele names as keys and the allele name (e.g. A*02:01:01) as value
#dictionary "readcount" is initialized with 0 for each allele
#dictionary "readspergroup" is initialized with 0 for each group (2digit, e.g. A*01)
#dictionary "allelesPerLocus" stores the number of alleles per locus.
def createRefDict(hlafasta,locus1,locus2,locus3):
	map={}
	allelesPerLocus[locus1]=0
	allelesPerLocus[locus2]=0
	allelesPerLocus[locus3]=0
	handle=open(hlafasta,'r')
	for record in SeqIO.parse(handle, "fasta") :
		l=record.description.split(' ')
		hlapseudoname=l[0]
		hlaallele=l[1]
		map[hlapseudoname]=hlaallele
		readcount[hlaallele]=0
		readspergroup[hlaallele.split(":")[0]]=0
		allelesPerLocus[hlaallele.split('*')[0]]+=1
	handle.close()
	return map

#Open sam-file and count mappings for each allele 
def readMapping(map,sam):
	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			l=line.split('\t')
	                hla=l[2]
                	readcount[map[hla]]+=1

#predict HLA class II type
def predictHLA2(sam,medians,output,medianflag):

	readspergroupDQAlist=[]
	readspergroupDQBlist=[]
	readspergroupDRBlist=[]
	readspergroupDQA={}
	readspergroupDQB={}
	readspergroupDRB={}
	maxallelepergroup={}
	
	#for each allele, to which at least 1 read map, find the allele which has the most reads within a group (2-digit-level, e.g. DQA*02") and save
	#i) the key of this allele as ambassador for the group => maxallelepergroup holds for each group the top-scoring allele
	#ii) the number of reads mapping to the top-scoring allele => readspergroup
	for key in readcount:
		if readcount[key] > 0:
			group=key.split(":")[0]
			if readspergroup[group] <=readcount[key]:
				readspergroup[group]=readcount[key]
				maxallelepergroup[group]=key

	#consider all DQA-,DQB-, and DRB-groups seperately
	#readspergroup<DQA|DQB|DRB>list = list of all reads mapping to the top-scoring groups minus the decision threshold (which is 0 in the first iteration)
	#readspergroup<DQA|DQB|DRB> = contains the same entries as the list, but the entries are uniquely accessible via the group-key (e.g. DQA*01)
	for key in readspergroup:	
		if key[0:3]=='DQA':
			readspergroupDQAlist.append(readspergroup[key]-float(medians[0]))
			readspergroupDQA[key]=readspergroup[key]-float(medians[0])
		if key[0:3]=='DQB':
			readspergroupDQBlist.append(readspergroup[key]-float(medians[1]))
			readspergroupDQB[key]=readspergroup[key]-float(medians[1])
		if key[0:3]=='DRB':
			readspergroupDRBlist.append(readspergroup[key]-float(medians[2]))
			readspergroupDRB[key]=readspergroup[key]-float(medians[2])
				
	#compute the decision threshold (median) for homozygosity vs. heterozygosity for the second iteration 
	medianDQA=numpy.median(readspergroupDQAlist)
	medianDQB=numpy.median(readspergroupDQBlist)
	medianDRB=numpy.median(readspergroupDRBlist)
	a = ""
	b = ""
	c = ""
	
	#Determine top-scoring group of the whole locus (DQA,DQBB,DRB) and store it
	#maxkey<DQA,DQB,DRB> = group (e.g. DQA*01) with the most reads
	#It can be that, e.g. in cancer cells a whole locus is lost. For that reason it is checked if the 
	#number of mapping reads of the top-scoring group maxkey<DQA,DQB,DRB> is > 0, otherwise "no" ist reported for this locus
	if len(readspergroupDQAlist)>0:
		maxkeyDQA=max(readspergroupDQA,key=lambda a:readspergroupDQA.get(a))
		if readspergroupDQA[maxkeyDQA] > 0:
			maxDQA=maxallelepergroup[maxkeyDQA]
			a = max(readspergroupDQA,key=lambda a:readspergroupDQA.get(a))
		else:
			a = "no"
	else:
		a = "no"

	if len(readspergroupDQBlist)>0:
		maxkeyDQB=max(readspergroupDQB,key=lambda a:readspergroupDQB.get(a))
		if readspergroupDQB[maxkeyDQB] > 0:
			maxDQB=maxallelepergroup[maxkeyDQB]
			b = max(readspergroupDQB,key=lambda a:readspergroupDQB.get(a))
		else:
			b = "no"
	else:
		b = "no"

	if len(readspergroupDRBlist)>0:
		maxkeyDRB=max(readspergroupDRB,key=lambda a:readspergroupDRB.get(a))
		if readspergroupDRB[maxkeyDRB] > 0:
			maxDRB=maxallelepergroup[maxkeyDRB]
			c = max(readspergroupDRB,key=lambda a:readspergroupDRB.get(a))
		else:
			c = "no"
	else:
		c = "no"
		
	DQA4digit={}
	DQB4digit={}
	DRB4digit={}

	for key in readcount:
		if key.split(":")[0]==maxkeyDQA:
			DQA4digit[key]=readcount[key] 
		if key.split(":")[0]==maxkeyDQB:
			DQB4digit[key]=readcount[key]
		if key.split(":")[0]==maxkeyDRB:
			DRB4digit[key]=readcount[key]
	
	DQA4digit=sorted(DQA4digit.items(), key=itemgetter(1),reverse=True)
	DQB4digit=sorted(DQB4digit.items(), key=itemgetter(1),reverse=True)
	DRB4digit=sorted(DRB4digit.items(), key=itemgetter(1),reverse=True)
	
	if medianflag:
		DQA4digithandle=open(sam+".4digitsDQA2","w")
		DQB4digithandle=open(sam+".4digitsDQB2","w")
		DRB4digithandle=open(sam+".4digitsDRB2","w")
		for key in DQA4digit:
			DQA4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in DQB4digit:
			DQB4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in DRB4digit:
			DRB4digithandle.write(key[0]+": "+str(key[1])+"\n")
		DQA4digithandle.close()
		DQB4digithandle.close()
		DRB4digithandle.close()

	else:
		DQA4digithandle=open(sam+".4digitsDQA1","w")
		DQB4digithandle=open(sam+".4digitsDQB1","w")
		DRB4digithandle=open(sam+".4digitsDRB1","w")
		for key in DQA4digit:
			DQA4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in DQB4digit:
			DQB4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in DRB4digit:
			DRB4digithandle.write(key[0]+": "+str(key[1])+"\n")
		DQA4digithandle.close()
		DQB4digithandle.close()
		DRB4digithandle.close()



	readspergroupDQA=sorted(readspergroupDQA.items(), key=itemgetter(1),reverse=True)
	readspergroupDQB=sorted(readspergroupDQB.items(), key=itemgetter(1),reverse=True)
	readspergroupDRB=sorted(readspergroupDRB.items(), key=itemgetter(1),reverse=True)

	dqavec=""
	dqbvec=""
	drbvec=""
	if medianflag:
		#in the 2nd iteration: 
		#1.) DO NOT remove the top-scoring group from the set <dqa,dqb,drb>vec as this is more strict when calculating the probability of the top scoring group being an outlier
		#The strings <dqa,dqb,drb>vec are used by the R-script to calculate the probability of the top scoring group being an outlier
		for key in maxallelepergroup:
			if key[0:3]=="DQA":
				dqavec=dqavec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0:3]=="DQB":
				dqbvec=dqbvec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0:3]=="DRB":
				drbvec=drbvec+str(readcount[maxallelepergroup[key]])+","
		#2.) Add the decision thresholds to the sets <dqa,dqb,drb>vec, as this enables measuring the distance of the top-scoring group to this distance
		if a=="no":
			dqavec+=str(medians[0])
		if b=="no":
			dqbvec+=str(medians[1])
		if c=="no":
			drbvec+=str(medians[2])
		#3.) In case of, e.g. a loss of whole HLA-locus (DQA/B,DRB), <dqa,dqb,drb>vec only contain medians[<0|1|2>].
		#To avoid errors in the R-Script, set to 0
		if dqavec==str(medians[0]):
			dqavec = "0"
			dqacount="0"
		else:
			dqacount=str(readcount[maxallelepergroup[maxkeyDQA]])
			if dqavec==str(readcount[maxallelepergroup[maxkeyDQA]])+",":
				dqavec=str(readcount[maxallelepergroup[maxkeyDQA]])+","+str(readcount[maxallelepergroup[maxkeyDQA]])+","
		if dqbvec==str(medians[1]):
			dqbvec = "0"
			dqbcount="0"
		else:
			dqbcount=str(readcount[maxallelepergroup[maxkeyDQB]])
			if dqbvec==str(readcount[maxallelepergroup[maxkeyDQB]])+",":
				dqbvec=str(readcount[maxallelepergroup[maxkeyDQB]])+","+str(readcount[maxallelepergroup[maxkeyDQB]])+","
		if drbvec==str(medians[2]):
			drbvec = "0"
			drbcount="0"
		else:
			drbcount=str(readcount[maxallelepergroup[maxkeyDRB]])
			if drbvec==str(readcount[maxallelepergroup[maxkeyDRB]])+",":
				drbvec=str(readcount[maxallelepergroup[maxkeyDRB]])+","+str(readcount[maxallelepergroup[maxkeyDRB]])+","
	else:
		#in the 1st iteration: remove the top-scoring group from the set <dqa,dqb,drb>vec as this increases the certainty when calculating the probability of the top scoring group being an outlier
		for key in maxallelepergroup:
			if key[0:3]=="DQA":
				if key!=maxkeyDQA:
					dqavec=dqavec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0:3]=="DQB":
				if key!=maxkeyDQB:
					dqbvec=dqbvec+str(readcount[maxallelepergroup[key]])+","	
		for key in maxallelepergroup:
			if key[0:3]=="DRB":
				if key!=maxkeyDRB:
					drbvec=drbvec+str(readcount[maxallelepergroup[key]])+","
		#2.) DO NOT add the decision thresholds to the sets <dqa,dqb,drb>vec
		if dqavec=="":
			dqacount="0"
			dqavec="0"
		else:
			dqacount=str(readcount[maxallelepergroup[maxkeyDQA]])
			if dqavec==str(readcount[maxallelepergroup[maxkeyDQA]])+",":
				dqavec=str(readcount[maxallelepergroup[maxkeyDQA]])+","+str(readcount[maxallelepergroup[maxkeyDQA]])+","
			if len(dqavec.split(","))==2:
				dqavec=dqavec+"0,"

		if dqbvec=="":
			dqbcount="0"
			dqbvec="0"
		else:
			dqbcount=str(readcount[maxallelepergroup[maxkeyDQB]])
			if dqbvec==str(readcount[maxallelepergroup[maxkeyDQB]])+",":
				dqbvec=str(readcount[maxallelepergroup[maxkeyDQB]])+","+str(readcount[maxallelepergroup[maxkeyDQB]])+","
			if len(dqbvec.split(","))==2:
				dqbvec=dqbvec+"0,"

		if drbvec=="":
			drbcount="0"
			drbvec="0"
		else:
			drbcount=str(readcount[maxallelepergroup[maxkeyDRB]])
			if drbvec==str(readcount[maxallelepergroup[maxkeyDRB]])+",":
				drbvec=str(readcount[maxallelepergroup[maxkeyDRB]])+","+str(readcount[maxallelepergroup[maxkeyDRB]])+","
			if len(drbvec.split(","))==2:
				drbvec=drbvec+"0,"

	#call R-script "commmand.R" to calculate the confidence of the top-scoring allele
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command.R --args "+dqacount+" "+dqavec+" "+maxkeyDQA+" "+dqbcount+" "+dqbvec+" "+maxkeyDQB+" "+drbcount+" "+drbvec+" "+maxkeyDRB).read()
	parseOutput=routput.split("\n")

	entries = []
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))

	if a=="no":
		if entries[1]!="NA":
			entries[1]=str(1-float(entries[1]))
	if b=="no":
		if entries[3]!="NA":
			entries[3]=str(1-float(entries[3]))
	if c=="no":
		if entries[5]!="NA":
			entries[5]=str(1-float(entries[5]))
	
	if medianflag:
                fourDigitString=""
                if not a=="no":
                        fourDigitString+=a+","+fourdigits.determine4digits_main(sam+".4digitsDQA2",dqavec,runName,2)+","
                else:
                        fourDigitString+=a+",,,"
                if not b=="no":
                        fourDigitString+=b+","+fourdigits.determine4digits_main(sam+".4digitsDQB2",dqbvec,runName,2)+","
                else:
                        fourDigitString+=b+",,,"
                if not c=="no":
                        fourDigitString+=c+","+fourdigits.determine4digits_main(sam+".4digitsDRB2",drbvec,runName,2)+","
                else:
                        fourDigitString+=c+",,,"
        else:
                fourDigitString=""
                if not a=="no":
                        fourDigitString+=a+","+fourdigits.determine4digits_main(sam+".4digitsDQA1",dqavec,runName,2)+","
                else:
                        fourDigitString+=a+",,,"
                if not b=="no":
                        fourDigitString+=b+","+fourdigits.determine4digits_main(sam+".4digitsDQB1",dqbvec,runName,2)+","
                else:
                        fourDigitString+=b+",,,"
                if not c=="no":
                        fourDigitString+=c+","+fourdigits.determine4digits_main(sam+".4digitsDRB1",drbvec,runName,2)+","
                else:
                        fourDigitString+=c+",,,"
	pred2File(entries,medianDQA,medianDQB,medianDRB,output,a,b,c,2)
	return fourDigitString

#predict HLA class I type 
def predictHLA1(sam,medians,output,medianflag): 
	
	readspergroupAlist=[]
	readspergroupBlist=[]
	readspergroupClist=[]
	readspergroupA={}
	readspergroupB={}
	readspergroupC={}
	maxallelepergroup={}
	
	#for each allele, to which at least 1 read map, find the allele which has the most reads within a group (2-digit-level, e.g. A*02") and save
	#i) the key of this allele as ambassador for the group => maxallelepergroup holds for each group the top-scoring allele
	#ii) the number of reads mapping to the top-scoring allele => readspergroup
	for key in readcount:
		if readcount[key] > 0:
			group=key.split(":")[0]
			if readspergroup[group] <=readcount[key]:
				readspergroup[group]=readcount[key]
				maxallelepergroup[group]=key

	
	#consider all A-,B-, and C-groups seperately
	#readspergroup<A|B|C>list = list of all reads mapping to the top-scoring groups minus the decision threshold (which is 0 in the first iteration)
	#readspergroup<A|B|C> = contains the same entries as the list, but the entries are uniquely accessible via the group-key (e.g. B*27)
	for key in readspergroup:	
		if key[0]=='A':
			readspergroupAlist.append(readspergroup[key]-float(medians[0]))
			readspergroupA[key]=readspergroup[key]-float(medians[0])
		if key[0]=='B':
			readspergroupBlist.append(readspergroup[key]-float(medians[1]))
			readspergroupB[key]=readspergroup[key]-float(medians[1])
		if key[0]=='C':
			readspergroupClist.append(readspergroup[key]-float(medians[2]))
			readspergroupC[key]=readspergroup[key]-float(medians[2])
	
	#compute the decision threshold (median) for homozygosity vs. heterozygosity for the second iteration 
	medianA=numpy.median(readspergroupAlist)
	medianB=numpy.median(readspergroupBlist)
	medianC=numpy.median(readspergroupClist)
	a = ""
	b = ""
	c = ""
	
	#Determine top-scoring group of the whole locus (A,B,C) and store it
	#maxkey<A,B,C> = group (e.g. A*02) with the most reads
	#It can be that, e.g. in cancer cells A whole locus is lost. For that reason it is checked if the 
	#number of mapping reads of the top-scoring group maxkey<A,B,C> is > 0, otherwise "no" ist reported for this locus
	if len(readspergroupAlist)>0:
		maxkeyA=max(readspergroupA,key=lambda a:readspergroupA.get(a))
		if readspergroupA[maxkeyA] > 0:
			maxA=maxallelepergroup[maxkeyA]
			a = max(readspergroupA,key=lambda a:readspergroupA.get(a))
		else:
			a = "no"
	else:
		a = "no"

	if len(readspergroupBlist)>0:
		maxkeyB=max(readspergroupB,key=lambda a:readspergroupB.get(a))
		if readspergroupB[maxkeyB] > 0:
			maxB=maxallelepergroup[maxkeyB]
			b = max(readspergroupB,key=lambda a:readspergroupB.get(a))
		else:
			b = "no"
	else:
		b = "no"

	if len(readspergroupClist)>0:
		maxkeyC=max(readspergroupC,key=lambda a:readspergroupC.get(a))
		if readspergroupC[maxkeyC] > 0:
			maxC=maxallelepergroup[maxkeyC]
			c = max(readspergroupC,key=lambda a:readspergroupC.get(a))
		else:
			c = "no"
	else:
		c = "no"
	
	A4digit={}
	B4digit={}
	C4digit={}

	for key in readcount:
		if key.split(":")[0]==maxkeyA:
			A4digit[key]=readcount[key] 
		if key.split(":")[0]==maxkeyB:
			B4digit[key]=readcount[key]
		if key.split(":")[0]==maxkeyC:
			C4digit[key]=readcount[key]
	
	A4digit=sorted(A4digit.items(), key=itemgetter(1),reverse=True)
	B4digit=sorted(B4digit.items(), key=itemgetter(1),reverse=True)
	C4digit=sorted(C4digit.items(), key=itemgetter(1),reverse=True)
	
	if medianflag:
		A4digithandle=open(sam+".4digitsA2","w")
		B4digithandle=open(sam+".4digitsB2","w")
		C4digithandle=open(sam+".4digitsC2","w")
		for key in A4digit:
			A4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in B4digit:
			B4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in C4digit:
			C4digithandle.write(key[0]+": "+str(key[1])+"\n")
		A4digithandle.close()
		B4digithandle.close()
		C4digithandle.close()

	else:
		A4digithandle=open(sam+".4digitsA1","w")
		B4digithandle=open(sam+".4digitsB1","w")
		C4digithandle=open(sam+".4digitsC1","w")
		for key in A4digit:
			A4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in B4digit:
			B4digithandle.write(key[0]+": "+str(key[1])+"\n")
		for key in C4digit:
			C4digithandle.write(key[0]+": "+str(key[1])+"\n")
		A4digithandle.close()
		B4digithandle.close()
		C4digithandle.close()

	
	
	readspergroupA=sorted(readspergroupA.items(), key=itemgetter(1),reverse=True)
	readspergroupB=sorted(readspergroupB.items(), key=itemgetter(1),reverse=True)
	readspergroupC=sorted(readspergroupC.items(), key=itemgetter(1),reverse=True)

	avec=""
	bvec=""
	cvec=""

	if medianflag:
		
		#in the 2nd iteration: 
		#1.) DO NOT remove the top-scoring group from the set <a,b,c>vec as this is more strict when calculating the probability of the top scoring group being an outlier
		#The strings <a,b,c>vec are used by the R-script to calculate the probability of the top scoring group being an outlier
		for key in maxallelepergroup:
			if key[0]=="A":
				avec=avec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0]=="B":
				bvec=bvec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0]=="C":
				cvec=cvec+str(readcount[maxallelepergroup[key]])+","
		#2.) Add the decision thresholds to the sets <a,b,c>vec, as this enables measuring the distance of the top-scoring group to this distance
		if a=="no":
			avec+=str(medians[0])
		if b=="no":		
			bvec+=str(medians[1])
		if c=="no":
			cvec+=str(medians[2])
		#3.) In case of, e.g. a loss of whole HLA-locus (A,B,C), <a,b,c>vec only contain median[<0|1|2>].
		#To avoid errors in the R-Script, set to 0
		if avec=="" or avec==medians[0]:
			avec = "0"
			acount="0"
		else:
			acount=str(readcount[maxallelepergroup[maxkeyA]])

		if bvec=="" or bvec==medians[1]:
			bvec = "0"
			bcount="0"
		else:
			bcount=str(readcount[maxallelepergroup[maxkeyB]])

		if cvec=="" or cvec==medians[2]:
			cvec = "0"
			ccount="0"
		else:
			ccount=str(readcount[maxallelepergroup[maxkeyC]])
	else:
		#in the 1st iteration: remove the top-scoring group from the set <a,b,c>vec as this increases the certainty when calculating the probability of the top scoring group being an outlier
		for key in maxallelepergroup:
			if key[0]=="A":
				if key!=maxkeyA:
					avec=avec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0]=="B":
				if key!=maxkeyB:
					bvec=bvec+str(readcount[maxallelepergroup[key]])+","
		for key in maxallelepergroup:
			if key[0]=="C":
				if key!=maxkeyC:
					cvec=cvec+str(readcount[maxallelepergroup[key]])+","
		#2.) DO NOT add the decision thresholds to the sets <a,b,c>vec
		if avec=="":
			acount="0"
			avec="0"
		else:
			acount=str(readcount[maxallelepergroup[maxkeyA]])
		if bvec=="":
			bcount="0"
			bvec="0"
		else:
			bcount=str(readcount[maxallelepergroup[maxkeyB]])
		if cvec=="":
			ccount="0"
			cvec="0"
		else:
			ccount=str(readcount[maxallelepergroup[maxkeyC]])

	#call R-script "commmand.R" to calculate the confidence of the top-scoring allele
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command.R --args "+acount+" "+avec+" "+maxkeyA+" "+bcount+" "+bvec+" "+maxkeyB+" "+ccount+" "+cvec+" "+maxkeyC).read()
	parseOutput=routput.split("\n")

	entries = []
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))

	if medianflag:
		fourDigitString=""
		if not a=="no":
			fourDigitString+=a+","+fourdigits.determine4digits_main(sam+".4digitsA2",avec,runName,1)+","
		else:
			fourDigitString+=a+",,,"
		if not b=="no":
			fourDigitString+=b+","+fourdigits.determine4digits_main(sam+".4digitsB2",bvec,runName,1)+","
		else:
			fourDigitString+=b+",,,"
		if not c=="no":
			fourDigitString+=c+","+fourdigits.determine4digits_main(sam+".4digitsC2",cvec,runName,1)+","
		else:
			fourDigitString+=c+",,,"
	else:
		fourDigitString=""
		if not a=="no":
			fourDigitString+=a+","+fourdigits.determine4digits_main(sam+".4digitsA1",avec,runName,1)+","
		else:
			fourDigitString+=a+",,,"
		if not b=="no":
			fourDigitString+=b+","+fourdigits.determine4digits_main(sam+".4digitsB1",bvec,runName,1)+","
		else:
			fourDigitString+=b+",,,"
		if not c=="no":
			fourDigitString+=c+","+fourdigits.determine4digits_main(sam+".4digitsC1",cvec,runName,1)+","
		else:
			fourDigitString+=c+",,,"
			
	if a=="no":
		if entries[1]!="NA":
			entries[1]=str(1-float(entries[1]))
	if b=="no":
		if entries[3]!="NA":
			entries[3]=str(1-float(entries[3]))
	if c=="no":
		if entries[5]!="NA":
			entries[5]=str(1-float(entries[5]))
	
	pred2File(entries,medianA,medianB,medianC,output,a,b,c,1)
	return fourDigitString
	
#write digital haplotype into file
def pred2File(entries,medianA,medianB,medianC,output,a,b,c,HLAclass):
	if HLAclass==1:
                locus1="A"
                locus2="B"
                locus3="C"
        else:
                locus1="DQA"
                locus2="DQB"
                locus3="DRB"
	out = open(output, 'w')
	out.write("HLA\tHLA1\tmedian-Value\talternative\tp-Value\n")
	out.write(locus1+"\t"+a+"\t")
	# write the median-Value for A
	out.write(str(medianA)+"\t")
	out.write(entries[0]+"\t"+entries[1]+"\n")

	out.write(locus2+"\t"+b+"\t")
	# write the median-Value for B
	out.write(str(medianB)+"\t")
	out.write(entries[2]+"\t"+entries[3]+"\n")

	out.write(locus3+"\t"+c+"\t")
	# write the median-Value for C
	out.write(str(medianC)+"\t")
	out.write(entries[4]+"\t"+entries[5]+"\n")
	out.close()
	print "The digital haplotype is written into "+output

#open mapping file and all read ids to the list "removeList", which map to one of the three groups in "alleles"
def createRemoveList(runName,map):
	removeList={}
	alleles = []
	sam=runName+"-iteration1.sam"
	alleles_in=runName+".digitalhaplotype1"
	alleles.append(linecache.getline(alleles_in, 2).split('\t', 2)[1])
	alleles.append(linecache.getline(alleles_in, 3).split('\t', 2)[1])
	alleles.append(linecache.getline(alleles_in, 4).split('\t', 2)[1])

	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0]
			hlapseudoname = line.split("\t")[2]
			if map[hlapseudoname].split(':')[0] in alleles:
				removeList[illuminaid]=1
	samhandle.close()
	return removeList

#open mapping file and all read ids to the list "removeList", which map to one of the three four-digit alleles in "alleles"
def createRemoveListFourDigits(runName,map,a1,b1,c1):
	removeList={}
	sam=runName+"-iteration1.sam"

	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0]
			hlapseudoname = line.split("\t")[2]
			fourdigits=map[hlapseudoname].split(':')[0]+":"+map[hlapseudoname].split(':')[1]
			if fourdigits==a1 or fourdigits==b1 or fourdigits==c1:
				removeList[illuminaid]=1
	samhandle.close()
	return removeList	
	
#Remove reads that mapped to the three top-scoring alleles and write the remaining reads into two new read files
def removeReads(runName,removeList):
	aligned1=runName+"_1.aligned"
	aligned2=runName+"_2.aligned"
	newReadFile1=runName+"-2nditeration_1.fq"
	newReadFile2=runName+"-2nditeration_2.fq"
	#r1 and r2, which are the input of bowtie in the 2nd iteration
	r1=open(newReadFile1,"w")
	r2=open(newReadFile2,"w")
	#open the 2 files, that contain the reads that mapped in the 1st iteration
	aligned_handle1=open(aligned1,"r")
	aligned_handle2=open(aligned2,"r")
	
	#One read entry consists of 4 lines: header, seq, "+", qualities.
	for record in SeqIO.parse(aligned_handle1, "fastq"):
		illuminaid=record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
		if not illuminaid in removeList:
			SeqIO.write(record, r1, "fastq")
	
	for record in SeqIO.parse(aligned_handle2, "fastq"):
		illuminaid=record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
		if not illuminaid in removeList:
			SeqIO.write(record, r2, "fastq")

#write the final prediction (both digital haplotypes) to file and stdout
def reportHLAgenotype(output1,output2,finaloutput,HLAclass):
	filehandle1=open(output1,'r').readlines()[1:4]
	filehandle2=open(output2,'r').readlines()[1:4]
	outfile=open(finaloutput, 'w')
	outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
	
	if HLAclass==1:
		a="A"
		b="B"
		c="C"
	else:
		a="DQA"
		b="DQB"
		c="DRB"
	
	for i in range(len(filehandle1)):
		filehandle1[i]=filehandle1[i][0:-1]
	for i in range(len(filehandle2)):
		filehandle2[i]=filehandle2[i][0:-1]

	a1 = filehandle1[0].split('\t',2)[1]
	a1score = filehandle1[0].split('\t')[4]
	b1 = filehandle1[1].split('\t',2)[1]
	b1score = filehandle1[1].split('\t')[4]
	c1 = filehandle1[2].split('\t',2)[1]
	c1score = filehandle1[2].split('\t')[4]

	a2 = filehandle2[0].split('\t',2)[1]
	if a2 == "no":
		a2 = "hoz("+filehandle2[0].split('\t')[3]+")"
	a2score = filehandle2[0].split('\t')[4]
	b2 = filehandle2[1].split('\t',2)[1]
	if b2 == "no":
		b2 = "hoz("+filehandle2[1].split('\t')[3]+")"
	b2score = filehandle2[1].split('\t')[4]
	c2 = filehandle2[2].split('\t', 2)[1]
	if c2 == "no":
		c2 = "hoz("+filehandle2[2].split('\t')[3]+")"
	c2score = filehandle2[2].split('\t')[4]
	#write complete HLA genotype to file
	outfile.write(a+"\t"+a1+"\t"+a1score+"\t"+a2+"\t"+a2score+"\n")
	outfile.write(b+"\t"+b1+"\t"+b1score+"\t"+b2+"\t"+b2score+"\n")
	outfile.write(c+"\t"+c1+"\t"+c1score+"\t"+c2+"\t"+c2score+"\n")
	outfile.close()
	#.. and print it to STDOUT
	print "-----------2 digit typing results-------------"
	print "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence"
	print a+"\t"+a1+"\t"+a1score+"\t"+a2+"\t"+a2score
	print b+"\t"+b1+"\t"+b1score+"\t"+b2+"\t"+b2score
	print c+"\t"+c1+"\t"+c1score+"\t"+c2+"\t"+c2score

def reportHLA4digitGenotype(final_a,final_b,final_c,finaloutput4digit,HLAclass):
	#write complete HLA genotype at four-digit-level to file
	outfile=open(finaloutput4digit, 'w')
	outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
	if HLAclass==1:
		a="A"
		b="B"
		c="C"
	else:
		a="DQA"
		b="DQB"
		c="DRB"
		
	outfile.write(a+"\t"+final_a.replace(",","\t")+"\n")
	outfile.write(b+"\t"+final_b.replace(",","\t")+"\n")
	outfile.write(c+"\t"+final_c.replace(",","\t")+"\n")
	outfile.close()
	#.. and print it to STDOUT
	print "-----------4 digit typing results-------------"
	print "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence"
	print "!"+a+"\t"+final_a.replace(",","\t")
	print "!"+b+"\t"+final_b.replace(",","\t")
	print "!"+c+"\t"+final_c.replace(",","\t")
	
#calculate locus-specific expression
def expression(locus1,locus2,locus3,length1,length2,length3,map,runName,readFile1):
	outfile=open(runName+".expression", 'w')
	aligned1=runName+"_1.aligned"
	sam=runName+"-iteration1.sam"
	logfile=runName+".bowtielog"
	print logfile
	totalreads=float(linecache.getline(logfile,1).split(':')[1])
	alleles_in=runName+".digitalhaplotype1"
	alleles=[]
	alleles.append(linecache.getline(alleles_in, 2).split('\t', 2)[1])
	alleles.append(linecache.getline(alleles_in, 2).split('\t')[3])
	alleles.append(linecache.getline(alleles_in, 3).split('\t', 2)[1])
	alleles.append(linecache.getline(alleles_in, 3).split('\t')[3])
	alleles.append(linecache.getline(alleles_in, 4).split('\t', 2)[1])
	alleles.append(linecache.getline(alleles_in, 4).split('\t')[3])
	
	#create read dictionary
	reads={}
	aligned_handle1=open(aligned1,"r")
	for record in SeqIO.parse(aligned_handle1, "fastq"):
		illuminaid=record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
		reads[illuminaid]={}
		reads[illuminaid][locus1]=0
		reads[illuminaid][locus2]=0
		reads[illuminaid][locus3]=0
	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0].split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
			hlapseudoname = line.split("\t")[2]
			if map[hlapseudoname].split(':')[0] in alleles:
				reads[illuminaid][map[hlapseudoname].split('*')[0]]+=1
	count={}
	count[locus1]=0
	count[locus2]=0
	count[locus3]=0
	for key in reads:
		n=0
		for locus in reads[key]:
			if reads[key][locus] > 0:
				n+=1
		for locus in reads[key]:
			if reads[key][locus] > 0:
				count[locus]+=float(1.0/float(n))
	
	#Calculate RPKM and print expression values for each locus to stdout
	for locus in count:
		if locus==locus1:
			print locus+": "+str(round(float((1000.0/length1))*float((1000000.0/totalreads))*count[locus],2))+" RPKM"
			outfile.write(locus+": "+str(round(float((1000.0/length1))*float((1000000.0/totalreads))*count[locus],2))+" RPKM\n")
		if locus==locus2:
			print locus+": "+str(round(float((1000.0/length2))*float((1000000.0/totalreads))*count[locus],2))+" RPKM"
			outfile.write(locus+": "+str(round(float((1000.0/length2))*float((1000000.0/totalreads))*count[locus],2))+" RPKM\n")
		if locus==locus3:
			print locus+": "+str(round(float((1000.0/length3))*float((1000000.0/totalreads))*count[locus],2))+" RPKM"
			outfile.write(locus+": "+str(round(float((1000.0/length3))*float((1000000.0/totalreads))*count[locus],2))+" RPKM\n")
        outfile.close()

#In case of ambiguous typings, the allele(s) with the best p-value (which is actually the smallest one, so the name of the function is misleading - sorry) is reported and thus the min(p) is returned
def getMaxP(file):
	p=[]
	for line in open(file,"r"):
		if not line[0]=="#":
			if line.split("\t")[1][0:-1]=="NA":
				return "NA"
			p.append(float(line.split("\t")[1][0:-1]))
	try:
		return min(p)
	except ValueError:
		return "NA"

#Return the p value of a prediction, which is stored in the intermediate textfile
def getP(locus,finaloutput):
	if locus=="A" or locus=="DQA":
		p=linecache.getline(finaloutput, 2).split('\t')[4]
	if locus=="B" or locus=="DQB":
		p=linecache.getline(finaloutput, 3).split('\t')[4]
	if locus=="C" or locus=="DRB":
		p=linecache.getline(finaloutput, 4).split('\t')[4]
	return p[0:-1]
	
if __name__ == '__main__':
	parser = OptionParser(usage="usage: %prog -1 readFile1 -2 readFile2 -r runName [-p <int>] [-3 <int>]", version="%prog 2.2")
	parser.add_option("-1",
			action="store", 
			dest="readFile1",
			help="File name of #1 mates (uncompressed or gzipped fastq)")
	parser.add_option("-2",
			action="store", 
			dest="readFile2",
			help="File name of #2 mates (uncompressed or gzipped fastq)")
	parser.add_option("-r", "--runName",
			action="store", 
			dest="runName",
			help="Name of this HLA typing run. Wil be used throughout this process as part of the name of the newly created files.")
	parser.add_option("-p", "--threads",
                        action="store",
                        dest="threads",
                        default="6",
                        help="Bowtie option: Launch <int> parallel search threads. Default (seq2HLA): 6")
	parser.add_option("-3", "--trim3",
			action="store",
			dest="trim3",
			default="0",
			help="Bowtie option: -3 <int> trims <int> bases from the low quality 3' end of each read. Default: 0")


	(options, args) = parser.parse_args()
	if not options.readFile1:   
		parser.error('File name #1 pair not given.')
	if not options.readFile2:   
		parser.error('File name #2 pair not given.')
	if not options.runName:   
		parser.error('Run name not given.')
	readFile1=options.readFile1
	readFile2=options.readFile2
	runName=options.runName
	bowtiebuildClassI=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/ClassIWithoutNQex2-3.plus75"
	bowtiebuildClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75"
	fastaClassI=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/ClassIWithoutNQex2-3.plus75.fasta"
	fastaClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75.fasta"
	trim3=str(options.trim3)
	threads=str(options.threads)
	main(runName,readFile1,readFile2,fastaClassI,fastaClassII,bowtiebuildClassI,bowtiebuildClassII,trim3,threads)


