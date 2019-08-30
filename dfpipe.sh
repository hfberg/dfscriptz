#!/bin/bash

# 16-01-2019 written by Hanna F. Berg fro Drug Farm
# A bash script to conveniently filter a fastq file for linker 1, linker 2 and poly A

in=$1 #Path to extracted sample R1
in2=$2 #Path to extracted sample R2
samp_name=$3 #subfolder with samp name in preprocessing folder, goes for both R1 and R2


#filter1

/home/drugfarm/Tools/bbmap/bbduk2.sh in1=$in in2=$in2 outm=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_filt1" outm2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_filt1" fliteral=CGACTCACTACAGGG k=15 skipr2=t hdist=3 -Xmx58g

#filter 2

/home/drugfarm/Tools/bbmap/bbduk2.sh in1=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_filt1" in2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_filt1" outm=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_filt2" outm2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_filt2" fliteral=TCGGTGACACGATCG k=15 skipr2=t hdist=3 -Xmx58g 

# poly A filter

/home/drugfarm/Tools/bbmap/bbduk2.sh in1=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_filt2" in2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_filt2"  outm=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_polyA" outm2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_polyA" fliteral=TTTTTTTTTTTT k=12 skipr2=t hdist=3 -Xmx58g


# Fastq to Sam
java -Xmx58g -jar /home/drugfarm/Tools/drop-seq-1.12/3rdParty/picard/picard.jar FastqToSam F1=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R1_polyA" F2=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_R2_polyA" O=/home/drugfarm/proj/preprocess/$samp_name/$samp_name"_FastqToSam" SAMPLE_NAME=$samp_name

# Explanation of what is done
# in=<file>           Main input. in=stdin.fq will pipe from stdin.
# in2=<file>          Input for 2nd read of pairs in a different file.

# outm=<file>         (outmatch) Write reads here that contain kmers matching the database.
# outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a different file.

# fliteral=<seq,seq>  Comma-delimited list of literal sequences for filtering.

# k=27                Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1.

# hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only). Memory use is proportional to (3*K)^hdist.

# skipr2=f            Don't do kmer-based operations on read 2.





 
