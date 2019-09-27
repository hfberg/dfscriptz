#!/bin/bash

# 16-01-2019 written by Hanna F. Berg for Drug Farm
# 20-08-2019 updated by Hanna F. Berg to increase the user friendliness

# A bash script to conveniently filter a fastq file for linker 1, linker 2 and poly A


progname=dfpipe
tools_path=`pwd`
files_to_delete=

function help () {
    cat >&2 <<EOF
USAGE: Start with using cd to go to the folder where your tools BBmap and Drop-seq are. The script is based on that the working directory is your tools folder. 

The $progname is pre-processing raw illumina fastq reads to prepare for the drop-seq pipeline. Filters on two linkers and the poly-A tag. Then extracts biolocigal information from R1 to R2 and merges into a bam file to be put into the drop-seq pipeline. $progname assumes BBMap and dropseq to be in the same folder as $progname.

-1 <Read 1>     : The full path to the gzipped fastq file of the paired Illumina read 1.  Required.
-2 <Read 2> 	: The full path to the gzipped fastq file of the paired Illumina read 2.  Required.
-n <name>	: Desired sample name. Required
-o <outputdir>  : Folder in which to write output. Required

EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    help
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

while getopts "1:2:n:o:h" OPTION
do
	case $OPTION in
		1)	in=$OPTARG;;
			#echo "You set flag -1";;

		2)	in2=$OPTARG;;
			#echo "The value of -2 is $OPTARG";;
			
		n)	samp_name=$OPTARG;;
			#echo "You set flag -n";;

		o)	output_dir=$OPTARG;;
			#echo "You set flag -o";;

		h ) 	help
          		exit 1;;

    		\? ) 	echo "You defined a flag which is not a valid option for $progname. write -h to get help on available options"
         		exit 1;;
	esac
done

shift $((OPTIND-1))


check_set "$in" "Sample Read 1" "-1"
check_set "$in2" "Sample Read 2" "-2"
check_set "$samp_name" "A sample name is required"  "-n"
check_set "$output_dir" "A path to an output directory has not been provided." "-o"

out=$output_dir/$samp_name
mkdir -p  $out

#filter1

`pwd`/bbmap/bbduk2.sh in1=$in in2=$in2 outm=$out/$samp_name"_R1_filt1" outm2=$out/$samp_name"_R2_filt1" fliteral=CGACTCACTACAGGG k=15 skipr2=t hdist=3 -Xmx58g

#filter 2

`pwd`/bbmap/bbduk2.sh in1=$out/$samp_name"_R1_filt1" in2=$out/$samp_name"_R2_filt1" outm=$out/$samp_name"_R1_filt2" outm2=$out/$samp_name"_R2_filt2" fliteral=TCGGTGACACGATCG k=15 skipr2=t hdist=3 -Xmx58g 

#remove files from filter 1.
rm $out/$samp_name"_R1_filt1"
rm $out/$samp_name"_R2_filt1"

echo " "
echo "****** Removed ******"
echo " "

echo "removed files" $out/$samp_name"_R1_filt1" "and" $out/$samp_name"_R2_filt1" 

echo " "
echo " "
echo " "

# poly A filter

`pwd`/bbmap/bbduk2.sh in1=$out/$samp_name"_R1_filt2" in2=$out/$samp_name"_R2_filt2"  outm=$out/$samp_name"_R1_polyA" outm2=$out/$samp_name"_R2_polyA" fliteral=TTTTTTTTTTTT k=12 skipr2=t hdist=3 -Xmx58g

#remove files from filter 2.
rm $out/$samp_name"_R1_filt2"
rm $out/$samp_name"_R2_filt2"

echo " "
echo "****** Removed ******"
echo " "

echo "removed files" $out/$samp_name"_R1_filt2" "and" $out/$samp_name"_R2_filt2" 

echo " "
echo " "
echo " "

# Fastq to Sam

java -Xmx58g -jar `pwd`/drop-seq-2.3.0/3rdParty/picard/picard.jar FastqToSam F1=$out/$samp_name"_R1_polyA" F2=$out/$samp_name"_R2_polyA" O=$out/$samp_name"_FastqToSam" SAMPLE_NAME=$samp_name

#remove file from poly A filter.
rm $out/$samp_name"_R1_polyA"
rm $out/$samp_name"_R2_polyA"

echo " "
echo "****** Removed ******"
echo " "

echo "removed files" $out/$samp_name"_R1_polyA" "and" $out/$samp_name"_R2_polyA" 

echo " "
echo " "
echo " "

# Explanation of what is done
# in=<file>           Main input. in=stdin.fq will pipe from stdin.
# in2=<file>          Input for 2nd read of pairs in a different file.

# outm=<file>         (outmatch) Write reads here that contain kmers matching the database.
# outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a different file.

# fliteral=<seq,seq>  Comma-delimited list of literal sequences for filtering.

# k=27                Kmer length used for finding contaminants.  Contaminants shorter than k will not be found. k must be at least 1.

# hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only). Memory use is proportional to (3*K)^hdist.

# skipr2=f            Don't do kmer-based operations on read 2.





 
