#!/bin/bash

######################################################################
# paMATRAX*.sh: A shell script for sequentially execute multiple 
# sequence alignment (mafft), trimming (trimAl) and tree-estimation 
#(raxml or fastree).
# If you use this script you should cite the actual
# programs in the dependencies list: 
# 
# The following programs are required and should be referenced
# in the  $PATH.
#
#  *gnu-parallel
#  *raxmlHPC 
#  *trimal 
#  *FastTree  
#  *mafft 
#  *Al2Phylo.py
######################################################################

#Program specific commands. User should modify this accordingly.
mafft_cmd="mafft --anysymbol --auto --quiet --thread 4"
trimal_cmd="ln -fs " #do not trim, just link the original files for next step
raxml_cmd="/beegfs/group_dv/software/source/standard-RAxML-8.2.9/raxmlHPC-PTHREADS-AVX2 -T 1 -f a -p 767 -x 23333 -#100 -m PROTGAMMAJTTX"
fasttree_cmd="fasttree -nt"
Al2Phylo_cmd="Al2Phylo.py -m 50 -p 0.30 -t 2"

export mafft_cmd
export trimal_cmd
export raxml_cmd
export Al2Phylo_cmd
export fasttree_cmd

#Initialize variables
EXT="fasta"
AFLAG=1
TFLAG=1
SFLAG=1
CFLAG=0
TinEXT='.fa'
TREE_BUILDER=0

usage() {
cat <<EOF

usage: $0 <options>

This script runs the phylogetic pipeline Align -> Trim-> Sanitize-> Tree
on all sequences in the cwd. It calls  gnu-parallel, mafft, trimAL, RAxML or
FastTree. If an output file is found it would be skipped.
Please cite the appropriate programs used on each step.

-h  |  Print this help
-e  |  Extension of the unaligned sequences (default: fasta)
-a  |  Stop after alignment
-t  |  Stop after trimming
-s  |  Stop after sanitation
-c  |  Sanitize trimmed alignments with Al2Phylo.py
-f  |  Use FastTree for building trees (default raxml).

These are the default parameters for each program. Modify accordingly:

    mafft:    $mafft_cmd
   trimAl:    $trimal_cmd
 Al2Phylo:    $Al2Phylo_cmd
    raxml:    $raxml_cmd
 fasttree:    $fasttree_cmd

EOF
}

while getopts "he:atscf" opt; do

    case "$opt" in
	h)
	    usage
	    exit 0
	    ;;
	e)
	    EXT=$OPTARG
	    ;;
	a) 
	    AFLAG=0
	    ;;
	t)
	    TFLAG=0
	    ;;
	s) SFLAg=0
	    ;;
	c) CFLAG=1
	    ;;
	f)
	    TREE_BUILDER=1
	    ;;
	?)
	    usage >&2
	    exit 1
	    ;;
    esac

done

main () {
    printf "\nStarting MSA"
    parallel --env mafft_cmd -j+0 'if [ ! -s {.}.al  ]; then $mafft_cmd {} > {.}.al 2>>mafft.log; fi' ::: *.$EXT;
    printf "\nAll alignemnets are completed. Alignments files writen with extension with extension .al"
    if [ $AFLAG -eq 0 ]
    then
	printf "\nPipeline stopped after alignement."
	exit 0
    else
	printf "\n\nStarting trimming."
	parallel --env trimal_cmd  -j+0 'if [ ! -s {.}.fa  ]; then $trimal_cmd {} {.}.fa; fi' ::: *.al;     
	printf "\nAll alignments were trimmed. Trimmed alignments written with extension .fa"
	if [ $TFLAG -eq 0 ]
	then
	    printf "\nPipeline stoped after trimming."
	    exit 0
	else
	
	    if [ $CFLAG -eq 1 ]
		
	    then
		printf "\n\nStarting cleaning"
		parallel --env Al2Phylo_cmd -j+0 'if [[ ! -s {.}_clean.fa && ! -e {= s:_clean\.fa::; =}_clean.fa ]]; then $Al2Phylo_cmd -in {} >> Al2Phylo.log; fi' ::: *.fa; 
		TinEXT='_clean.fa'
		printf '\nAll alignments were cleaned. Cleaned alignments end with _clean.fa'
	    fi    
	    
	    if [ $SFLAG -eq 0 ]
	    then
		printf "\nPipeline stopped after sanitation."
		exit 0
	    else		
		if [ $TREE_BUILDER -eq 0 ]
		then
		    printf "\n\nStarting tree estimation with raxml"
		    parallel --env raxml_cmd -j+0 'if [ ! -s RAxML_info.{.}.out  ]; then  $raxml_cmd -s {} -n {.}.out 2>> raxml.log; fi' ::: *$TinEXT;
		else
		    printf "\n\nStarting tree estimation with FastTree"
		    parallel --env fasttree_cmd -j+0 'if [ ! -s {.}.tre  ]; then  $fasttree_cmd  {} > {.}.tre  2>> fasttree.log; fi' ::: *$TinEXT;
		fi
	    fi
	fi	
    fi
}
	    

shift $((OPTIND-1))
als=`ls -1 *.$EXT 2>/dev/null | wc -l` 
#printf $als
if [ $als -lt 1 ]
then
    printf "ERROR: No input files found in the current directory"
else
    find . -empty -delete
    printf "\n%s%d%s" "There are " $als " files found in the current directory." 
    main
    printf "\n\nAll files in the directory have been processed. We are done with paMATRAX+.\n"
fi
exit

