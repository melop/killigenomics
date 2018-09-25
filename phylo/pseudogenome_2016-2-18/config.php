<?php
/*
#
#  Configuration script for Cui et al. 2013
#  Alignment production
#  Not all configurations needed in this file, but do not delete
#
*/

// Specify these variables
define("FASTAHACK_PATH" , "fastahack"); // path to fastahack
define("FASTAHACK_BUFFERSIZE" , 1000*1000); // buffer size in bytes
define("BCFTOOLS_PATH" , "bcftools1.2"); //path to bcftools
define("BCFTOOLS_BUFFERSIZE" , 5000*1000); // buffer size in bytes
define("BCFTOOLS_QUALFILTER" , 20); //default quality filter - OBSOLETE! always use -q to explicitly specify this parameter instead.
/*
# End of configuration items needed for getTranscript*.php scripts.
# Ignore any configurations below this line.
*/



define("TEMP_DIR","/data/melop/autoblast/temp/");
define("RAXML_PATH","/media/Data/apps/RAxML-7.2.8-ALPHA/raxmlHPC-SSE3");
define("CONSEL_PATH","/media/Data/apps/consel/src/");
define("AU_TEMP_PATH","/media/Data/apps/genomewide_AU/TMP/");

define("MASK_AVAILABLE", 0);
define("MASK_OCCUPIED", pow(2,0));
define("MASK_BLASTN", pow(2,1));
define("MASK_BLASTX", pow(2,5));
define("MASK_DB_NT", pow(2,2)); //Used NT database
define("MASK_DCMEGABLAST", pow(2,3));
define("MASK_DB_NR", pow(2,4)); //Used NR database
define("BLAST_DB", "nr"); // Use NT database

define("PBS_TEMPLATE","autoblast.pbs.template");
define("Ext_Writing", ".wri");
define("Ext_Fasta", ".fas");
define("Ext_Output", ".xml");
define("Ext_Terminate", ".ter");

define("Previous_Instance_ID", 0); //7750785

//define("BLAST_EXE", "/apps/user/blast/blastn.old"); // 2.2.24 compiled by ray
//define("BLAST_EXE", "/apps/ncbi-blast/2.2.25+/GCC412-Release64/bin/blastn"); //compiled by henrik
define("BLAST_EXE", "/apps/user/ncbi-blast-2.2.25+/bin/blastx"); //precompiled

define("BLAST_CMD", "-db ".BLAST_DB." -query %FASTA_INPUT% -out %OUTPUT_FILE% -num_threads 8 -max_target_seqs 20 -outfmt 5");

define("WALL_TIME", "1:00:00");
define("NUM_PBS",1);
define("NUM_NODES",1); //nodes PER PBS!
define("NUM_PROCESSES_PER_NODE", 2);
define("CPU_TYPE", "shanghai");
define("PBS_QUEUE_NAME", "short");
define("QSUB_CMD","qsub -l nodes=".NUM_NODES.":ppn=8:".CPU_TYPE." %PBS_FILE%");

$TITLE_KEYWORD_EXCLUDE = Array(""); // use lowercase only!

set_time_limit(0);

define("DATABASE_SERVER","");
define("DATABASE_USER", "");
define("DATABASE_PWD", "");
define("DATABASE_LIMIT_TABLE", "");
define("DATABASE_LIMIT_TABLE_ID", "");
define("DATABASE_LIMIT_TABLE_BLASTFLAG", "");
define("DATABASE_TABLE", "");
define("BLAST_RET_TABLE", "");

define("CLUSTAL_PATH", "/home/gilnet/apps/clustal/clustalo");
define("CLUSTALW_PATH", "/home/gilnet/apps/clustal/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2");
define("CLUSTAL_TMP", "/home/gilnet/apps/dNdS/tmp");
define("CODEML_PATH", "/media/Data/apps/paml44/bin/codeml");
define("CHI2_PATH", "/media/Data/apps/paml44/bin/chi2");
define("CODEML_TMP", "/media/Data/apps/genomewide_dnds/codeml_tmp/");
define("CODEML_CTRL_TEMPLATE" , "/media/Data/apps/genomewide_dnds/codeml.ctl.template");

?>
