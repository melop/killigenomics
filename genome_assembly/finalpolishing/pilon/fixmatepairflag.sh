mkdir fixedbam

#Pipe each bam to fixmatepairflag.php to get the flags fixed for pilon.
#hhvm can be replaced by "php"

samtools1.2 view -h ./mapped/3kbMP.lib1.merged.bam | hhvm fixmatepairflag.php | samtools1.2 view -Sb - > fixedbam/3kbMP.lib1.merged.fixflag.bam &

samtools1.2 view -h ./mapped/8kbMP.lib2.merged.bam | hhvm fixmatepairflag.php | samtools1.2 view -Sb - > fixedbam/8kbMP.lib2.merged.fixflag.bam &

samtools1.2 view -h ./mapped/12kbMP.lib1.merged.bam | hhvm fixmatepairflag.php | samtools1.2 view -Sb - > fixedbam/12kbMP.lib1.merged.fixflag.bam & # this is only 3kb

samtools1.2 view -h ./mapped/12kbMP.lib2.merged.bam | hhvm fixmatepairflag.php | samtools1.2 view -Sb - > fixedbam/12kbMP.lib2.merged.fixflag.bam &

wait

#Index
samtools1.2 index fixedbam/3kbMP.lib1.merged.fixflag.bam &

samtools1.2 index fixedbam/12kbMP.lib1.merged.fixflag.bam &

samtools1.2 index fixedbam/8kbMP.lib2.merged.fixflag.bam &

samtools1.2 index fixedbam/12kbMP.lib2.merged.fixflag.bam &


wait
