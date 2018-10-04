if [ -e data/data/GRZ.brain.12.wk.JG10.fq.gz.done ]; then
 echo data/data/GRZ.brain.12.wk.JG10.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036126 | gzip -c > data/GRZ.brain.12.wk.JG10.fq.gz && touch data/data/GRZ.brain.12.wk.JG10.fq.gz.done &
 fi 
if [ -e data/data/GRZ.brain.12.wk.JG11.fq.gz.done ]; then
 echo data/data/GRZ.brain.12.wk.JG11.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036127 | gzip -c > data/GRZ.brain.12.wk.JG11.fq.gz && touch data/data/GRZ.brain.12.wk.JG11.fq.gz.done &
 fi 
if [ -e data/data/GRZ.brain.12.wk.JG13.fq.gz.done ]; then
 echo data/data/GRZ.brain.12.wk.JG13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036128 | gzip -c > data/GRZ.brain.12.wk.JG13.fq.gz && touch data/data/GRZ.brain.12.wk.JG13.fq.gz.done &
 fi 
if [ -e data/data/GRZ.brain.12.wk.JG14.fq.gz.done ]; then
 echo data/data/GRZ.brain.12.wk.JG14.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036129 | gzip -c > data/GRZ.brain.12.wk.JG14.fq.gz && touch data/data/GRZ.brain.12.wk.JG14.fq.gz.done &
 fi 
if [ -e data/data/GRZ.brain.12.wk.JG15.fq.gz.done ]; then
 echo data/data/GRZ.brain.12.wk.JG15.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036125 | gzip -c > data/GRZ.brain.12.wk.JG15.fq.gz && touch data/data/GRZ.brain.12.wk.JG15.fq.gz.done &
 fi 
if [ -e data/data/GRZ.liver.12.wk.JG10.fq.gz.done ]; then
 echo data/data/GRZ.liver.12.wk.JG10.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036131 | gzip -c > data/GRZ.liver.12.wk.JG10.fq.gz && touch data/data/GRZ.liver.12.wk.JG10.fq.gz.done &
 fi 
if [ -e data/data/GRZ.liver.12.wk.JG11.fq.gz.done ]; then
 echo data/data/GRZ.liver.12.wk.JG11.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036132 | gzip -c > data/GRZ.liver.12.wk.JG11.fq.gz && touch data/data/GRZ.liver.12.wk.JG11.fq.gz.done &
 fi 
if [ -e data/data/GRZ.liver.12.wk.JG13.fq.gz.done ]; then
 echo data/data/GRZ.liver.12.wk.JG13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036133 | gzip -c > data/GRZ.liver.12.wk.JG13.fq.gz && touch data/data/GRZ.liver.12.wk.JG13.fq.gz.done &
 fi 
if [ -e data/data/GRZ.liver.12.wk.JG14.fq.gz.done ]; then
 echo data/data/GRZ.liver.12.wk.JG14.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036134 | gzip -c > data/GRZ.liver.12.wk.JG14.fq.gz && touch data/data/GRZ.liver.12.wk.JG14.fq.gz.done &
 fi 
if [ -e data/data/GRZ.liver.12.wk.JG15.fq.gz.done ]; then
 echo data/data/GRZ.liver.12.wk.JG15.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036130 | gzip -c > data/GRZ.liver.12.wk.JG15.fq.gz && touch data/data/GRZ.liver.12.wk.JG15.fq.gz.done &
 fi 
if [ -e data/data/GRZ.skin.12.wk.JG13.fq.gz.done ]; then
 echo data/data/GRZ.skin.12.wk.JG13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036136 | gzip -c > data/GRZ.skin.12.wk.JG13.fq.gz && touch data/data/GRZ.skin.12.wk.JG13.fq.gz.done &
 fi 
if [ -e data/data/GRZ.skin.12.wk.JG14.fq.gz.done ]; then
 echo data/data/GRZ.skin.12.wk.JG14.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036139 | gzip -c > data/GRZ.skin.12.wk.JG14.fq.gz && touch data/data/GRZ.skin.12.wk.JG14.fq.gz.done &
 fi 
if [ -e data/data/GRZ.skin.12.wk.JG5.fq.gz.done ]; then
 echo data/data/GRZ.skin.12.wk.JG5.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036138 | gzip -c > data/GRZ.skin.12.wk.JG5.fq.gz && touch data/data/GRZ.skin.12.wk.JG5.fq.gz.done &
 fi 
if [ -e data/data/GRZ.skin.12.wk.JG6.fq.gz.done ]; then
 echo data/data/GRZ.skin.12.wk.JG6.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036135 | gzip -c > data/GRZ.skin.12.wk.JG6.fq.gz && touch data/data/GRZ.skin.12.wk.JG6.fq.gz.done &
 fi 
if [ -e data/data/GRZ.skin.12.wk.JG7.fq.gz.done ]; then
 echo data/data/GRZ.skin.12.wk.JG7.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2036137 | gzip -c > data/GRZ.skin.12.wk.JG7.fq.gz && touch data/data/GRZ.skin.12.wk.JG7.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.5.wk.JM100.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.5.wk.JM100.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030817 | gzip -c > data/MZM-0410.brain.5.wk.JM100.fq.gz && touch data/data/MZM-0410.brain.5.wk.JM100.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.5.wk.JM101.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.5.wk.JM101.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030818 | gzip -c > data/MZM-0410.brain.5.wk.JM101.fq.gz && touch data/data/MZM-0410.brain.5.wk.JM101.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.5.wk.JM103.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.5.wk.JM103.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030819 | gzip -c > data/MZM-0410.brain.5.wk.JM103.fq.gz && touch data/data/MZM-0410.brain.5.wk.JM103.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.5.wk.JM96.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.5.wk.JM96.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030820 | gzip -c > data/MZM-0410.brain.5.wk.JM96.fq.gz && touch data/data/MZM-0410.brain.5.wk.JM96.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.5.wk.JM97.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.5.wk.JM97.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030821 | gzip -c > data/MZM-0410.brain.5.wk.JM97.fq.gz && touch data/data/MZM-0410.brain.5.wk.JM97.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.5.wk.JM100.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.5.wk.JM100.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873320 | gzip -c > data/MZM-0410.liver.5.wk.JM100.fq.gz && touch data/data/MZM-0410.liver.5.wk.JM100.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.5.wk.JM101.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.5.wk.JM101.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873321 | gzip -c > data/MZM-0410.liver.5.wk.JM101.fq.gz && touch data/data/MZM-0410.liver.5.wk.JM101.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.5.wk.JM103.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.5.wk.JM103.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873322 | gzip -c > data/MZM-0410.liver.5.wk.JM103.fq.gz && touch data/data/MZM-0410.liver.5.wk.JM103.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.5.wk.JM97.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.5.wk.JM97.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873323 | gzip -c > data/MZM-0410.liver.5.wk.JM97.fq.gz && touch data/data/MZM-0410.liver.5.wk.JM97.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.5.wk.JM98.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.5.wk.JM98.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873324 | gzip -c > data/MZM-0410.liver.5.wk.JM98.fq.gz && touch data/data/MZM-0410.liver.5.wk.JM98.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.5.wk.JM100.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.5.wk.JM100.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873345 | gzip -c > data/MZM-0410.skin.5.wk.JM100.fq.gz && touch data/data/MZM-0410.skin.5.wk.JM100.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.5.wk.JM102.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.5.wk.JM102.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873346 | gzip -c > data/MZM-0410.skin.5.wk.JM102.fq.gz && touch data/data/MZM-0410.skin.5.wk.JM102.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.5.wk.JM103.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.5.wk.JM103.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873347 | gzip -c > data/MZM-0410.skin.5.wk.JM103.fq.gz && touch data/data/MZM-0410.skin.5.wk.JM103.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.5.wk.JM94.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.5.wk.JM94.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873348 | gzip -c > data/MZM-0410.skin.5.wk.JM94.fq.gz && touch data/data/MZM-0410.skin.5.wk.JM94.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.5.wk.JM95.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.5.wk.JM95.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873349 | gzip -c > data/MZM-0410.skin.5.wk.JM95.fq.gz && touch data/data/MZM-0410.skin.5.wk.JM95.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.12.wk.JM104.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.12.wk.JM104.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030822 | gzip -c > data/MZM-0410.brain.12.wk.JM104.fq.gz && touch data/data/MZM-0410.brain.12.wk.JM104.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.12.wk.JM105.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.12.wk.JM105.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030823 | gzip -c > data/MZM-0410.brain.12.wk.JM105.fq.gz && touch data/data/MZM-0410.brain.12.wk.JM105.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.12.wk.JM107.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.12.wk.JM107.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030824 | gzip -c > data/MZM-0410.brain.12.wk.JM107.fq.gz && touch data/data/MZM-0410.brain.12.wk.JM107.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.12.wk.JM109.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.12.wk.JM109.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030825 | gzip -c > data/MZM-0410.brain.12.wk.JM109.fq.gz && touch data/data/MZM-0410.brain.12.wk.JM109.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.12.wk.JM113.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.12.wk.JM113.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030826 | gzip -c > data/MZM-0410.brain.12.wk.JM113.fq.gz && touch data/data/MZM-0410.brain.12.wk.JM113.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.12.wk.JM104.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.12.wk.JM104.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873325 | gzip -c > data/MZM-0410.liver.12.wk.JM104.fq.gz && touch data/data/MZM-0410.liver.12.wk.JM104.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.12.wk.JM105.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.12.wk.JM105.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873326 | gzip -c > data/MZM-0410.liver.12.wk.JM105.fq.gz && touch data/data/MZM-0410.liver.12.wk.JM105.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.12.wk.JM106.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.12.wk.JM106.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873327 | gzip -c > data/MZM-0410.liver.12.wk.JM106.fq.gz && touch data/data/MZM-0410.liver.12.wk.JM106.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.12.wk.JM107.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.12.wk.JM107.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873328 | gzip -c > data/MZM-0410.liver.12.wk.JM107.fq.gz && touch data/data/MZM-0410.liver.12.wk.JM107.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.12.wk.JM115.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.12.wk.JM115.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873329 | gzip -c > data/MZM-0410.liver.12.wk.JM115.fq.gz && touch data/data/MZM-0410.liver.12.wk.JM115.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.12.wk.JM104.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.12.wk.JM104.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873350 | gzip -c > data/MZM-0410.skin.12.wk.JM104.fq.gz && touch data/data/MZM-0410.skin.12.wk.JM104.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.12.wk.JM109.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.12.wk.JM109.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873351 | gzip -c > data/MZM-0410.skin.12.wk.JM109.fq.gz && touch data/data/MZM-0410.skin.12.wk.JM109.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.12.wk.JM110.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.12.wk.JM110.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873352 | gzip -c > data/MZM-0410.skin.12.wk.JM110.fq.gz && touch data/data/MZM-0410.skin.12.wk.JM110.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.12.wk.JM112.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.12.wk.JM112.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873353 | gzip -c > data/MZM-0410.skin.12.wk.JM112.fq.gz && touch data/data/MZM-0410.skin.12.wk.JM112.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.12.wk.JM116.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.12.wk.JM116.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873354 | gzip -c > data/MZM-0410.skin.12.wk.JM116.fq.gz && touch data/data/MZM-0410.skin.12.wk.JM116.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.20.wk.JM75.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.20.wk.JM75.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030827 | gzip -c > data/MZM-0410.brain.20.wk.JM75.fq.gz && touch data/data/MZM-0410.brain.20.wk.JM75.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.20.wk.JM79.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.20.wk.JM79.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030828 | gzip -c > data/MZM-0410.brain.20.wk.JM79.fq.gz && touch data/data/MZM-0410.brain.20.wk.JM79.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.20.wk.JM82.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.20.wk.JM82.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030829 | gzip -c > data/MZM-0410.brain.20.wk.JM82.fq.gz && touch data/data/MZM-0410.brain.20.wk.JM82.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.20.wk.JM84.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.20.wk.JM84.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030830 | gzip -c > data/MZM-0410.brain.20.wk.JM84.fq.gz && touch data/data/MZM-0410.brain.20.wk.JM84.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.20.wk.JM86.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.20.wk.JM86.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030831 | gzip -c > data/MZM-0410.brain.20.wk.JM86.fq.gz && touch data/data/MZM-0410.brain.20.wk.JM86.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.20.wk.JM75.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.20.wk.JM75.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873330 | gzip -c > data/MZM-0410.liver.20.wk.JM75.fq.gz && touch data/data/MZM-0410.liver.20.wk.JM75.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.20.wk.JM79.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.20.wk.JM79.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873331 | gzip -c > data/MZM-0410.liver.20.wk.JM79.fq.gz && touch data/data/MZM-0410.liver.20.wk.JM79.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.20.wk.JM82.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.20.wk.JM82.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873332 | gzip -c > data/MZM-0410.liver.20.wk.JM82.fq.gz && touch data/data/MZM-0410.liver.20.wk.JM82.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.20.wk.JM84.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.20.wk.JM84.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873333 | gzip -c > data/MZM-0410.liver.20.wk.JM84.fq.gz && touch data/data/MZM-0410.liver.20.wk.JM84.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.20.wk.JM85.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.20.wk.JM85.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873334 | gzip -c > data/MZM-0410.liver.20.wk.JM85.fq.gz && touch data/data/MZM-0410.liver.20.wk.JM85.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.20.wk.JM79.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.20.wk.JM79.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873355 | gzip -c > data/MZM-0410.skin.20.wk.JM79.fq.gz && touch data/data/MZM-0410.skin.20.wk.JM79.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.20.wk.JM84.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.20.wk.JM84.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873356 | gzip -c > data/MZM-0410.skin.20.wk.JM84.fq.gz && touch data/data/MZM-0410.skin.20.wk.JM84.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.20.wk.JM87.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.20.wk.JM87.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873357 | gzip -c > data/MZM-0410.skin.20.wk.JM87.fq.gz && touch data/data/MZM-0410.skin.20.wk.JM87.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.20.wk.JM88.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.20.wk.JM88.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873358 | gzip -c > data/MZM-0410.skin.20.wk.JM88.fq.gz && touch data/data/MZM-0410.skin.20.wk.JM88.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.20.wk.JM91.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.20.wk.JM91.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873359 | gzip -c > data/MZM-0410.skin.20.wk.JM91.fq.gz && touch data/data/MZM-0410.skin.20.wk.JM91.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.27.wk.JM52.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.27.wk.JM52.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030832 | gzip -c > data/MZM-0410.brain.27.wk.JM52.fq.gz && touch data/data/MZM-0410.brain.27.wk.JM52.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.27.wk.JM62.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.27.wk.JM62.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030833 | gzip -c > data/MZM-0410.brain.27.wk.JM62.fq.gz && touch data/data/MZM-0410.brain.27.wk.JM62.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.27.wk.JM66.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.27.wk.JM66.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030834 | gzip -c > data/MZM-0410.brain.27.wk.JM66.fq.gz && touch data/data/MZM-0410.brain.27.wk.JM66.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.27.wk.JM69.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.27.wk.JM69.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030835 | gzip -c > data/MZM-0410.brain.27.wk.JM69.fq.gz && touch data/data/MZM-0410.brain.27.wk.JM69.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.27.wk.JM72.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.27.wk.JM72.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030836 | gzip -c > data/MZM-0410.brain.27.wk.JM72.fq.gz && touch data/data/MZM-0410.brain.27.wk.JM72.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.27.wk.JM61.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.27.wk.JM61.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873335 | gzip -c > data/MZM-0410.liver.27.wk.JM61.fq.gz && touch data/data/MZM-0410.liver.27.wk.JM61.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.27.wk.JM69.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.27.wk.JM69.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873336 | gzip -c > data/MZM-0410.liver.27.wk.JM69.fq.gz && touch data/data/MZM-0410.liver.27.wk.JM69.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.27.wk.JM72.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.27.wk.JM72.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873337 | gzip -c > data/MZM-0410.liver.27.wk.JM72.fq.gz && touch data/data/MZM-0410.liver.27.wk.JM72.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.27.wk.JM74.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.27.wk.JM74.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873338 | gzip -c > data/MZM-0410.liver.27.wk.JM74.fq.gz && touch data/data/MZM-0410.liver.27.wk.JM74.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.27.wk.JM83.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.27.wk.JM83.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873339 | gzip -c > data/MZM-0410.liver.27.wk.JM83.fq.gz && touch data/data/MZM-0410.liver.27.wk.JM83.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.27.wk.JM55.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.27.wk.JM55.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873360 | gzip -c > data/MZM-0410.skin.27.wk.JM55.fq.gz && touch data/data/MZM-0410.skin.27.wk.JM55.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.27.wk.JM62.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.27.wk.JM62.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873361 | gzip -c > data/MZM-0410.skin.27.wk.JM62.fq.gz && touch data/data/MZM-0410.skin.27.wk.JM62.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.27.wk.JM72.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.27.wk.JM72.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873362 | gzip -c > data/MZM-0410.skin.27.wk.JM72.fq.gz && touch data/data/MZM-0410.skin.27.wk.JM72.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.27.wk.JM74.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.27.wk.JM74.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873363 | gzip -c > data/MZM-0410.skin.27.wk.JM74.fq.gz && touch data/data/MZM-0410.skin.27.wk.JM74.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.27.wk.JM83.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.27.wk.JM83.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873364 | gzip -c > data/MZM-0410.skin.27.wk.JM83.fq.gz && touch data/data/MZM-0410.skin.27.wk.JM83.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM1.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM1.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030837 | gzip -c > data/MZM-0410.brain.39.wk.JM1.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM1.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM12.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM12.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030838 | gzip -c > data/MZM-0410.brain.39.wk.JM12.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM12.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM13.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030839 | gzip -c > data/MZM-0410.brain.39.wk.JM13.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM13.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM15.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM15.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030840 | gzip -c > data/MZM-0410.brain.39.wk.JM15.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM15.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM18.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM18.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1030841 | gzip -c > data/MZM-0410.brain.39.wk.JM18.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM18.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM24.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM24.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225623 | gzip -c > data/MZM-0410.brain.39.wk.JM24.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM24.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM38.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM38.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225624 | gzip -c > data/MZM-0410.brain.39.wk.JM38.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM38.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM41.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM41.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225626 | gzip -c > data/MZM-0410.brain.39.wk.JM41.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM41.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM46.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM46.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225628 | gzip -c > data/MZM-0410.brain.39.wk.JM46.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM46.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM6.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM6.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225630 | gzip -c > data/MZM-0410.brain.39.wk.JM6.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM6.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.brain.39.wk.JM8.fq.gz.done ]; then
 echo data/data/MZM-0410.brain.39.wk.JM8.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225632 | gzip -c > data/MZM-0410.brain.39.wk.JM8.fq.gz && touch data/data/MZM-0410.brain.39.wk.JM8.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.39.wk.JM1.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.39.wk.JM1.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873340 | gzip -c > data/MZM-0410.liver.39.wk.JM1.fq.gz && touch data/data/MZM-0410.liver.39.wk.JM1.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.39.wk.JM12.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.39.wk.JM12.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873341 | gzip -c > data/MZM-0410.liver.39.wk.JM12.fq.gz && touch data/data/MZM-0410.liver.39.wk.JM12.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.39.wk.JM13.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.39.wk.JM13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873342 | gzip -c > data/MZM-0410.liver.39.wk.JM13.fq.gz && touch data/data/MZM-0410.liver.39.wk.JM13.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.39.wk.JM15.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.39.wk.JM15.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873343 | gzip -c > data/MZM-0410.liver.39.wk.JM15.fq.gz && touch data/data/MZM-0410.liver.39.wk.JM15.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.liver.39.wk.JM18.fq.gz.done ]; then
 echo data/data/MZM-0410.liver.39.wk.JM18.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873344 | gzip -c > data/MZM-0410.liver.39.wk.JM18.fq.gz && touch data/data/MZM-0410.liver.39.wk.JM18.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM1.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM1.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873365 | gzip -c > data/MZM-0410.skin.39.wk.JM1.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM1.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM13.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM13.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873366 | gzip -c > data/MZM-0410.skin.39.wk.JM13.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM13.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM15.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM15.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873367 | gzip -c > data/MZM-0410.skin.39.wk.JM15.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM15.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM18.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM18.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873368 | gzip -c > data/MZM-0410.skin.39.wk.JM18.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM18.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM24.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM24.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR1873369 | gzip -c > data/MZM-0410.skin.39.wk.JM24.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM24.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM12.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM12.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225622 | gzip -c > data/MZM-0410.skin.39.wk.JM12.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM12.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM38.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM38.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225625 | gzip -c > data/MZM-0410.skin.39.wk.JM38.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM38.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM41.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM41.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225627 | gzip -c > data/MZM-0410.skin.39.wk.JM41.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM41.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM46.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM46.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225629 | gzip -c > data/MZM-0410.skin.39.wk.JM46.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM46.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM6.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM6.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225631 | gzip -c > data/MZM-0410.skin.39.wk.JM6.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM6.fq.gz.done &
 fi 
if [ -e data/data/MZM-0410.skin.39.wk.JM8.fq.gz.done ]; then
 echo data/data/MZM-0410.skin.39.wk.JM8.fq.gz.done finished; 
else
 fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@$sn' SRR2225633 | gzip -c > data/MZM-0410.skin.39.wk.JM8.fq.gz && touch data/data/MZM-0410.skin.39.wk.JM8.fq.gz.done &
 fi 

wait
