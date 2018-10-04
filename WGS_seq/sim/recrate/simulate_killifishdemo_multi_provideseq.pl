#perl! -w                                                                                                                                                                           
#script to simulate family data and a population sample from macs and format the data                                                                         
if(@ARGV < 5){

    print "perl simulate_SNP_data.pl fasta_file popname num_indivs avg_COV stdev_COV Locus_len\n"; exit;

}

my $read1len = 142;

my $read2len = 151;

my $fragsize = 336 - 9; # average insert size based on bioanalyzer results sent by ngxbio NGX Bio Order 15204. 9 is the number of bases in the R1 barcode.

my $file=shift(@ARGV); chomp $file;

my $popname=shift(@ARGV); chomp $popname;

my $total_indiv=shift(@ARGV); chomp $total_indiv;

my $avg_cov=shift(@ARGV); chomp $avg_cov;

my $stdev_cov=shift(@ARGV); chomp $stdev_cov;

my $locuslen=shift(@ARGV); chomp $locuslen;

my $avg_reads = int(($avg_cov * $locuslen) / ($read1len + $read2len));

my $stdev_reads= int(($stdev_cov * $locuslen) / ($read1len + $read2len));

my $number=int(rand(1000000)); chomp $number;

my $seqerror = 0.0033; # from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y

my $errindelprop = 0.001235; # from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y

print "locus len: $locuslen avg reads $avg_reads, stdev reads $stdev_reads\n";
my $folder="pop1_and_pop2_reads.fq_parsed_"."$popname";
if(-d $folder){
    print "removing previous folder\n";
system("rm -r $folder");
}#remove if exists

system("mkdir -p $folder");

my $hap1_num=1; my $hap2_num=2;
for my $j (1..$total_indiv){
    my $indiv="./$folder/indiv"."$j".".fa";

    my $hap1_name="$popname".".hap"."$hap1_num";
    my $hap2_name="$popname".".hap"."$hap2_num";

    print "$hap1_name\t$hap2_name\n";

    system("grep -wA 1 $hap1_name $file > $indiv");
    system("grep -wA 1 $hap2_name $file >> $indiv");

    if ($j == 1){
        system("mkdir $folder");
    }#output folder                                                                                                                                                                                                 
    my $fastq1="$popname"."$j"."_read1.fastq";
    my $fastq2="$popname"."$j"."_read2.fastq";

    my $num_reads=qx(Rscript getreads.R $avg_reads $stdev_reads | perl -pi -e 's/ +/\t/g' | cut -f 2); chomp $num_reads;

    my $seed=rand(1000);

    print "seed for simulation is $seed\n";
    system("./wgsim -N$num_reads -1$read1len -2$read2len -d$fragsize -S$seed -e$seqerror -r0 -R$errindelprop $indiv ./$folder/$fastq1 ./$folder/$fastq2");

    system("perl -pi -e 's/^(:{2,})/'A' x length(\$1)/e' ./$folder/$fastq1");

    system("perl -pi -e 's/^(:{2,})/'A' x length(\$1)/e' ./$folder/$fastq2");

    #system("rm $indiv");

    system("gzip ./$folder/$fastq1");
    system("gzip ./$folder/$fastq2");

    $hap2_num=$hap1_num+3;
    $hap1_num=$hap1_num+2;

}

