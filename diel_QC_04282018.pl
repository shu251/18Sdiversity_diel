#! /usr/bin/perl -w
open INFILE, "prefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;

# SHu - 04-28-2018

for $i(@prefix)
{
# Merge raw reads
print "join_paired_ends.py -f /galadriel/sarah/Diel_18S_2/raw_data/",$i,"_L001_R1_001.fastq -r /galadriel/sarah/Diel_18S_2/raw_data/",$i,"_L001_R2_001.fastq -o excess_",$i," -j 20\n"; 
# Retrieve merged sequences
print "mv excess_",$i,"/fastqjoin.join.fastq ",$i,"_merged.fastq\n";

# Split - library -QC
print "split_libraries_fastq.py -i ",$i,"_merged.fastq -m /galadriel/sarah/Diel_18S_2/map_dir/",$i,"_map.txt --barcode_type 'not-barcoded' --sample_ids ",$i," -q 29 -n 0 -o excess_",$i,"/split_",$i,"\n";
print "mv excess_",$i,"/split_",$i,"/seqs.fna ",$i,"_merged_QC.fasta\n";

# Remove primers
# Remove forward primers, report in excess dir
print "cutadapt -g CCAGCASCYGCGGTAATTCC ",$i,"_merged_QC.fasta > tmpFOR_",$i,".fasta 2> excess_",$i,"/primerreportFOR_",$i,".txt\n";
# Remove reverse primers, report in excess dir
print "cutadapt -a TYRATCAAGAACGAAAGT tmpFOR_",$i,".fasta > ",$i,"_merged_QC_trim.fasta 2> excess_",$i,"/primerreportREV_",$i,".txt\n";
# primer "rm tmpFOR.fasta"

# Remove sequences by length
print "/beleriand/python_fun/seqlength_cutoff.pl ",$i,"_merged_QC_trim.fasta 150 500 ",$i,"_merged_QC_trim_len.fasta\n";

# Get stats for QC process
print "count_seqs.py -i /galadriel/sarah/Diel_18S_2/raw_data/",$i,"_L001_R1_001.fastq,",$i,"_merged.fastq,",$i,"_merged_QC.fasta,",$i,"_merged_QC_trim.fasta,",$i,"_merged_QC_trim_len.fasta,",$i,"_merged_QC_trim_len_nc.fasta > stats_",$i,".txt\n";

# Clean up directory
print "mv ",$i,"_merged.fastq excess_",$i,"\n";
print "mv ",$i,"_merged_QC.fasta excess_",$i,"\n";
print "mv ",$i,"_merged_QC_trim.fasta excess_",$i,"\n";
print "mv ",$i,"_merged_QC_trim_len.fasta excess_",$i,"\n";
} 
# Combine quality checked sequences
print "cat ",$i,"*nc.fasta >> diel_allseqs.fasta\n";

# Chimera check
print "vsearch --uchime_ref diel_allseqs.fasta --nonchimeras alldielseqs_refNC.fasta --db $PWD/pr2_version_4.72_mothur.fasta\n";

# OTU clustering 
# open ref with ref-based chimera checking:
print "pick_open_reference_otus.py -i alldielseqs_refNC.fasta -o pick_open_refNC -m uclust -r $PWD/pr2_version_4.72_mothur.fasta --suppress_step4 --suppress_taxonomy_assignment\n";
# Produced an error, but this is because I suppressed those tax assignment steps

# Assign taxonomy and generate OTU table
print "assign_taxonomy.py -i rep_set.fna -t $PWD/pr2_version_4.72_mothur.tax -r $PWD/pr2_version_4.72_mothur.fasta -m uclust -o uclust_tax --similarity 0.90 --uclust_max_accepts 3 --min_consensus_fraction 0.51\n";
print "make_otu_table.py -i final_otu_map.txt -o V4_tagseq_diel.biom -t $PWD/pr2_version_4.72_mothur.tax\n";
