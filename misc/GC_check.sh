#!/bin/bash
# This script is to calculate GC content for coding and non-coding regions of reference genome

# Reference genome: /local/storage/CV30_reference_genome/CV30_masked_nomtDNA.fasta
# No haplotigs CDS: gff_cds.sorted.merged.nohaplotigs.bed

bedtools nuc -fi /local/storage/CV30_reference_genome/CV30_masked_nomtDNA.fasta -bed gff_cds.sorted.merged.nohaplotigs.bed > gff_cds.sorted.merged.nohaplotigs_gc_content.txt

# Column 5 = pct_gc (proportion GC)
# Column 7: num_C
# Column 8: num_G
# Column 12 = seq_len

# Then get global GC content: 
awk 'BEGIN{total_gc=0; total_len=0} 
     !/^#/ {
        total_gc += ($7 + $8);
        total_len += $12;
     }
     END{
        printf("Overall CDS GC content: %.4f%%\n", (total_gc / total_len) * 100);
     }' gff_cds.sorted.merged.nohaplotigs_gc_content.txt

# No haplotigs non-CDS: 
cut -f1,2 /local/storage/CV30_reference_genome/CV30_masked_nomtDNA.fasta.fai > genome_file.txt
bedtools complement -i gff_cds.sorted.merged.nohaplotigs.bed -g genome_file.txt > complement_gff_cds.sorted.merged.nohaplotigs.bed

bedtools nuc -fi /local/storage/CV30_reference_genome/CV30_masked_nomtDNA.fasta -bed complement_gff_cds.sorted.merged.nohaplotigs.bed > complement_gff_cds.sorted.merged.nohaplotigs.bed_gc_content.txt

# Then get global GC content: 
awk 'BEGIN{total_gc=0; total_len=0} 
     !/^#/ {
        total_gc += ($7 + $8);
        total_len += $12;
     }
     END{
        printf("Overall CDS GC content: %.4f%%\n", (total_gc / total_len) * 100);
     }' complement_gff_cds.sorted.merged.nohaplotigs.bed_gc_content.txt
