#!/bin/bash

fexons(){
## consider all exons belonging to a gene
## bed12 input file should contain gene_id at 4th column
        #mycat /main/hmtools/data/hg19_ensGene_coding.bed.gz \
        bed12_to_exon $1 \
        | awk -v OFS="\t" '{ print $1"@"$4"@"$6, $2,$3;}' \
        | sort -k1,1 -k2,3n -u | bed_flat -\
        | tr "@" "\t" |  awk -v OFS="\t" '{ print $1,$4,$5,$2,0,$3;}' \
        | sort -k1,1 -k2,3n 
}

genename(){
        cat $1 | perl -e 'use strict;
                my %eg=();
                my $file=$ARGV[0];
                open(F,$file) or die;
                while(<F>){ my ($k,$v)=split /\s/,$_; $eg{$k}=$v; }
                close(F);
                while(<STDIN>){ chomp;
                        if($_=~/(ENSG\d+)/){
                                my $tmp=$eg{$1};
                                if(defined $tmp){
                                        $_=~s/ENSG\d+/$tmp/g;
                                }
                        }
                        print $_,"\n";
                }
        ' $HMHOME/data/hg19/ensgToGenename.txt
}

