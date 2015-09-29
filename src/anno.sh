#!/bin/bash


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

