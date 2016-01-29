#NPROC=8;
#BOWTIE_IDX=/storage01/home/kimji/Hm/BI2/hg19
#BOWTIE2=/storage01/home/kimji/Hm/bin/bowtie2
#TRANSCRIPTOME_IDX=/storage01/home/kimji/Hm/BI2/Homo_sapiens.Ensembl.GRCh37.65
#DATA=(
#../Fqs/ALK_ATI/SRR2040078_1.fastq.gz
#../Fqs/ALK_ATI/SRR2040078_2.fastq.gz
#mm15_rnaseq
#
#../Fqs/ALK_ATI/SRR2073855_1.fastq.gz
#../Fqs/ALK_ATI/SRR2073855_2.fastq.gz
#atc28_rnaseq
#)

tophat2.run_pair(){
	RUN=${1:-bash}; ## put run command
	TOPHAT2=${TOPHAT2:-tophat2}
	if [ -z $NPROC  ] || \
	   [ -z $BOWTIE2_IDX ] || \
	   [ ! -x $TOPHAT2 ] || \
	   [ -z $TRANSCRIPTOME_IDX ]; then 
		echo "see usage"; return;
	fi
	for (( i=0; i<${#DATA[@]};i+=3 ));do
		FQ1=${DATA[$i]};
		FQ2=${DATA[$i+1]};
		OUT=${DATA[$i+2]};

		echo "$FQ1 $FQ2 => $OUT/a.bam"
		if [[ -f $FQ1 && -f $FQ2 ]];then
			echo "run .. $RUN"
			cmd='
				#!/bin/sh
				#BSUB -n $NPROC
				#BSUB -o bsub.out.%J
				#BSUB -e bsub.err.%J
				#BSUB -J tophat2_${n}
				#mkdir -p $OUT
				'$TOPHAT2' --no-coverage-search  --no-novel-juncs --no-novel-indels \
					--transcriptome-index '$TRANSCRIPTOME_IDX' \
					-p '$NPROC' \
					-o '$OUT' '$BOWTIE2_IDX' '$FQ1' '$FQ2'
			'
			echo "$cmd" | eval "$RUN"
		fi
	done
}

