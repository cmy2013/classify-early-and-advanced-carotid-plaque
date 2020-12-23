STAR=~/Software/STAR-2.5.3a/bin/Linux_x86_64/STAR
data_home=fq
gtf=~/common_data/Homo_sapiens.GRCh38.99.gtf

while read line; do
     file1=$line
     echo $file1
     read line
     file2=$line
     echo $file2
	 output_name_tmp=$(cut -d'/' -f7 <<< "$line")	 
	 output_name=$(cut -d'_' -f1 <<< "$output_name_tmp")	 
	 echo $output_name


	 # $STAR --runThreadN 20 --genomeDir STAR_index/Ensembl38 --readFilesIn ${data_home}/${file1} ${data_home}/${file2} --readFilesCommand - --sjdbGTFfile $gtf --outFileNamePrefix output/${output_name} --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --limitBAMsortRAM 321098844080 --outSAMunmapped Within &> $output_name.log.o 2> $output_name.log.e
	 $STAR --runThreadN 20 --genomeDir STAR_index/Ensembl38 --readFilesIn ${data_home}/${file1} ${data_home}/${file2} --readFilesCommand - --sjdbGTFfile $gtf --outFileNamePrefix star_output/${output_name} --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --limitBAMsortRAM 321098844080 --outSAMunmapped Within &> $output_name.log.o 2> $output_name.log.e
done < sample_list.txt


# ############################################################################
# ##mapping rate
# file_home="/data1_2/rhu1/07022020_TUJIAN/Analysis/output"
# for file in $(ls $file_home/*.progress.out);
# do
# echo $(basename $file)
# tail -n 2 $file
# done


rsem=~/Software/RSEM-1.3.0/rsem-calculate-expression
data_folder=star_output
rsem_index=RSEM_ref/Ensemble38/human_ensembl

$rsem --bam ${data_folder}/SRR6059634Aligned.toTranscriptome.out.bam $rsem_index SRR6059634.rsem --paired-end -p 6 --no-bam-output &> SRR6059634.log.o.rsem
$rsem --bam ${data_folder}/SRR6059635Aligned.toTranscriptome.out.bam $rsem_index SRR6059635.rsem --paired-end -p 6 --no-bam-output &> SRR6059635.log.o.rsem
$rsem --bam ${data_folder}/SRR6059636Aligned.toTranscriptome.out.bam $rsem_index SRR6059636.rsem --paired-end -p 6 --no-bam-output &> SRR6059636.log.o.rsem
$rsem --bam ${data_folder}/SRR6059637Aligned.toTranscriptome.out.bam $rsem_index SRR6059637.rsem --paired-end -p 6 --no-bam-output &> SRR6059637.log.o.rsem
$rsem --bam ${data_folder}/SRR6059638Aligned.toTranscriptome.out.bam $rsem_index SRR6059638.rsem --paired-end -p 6 --no-bam-output &> SRR6059638.log.o.rsem
$rsem --bam ${data_folder}/SRR6059639Aligned.toTranscriptome.out.bam $rsem_index SRR6059639.rsem --paired-end -p 6 --no-bam-output &> SRR6059639.log.o.rsem
$rsem --bam ${data_folder}/SRR6059640Aligned.toTranscriptome.out.bam $rsem_index SRR6059640.rsem --paired-end -p 6 --no-bam-output &> SRR6059640.log.o.rsem
$rsem --bam ${data_folder}/SRR6059641Aligned.toTranscriptome.out.bam $rsem_index SRR6059641.rsem --paired-end -p 6 --no-bam-output &> SRR6059641.log.o.rsem
$rsem --bam ${data_folder}/SRR6059642Aligned.toTranscriptome.out.bam $rsem_index SRR6059642.rsem --paired-end -p 6 --no-bam-output &> SRR6059642.log.o.rsem
$rsem --bam ${data_folder}/SRR6059643Aligned.toTranscriptome.out.bam $rsem_index SRR6059643.rsem --paired-end -p 6 --no-bam-output &> SRR6059643.log.o.rsem
$rsem --bam ${data_folder}/SRR6059644Aligned.toTranscriptome.out.bam $rsem_index SRR6059644.rsem --paired-end -p 6 --no-bam-output &> SRR6059644.log.o.rsem
$rsem --bam ${data_folder}/SRR6059645Aligned.toTranscriptome.out.bam $rsem_index SRR6059645.rsem --paired-end -p 6 --no-bam-output &> SRR6059645.log.o.rsem
$rsem --bam ${data_folder}/SRR6059646Aligned.toTranscriptome.out.bam $rsem_index SRR6059646.rsem --paired-end -p 6 --no-bam-output &> SRR6059646.log.o.rsem
$rsem --bam ${data_folder}/SRR6059647Aligned.toTranscriptome.out.bam $rsem_index SRR6059647.rsem --paired-end -p 6 --no-bam-output &> SRR6059647.log.o.rsem
$rsem --bam ${data_folder}/SRR6059648Aligned.toTranscriptome.out.bam $rsem_index SRR6059648.rsem --paired-end -p 6 --no-bam-output &> SRR6059648.log.o.rsem
$rsem --bam ${data_folder}/SRR6059649Aligned.toTranscriptome.out.bam $rsem_index SRR6059649.rsem --paired-end -p 6 --no-bam-output &> SRR6059649.log.o.rsem
$rsem --bam ${data_folder}/SRR6059650Aligned.toTranscriptome.out.bam $rsem_index SRR6059650.rsem --paired-end -p 6 --no-bam-output &> SRR6059650.log.o.rsem
$rsem --bam ${data_folder}/SRR6059651Aligned.toTranscriptome.out.bam $rsem_index SRR6059651.rsem --paired-end -p 6 --no-bam-output &> SRR6059651.log.o.rsem
$rsem --bam ${data_folder}/SRR6059652Aligned.toTranscriptome.out.bam $rsem_index SRR6059652.rsem --paired-end -p 6 --no-bam-output &> SRR6059652.log.o.rsem
$rsem --bam ${data_folder}/SRR6059653Aligned.toTranscriptome.out.bam $rsem_index SRR6059653.rsem --paired-end -p 6 --no-bam-output &> SRR6059653.log.o.rsem
$rsem --bam ${data_folder}/SRR6059654Aligned.toTranscriptome.out.bam $rsem_index SRR6059654.rsem --paired-end -p 6 --no-bam-output &> SRR6059654.log.o.rsem
$rsem --bam ${data_folder}/SRR6059655Aligned.toTranscriptome.out.bam $rsem_index SRR6059655.rsem --paired-end -p 6 --no-bam-output &> SRR6059655.log.o.rsem
$rsem --bam ${data_folder}/SRR6059656Aligned.toTranscriptome.out.bam $rsem_index SRR6059656.rsem --paired-end -p 6 --no-bam-output &> SRR6059656.log.o.rsem
$rsem --bam ${data_folder}/SRR6059657Aligned.toTranscriptome.out.bam $rsem_index SRR6059657.rsem --paired-end -p 6 --no-bam-output &> SRR6059657.log.o.rsem
$rsem --bam ${data_folder}/SRR6059658Aligned.toTranscriptome.out.bam $rsem_index SRR6059658.rsem --paired-end -p 6 --no-bam-output &> SRR6059658.log.o.rsem
$rsem --bam ${data_folder}/SRR6059659Aligned.toTranscriptome.out.bam $rsem_index SRR6059659.rsem --paired-end -p 6 --no-bam-output &> SRR6059659.log.o.rsem
$rsem --bam ${data_folder}/SRR6059660Aligned.toTranscriptome.out.bam $rsem_index SRR6059660.rsem --paired-end -p 6 --no-bam-output &> SRR6059660.log.o.rsem
$rsem --bam ${data_folder}/SRR6059661Aligned.toTranscriptome.out.bam $rsem_index SRR6059661.rsem --paired-end -p 6 --no-bam-output &> SRR6059661.log.o.rsem
$rsem --bam ${data_folder}/SRR6059662Aligned.toTranscriptome.out.bam $rsem_index SRR6059662.rsem --paired-end -p 6 --no-bam-output &> SRR6059662.log.o.rsem
$rsem --bam ${data_folder}/SRR6059663Aligned.toTranscriptome.out.bam $rsem_index SRR6059663.rsem --paired-end -p 6 --no-bam-output &> SRR6059663.log.o.rsem
$rsem --bam ${data_folder}/SRR6059664Aligned.toTranscriptome.out.bam $rsem_index SRR6059664.rsem --paired-end -p 6 --no-bam-output &> SRR6059664.log.o.rsem
$rsem --bam ${data_folder}/SRR6059665Aligned.toTranscriptome.out.bam $rsem_index SRR6059665.rsem --paired-end -p 6 --no-bam-output &> SRR6059665.log.o.rsem

