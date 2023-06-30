
counter=1
ls /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/*.dupmarked.realigned.recalibrated.MT.bam | while read -r bam1; do

    bai1=$bam1.bai
    read -r bam2
    bai2=$bam2.bai
    read -r bam3
    bai3=$bam3.bai
    read -r bam4
    bai4=$bam4.bai
    read -r bam5
    bai5=$bam5.bai
    read -r bam6
    bai6=$bam6.bai
    read -r bam7
    bai7=$bam7.bai
    read -r bam8
    bai8=$bam8.bai
    read -r bam9
    bai9=$bam9.bai
    read -r bam10
    bai10=$bam10.bai
    read -r bam11
    bai11=$bam11.bai
    read -r bam12
    bai12=$bam12.bai
    read -r bam13
    bai13=$bam13.bai
    read -r bam14
    bai14=$bam14.bai
    read -r bam15
    bai15=$bam15.bai
    
    echo "

    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam1* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam2* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam3* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam4* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam5* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam6* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam7* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam8* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam9* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam10* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam11* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam12* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam13* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam14* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 
    rsync --progress -rhL /g/data3/vj26/results/phase2/hs37d5x/bwa_IR_BQSR_mt/$bam15* /g/data2/gd7/research/claput/MT_MGRB/MGRB_bams 

    " > /g/data2/gd7/research/claput/MT_MGRB/cp_MGRB_bams_sh/run$counter.sh

done

for inbam in /g/data2/gd7/research/claput/MT_MGRB/test_bams/*.bam
do echo $inbam
sample_name=$(basename $inbam | cut -d . -f1 ) 
echo $sample_name
shTem=/g/data2/gd7/research/claput/MT_MGRB/mity_call_sh/$sample_name

-q copyq -l mem=500GB -e $shTem.e -o $shTem.o -l wd $shTem.sh

done


# expressbw for testing 