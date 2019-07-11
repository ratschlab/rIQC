#BSUB -n 1
#BSUB -R "rusage[mem=10000]"
#BSUB -W 1:00
#BSUB -J "degArr[1-130]"

set -e

outdir="/cluster/work/grlab/projects/m53"
mkdir -p $outdir
m53=/cluster/home/akahles/git/software/m53/checkBias_2.0.py
anno=${outdir}/annotation/gencode.v19.annotation.hs37d5_chr.gtf 
anno_tmp=${outdir}/annotation/gencode.v19.annotation.hs37d5_chr.deg.tmp
genome=${outdir}/genome/genome.fa.gz

### chunk size
CS=1

### list of all files
filelist=${outdir}/annotation/filelist_gtex_bams.txt
filenum=$(wc -l $filelist | cut -f 1 -d ' ')

### chunking logix
I=$(($LSB_JOBINDEX - 1))
L=$(($I * $CS))
if [ "$(($L + $CS))" -gt "$filenum" ]
then
    U=$filenum
else
    U=$((($I + 1) * $CS))
fi

mkdir -p ${outdir}/bam_based

### get files
for f in $(head -n $U $filelist | tail -n $(($U - $L)) )
do
    echo $running $f
    if [ ! -f "${outdir}/bam_based/$(basename $f)_sample_a_ratio_uq.tsv" ]
    then
        python $m53 -a $anno -m $anno_tmp -n ${f} -o ${outdir}/bam_based/$(basename $f)
    fi
done
