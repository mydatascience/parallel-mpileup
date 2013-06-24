#USAGE ./launch.sh ref.fasta aln.bam threads 

ref=$1
bam=$2
threads=$3

dir=$(dirname $(readlink -e $0))

$dir/group_divider $ref $threads
ls *.bed | xargs -I {} -n 1 -P 16 sh -c "samtools mpileup -uf $ref -l '{}' $bam | bcftools view -cvg - > tmp.{}.vcf"
vcf-concat tmp.*.vcf | vcf-sort > res.vcf
