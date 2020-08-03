rm "${name}.R1.trimmed.fq" "${name}.R2.trimmed.fq" "${name}.sam" "${name}.bam"
ram=4
threads=2
# read1=TB-B265-054_S32_L001_R1_001.fastq.gz
# read2=TB-B265-054_S32_L001_R2_001.fastq.gz
read1=fake.1.txt
read2=fake.2.txt
name=test1
adapters_file=illumina_adapters_all.fasta
minL=50
ref=Rv0678.fasta
minVarCall=5

# trim read: bbduk -------------------------
/Users/SemiQuant/bioinfomatics/bbmap/bbduk.sh -Xmx"$ram"g threads=$threads \
  in="$read1" in2="$read2" \
  out="${name}.R1.trimmed.fq" out2="${name}.R2.trimmed.fq" \
  minlen="$minL" \
  stats="${name}_trim_stats.txt" \
  ordered=t \
  hdist=1 \
  qtrim=r \
  trimq=20 \
  ref="$adapters_file" \
  ktrim=r \
  mink=9



read1="${name}.R1.trimmed.fq"
read2="${name}.R2.trimmed.fq"




# align reads: bowtie2 -------------------------
# bowtie2-build --threads "$threads" "$ref" "${ref/.f*/}"

/Users/SemiQuant/bioinfomatics/bowtie2-2.4.1-macos-x86_64/bowtie2 --local --seed 1987 \
  -1 "$read1" -2 "$read2" \
  --rg-id "${name}" --rg "SM:${name}" \
  -p $threads -x "${ref/.f*/}" \
  -S "${name}.sam"

samtools view -b -h -@ $threads --reference "$ref" "${name}.sam" |  samtools sort - -o "${name}.bam"



freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 \
  --haplotype-length 600
  








# call variants: freebayes -------------------------
# /Users/SemiQuant/bioinfomatics/freebayes-v1.3.1 -f "$ref" --gvcf -g 10000000 -C $minVarCall --ploidy 1 "${name}.bam" > var.vcf
freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 > var.vcf





samtools phase -F -k 500 -b "phase" -q 5 -Q 20 "${name}.bam"
# Prefix of BAM output. When this option is in use, phase-0 reads will be saved in file STR.0.bam and phase-1 reads in STR.1.bam. Phase unknown reads will be randomly allocated to one of the two files. Chimeric reads with switch errors will be saved in STR.chimeric.bam. [null]


samtools phase -F -k 5 -b "phase" -q 1 -Q 5 "${name}.bam"


/Users/SemiQuant/Bioinformatics/Projects/Hasan/test/whatshap-env/bin/whatshap phase -o phased.vcf --reference="$ref" var.vcf "${name}.bam"


https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller




python2.7 /Users/SemiQuant/Downloads/BHap_v1/run_BHap.py \
  -d ./tmpName \
  -r "$ref" \
  -t fastq \
  -1 $read1 \
  -2 $read2 \
  -l 100 \
  -c 10 \
  -g 598




shorah "${name}.bam"

 amplicon [-h] [-v] -b BAM -f REF [-a FLOAT] [-r chrm:start-stop] [-R INT] [-x INT] [-S FLOAT] [-I] [-of {csv,vcf} [{csv,vcf} ...]] [-c INT] [-d] [-m FLOAT]
shorah amplicon -b "${name}.bam" -f "$ref"



shorah amplicon --bam "test1.bam" -r "RV0678_-100toend:1-590" -f "Rv0678.fasta"
shorah shotgun --bam "${name}.bam" -r "RV0678_-100toend:1-598" -f "$ref"


shorah amplicon --bam "/Users/SemiQuant/Downloads/shorah-1.99.0/examples/amplicon_test/ampli_sorted.bam" -r "reference:1-73" -f "/Users/SemiQuant/Downloads/shorah-1.99.0/examples/amplicon_test/reference.fasta"


singularity shell ../del/shorah_1.99.0--py37h41a77f4_0.sif


cp /Users/SemiQuant/Downloads/bdq_smorTrails/reference.fasta /Users/SemiQuant/Downloads/bdq_smorTrails/clean-bowtie2.bam ./
samtools sort clean-bowtie2.bam -o sorted.bam
shorah amplicon --bam "sorted.bam" -r "Rv0678_amplicon4:5-50" -f "reference.fasta"




singularity shell ../del/shorah_1.99.0--py37h41a77f4_0.sif
rm /Users/SemiQuant/Bioinformatics/Projects/Hasan/test/shora_trial/shorah.log
shorah amplicon --bam "sorted.bam" -r "RV0678_-100toend:1-150" -f "reference.fasta"




# gatk CreateSequenceDictionary -R "$ref"

 gatk HaplotypeCaller \
   -R "$ref" \
   -I "${name}.bam" \
   -O gatk.output.vcf.gz \
   -bamout bamout.bam \
   --base-quality-score-threshold 30 \
   --min-base-quality-score 30 \
   -ploidy 20

gunzip gatk.output.vcf.gz
cat gatk.output.vcf



freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 \
  --max-complex-gap 498 > v1.vcf



freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 \
  --haplotype-length 498 > v2.vcf



freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 \
  --max-complex-gap 498 \
  --haplotype-length 498 > v3.vcf