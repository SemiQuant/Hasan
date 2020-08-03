freebayes does phasing but only on a single read, i.e, it wont join a mutation from read 2 with read 1.

Try out samtools phase and gatk haploype caller

One you have the phase information, use bayse to predict genotypes?





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


# error reestimation and filtering?: https://benjjneb.github.io/dada2/tutorial.html -------------------------




# correct consensus overlapping: BBMerge -------------------------

make this more stringent by masking out those that dont agree


/Users/SemiQuant/bioinfomatics/bbmap/bbmerge.sh  -Xmx"$ram"g threads=$threads \
  in="$read1" in2="$read2" \
  merge=f \
  out="${name}.R1.trimmed.overlappCor.fq" out2="${name}.R2.trimmed.overlappCor.fq"

# Allowing perfect overlaps only: pfilter=1

read1="${name}.R1.trimmed.overlappCor.fq"
read2="${name}.R2.trimmed.overlappCor.fq"


# align reads: bowtie2 -------------------------
# bowtie2-build --threads "$threads" "$ref" "${ref/.f*/}"

/Users/SemiQuant/bioinfomatics/bowtie2-2.4.1-macos-x86_64/bowtie2 --local --seed 1987 \
  -1 "$read1" -2 "$read2" \
  --rg-id "${name}" --rg "SM:${name}" \
  -p $threads -x "${ref/.f*/}" \
  -S "${name}.sam"

samtools view -b -h -@ $threads --reference "$ref" "${name}.sam" |  samtools sort - -o "${name}.bam"

# remove primers -------------------------
pysam


# call variants: freebayes -------------------------
# /Users/SemiQuant/bioinfomatics/freebayes-v1.3.1 -f "$ref" --gvcf -g 10000000 -C $minVarCall --ploidy 1 "${name}.bam" > var.vcf
freebayes -f "$ref" --gvcf -g 10000000 \
  -C $minVarCall -F 0.001 --ploidy 1 "${name}.bam" \
  --min-alternate-count 5 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous \
  --min-base-quality 30 > var.vcf

### notes

--report-monomorphic

--use-reference-allele


Generate frequency-based calls for all variants passing input thresholds. You'd do this in the case that you didn't know the number of samples in the pool.
freebayes -f ref.fa -F 0.01 -C 1 --pooled-continuous aln.bam >var.vcf


Generate long haplotype calls over known variants:

freebayes -f ref.fa --haplotype-basis-alleles in.vcf.gz \
                    --haplotype-length 50 aln.bam


Parallel operation (use 36 cores in this case):

freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 \
    -f ref.fa aln.bam >var.vcf




Calling variants in a population
freebayes is designed to be run on many individuals from the same population (e.g. many human individuals) simultaneously. The algorithm exploits a neutral model of allele diffusion to impute most-confident genotypings across the entire population. In practice, the discriminant power of the method will improve if you run multiple samples simultaneously. In other words, if your study has multiple individuals, you should run freebayes against them at the same time. This also ensures consistent reporting of information about evidence for all samples at any locus where any are apparently polymorphic.

To call variants in a population of samples, each alignment must have a read group identifier attached to it (RG tag), and the header of the BAM file in which it resides must map the RG tags to sample names (SM). Furthermore, read group IDs must be unique across all the files used in the analysis. One read group cannot map to multiple samples. The reason this is required is that freebayes operates on a virtually merged BAM stream provided by the BamTools API. If merging the files in your analysis using bamtools merge would generate a file in which multiple samples map to the same RG, the files are not suitable for use in population calling, and they must be modified.

Users may add RG tags to BAM files which were generated without this information by using (as mentioned in "Calling variants" above) bamaddrg. If you have many files corresponding to many individuals, add a unique read group and sample name to each, and then open them all simultaneously with freebayes. The VCF output will have one column per sample in the input.


# filter variants -------------------------










