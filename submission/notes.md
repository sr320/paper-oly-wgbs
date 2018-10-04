
```
Samples
1NF11
1NF15
1NF16
1NF17
2NF5
2NF6
2NF7
2NF8
```

Oyster weights
https://docs.google.com/spreadsheets/d/1aFyr7--nwnetk_AtYEkpCR8O6C6Nw4DnNiGUh8DUvh8/edit#gid=0

---




Refs

PBJelly - https://doi.org/10.1371/journal.pone.0047768

SOAPdenovo - https://doi.org/10.1093/bioinformatics/btn025

methylKit - https://doi.org/10.1186/gb-2012-13-10-r87

QUAST - https://doi.org/10.1093/bioinformatics/btt086



## SUPPLEMENTAL FILE #1

---

#### PBJelly Protocol.xml
```
<jellyProtocol>
    <reference>/home/sam/data/oly_BGI_scaffolds.fasta</reference>  
    <outputDir>/home/sam/analyses/20171130_oly_pbjelly</outputDir>
    <blasr>-minMatch 8 -minPctIdentity 70 -bestn 1 -nCandidates 20 -maxScore -500 -nproc 48 -noSplitSubreads</blasr>
    <input baseDir="/home/sam/data/">
        <job>m130619_081336_42134_c100525122550000001823081109281326_s1_p0.fastq</job>
        <job>m170211_224036_42134_c101073082550000001823236402101737_s1_X0_filtered_subreads.fastq</job>
	<job>m170301_100013_42134_c101174162550000001823269408211761_s1_p0_filtered_subreads.fastq</job>
	<job>m170301_162825_42134_c101174162550000001823269408211762_s1_p0_filtered_subreads.fastq</job>
	<job>m170301_225711_42134_c101174162550000001823269408211763_s1_p0_filtered_subreads.fastq</job>
	<job>m170308_163922_42134_c101174252550000001823269408211742_s1_p0_filtered_subreads.fastq</job>
	<job>m170308_230815_42134_c101174252550000001823269408211743_s1_p0_filtered_subreads.fastq</job>
	<job>m170315_001112_42134_c101169372550000001823273008151717_s1_p0_filtered_subreads.fastq</job>
	<job>m170315_063041_42134_c101169382550000001823273008151700_s1_p0_filtered_subreads.fastq</job>
	<job>m170315_124938_42134_c101169382550000001823273008151701_s1_p0_filtered_subreads.fastq</job>
	<job>m170315_190851_42134_c101169382550000001823273008151702_s1_p0_filtered_subreads.fastq</job>
    </input>
</jellyProtocol>
```

---

#### PBJelly code

source /home/shared/PBSuite_15.8.24/setup.sh
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py setup /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py mapping /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py support /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py extraction /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py assembly /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml
time python /home/shared/PBSuite_15.8.24/bin/Jelly.py output /home/sam/analyses/20171130_oly_pbjelly/Protocol.xml



---


#### Bismark code

# Prep genome
/home/shared/Bismark-0.19.1/bismark_genome_preparation \
--path_to_bowtie /home/shared/bowtie2-2.3.4.1-linux-x86_64/ \
--verbose /home/sam/data/oly_methylseq/oly_genome/ \
2> 20180507_bismark_genome_prep.err

# Directories and programs

bismark_dir="/gscratch/srlab/programs/Bismark-0.19.0"
trimmed="/gscratch/scrubbed/samwhite/data/O_lurida/BSseq/whole_genome_BSseq_reads/20180830_trimgalore"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
genome="/gscratch/scrubbed/samwhite/data/O_lurida/BSseq/20180503_oly_genome_pbjelly_sjw_01_bismark/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"

${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
--genome ${genome} \
--score_min L,0,-0.6 \
-p 28 \
--non_directional \
${trimmed}/1_ATCACG_L001_R1_001_trimmed.fq.gz \
${trimmed}/2_CGATGT_L001_R1_001_trimmed.fq.gz \
${trimmed}/3_TTAGGC_L001_R1_001_trimmed.fq.gz \
${trimmed}/4_TGACCA_L001_R1_001_trimmed.fq.gz \
${trimmed}/5_ACAGTG_L001_R1_001_trimmed.fq.gz \
${trimmed}/6_GCCAAT_L001_R1_001_trimmed.fq.gz \
${trimmed}/7_CAGATC_L001_R1_001_trimmed.fq.gz \
${trimmed}/8_ACTTGA_L001_R1_001_trimmed.fq.gz

# Deduplicate bam files

${bismark_dir}/deduplicate_bismark \
--bam \
--single \
*.bam

# Methylation extraction

${bismark_dir}/bismark_methylation_extractor \
--bedgraph \
--counts \
--scaffolds \
--remove_spaces \
--multicore 28 \
--buffer_size 75% \
*deduplicated.bam

# Bismark processing report

${bismark_dir}/bismark2report

#Bismark summary report

${bismark_dir}/bismark2summary

# Sort files for methylkit and IGV

find *deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} ${samtools} \
sort --threads 28 {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ 16" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ 28 {}.sorted.bam

---

Refs

PBJelly - https://doi.org/10.1371/journal.pone.0047768

SOAPdenovo - https://doi.org/10.1093/bioinformatics/btn025

## SUPPLEMENTAL FILE #1


`https://raw.githubusercontent.com/sr320/nb-2018/master/O_lurida/analyses/0919_igv_session.xml`
