
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
A draft genome for the Olympia oyster was created using a combination of short read sequence data (Illumina HiSeq4000) combined with long read sequence data (PacBio RSII) using PBJelly (English et al, 2012).

Illumina short reads (NCBI SRA: SRP072461) were assembled using SOAPdenovo (Li et al, 2008). The scaffolds (n=765,755) from this assembly were combined with the PacBio long read data (NCBI SRA: SRR5809355) using PBJelly (English et al, 2012). Assembly with PBJelly was performed using the default settings. See Supplemental File #1 for code used to run assembly.


Refs

PBJelly - https://doi.org/10.1371/journal.pone.0047768

SOAPdenovo - https://doi.org/10.1093/bioinformatics/btn025

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