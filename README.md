# FuMa #

	                        7
	                     .:OMNZ7Z$I,..78
	                    788:.:,.....DMD:
	                   8DO$,. .~,...I8,
	                  $DZI,..:ZO$?$D$$
	                .88$..=7..=D=:,DIO?
	                8DZ,..?NMO...$?.DD
	    .         ~8Z......ZM+..87O=88
	   .7$7      ONM~....,.,=.?MMNM~.:ZZ
	 .:NO+II.   ,NNN7.....O$..:DMN,..,O8
	  OMI.:I.  .ONNNN....NN8OOO$..... NN
	 .$NI..:.  .$8Z$Z,.+DDO.7DN~..Z8.~D=
	 ..Z=...~.  ONDZ..:7OO,.+DD,..8MMM$
	  .7+...Z,..I8O?+?++I?..78,..,ONMM
	  ..O=..Z7?==,..7I++?...,O$..ZNMZ
	   . 8Z7..:++.....~Z7:=..$NDZ?I:
	    .?D~...++,..~:IO7,??.OMM
	     .,...=?,....,7$...7ON8
	      ,?.. 7I,...:I7=,IDMD
	     .=?.  Z$:...:+=+?8MO
	   .?ZI,..=+,... ..IODMN
	   .7Z$?...   IZ$IZDMMN
	   .$I,.    ,ZZ8DNNNMI
	  .=$=. .. ON8+~7?..=
	  .==,. ,?DMD:.:$$~+?
	   Z~.. .ONN=..+$$I+?
	  .Z:.  =DD:..?Z??7I
	  .Z:,?8MM. .+8OOOO:
	  .$DDNMM   7DMMND
	  ...~~..    .ID

### Introduction ###

FuMa (Fusion Matcher) matches predicted fusion events (both genomic and transcriptomic) according to chromosomal location or assocatiated gene annotation(s) where the latter should be genome build inspecific.
The FuMa project currently supports input files from:

+	ChimeraScan<sup>[1]</sup>
+	DeFuse<sup>[2]</sup>
+	Illumina HiSeq<sup>[3]</sup>
+	Complete Genomics<sup>[4]</sup>
+	Tophat Fusion (both the candidates as the filtered results)<sup>[5]</sup>
+	(Trinity -->) GMAP<sup>[5]</sup>

Because RNA-Sequencing deals with samples that may have undergrond splicing, reads may split up because of biological processes. If a fusion event takes place, the same thing may happen. Therefore we hypothesize that using spanning read distances may be unreliable, because there are known introns of > 100kb. Therefore, FuMa translates the breakpoint to gene names, and only overlaps breakpoints with the same genename(s).


### Installation ###
Currently no installation is required; make sure the python (2) files can be called from the commandline.

### Usage ###
The commandline usage of FuMa is:

	usage: overlay_fusions.py [-h]
	                          [-a [ADD_GENE_ANNOTATION [ADD_GENE_ANNOTATION ...]]]
	                          -s ADD_SAMPLE [ADD_SAMPLE ...]
	                          [-l [LINK_SAMPLE_TO_ANNOTATION [LINK_SAMPLE_TO_ANNOTATION ...]]]
	                          [-f {summary,extensive}] [-o OUTPUT]

##### Example #####
Given a working directory with the following _tree_ structure:

	├── fuma -> ../../../fuma_working_directory
	├── input
	│   ├── dna
	│   │   └── complete_genomics
	│   │       ├── highConfidenceJunctionsBeta-GS000012345-ASM-N1.tsv
	│   │       └── highConfidenceJunctionsBeta-GS000012345-ASM-T1.tsv
	│   ├── references -> ../../../references/
	│   │   ├── refseq_genes_hg18.bed
	│   │   └── refseq_genes_hg19.bed
	│   └── rna
	│       ├── chimerascan
	│       │   └── chimeras.bedpe
	│       ├── defuse
	│       │   ├── results.filtered.tsv
	│       ├── tophat_fusion_post
	│       │   ├── potential_fusion.txt
	│       └── trinity_gmap
	│           └── Trinity.transloc
	└── results

We call fuma as follows:

	./fuma/fuma.py \
	    -a  "hg18:input/references/refseq_genes_hg18.bed" \
	        "hg19:input/references/refseq_genes_hg19.bed" \
	    \
	    -s  "tophat_fusion_post:tophatfusionpost:input/rna/tophat_fusion_post/potential_fusion.txt" \
	    \
	        "chimerascan:chimerascan:input/rna/chimerascan/chimeras.bedpe" \
	        \
	        "defuse:defuse:input/rna/defuse/results.tsv" \
	        \
	        "complete_genomics_normal:CompleteGenomics:input/dna/complete_genomics/allJunctionsBeta-GS000012345-ASM-N1.tsv" \
	        "complete_genomics_tumor:CompleteGenomics:input/dna/complete_genomics/allJunctionsBeta-GS000012345-ASM-T1.tsv" \
	        \
	        "trinity:trinitygmap:input/rna/trinity_gmap/Trinity.transloc" \
	    \
	    -l  "tophat_fusion_post:hg18" \
	        "chimerascan:hg18" \
	        "defuse:hg18" \
	        "complete_genomics_normal:hg18" \
	        "complete_genomics_tumor:hg18" 
	        "trinity:hg19" \
	    \
	    -o  "results/"

##### Output (summary) #####

### References ###
<sup>[1]</sup> ...

<sup>[2]</sup> ....

<sup>[3]</sup> ...

<sup>[4]</sup> ....

<sup>[5]</sup> ...

<sup>[6]</sup> ...

