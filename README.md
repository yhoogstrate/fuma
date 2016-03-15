# FuMa (Fusion Matcher) [![Build Status](https://travis-ci.org/yhoogstrate/fuma.svg?branch=master)](https://travis-ci.org/yhoogstrate/fuma) #

A manuscript describing how FuMa works has been published in Oxford Bioinformatics at [http://dx.doi.org/10.1093/bioinformatics/btv721](http://dx.doi.org/10.1093/bioinformatics/btv721). Reference:

```
Youri Hoogstrate, René Böttcher, Saskia Hiltemann, Peter J. van der Spek, Guido Jenster, and Andrew P. Stubbs
FuMa: reporting overlap in RNA-seq detected fusion genes
Bioinformatics first published online December 10, 2015
```

- [Introduction](#introduction)
- [Technical Implementation](#technical-implementation)
    - [Overlap-based matching](#overlap-based-matching)
    - [Subset-based matching](#subset-based-matching)
    - [Exact gene set matching (EGM)](#exact-gene-set-matching-egm)
    - [Differences between matching types](#differences-between-matching-types)
        - [Example 1: long genes](#example-1-long-genes)
        - [Example 2: set expansion and shrinkage](#example-2-set-expansion-and-shrinkage)
- [Installation](#installation)
    - [Ubuntu](#ubuntu)
    - [Galaxy](#galaxy)
- [Usage](#usage)
    - [Command line](#command-line)
         - [-a ADD_GENE_ANNOTATION](#-a-add_gene_annotation)
         - [Obtain BED file -> fuma-gencode-gtf-to-bed](#obtain-bed-file---fuma-gencode-gtf-to-bed)
         - [-s ADD_SAMPLE](#-s-add_sample)
         - [-l LINK_SAMPLE_TO_ANNOTATION](#-l-link_sample_to_annotation)
         - [-m MATCHING_METHOD](#-m-matching_method)
         - [-f OUTPUT_FORMAT](#-f-output_format)
         - [--strand-specific-matching](#--strand-specific-matching)
         - [--acceptor-donor-order-specific-matching](#--acceptor-donor-order-specific-matching)
         - [Input formats](#input-formats)
         - [--verbose](#--verbose)
    - [Galaxy](#galaxy-1)
- [Examples](#examples)
    - [Example 01: one sample, two tools](#example-01-one-sample-two-tools)
    - [Example 02: one sample, one tool, different reference genomes](#example-02-one-sample-one-tool-different-reference-genomes)
    - [Example 03: Edgren dataset as part of Chimera supplement](#example-03-edgren-dataset-as-part-of-chimera-supplement)
- [References](#references)

## Introduction ##
This is the Manual as part of the Supplementary Material that belongs to the manuscript: *FuMa: reporting overlap in RNA-seq detected fusion genes* (under submission). FuMa (Fusion Matcher) matches predicted fusion events (both genomic and transcriptomic) according to chromosomal location and corresponding annotated genes. It is the organisation of the transcriptome (provided by the user) that forms the basis for FuMa to consider fusion genes to be identical or not. The provided gene annotation can be adjusted to define the biological question. For example, if it is desired to only consider fusion events that occur within exons, FuMa can be provided a list of such regions instead of entire genes. Currently FuMa supports input files from:

+	Chimera (Beccuti et al., 2014)
+	ChimeraScan (Iyer et al., 2011)
+	CompleteGenomics (Carnevali et al., 2012)
+	DeFuse (McPherson et al., 2011)
+	EricScript (Benelli et al., 2012)
+	FusionCatcher (Nicorici et al., 2014)
+	FusionMap (Ge et al., 2011)
+	GMAP (Wu and Watanabe, 2005)
+	JAFFA (Davidson et al., 2015)
+	STAR (Dobin et al., 2013)
+	STAR Fusion ([https://github.com/STAR-Fusion/STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion))
+	TopHat-Fusion (Kim and Salzberg, 2011)

## Technical Implementation ##
Matching fusion genes based on the genomic location shows limited accuracy. Therefore it is more convenient to use the gene names overlapping the breakpoints instead. Since ∼10% of the annotated human genes are overlapping (Sanna et al., 2008), and more genes and transcripts are being discovered by the RNA-Seq technology, breakpoints frequently span multiple genes. This complicates matching based on gene names and to account for that, matching two fusion genes in FuMa is achieved using set-theory based matching (overlap or subset). First both genomic partners of a fusion event are annotated with overlapping gene(s). Then fusion genes will be matches based on the opverlapping genes (using set-theory). The overlap- and subset matching approach have the advantage over the more stringent exact gene matching (EGM) approach that a certain level of overlapping genes are considered as acceptable. They behave quite similar but have features that require a more detailed explanation.

### Overlap-based matching ###
Overlap based matching is the default matching scheme of FuMa. It considers two fusion genes identical if both the genes sets, the left and the right, have at least one overlapping gene in common. We provide a more detailed description ([Fig. S1a: Example of overlap matching](#fig-s1a-example-of-overlap-matching)) and outline a corresponding truth-table. This scheme is less stringent than matching using subset based matching and has a few noteworthy characteristics:

 - **Long genes.** Long genes may span more other genes by chance. Therefore, two distant fusion genes that, by chance, also fall in the same long gene, may be matched only because they both overlap this same long gene. (See [example 1](#example-1-long-genes))
 - **Set expansion or set shrinkage.** When two (input) fusion genes match, the matched fusion gene has to have annotated genes based on the gene sets of the two (input) fusion genes. There are two sets that make sense to return; the intersect (all genes that must be present in both fusion genes) or the union (all genes, that must be present in at least one of them). When we use the union, those genes that are present in only one or in both gene sets, we introduce a problem we refer to as *set expansion*, which will result in an outcome that is dependent on the order of matching and on the iteration depth. This is very undesirable behavior and therefore FuMa returns the intersect instead. But the intersect of two gene sets may create a gene sets that are smaller than both gene sets initially used for matching. We refer to this as *set shrinkage*. For example, if set (*GREEN*, *BLUE*) is being matched with (*BLUE*, *RED*), the set of overlapping genes will be (*BLUE*). This is different from the subset method, because there the smallest initial gene set is being returned, since that's the set shared by both fusion genes. Therefore the gene sets in the subset method will never become smaller than the gene set of the smallest input gene set, while for the overlap based method the matched subset is not neccesairily equal to any set observed at the breakpoints. (See [example 2](#2-set-expansion-and-shrinkage))

![Fig. S1: Subset matching methodology](https://github.com/yhoogstrate/fuma/raw/master/share/Fig_S1.png)
###### Fig. S1a: Example of overlap matching ######
*Both scenario's (left and right) illustrate two predicted fusion genes, Fusion #1 and Fusion #2. Both have the same right location (red dashed line through the yellow gene), located in one single gene annotation, the yellow gene. Fusion #1 has two annotated genes on its left location: the green and the blue gene. In the right scenario, Fusion #2 is located in the blue and purple gene while in the left scenario it is only located within the blue gene. In the left scenario, the two fusions are considered identical because the left gene set of Fusion #2 (blue) overlaps the left gene set of Fusion #1 (blue and green). Also in the right scenario, the left gene sets (purple, blue) and (green, blue) are overlaping and the fusion genes are therefore considered to be identical, but the set is reduced to (blue) since that's the part that overlaps. The corresponding table of FuMa's overlap based matching strategy is given below [Table S1a: Overlap-based truth table](#table-s1a-overlap-based-truth-table). Depending on the genes spanning the breakpoints (first four columns), FuMa determines whether the fusion genes match (fifth column). The first four columns represent the gene sets (delimited with a comma) spanning the left and right locations. These gene names correspond to the colors used in figure above. The 5th column indicates whether FuMa considers the two fusions a match or not. The 6th and 7th columns represent the gene sets of the merged fusion gene as result of matching Fusion #1 and #2. The first examples matches because (blue) overlaps (blue, green), the second example matches because (blue, purple) and (blue, green) have blue in common.*

###### Table S1a: Overlap-based truth table ######

| Fusion #1    |        | Fusion #2   |        | Returning Match |      |        |
|--------------|--------|-------------|--------|-----------------|------|--------|
| Left         | Right  | Left        | Right  | Match           | Left | Right  |
| Blue         | Yellow | Blue, Green | Yellow | True            | Blue | Yellow |
| Blue, Purple | Yellow | Blue, Green | Yellow | True            | Blue | Yellow |

### Subset-based matching ###
The subset matching approach of FuMa considers two fusion genes identical if one of the left gene sets is a subset of the other left gene set, and one of the right gene sets is a subset of the other right gene set. Consequently for both the left and the right gene set, the intersect (subset) will be returned to the matched fusion gene. To illustrate how the subset matching methodology works, we give an example ([Fig. S1b: Example of the subset matching methodology](#fig-s1b-example-of-the-subset-matching-methodology)) and outline the corresponding truth-table.

![Fig. S1: Subset matching methodology](https://github.com/yhoogstrate/fuma/raw/master/share/Fig_S1.png)
###### Fig. S1b: Example of the subset matching methodology ######
*Both scenario's (left and right) illustrate two predicted fusion genes. In addition, both fusion genes have the same right location (red dashed line), located in one single gene annotation, the yellow gene. Also, Fusion #1 has two annotated genes on its left location: the green- and the blue gene. In the right scenario, Fusion #2 is located in the blue- and purple gene while in the left scenario it is only located within the blue gene. Therefore, in the left scenario, the two fusions are considered identical because the left gene set of Fusion #2 (blue) is a subset of the left gene set of Fusion #1 (blue and green). In the right scenario, the left gene sets (purple, blue) and (green, blue) are no subsets of each other and the fusion genes are therefore considered as distinct fusion genes. The corresponding table of FuMa's subset based matching strategy is given below [Table S1b: Subset-based truth table](#table-s1b-subset-based-truth-table). Depending on the genes spanning the breakpoints (first four columns), FuMa determines whether the fusion genes match (fifth column). The first four columns represent the gene sets (delimited with a comma) spanning the left and right locations. These gene names correspond to the colors used in the figure above. The 5th column indicates whether FuMa considers the two fusions a match or not. The 6th and 7th columns represent the gene sets of the merged fusion gene as result of matching Fusion #1 and #2. The first examples matches because (blue) is a valid subset of (blue, green) while the second example does not match because the left gene sets contain either (purple) or (green) which are mutually exclusive.*

###### Table S1b: Subset-based truth table ######

| Fusion #1    |        | Fusion #2   |        | Returning Match |      |        |
|--------------|--------|-------------|--------|-----------------|------|--------|
| Left         | Right  | Left        | Right  | Match           | Left | Right  |
| Blue         | Yellow | Blue, Green | Yellow | True            | Blue | Yellow |
| Blue, Purple | Yellow | Blue, Green | Yellow | False           |      |        |

### Exact gene set matching (EGM) ###

EGM consider fusion genes to be identical if their left and right gene sets are exactly identical. This is the most stringt matching scheme.

### Differences between matching types ###

The matching schemes have different noteworthy characteristics outlined in the following sections.

#### Example 1: long genes ####

	    f1                     f2
	    |                      |
	[ gene-A ]             [ gene-B ]
	[---------- long gene ----------]

In the illustrated example situation above, fusion genes *f1* and *f2* shall be matched using the overlap approach, since they both overlap *long gene*. In the case long gene is a really huge gene, it may span many other genes. Any fusion annotated upon this very long gene will in the overlap based matching be considered a match with any other fusion gene annotated within the long gene. When the subset matching was used, they would not have been considered a match, since (*gene-A*, *long gene*) is not a subset of (*gene-B*, *long gene*).

#### Example 2: set expansion and shrinkage ####

When the overlap based matching is used and consideres two fusion genes a match, a consensus left- and right gene set has to be returned for the merged fusion gene. There are two sets that can practically be returned, but both have some characteristics that are worthwile to mention.

#### Set shrinkage ####

The priciple of *set shrinkage* occurs when the returning gene set contain is the intersect of the two sets; contains only those genes that overlap. Consider two example fusion genes that have the following gene sets:

	Fusion1: GeneA, GeneB, GeneC
	                  |      |
	Fusion2:        GeneB, GeneC, GeneD, GeneE

The fusion genes are considered to be a match and the merged fusion gene should contain a new gene set. The intersect of the gene sets of *Fusion1* and *Fusion2* is (*GeneB*, *GeneC*). Hence, genes *GeneA*, *GeneD* and *GeneE* are taken out of the merged fusion.

When we continue matching with e.g. *Fusion3*:

	Fusion1,2*:     GeneB, GeneC
	                  |
	Fusion3:        GeneB,        GeneD

Both fusion genes have only *GeneB* in common, and the merged fusion gene will thus only contain *GeneB*. So *GeneC* is now also lost, although it was present in *Fusion1* and *Fusion2*. *GeneB* is the only gene shared in all three fusion genes, but it may be important to know that *GeneC* was shared in two of fusion genes. This information is lost because of the nature of the overlap matching approach in combination with returning the intersect. We refer to this as the set shrinkage issue. Note that the intersect is the implemented method for overlap based matching. When the subset approach was used instead, Fusion3 would not have been considered a match with the merged fusion gene *Fusion1,2\**.

#### Set expansion ####

**This section illustratates a methodology and this is not actually implemented in FuMa.**

When a merged fusion gene would contain the union of the genes, we would encounter a so called set expansion which will introduce order and iteration depentent results. To illustrate the problem of set expansion, imagine the following breakpoints:

 1. b1 = ```(A,A')```
 2. b2 = ```(A,A'')```
 3. b3 = ```(A'',B)```

To visualize such situation, we are most likely dealing with an annotation similar to this:

	       b1      b2       b3     
	       |       |        |      
	[---A'---]     |        |      
	     [-----A-----]      |      
	             [-----A''-----]   
	                      [---B---]

When we match these three breakpoints using the overlap-based method that returns any of the genes involved in any fusion gene, the results will become dependent on the order of matching and on the iteration depth. For this example we denote the following possible orders of matching:

 1. ```(b1 & b2) & b3```
 2. ```(b1 & b3) & b2```
 3. ```(b2 & b3) & b1```
 
When we match in **order 1**, we observe the following:

 1. Iteration 1: 
    - ```(A,A') & (A,A'') -> (b1 & b2) = (A,A',A'')*```
 2. Iteration 2:
    - ```(A,A',A'')* & (A'',B) -> (b1 & b2 & b3) = (A,A',A'',B)```

When we match in **order 2**, we observe the following:

 1. Iteration 1:
    - ```(A,A') & (A'',B) -> ``` no match; b1 and b3 are not considered to be identical
 
When we match in **order 3**, we observe the following:

 1. Iteration 1:
    - ```(A,A'') & (A'',B) -> (b2 & b3) = (A,A'',B)*```
 2. Iteration 2:
    - ```(A,A'',B)* & (A,A') -> (b1 & b2 & b3) = (A,A',A'',B)```

This illustrates that *b1* and *b3* are considered identical in *order 1* and *order 3*, but not in *order 2*. 

The second problem we encounter is that the gene sets have become larger. Before matching, the gene sets all had a size of 2 genes, after the first iteration the size of the matches were 3 genes and after the second iteration the size of the genes sets have become 4 genes. Therefore, the merged fusion gene can be matched with more fusion genes than each of the input fusion genes themselves. Therefore it is not a convenient strategy to return the entire set of genes.

## Installation ##
### Debian, Ubuntu and derivatives ###
FuMa requires Python 2.7, depends on HTSeq and can be obtained via git. We recommand the following commands to install FuMa (on Ubuntu and Debian derivate systems):

	sudo apt-get install build-essential python-dev git python-pip
	sudo pip uninstall fuma
	
	git clone https://github.com/yhoogstrate/fuma.git
	
	cd fuma
	
	python setup.py build
	python setup.py test
	sudo python setup.py install
	
	fuma --version

### Galaxy ###
Because usage of FuMa via the command line can be experienced as complicated, we also provide FuMa as Galaxy tool (Goecks
et al., 2010; Blankenberg et al., 2010; Giardine et al., 2005). The toolshed repository is in which FuMa is available is:

[https://toolshed.g2.bx.psu.edu/view/yhoogstrate/fuma](https://toolshed.g2.bx.psu.edu/view/yhoogstrate/fuma)

To install FuMa via Galaxy, you have to make sure you have the main toolshed [https://toolshed.g2.bx.psu.edu/](https://toolshed.g2.bx.psu.edu/) is configured in the servers tool_sheds_conf.xml. To install FuMa within galaxy, follow the procedure via the galaxy admin panel. We have made FuMa publicly available at the following galaxy instance:

[https://bioinf-galaxian.erasmusmc.nl/galaxy/](https://bioinf-galaxian.erasmusmc.nl/galaxy/)

We have made the example data available as shared data library at the following url:

[https://bioinf-galaxian.erasmusmc.nl/galaxy/library/list#folders/F313c46a90355d6dd](https://bioinf-galaxian.erasmusmc.nl/galaxy/library/list#folders/F313c46a90355d6dd)

## Usage ##
### Command line ###
To run FuMa via the command line, each dataset should be given as a separate file. Similarly, the corresponding gene annotation has to be linked to each dataset. Similarly, the file format has to be specified for each input dataset. This is a rather complex information structure and therefore, unfortunately, the command line arguments may be experienced as complicated. The command line usage of FuMa is:

	usage: fuma [-h] [-V] [--formats] [-m {overlap,subset,egm}]
	            [--strand-specific-matching] [--verbose]
	            [-a [ADD_GENE_ANNOTATION [ADD_GENE_ANNOTATION ...]]] -s ADD_SAMPLE
	            [ADD_SAMPLE ...]
	            [-l [LINK_SAMPLE_TO_ANNOTATION [LINK_SAMPLE_TO_ANNOTATION ...]]]
	            [-f {summary,list,extensive}] [-o OUTPUT]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -V, --version         show program's version number and exit
	  --formats             show accepted dataset formats
	  -m {overlap,subset,egm}, --matching-method {overlap,subset,egm}
	                        The used method to match two gene sets. Overlap
	                        matches when two gene set have one or more genes
	                        overlapping. Subset matches when one gene set is a
	                        subset of the other. EGM is exact gene matching; all
	                        genes in both sets need to be identical to match.
	  --strand-specific-matching
	                        Take strand specificness into account (5' -> 3' ? 3'
	                        -> 5')
	  --verbose             increase output verbosity
	  -a [ADD_GENE_ANNOTATION [ADD_GENE_ANNOTATION ...]], --add-gene-annotation [ADD_GENE_ANNOTATION [ADD_GENE_ANNOTATION ...]]
	                        annotation_alias:filename * file in BED format
	  -s ADD_SAMPLE [ADD_SAMPLE ...], --add-sample ADD_SAMPLE [ADD_SAMPLE ...]
	                        sample_alias:format:filename (available formats: fuma
	                        --formats)
	  -l [LINK_SAMPLE_TO_ANNOTATION [LINK_SAMPLE_TO_ANNOTATION ...]], --link-sample-to-annotation [LINK_SAMPLE_TO_ANNOTATION [LINK_SAMPLE_TO_ANNOTATION ...]]
	                        sample_alias:annotation_alias
	  -f {summary,list,extensive}, --format {summary,list,extensive}
	                        Output-format
	  -o OUTPUT, --output OUTPUT
	                        output filename; '-' for stdout
	
	
	For more info please visit:
	<https://github.com/yhoogstrate/fuma>


#### -a ADD_GENE_ANNOTATION ####
Gene annotations have to be provided in a tab-delimited file, with the first column containing the genes chromosome, the second and the third column the (1-based) start and end position, and the fourth column the (unique) gene identifier or name, as shown in the example below:

	chr1   100000000  120000000  GeneNameA
	chr2   100000000  120000000  GeneNameB
	chr21  100000000  120000000  GeneNameC
	chr22  100000000  120000000  GeneNameD
	chrX   140000000  160000000  GeneNameX
	chrY   140000000  160000000  GeneNameY

This format is compatible with the BED format [https://genome.ucsc.edu/FAQ/FAQformat.html#format1](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), but requires that the 4th column is present and requires it to contain unique gene names. Additional columns are allowed, but are nowhere taken into account. **Do not provide BED files that describe one exon per line** because this will exclude the introns, but provide BED files that describe one gene per line instead. For files with one exon per line, we can not merge  exons into genes because when they are merged on the basis of the gene names, duplicates on the same chromosome that span a large distance may introduce overlap and large uncertainty.

In contrast, if you explicitly want to match only in exon regions, you should use BED files with one exon per line. In that case is advised to provide non-unique gene names, like the following example:

	chr1  100000000  100001000  GeneNameA
	chr1  100002000  100003000  GeneNameA
	chr1  100005000  100006000  GeneNameA
	chr2  100000000  100100000  GeneNameB
	chr2  100101000  100103000  GeneNameB

In FuMa the gene annotation argument is provided as unique alias followed by the filename, separated with a colon:

	-a "hg19:somefile.bed"

In this case the alias of the BED-file, hg19, will later be used to link it to datasets. In case you want multiple references, you can provide arguments delimited with whitespaces:

	-a "hg18:somefile_hg18.bed" "hg19:somefile_hg19.bed"

#### Obtain BED file -> fuma-gencode-gtf-to-bed ####

Because obtaining such files turns out to more difficult than expected, we have provided an extra utility named `fuma-gencode-gtf-to-bed`.
The user should start with download a GTF file from (at least tested with) GenCode. Then user should proceed with running the following command:

	fuma-gencode-gtf-to-bed -o converted.bed input.gtf

The utility will use all annotations in the GTF file and will aggregate all exons per `transcript_id`, while it will use the gene_id as unique identifier in the BED file. The reason for this is that if transcripts that belong to the same gene while they are quite distant to each other (or homologues using the same name, which happens), they will be annotated per transcript such that the long distance between the transcripts will not unneccesairily be marked as part of that gene. In case multiple transcripts from the same gene are annotated upon each other, FuMa will treat them as the same gene as long as their identifier is the same, which is the case since the `gene_id` is being used for this.

This tool should work for all GTF files for which all entries have a proper and uniquely wise correct definition of the `gene_id` and `transcript_id`.

#### -s ADD_SAMPLE  ####
To provide FuMa a fusion gene detection experiment, it should be provided with the "-s" argument which should follow the following syntax:

*sample_alias*:*format*:*filename*

The *sample_alias* will be used for two things: (1) as column header and alias in the final output and (2) to link the references to the samples. The format is the file format in which the fusion genes are described. Note that some tools have multiple output formats. These are usually the file formats for interim output files.

#### -l LINK_SAMPLE_TO_ANNOTATION ####

Each dataset must be annotated with only one gene annotation. This can be achieved using the following argument syntax:

*sample_alias*:*annotation_alias*

In case you have a particular same *s* and a reference *ref*, you can link *s* to *ref* as follows:

	-l "s:ref"

In case you have two samples, one on *ref1* and one on *ref2*, you can provide it as follows:

	-l "defuse_hg18:hg18" "chimerascan_hg19:hg19"

#### -l MATCHING_METHOD ####

FuMa has the option to use three methods to match fusion genes; 'overlap', 'subset' and 'egm' (default is 'overlap'). These method can be selected with the ```-m``` or ```--matching-method```, argument as follows:

	fuma -m egm [ ... ]

	fuma --matching-method subset [ ... ]

#### -f OUTPUT_FORMAT ####

FuMa has the built-in option for several output formats. The most straight-forward format is the '*list*' output format which contains per (matched) fusion gene, for each matching tool, the genomic locations and identifier(s) or an empty column if the tool didn't pick it up. In the following example we have three fusion genes; one detected by TopHat fusion, one by STAR and one by both. The corresponding output in '*list*' format would be something like:

| Left Genes | Right Genes | STAR             | TopHat Fusion
|:-----------|:------------|:-----------------|:-------------
| FOO1       | BAR1        | UID_A=chr1:12-34 |
| FOO2       | BAR2        |                  | TID_A=chr4:66-77
| DOX1       | BOX5        | UID_B=chr5:85-95 | TID_B=chr5:88-99

Occasionally tools predict multiple fusion events within the same left- and right genes, which FuMa will consider as duplicates. In case we observe a duplicate, we simply provide both identifiers delimited with a comma into one cell, such that duplicate entries can always be traced back in the output:

| Left Genes | Right Genes | FusionMap
|:-----------|:------------|:---------
| FOO1       | BAR1        | UID_A=chr1:12-34,UID_B=chr1:12-34

When a breakpoint location spans multiple gene annotations, the genes in the column are delimited with a colon:

| Left Genes | Right Genes | OncoFuse
|:-----------|:------------|:---------
| FOO1:FOO2  | BAR1        | UID_A=chr1:12-34

The Galaxy wrapper has the option to replace the columns to TRUE or FALSE depending on whether a match was found or not.

The output format '*extensive*' is file format similar to the format Complete Genomics provides (http://www.completegenomics.com/documents/DataFileFormats_Cancer_Pipeline_2.4.pdf from p135) and that only contains those fusion genes that have at least one match. This format is in particular useful if the output of one run needs to be (re-)used for another run.

The output format '*summary*' is a set of tables that contains the numbers of detected matches per dataset combination, useful for creating Venn diagrams.

#### --strand-specific-matching ####

FuMa has the built-in option to separate fusion genes based on the predicted strand of the acceptor or donor. In the following example we have fusion genes #1 and #2, with exactly the same breakpoints, but the transcripts of the second gene are predicted to have different strands.

	#1:
	        b1 (+) ->          <- (-) b2
	        |                         |
	[ --- Gene A --- ]        [ --- Gene B --- ]
	
	#2:
	        b1 (+) ->                 b2 (+) ->
	        |                         |
	[ --- Gene A --- ]        [ --- Gene B --- ]

To let FuMa consider these fusion as distinct fusion genes because of the different strands, the user has to enable strand specific matching by including the ```--strand-specific-matching``` argument:

	fuma \
	    --strand-specific-matching \
	    --acceptor-donor-order-specific-matching \
	    -a  "hg19:genes_hg19.bed" \
	    \
	    -s  "chimerascan:chimerascan:FOO_chimerascan/chimeras.bedpe" \
	        "defuse:defuse:FOO_defuse/results.tsv" \
	    -l  "chimerascan:hg19" \
	        "defuse:hg19" \
	    -f  "list" \
	    -o  "chimerascan_defuse_overlap.txt"

It is recommended to use this option together with the  **```--acceptor-donor-order-specific-matching```** option.

#### --acceptor-donor-order-specific-matching ####
The order in which the acceptor and donor gene are denoted is for certain tools determinant where the transcript started. This information may be crucial to explain the function and biological role of a fusion gene. For example, TMPRSS2-ERG, a fusion gene found in about 50% of all screened prostate cancers, uses regulatory elements from the androgen driven gene TMPRSS2, fused to the gene ERG that has an oncogenic role in human prostate cancer (Tomlins et. al, 2008). These principles would not apply if the order of these genes would be vice versa.

FuMa has the built-in option to separate fusion genes based on the order of the denotation of the acceptor or donor. In the following example we have fusion genes #1 and #2, with exactly the same breakpoints, but the order of the acceptor and donor gene has changed.

	#1:
	        break1                    break2
	        |                         |
	[ --- Gene A --- ]        [ --- Gene B --- ]
	
	#2:
	        break1                    break2
	        |                         |
	[ --- Gene B --- ]        [ --- Gene A --- ]

To let FuMa consider these fusion as distinct fusion genes because of the different order of the donor and acceptor, the user has to enable strand specific matching by including the ```--acceptor-donor-order-specific-matching``` argument:

	fuma \
	    --acceptor-donor-order-specific-matching \
	    -a  "hg19:genes_hg19.bed" \
	    \
	    -s  "chimerascan:chimerascan:FOO_chimerascan/chimeras.bedpe" \
	        "defuse:defuse:FOO_defuse/results.tsv" \
	    -l  "chimerascan:hg19" \
	        "defuse:hg19" \
	    -f  "list" \
	    -o  "chimerascan_defuse_overlap.txt"

**It is important to state that some file formats (interim output and discordant reads) do not take this information into account.**

#### Input formats ####

FuMa supports the following file formats:

| Tools              | File                  | Format string
|:-------------------|:----------------------|:-------------
| Chimera            | prettyPrint() output  | chimera
| ChimeraScan        | chimeras.bedpe        | chimerascan
| Complete Genomics  | highConfidenceJu*.tsv | complete-genomics
| Complete Genomics  | allJunctionsBeta*.tsv | complete-genomics
| DeFuse             | results.txt           | defuse
| DeFuse             | results.classify.txt  | defuse
| DeFuse             | results.filtered.txt  | defuse
| EricScript         | .results.total.txt    | ericscript *************
| Fusion Catcher     | final-list_cand*.txt  | fusion-catcher_final
| FusionMap          |                       | fusionmap
| JAFFA              | jaffa_results.cvs     | jaffa
| Trinity + GMAP     |                       | trinity-gmap
| OncoFuse           |                       | oncofuse
| RNA STAR           | Chimeric.out.junction | rna-star_chimeric
| SOAPFuse           | final.*.for.genes.txt | soapfuse-final-gene
| SOAPFuse           | final.*.for.trans.txt | soapfuse-final-transcript
| STAR Fusion        | _candidates.final     | star-fusion_final
| TopHat Fusion pre  | fusions.out           | tophat-fusion_pre
| TopHat Fusion post | potential_fusion.txt  | tophat-fusion_post_potential_fusion
| TopHat Fusion post | result.txt            | tophat-fusion_post_result
| TopHat Fusion post | result.html           | tophat-fusion_post_result_html

************* EricScript often contains entries with unknown breakpoints.
Because no genomic coordinates are given those fusion genes can not be
imported into FuMa and only those with breakpoints will be taken into account.

Or run the following command line argument to get an overview of the versions at the command line:

	fuma --formats

#### --verbose ####

If you would like to see additional statistics during runtime (or post-runtime
if you store the output) you should run FuMa with the `--verbose` argument:

	fuma \
	    -a  "hg19:genes_hg19.bed" \
	    \
	    -s  "chimerascan:chimerascan:FOO_chimerascan/chimeras.bedpe" \
	        "defuse:defuse:FOO_defuse/results.tsv" \
	    -l  "chimerascan:hg19" \
	        "defuse:hg19" \
	    -f  "list" \
	    -o  "chimerascan_defuse_overlap.txt" \
	    --verbose

This allows the user to inspect the numbers of duplicate fusions, the
number of parsed genes from the gene set and showing which datasets
are being compared at run time.

* Note: As of 2.12.1 this argument is required, in preliminary versions
this was by default enabled.

### Galaxy ###

After having FuMa installed in Galaxy via the toolshed, it can be opened by typing '*fuma*' in the '*search tools*' field on the left panel in galaxy. When it has opened, the interface should be similar to [Fig. S2: FuMa in Galaxy](#fig-s2-fuma-in-galaxy). The main input of the Galaxy wrapper is a set of datasets. You can as add many datasets as the server can handle in terms of resources. For each dataset the user needs to specify (1) the history item in galaxy that contains the output file of the fusion gene detection experiment, (2) the corresponding file format and name of the tool that corresponds to the history item and (3) a corresponding gene annotation file (in BED format). Lastly, the user can specify the desired output format and proceed with the analysis.

![Fig. S2: FuMa in Galaxy](https://github.com/yhoogstrate/fuma/raw/master/share/Fig_S2.png)

###### Fig. S2: FuMa in Galaxy #######

## Examples ##

### Example 01: one sample, two tools ###

Imagine we have run sample FOO with Defuse and ChimeraScan, on the same reference genome (hg19). The corresponding gene annotation on hg19 is genes_hg19.bed and the output should be stored in chimerascan_defuse_overlap.txt. The command line argument to run this analysis would be:

	fuma \
	    -a  "hg19:genes_hg19.bed" \
	    \
	    -s  "chimerascan:chimerascan:FOO_chimerascan/chimeras.bedpe" \
	        "defuse:defuse:FOO_defuse/results.tsv" \
	    -l  "chimerascan:hg19" \
	        "defuse:hg19" \
	    -f  "list" \
	    -o  "chimerascan_defuse_overlap.txt"

### Example 02: one sample, one tool, different reference genomes ###

When want to compare the differences between runs on different genome builds, we can add each runs and define a different gene annotation for each run. Imagine we have run a sample with TopHat-Fusion on reference genomes hg18 and hg19, we can run FuMa as follows:

	fuma \
	    -a  "hg18:genes_hg18.bed" \
	        "hg19:genes_hg19.bed" \
	    \
	    -s  "thf_hg18:Tophat-Fusion Post result:thf_hg18/result.txt" \
	        "thf_hg19:Tophat-Fusion Post result:thf_hg19/result.txt" \
	    -l  "thf_hg18:hg18" \
	        "thf_hg19:hg19" \
	    -f  "list" \
	    -o  "thf_hg18_hg19_overlap.txt"

It is important that the gene annotations genes_hg18.bed and genes_hg19.bed contain similar gene names, since matching is based on these names. Therefore it is recommanded to remove gene names that are specific per annotation; the latest genes only available in hg19 will never be matched with hg18 simply because they do not exist in hg18.

### Example 03: Edgren dataset as part of Chimera supplement ###

The publicly available data from the Edgren dataset has been performed on FusionMap, ChimeraScan and DeFuse as proof of concept data for the Chimera package (Edgren et al., 2011; Beccuti et al., 2014). To obtain these result you should run the following command:

	wget http://www.bioconductor.org/packages/release/bioc/src/contrib/chimera_1.10.0.tar.gz
	tar -xzf chimera_1.10.0.tar.gz

Within the source of the chimera package, you can find the files with the following command line command:

	find . -type f | grep -i -E "Edgr[e]{1,2}n"

Please check whether the output is identical to:

	./chimera/inst/examples/Edgreen_fm.txt
	./chimera/inst/examples/edgren.stat.detection.txt
	./chimera/inst/examples/Edgren_df.tsv
	./chimera/inst/examples/Edgren_cs.txt
	./chimera/inst/examples/Edgren_true.positives.txt

To get a gene reference and the True positivies with genomic coordinates, run at the command line:

	wget https://testtoolshed.g2.bx.psu.edu/repos/yhoogstrate/fuma/raw-file/tip/test-data/refseq_genes_hg19.bed
	wget https://testtoolshed.g2.bx.psu.edu/repos/yhoogstrate/fuma/raw-file/tip/test-data/edgren_tp.txt

We can proceed with FuMa by running at the command line:

	edir="./chimera/inst/examples/"
	fuma \
	    -a  "hg19:refseq_genes_hg19.bed" \
	    \
	    -s  "chimerascan:chimerascan:"$edir"Edgren_cs.txt" \
	        "defuse:defuse:"$edir"Edgren_df.tsv" \
	        "fusionmap:fusionmap:"$edir"Edgreen_fm.txt" \
	        "edgren_TP:fusionmap:edgren_tp.txt" \
	    -l  "fusionmap:hg19" \
	        "defuse:hg19" \
	        "chimerascan:hg19" \
	        "edgren_TP:hg19" \
	    -f  "list" \
	    -o  "edgren_fuma_list.txt"

To convert the columns to boolean values, we run:

	fuma-list-to-boolean-list \
	-o "edgren_fuma_booleanlist.txt" \
	   "edgren_fuma_list.txt"

To find all fusion genes present in 3 or more datasets, run at the command:

	grep -E "Left-genes|TRUE.*?TRUE.*?TRUE.*?" "edgren_fuma_booleanlist.txt"

This will return the following list of 20 fusion genes:

	NM_018837:NM_198596:NM_001161841	NM_006420	TRUE	TRUE	TRUE	TRUE

	Left-genes	Right-genes	chimerascan	defuse	fusionmap	edgren_TP
	TEKT4P2	BRWD1	TRUE	TRUE	TRUE	FALSE
	MED1	ACSF2	TRUE	TRUE	TRUE	FALSE
	BCAS3	MED13	TRUE	TRUE	TRUE	FALSE
	SUMF1	LRRFIP2	TRUE	TRUE	FALSE	TRUE
	CMTM7	GLB1	TRUE	TRUE	FALSE	TRUE
	NUP214	NOTCH1	TRUE	TRUE	FALSE	TRUE
	EIF3H	CYTH1	TRUE	TRUE	FALSE	TRUE
	SNF8	RPS6KB1	TRUE	TRUE	FALSE	TRUE
	BCAS3	BCAS4	TRUE	TRUE	FALSE	TRUE
	IKZF3	VAPB	TRUE	TRUE	FALSE	TRUE
	CEP250	ZMYND8	TRUE	TRUE	FALSE	TRUE
	TTI1	DIDO1	TRUE	TRUE	FALSE	TRUE
	BSG	NFIX	TRUE	TRUE	FALSE	TRUE
	MYO9B	RAB22A	TRUE	TRUE	FALSE	TRUE
	ANKHD1-EIF4EBP3:ANKHD1	PCDH1	TRUE	TRUE	TRUE	TRUE
	ACACA	STAC2	TRUE	TRUE	TRUE	TRUE
	MYO19	SKA2	TRUE	TRUE	TRUE	TRUE
	SULF2	ARFGEF2	TRUE	TRUE	TRUE	TRUE
	TATDN1	GSDMB	TRUE	TRUE	TRUE	TRUE
	PKIA	RARA	TRUE	TRUE	TRUE	TRUE

## References ##
- Beccuti, M., Carrara, M., Cordero, F., Lazzarato, F., Donatelli, S., Nadalin, F., Policriti, A., and Calogero, R. A. (2014). Chimera: a Bioconductor package for secondary analysis of fusion products. Bioinformatics (Oxford, England), 30(24), 3556--7.
- Benelli M, Pescucci C, Marseglia G, Severgnini M, Torricelli F, Magi A. Discovering chimeric transcripts in paired-end RNA-seq data by using EricScript. Bioinformatics. 2012; 28(24): 3232-3239.
- Blankenberg, D., Kuster, G. V., Coraor, N., Ananda, G., Lazarus, R., Mangan, M., Nekrutenko, A., and Taylor, J. (2010). Galaxy: A web-based genome analysis tool for experimentalists. Current protocols in molecular biology, pages 19--10.
- Carnevali, P., Baccash, J., Halpern, A. L., Nazarenko, I., Nilsen, G. B., Pant, K. P., Ebert, J. C., Brownley, A., Morenzoni, M., Karpinchyk, V., Martin, B., Ballinger, D. G., and Drmanac, R. (2012). Computational Techniques for Human Genome Resequencing Using Mated Gapped Reads.
- Davidson, M., Majewski, I., Oshlack, A. (2015). JAFFA: High sensitivity transcriptome-focused fusion gene detection. Genome Medicine, 7(1), 1-12.
- Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., and Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15--21.
- Edgren, H., Murumagi, A., Kangaspeska, S., Nicorici, D., Hongisto, V., Kleivi, K., Rye, I. H., Nyberg, S., Wolf, M., Borresen-Dale, A.-L., and Kallioniemi, O. (2011). Identification of fusion genes in breast cancer by paired-end RNA-sequencing. Genome biology, 12(1), R6.
- Ge, H., Liu, K., Juan, T., Fang, F., Newman, M., and Hoeck, W. (2011). Fusionmap: detecting fusion genes from next-generation sequencing data at base-pair resolution. Bioinformatics.
- Giardine, B., Riemer, C., Hardison, R. C., Burhans, R., Elnitski, L., Shah, P., Zhang, Y., Blankenberg, D., Albert, I., Taylor, J., Miller, W. C., Kent, W. J., and Nekrutenko, A. (2005). Galaxy: a platform for interactive large-scale genome analysis. Genome research, 15(10), 1451--1455.
- Goecks, J., Nekrutenko, A., Taylor, J., and Team, T. G. (2010). Galaxy: a comprehensive approach for supporting accessible, reproducible, and transparent computational research in the life sciences. Genome Biol, 11(8), R86.
- Iyer, M. K., Chinnaiyan, A. M., and Maher, C. A. (2011). Chimerascan: a tool for identifying chimeric transcription in sequencing data. Bioinformatics, 27(20), 2903--2904.
- Kim, D. and Salzberg, S. L. (2011). TopHat-Fusion: an algorithm for discovery of novel fusion transcripts. Genome biology, 12(8), R72.
- McPherson, A., Hormozdiari, F., Zayed, A., Giuliany, R., Ha, G., Sun, M. G. F., Griffith, M., Moussavi, A., Senz, J., Melnyk, N., Pacheco, M., Marra, M. A., Hirst, M., Nielsen, T. O., Sahinalp, S. C., Huntsman, D., and Shah, S. P. (2011). Defuse: An algorithm for gene fusion discovery in tumor rna-seq data. PLoS Computational Biology, 7.
- Nicorici, D., Satalan, M., Edgren, H., Kangaspeska, S., Murumagi, A., Kallioniemi, O., Virtanen, S., and Kilkku, O. (2014). Fusioncatcher - a tool for finding somatic fusion genes in paired-end rna-sequencing data. Technical report.
- Sanna, C. R., Li, W.-H., and Zhang, L. (2008). Overlapping genes in the human and mouse genomes. BMC genomics, 9, 169.
- Tomlins, S. A., Laxman, B., Varambally, S., Cao, X., Yu, J., Helgeson, B. E., Cao, Q., Prensner, J. R., Rubin, M. A., Shah, R. B., Mehra, R., and Chinnaiyan, A. M. (2008). Role of
the tmprss2-erg gene fusion in prostate cancer. Neoplasia, 10(2), 177--188.
- Wu, T. D. and Watanabe, C. K. (2005). GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics (Oxford, England), 21(9), 1859--75.
