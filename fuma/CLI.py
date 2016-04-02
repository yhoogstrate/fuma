#!/usr/bin/env python

"""[License: GNU General Public License v3 (GPLv3)]
 
 This file is part of FuMa.
 
 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 FuMa is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 Documentation as defined by:
 <http://epydoc.sourceforge.net/manual-fields.html#fields-synonyms>
"""

import fuma
import sys,argparse,textwrap,datetime


def show_formats():
	return """FuMa supports the following file formats:

Tools              | File                  | Format string
----------------------------------------------------------
Chimera            | prettyPrint() output  | chimera
ChimeraScan        | chimeras.bedpe        | chimerascan
Complete Genomics  | highConfidenceJu*.tsv | complete-genomics
Complete Genomics  | allJunctionsBeta*.tsv | complete-genomics
DeFuse             | results.txt           | defuse
DeFuse             | results.classify.txt  | defuse
DeFuse             | results.filtered.txt  | defuse
EricScript         | .results.total.txt    | ericscript *************
Fusion Catcher     | final-list_cand*.txt  | fusion-catcher_final
FusionMap          |                       | fusionmap
JAFFA              | jaffa_results.cvs     | jaffa
Trinity + GMAP     |                       | trinity-gmap
OncoFuse           |                       | oncofuse
RNA STAR           | Chimeric.out.junction | rna-star_chimeric
SOAPFuse           | final.*.for.genes.txt | soapfuse-final-gene
SOAPFuse           | final.*.for.trans.txt | soapfuse-final-transcript
STAR Fusion        | _candidates.final     | star-fusion_final
TopHat Fusion pre  | fusions.out           | tophat-fusion_pre
TopHat Fusion post | potential_fusion.txt  | tophat-fusion_post_potential_fusion
TopHat Fusion post | result.txt            | tophat-fusion_post_result
TopHat Fusion post | result.html           | tophat-fusion_post_result_html
 
 The file formats that are supported in the direction (5' -> 3')
 specific mode are:
 
- chimerascan
- defuse
- fusion-catcher_final
- tophat-fusion_pre
- tophat-fusion_post_potential_fusion
- rna-star_chimeric
- soapfuse-final-gene
- soapfuse-final-transcript


************* EricScript often contains entries with unknown breakpoints.
Because no genomic coordinates are given those fusion genes can not be
imported into FuMa and only those with breakpoints will be taken into account.
"""

def CLI(argv=None):
	"""Command Line Interface
	
	base command line interface of FuMa
	"""
	parser = argparse.ArgumentParser()
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/fuma>")
	parser.add_argument('-V','--version', action='version', version=textwrap.dedent("%(prog)s "+fuma.__version__+"\n\nCopyright (C) 2013-"+str(datetime.datetime.now().year)+" Youri Hoogstrate.\n\nLicense GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"))
	parser.add_argument('--formats', action='version', version=show_formats(), help="show accepted dataset formats")
	
	parser.add_argument("-m",'--matching-method',choices=['overlap','subset','egm'],default='subset',help='The used method to match two gene sets. Overlap matches when two gene set have one or more genes overlapping. Subset matches when one gene set is a subset of the other. EGM is exact gene matching; all genes in both sets need to be identical to match.')
	
	parser.add_argument(   '--strand-specific-matching', dest='strand_specific_matching', action='store_true',  help='Consider fusion genes distinct  when the breakpoints have different strands: (A<-,B<-) != (->A,B<-); default')
	parser.add_argument('--no-strand-specific-matching', dest='strand_specific_matching', action='store_false', help='Consider fusion genes identical when the breakpoints have different strands: (A<-,B<-) == (->A,B<-)')
	parser.set_defaults(strand_specific_matching=True)
	
	parser.add_argument(   '--acceptor-donor-order-specific-matching', dest='acceptor_donor_order_specific_matching', action='store_true',  help='Consider fusion genes distinct  when the donor and acceptor sites are swapped: (A,B) != (B,A)')
	parser.add_argument('--no-acceptor-donor-order-specific-matching', dest='acceptor_donor_order_specific_matching', action='store_false', help='Consider fusion genes identical when the donor and acceptor sites are swapped: (A,B) == (B,A); default')
	parser.set_defaults(acceptor_donor_order_specific_matching=False)
	
	parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
	
	parser.add_argument("-a","--add-gene-annotation",help="annotation_alias:filename  * file in BED format",nargs="*")
	
	parser.add_argument("-s","--add-sample",nargs="+",required=True,help="sample_alias:format:filename (available formats: %(prog)s --formats)")
	parser.add_argument("-l","--link-sample-to-annotation",help="sample_alias:annotation_alias",nargs="*")
	
	parser.add_argument("-f","--format",default="list",choices=["summary","list","extensive"],help="Output-format")
	
	parser.add_argument("-g","--long-gene-size",default=200000,type=int,help="Gene-name based matching is more sensitive to long genes. This is the gene size used to mark fusion genes spanning a 'long gene' as reported the output. Use 0 to disable this feature.")
	
	parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="output_fuma.txt")
	
	if(argv == None):
		return parser.parse_args()
	else:
		# Argumented parameters are used in the unit tests.
		return parser.parse_args(argv)


def CLI_ensmble_gtf_to_bed_converter(argv=None):
	"""
		CLI for gtf to bed converter
	"""
	
	parser = argparse.ArgumentParser()
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/fuma>")
	parser.add_argument('-V','--version', action='version', version=textwrap.dedent("%(prog)s "+fuma.__version__+"\n\nCopyright (C) 2013-"+str(datetime.datetime.now().year)+" Youri Hoogstrate.\n\nLicense GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"))
	
	parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="-")
	parser.add_argument("genecode_gtf_file",nargs=1,help="Input GTF file, e.g. 'gencode_gtf_file.gtf' - not as .gz")
	
	if(argv == None):
		return parser.parse_args()
	else:
		# Argumented parameters are used in the unit tests.
		return parser.parse_args(argv)
