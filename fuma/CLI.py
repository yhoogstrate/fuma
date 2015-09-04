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
Fusion Catcher     | final-list_cand*.txt  | fusion-catcher_final
FusionMap          |                       | fusionmap
Trinity + GMAP     |                       | trinity-gmap
OncoFuse           |                       | oncofuse
RNA STAR           | Chimeric.out.junction | rna-star_chimeric
STAR Fusion        | _candidates.final     | star-fusion_final
TopHat Fusion pre  | fusions.out           | tophat-fusion_pre
TopHat Fusion post | potential_fusion.txt  | tophat-fusion_post_potential_fusion
TopHat Fusion post | result.txt            | tophat-fusion_post_result
 
 The file formats that are supported in the direction (5' -> 3')
 specific mode are:
 
- chimerascan
- defuse
- fusion-catcher_final
- tophat-fusion_pre
- tophat-fusion_post_potential_fusion
- rna-star_chimeric
"""

def CLI():
	"""Command Line Interface
	"""
	parser = argparse.ArgumentParser()
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/fuma>")
	parser.add_argument('-V','--version', action='version', version=textwrap.dedent("%(prog)s "+fuma.__version__+"\n\nCopyright (C) 2013-"+str(datetime.datetime.now().year)+" Youri Hoogstrate.\n\nLicense GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"))
	parser.add_argument('--formats', action='version', version=show_formats(), help="show accepted dataset formats")
	parser.add_argument('--egm', action='store_const', const=True, default=False, help='Exact Gene-list Matching approach (not recommended)')
	parser.add_argument('--overlap-based-matching', action='store_const', const=True, default=False, help='Overlap based matching (experimental)')
	
	parser.add_argument('--strand-specific-matching', action='store_const', const=True, default=False, help='Take strand specificness into account (5\' -> 3\' ? 3\' -> 5\')')
	
	parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
	
	parser.add_argument("-a","--add-gene-annotation",help="annotation_alias:filename  * file in BED format",nargs="*")
	
	parser.add_argument("-s","--add-sample",nargs="+",required=True,help="sample_alias:format:filename (available formats: %(prog)s --formats)")
	parser.add_argument("-l","--link-sample-to-annotation",help="sample_alias:annotation_alias",nargs="*")
	
	parser.add_argument("-f","--format",default="list",choices=["summary","list","extensive"],help="Output-format")
	
	parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="overlap/")
	
	return parser.parse_args()
