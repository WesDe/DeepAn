# SVAn
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Description
SVAn aims at annotating insertion variants in a given vcf file. Insertion calls must be sequence-resolved (precise location and resolved inserted sequence). Insertion type classification is mainly based on Dbvar and Sequence Ontology insertion descriptions.
Five insertion types can be annotated : novel sequence, mobile element, tandem repeat, tandem duplication and dispersed duplication. Insertions that can not be associated to a type are labelled as unassigned.
It also annotates for each insertion the repeat context of the insertion site and the junctional homology size at the breakpoint.
## Requirements 
- Python3
- biopython, numpy python module
- Blat (for linux : http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat )
- TRF (https://tandem.bu.edu/trf/trf.download.html)
- Dfam and the hmm model associated to the specie studied (https://www.dfam.org ,see dfamm hmm model)
- hmmer (see http://hmmer.org/documentation.html)

## Installation
    git clone https://github.com/WesDe/DeepAn.git

## How to use :
 
### VCF annotation :
    ppython3 SVAn.py -v your_vcf.vcf -b /path/to/blat_executable -t /path/to/trf_executable -m /path/to/dfamscan.pl  -h /path/to/homo_sapiens_dfam.hmm -r /path/to/reference_genome.fa -h threshold -c NUM_CPU
