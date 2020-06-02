# DeepAn
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Description
DeepAn aims to annotate insertion in vcf file which are sequence resolved. It is mainly based on Dbvar and Sequence Ontology insertion description.
Five insertion types can be annotated : novel sequence, mobile element, tandem repeat, tandem duplication and dispersed duplication. Insertions that can't be associated to a type are identified as unassigned.
It allows to annotate repeat location of the insertion and the junctional homolgy size observed in every insertion.

## Requirements 
- Python3
- biopython module (pip3 install biopython)
- Blat (for linux : http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat )
- TRF (https://tandem.bu.edu/trf/trf.download.html)
- Dfam and the hmm model associated to the specie studied (https://www.dfam.org ,see dfamm hmm model) 

## Installation
    git clone https://github.com/WesDe/DeepAn.git

## How to use :

### Insertion annotation, step 1 : 
Verify that the vcf to annotate contains sequence resolved insertion and are described in ALT or info section

### Insertion annotation, step 2 :
Transform the vcf file in FASTA format with :

    python3 Conversion_vcf_fasta/Vcf_to_fa.py vcf_file.vcf outputname.fa

### Insertion annotation, step 3 :
Detection of potential tandem repeat :

    python3 TRF/TRF_ALT.py insertion_file.fa path_to_trf Potential_TRF.csv

Detection of potential mobile element :

sh Mobile_element/dfam.sh path_to_dfamm_executable vcf_file.fa hmm_model Potential_mobile_element.csv

Detection of potential duplication :

    python3 Blat_WG.py Query_Blat_inser_WG.py reference_genome.fa insertion_file.fa path_to_blat Potential_duplication.csv

Be careful this step may take multiple days to detect the ensemble of potential duplication in the whole genome.


### Junctional homology detection :
Detection of potential large homology :

    python3 homology/Query_blat_large_homology.py reference_genome.fasta vcf_file.vcf path_to_blat Potential_large_microhomology.psl

Detection of potential small homology :

    python3 homology/Annotation_microhomology.py reference_genome.fasta vcf_file.vcf Potential_small_microhomology.csv

### VCF annotation :
    python3 Annotation_vcf.py -v vcf_to_annotate.vcf -a Potential_duplication.psl -m Potential_mobile_element.csv -t Potential_TRF.csv -u Potential_small_microhomology.csv -l Potential_large_microhomology.psl -o output_name_vcf_annotated.vcf -h threshold