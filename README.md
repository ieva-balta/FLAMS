# Course Project Notice

This repository is a **modified fork** of [FLAMS](https://github.com/hannelorelongin/FLAMS), created for coursework as part of the course **Integrated Bioinformatics Project I0U20a** of the **Master of Bioinformatics** programme at **KU Leuven**.

The goal of the project was to adapt FLAMS to function with **[UniProt](https://www.uniprot.org/) as primary PTM information source**. We label this fork as FLAMS v1.2.0, since the position-based search logic remains unchanged and only the database source was modified. This project was done by *Ieva Baltā* (ieva-balta), *Maria Jose Caceres Valdiviezo* (majocava), and *Natasya Limanta* (naaattella), under the supervision of *Hannelore Login* (hannelorelongin).

To acknowledge contributions clearly, any scripts or modules we modified include our names as additional authors alongside the original authors.

# FLAMS: Find Lysine Acylations & other Modification Sites

A bioinformatics tool to analyze the conservation of post-translational modifications (PTMs), by means of a position-based search against UniProt entries that contain PTM information and are a part of [entries with experimental evidence at protein level](https://www.uniprot.org/help/protein_existence) and where the feature [Evidence Code Ontology (ECO) identifiers](https://www.uniprot.org/help/evidences) are 0000269, 0000314, 0007744, or 0007829.

# Table of contents

1.  [Introduction](#introduction)
2.  [System requirements](#system-requirements)
    1.  [General dependencies](#general-dependencies)
    2.  [Third-party dependencies](#third-party-dependencies)
3.  [Installation](#installation)
4.  [Usage](#usage)
    1. [Example use case](#example-use-case)
5.  [Output](#output)
6.  [Supported PTMs](#supported-ptms)
    1. [Supported PTM databases](#supported-ptm-databases)
    2. [Supported PTM types](#supported-ptm-types)
7.  [Contact](#contact)
8.  [References](#references)
9.  [License](#license)

## Introduction

FLAMS is a bioinformatics tool to analyze the conservation of post-translational modifications, by means of a position-based search against the UniProt database(The UniProt Consortium , UniProt: the Universal Protein Knowledgebase in 2025, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D609–D617). FLAMS can be used (i) to quickly verify whether modifications in a specific protein have been reported before, (ii) to assess whether findings in one species might translate to other species, and (iii) to systematically assess the novelty and conservation of reported  modification sites.

The tool takes as input a protein (identifier or sequence) and the position of an amino acid. This repository contains the command-line tool `FLAMS`, which obtains an overview of the previously reported post-translational modifications matching your query, by using the following scripts:

* *input.py*: processing the user-provided input
* *uniprot.py* and *setup.py*: downloading and preparing the modification-specific databases
* *run_blast.py*: searching your query against the databases of proteins with post-translational modifications
* *display.py*: formatting the list of conserved post-translational modifications to a tab delimited output file
* *utils.py*: dealing with OS-dependent directory systems

FLAMS is also available as a web service at https://www.biw.kuleuven.be/m2s/cmpg/research/CSB/tools/flams/ (not updated with the UniProt database).

## System requirements

Linux 64-bit, Windows and Mac OS supported.

### General dependencies

* Python3 (>=3.10, <3.12)

### Third-party dependencies

* [BLAST+ (>=2.13, tested until 2.16)](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)

## Installation

The UniProt version is not integrated into conda or pip. To install clone this repository, create a conda environment with all the dependencies (can use the conda environment of the original FLAMS) and run `python -m pip install ./PathToLocalFLAMS`.

Make sure to have BLAST+ installed locally and available in PATH. For more information on how to install BLAST+ on Windows, click [here](https://www.ncbi.nlm.nih.gov/books/NBK52637/) .

## Usage

Run the tool:

`FLAMS [-h] (--in inputFilePath | --id UniProtID | --batch batchFilePath) [-p position] [--range errorRange] [-o outputFilePath] [-d dataDir] [-t threadsBLAST] [-e evalueBLAST]
            [-m modification [modification ...]] `

Required argument:
* one of:
  * `inputFilePath`, used with `--in`, is the path to a .fasta file with the protein you wish to query against. (has to contain only 1 protein)
  * `UniProtID`, used with `--id`, is the UniProt ID of the protein you wish to query against.
  * `batchFilePath`, used with `--batch`, is the path to a tab seperated file for batch runs. The file should contain 1 entry per line, with UniProt IDs in the 1st column, and positions in the 2nd column.

Required argument when running FLAMS with --id/--in:
* `position` is the position of a (modified) residue in the protein, which you want to query against.

Optional arguments:
* `errorRange` is an number of positions before and after `position` to also search for modifications. [default: 0]
* `outputFilePath` is the path to where the result will be saved (in a .tsv file format). [default: out.tsv] If FLAMS is run with --batch, the specified -o/--output is used as preposition, followed by '\_$UniProtID\_$position.tsv'. [default: '']
* `dataDir` is the path to the directory where intermediate files (the UniProt sequence files) are stored. [default: $PWD/data]
* `threadsBLAST` is a BLAST parameter, allows you to speed up the search by multithreading. [default: 1]
* `evalueBLAST` is a BLAST parameter, allows you to filter out low quality BLAST hits. [default: 0.01]
* `modification` is a space-separated list of modifications (all lower case) to search for at the given position. We also provide aggregated combinations for each amino acid ($AA-All), and Acylations and Ubs. For a full list of all supported PTMs, and how they are named, see the [Supported PTM types](#supported-ptm-types) section of the README. In general, PTMs are written all lowercase, and spaces within a PTM name are replaced by underscores. [default: K-All]

### Example use case

We provide two example use cases for FLAMS:

With the following command, you search whether the TatA (UniProt ID: A0A916NWA0) acetylation on K66 in *Dehalococcoide mccartyi* strain CBDB1, as described by [Greiner-Haas (2021)](https://doi.org/10.3390/microorganisms9020365), had been previously detected.

`FLAMS --in A0A916NWA0.fa -p 66 -m acetylation -o tatA.tsv`

With the following command, you search whether the *Mycobabcterium smegmatis*' FadD2 (UniProt ID: A0QQ22) K537 is known to carry any modifications of the 'acylations' category, similar to what was reported by [Xu (2020)](https://doi.org/10.1128/mSystems.00424-19).

`FLAMS --id A0QQ22 -p 537 -m Acylations -o FadD2.tsv`

With the following commands, you search what modification could the Q at the 6th position of the *Homo sapiens* Histone H3.1 (UniProt ID: P68431) have. 

`FLAMS --id P68431 -p 6 -m Q-All -o Q-All_P68431.tsv`

You can find the example input and output data in the folder `test_data`. The output data is organized in folders reflecting the FLAMS version used to generate it, as the output can vary depending on the exact FLAMS version (due to FLAMS database updates).

For more example use cases, see the Supplementary information of the paper.

## Output

The output file is a .tsv containing one row per modification that matched the query, i.e., a modification aligning (within the user-specified range) to the query position, in a protein similar to the query protein. In case of batch jobs (ran with --batch), one output file per query (= a single line in the batch job file) will be generated.

The output file contains 13 columns:

* UniProt ID: UniProt identifier of the matched protein
* Protein name: protein name of the matched protein
* Modification: the type of modification found in the matched protein
* $AA location: the location of this matched modification in the matched protein
* $AA window: the local sequence containing the conserved modification (window of five amino acids before and after°)
* Species: the textual description of the species of the matched protein
* BLAST E Value: E value of BLASTp search of the matched protein against your query protein
* BLAST identity: % identity of BLASTp search of the matched protein against your query protein
* BLAST coverage: % coverage of BLASTp search of the matched protein against your query protein
* Database: UniProt subdatabase (Swiss-Prot or TrEMBL)
* ECO codes: evidence codes of matched protein modification
* Sources: source type and evidence links of matched protein modification, format type:ID, for large-scale study data only PubMed IDs are displayed
* LSS Database: large-scale study database of matched protein modification (PTMXchange or PRIDE)°°
* LSS IDs: large-scale study IDs of matched protein modification°°
* LSS Confidence scores: large-scale study confidence scores of matched protein modification (Gold, Silver, Bronze)°°

°: window can be smaller than the [-5;+5] window if the sequence alignment ends sooner, which can happen for modified sites near the start/end of the protein.

°°: only displayed for the PTMs inferred from large-scale study data fetched using EBI Proteomics API; only for acetylation, phosphorylation, sumoylation and ubiquitination.

## Supported PTMs

### Supported PTM databases

FLAMS updates its search databases regularly. To get an overview of the supported databases, see the table below.

|FLAMS version|CPLM version|dbPTM version|database available for download|UniProt release|
|:----|:----|:----|:----|:----|
|v1.2.0| | |no|2025_04|
|v1.1.6|v4 (Feb '25 update)|2025_July|[yes](https://doi.org/10.5281/zenodo.16737546)|2025_03|
|v1.1.5|v4|2025_January|[yes](https://doi.org/10.5281/zenodo.14616210)|2024_06|
|v1.1.4|v4|2024_April|[yes](https://doi.org/10.5281/zenodo.10958721)|2024_02|
|v1.1.0-3|v4|2023_November|[yes](https://doi.org/10.5281/zenodo.10171879)|2023_05|
|v1.0|v4| |[yes](https://cplm.biocuckoo.cn/Download.php)|NA|

Please note that the software doesn't store all UniProt entries. Only the entries with [experimental evidence at protein level](https://www.uniprot.org/help/protein_existence) (existence:1) and features with [Evidence Code Ontology (ECO) identifiers](https://www.uniprot.org/help/evidences) 0000269, 0000314, 0007744, or 0007829 are stored. 

**Modification of database:**
If you wish to use other filters, you can modify the *uniprot.py* script (valid_ECO_codes for the ECO code filter and/or the url in get_uniprot_records() for existence code). You can make changes in the MODIFICATION dictionary in *setup.py* to add/remove modification types, edit amino acid lists, or the regular expressions used to group entries. All changes should be paired with a change of the version number in *setup.py* to download a new version of the database.

### Supported PTM types

An overview of the PTM types, how to call them in FLAMS, and on which amino acid they can be found° is given in the table below. This table can also be found as a tab seperated file named FLAMS_supported_ptms_v20.txt .

|FLAMS PTM name|A (Ala)|C (Cys)|D (Asp)|E (Glu)|F (Phe)|G (Gly)|H (His)|I (Ile)|K (Lys)|L (Leu)|M (Met)|N (Asn)|P (Pro)|Q (Gln)|R (Arg)|S (Ser)|T (Thr)|V (Val)|W (Trp)|Y (Tyr)|Acylations|Ubs|
|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
|acetylation|X|X|X|X| |X| | |X| |X| |X| |X|X|X|X| |X|X| |
|adp-ribosylation| |X|X|X| |X|X| |X| | |X| | |X|X| | | |X| | |
|adp-riboxanation| | |X| | | |X| | | | |X| |X|X|X| | | | | | |
|amidation|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | |
|ampylation| | | | | | | | | | | | | | | |X|X| | |X| | |
|benzoylation| | | | | | | | |X| | | | | | | | | | | |X| |
|beta-hydroxybutyrylation| | | | | | | | |X| | | | | | | | | | | |X| |
|biotinylation| | | | | | | | |X| | | | | | | | | | | | | |
|blocked_amino_end|X|X|X|X| |X|X|X| |X|X|X|X|X|X|X|X|X| | | | |
|bromination|X|X|X| |X|X|X|X|X|X| | |X|X|X| | |X|X|X| | |
|butyrylation| | | | | | | | |X| | | | | | | | | | | |X| |
|carboxyethylation| | | | | | | | |X| | | | | | | | | | | | | |
|carboxylation| | | | | | | | |X| | | | | | | | | | | | | |
|cholesterol_ester| | | | | |X| | | | | | | | | | | | | | | | |
|citrullination| | | | | | | | | | | | | | |X| | | | | | | |
|crotonylation| | | | | | | | |X| | | | | | | | | | | |X| |
|cyclopeptide|X|X|X|X|X|X| |X|X|X|X|X|X|X|X|X| | |X|X| | |
|cysteinylation|X|X| |X| |X|X|X|X|X| |X|X|X|X|X|X|X|X|X| | |
|c-linked_glycosylation| | | | | | | | | | | | | | | | | | |X| | | |
|deamidation| | | | | | | | | | | |X| |X| | | | | | | | |
|deamination| | | | | | | | |X| | | | | | | | | | | | | |
|decanoylation| | | | | | | | | | | | | | | |X|X| | | | | |
|decarboxylation| | |X| | | | | | | | | | | | | |X| | | | | |
|dehydration|X|X|X| | |X|X| |X|X| |X| | |X|X|X|X| |X| | |
|dephosphorylation| | | | | | | | | | | | | | | |X|X| | |X| | |
|disulfide_bond| |X| | | | | | | | | | | | | | | | | | | | |
|d-glucuronoylation| | | | | |X| | | | | | | | | | | | | | | | |
|fadylation|X|X| | | |X|X|X|X|X| |X| | |X|X|X|X| |X| | |
|farnesylation| |X| | | | | | | | | | | | | | | | | | | | |
|formation_of_an_isopeptide_bond| | | |X| | | | | | | | | |X| | | | | | | | |
|formylation| | | | | |X| | |X| |X| | | | | | | | | |X| |
|genarylation| |X| | | |X| |X| |X| |X| |X|X| | | |X|X| | |
|geranylgeranylation| |X| | | | | | | | | | | | | | | | | | | | |
|glutamylation|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | |
|glutarylation| | | | | | | | |X| | | | | | | | | | | |X| |
|glutathionylation| |X| | | | | | | | | | | | | | | | | | | | |
|glycation| | | | | | | | |X| | | | | | | | | | | | | |
|gmpylation| |X| | | | | | | | | | | | | | | | | | | | |
|gpi-anchor|X|X|X| | |X| | | | | |X| | | |X|X| | | | | |
|histidylation|X| | | | |X|X|X|X|X| |X| | | |X| |X| | | | |
|hydroxylation| |X|X|X|X| |X|X|X|X| |X|X| |X|X|X|X|X|X| | |
|hypusine| | |X|X| | | | |X| | |X| | | |X| | | | | | |
|imidazolation|X|X| | | | |X| |X|X|X|X| |X|X|X|X|X| | | | |
|iodination| | | | | | | | | | | | | | | | | | | |X| | |
|isomerization|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | |
|lactoylation| | | | | | | | |X| | | | | | | | | | | | | |
|lipoylation| | | | | | | | |X| | | | | | | | | | | | | |
|malonylation| | | | | | | | |X| | | | | | | | | | | |X| |
|methylation|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| |X| | |
|myristoylation| |X| | | |X| | |X| | | | | | | | | | | | | |
|neddylation| | | | | | | | |X| | | | | | | | | | | | |X|
|nitration| | | | | | | | | | | | | | | | | | | |X| | |
|n-carbamoylation|X| | | | | | | | | | | | | | | | | | | | | |
|n-linked_glycosylation| | |X| | | | |X|X| | |X| | |X|X|X|X|X| | | |
|n-palmitoylation| |X| | | |X| | |X| | | | | | | | | | | | | |
|octanoylation| | | | | | | | | | | | | | | |X|X| | | | | |
|oxidation| |X| | | | | | | |X|X| | | | |X| | |X| | | |
|o-linked_glycosylation| | | | | | | | |X| | | |X| | |X|X| | |X| | |
|o-palmitoleoylation| | | | | | | | | | | | | | | |X| | | | | | |
|o-palmitoylation| | | | | | | | | | | | | | | |X|X| | | | | |
|palmitoylation|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | |
|phosphatidylethanolamine_amidation| | | | | |X| | | | | | | | | | | | | | | | |
|phosphoglycerylation| | | | | | | | |X| | | | | | | | | | | | | |
|phosphorylation|X|X|X|X|X|X|X|X|X|X| |X|X|X|X|X|X|X|X|X| | |
|prenylation| | | |X| | | | | | | | | | | | | | |X| | | |
|propionylation| | | | | | | | |X| | | | | | | | | | | |X| |
|pupylation| | | | | | | | |X| | | | | | | | | | | | |X|
|pyridoxal_phosphate_addition|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | |
|pyrrolidone_carboxylic_acid| | | |X| | | | | | | | | |X| | | | | | | | |
|pyrrolylation| |X| | | | | | | | | | | | | | | | | | | | |
|pyruvate| |X| | | | | | | | | | | | | |X| | | | | | |
|quinone| | | | | | | | | | | | | | | | | | | | | | |
|serotonylation| | | | | | | | | | | | | |X| | | | | | | | |
|stearoylation| |X| | | | | | | | | | | | | | | | | | | | |
|succinylation| |X| | | | | | |X| | | | | | | | | |X| |X| |
|sulfation| |X| | | | | | | | | | | | | |X|X| | |X| | |
|sulfhydration| |X| | | | | | | | | | | | | | | | | | | | |
|sulfilimine_crosslink| | | | | | | | |X| |X| | | | | |X| | | | | |
|sulfoxidation| | | | | | | | | | |X| | | | | | | | | | | |
|sumoylation| | | | | | | | |X| | | | | | | | | | | | |X|
|s-archaeol| |X| | | | | | | | | | | | | | | | | | | | |
|s-carbamoylation| |X| | | | | | | | | | | | | | | | | | | | |
|s-cyanation| |X| | | | | | | | | | | | | | | | | | | | |
|s-cysteinylation| |X| | | | | | | | | | | | | | | | | | | | |
|s-diacylglycerol| |X| | | | | | | | | | | | | | | | | | | | |
|s-linked_glycosylation| |X| | | | | | | | | | | | | | | | | | | | |
|s-nitrosylation| |X| | | | | | | | | | | | | | | | | | | | |
|s-palmitoylation| |X| | | | | | | | | | | | | | | | | | | | |
|thiocarboxylation| | | | | |X| | | | | | | | | | | | | | | | |
|thioester_crosslink| | | | | | | | | | | | | | | | | | | | | | |
|ubiquitination| |X| | | | | | |X| | | | | |X|X| | | | | |X|
|umpylation| | | | | | | | | | | | | | | |X|X| | |X| | |
|2-hydroxyisobutyrylation| | | | | | | | |X| | | | | | | | | | | |X| |

°: for PTM types that exist in the CPLM/dbPTM version of FLAMS the list of amino acids that could carry the modification type was left the same, however, for the newly added types the list was inferred from the downloaded entries that were sorted into that type and might not be biologically relevant. You can update the list in MODIFICATIONS in *setup.py*.

## Contact

Laboratory of Computational Systems Biology, KU Leuven.

## References

If you use FLAMS in your work, please cite:

Longin, H. *et al* (2024) "FLAMS: Find Lysine Acylations and other Modification Sites." Bioinformatics. 40(1):btae005.

In addition, FLAMS relies on third-party software & databases:

Altschul, S.F. *et al* (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

The UniProt Consortium , UniProt: the Universal Protein Knowledgebase in 2025, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D609–D617

## License

FLAMS is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
