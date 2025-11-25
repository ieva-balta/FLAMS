#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen, ieva-balta, majocava, naaattella
"""

import subprocess
import logging
import requests
import sys
from io import BytesIO
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq

# from ..databases import cplmv4
# from ..databases import dbptm
from . import uniprot
from ..utils import get_data_dir

""" setup
This script contains our modification database, and all associated functions to generate and maintain BLAST databases for each modification type.
"""

@dataclass
class ModificationDatabase:
    """
    This dataclass contains a list of tuples, containing the information necessary to get the FASTA file of the modification.

    Parameters
    ----------
    module: Any
        Refers to database module, necessary to retrieve the fasta files on the modification described by Descriptor
    descriptor:
        List of RegEx expressions

    """
    module: Any
    descriptor: List[str]


@dataclass
class ModificationType:
    """
    This dataclass contains the different types of modifications, which can be identied by type and database version, and contains their modification databases.

    Parameters
    ----------
    type: str
        Label to identify modification
    version: float
        Label to identify database version
    dbs: List[ModificationDatabase]
        List of modification database

    """
    type: str
    version: float
    dbs: List[ModificationDatabase]
    aas: List[str]


# # Here we store a dict of Zenodo URLs that can be queried for.
# version_urls = {
#     1.0: "https://zenodo.org/records/10143464/files/{0}-{1}.zip?download=1",
#     1.1: "https://zenodo.org/records/10171879/files/{0}-{1}.zip?download=1",
#     1.2: "https://zenodo.org/records/10958721/files/{0}-{1}.fasta.zip?download=1",
#     1.3: "https://zenodo.org/records/14616210/files/{0}-{1}.zip?download=1",
#     1.4: "https://zenodo.org/records/16737546/files/{0}-{1}.zip?download=1"
#     }

# version for UniProt download, starting with 2.0
version = 2.0

# Here we store a dict of modifications that can be queried for.
    # sorted alphabetically
    # allows duplication (same record can be associated with multiple modification types)
    # new types are tagged with comment NEW TYPE
    # for new types the amino acid list is taken from fasta files, it is NOT based on literature
    # for old types the amino acid list is left the same
    # types not found in UniProt are commented out
    # stores the RegEx list under descriptor
MODIFICATIONS = {
    "acetylation": ModificationType(
        "acetylation", version,
        [
        # ModificationDatabase(cplmv4, "Acetylation"), ModificationDatabase(dbptm,"Acetylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*acetyl[-\w_]*"])
        ],
        ["A","C","D","E","G","K","M","P","R","S","T","V","Y"]
    ),
    "adp-ribosylation": ModificationType(
        "adp-ribosylation", version,
        [
        # ModificationDatabase(dbptm, "ADP-ribosylation"),
             ModificationDatabase(uniprot, [r"ADP[-_\s]*ribosyl[-\w_]*"])
        ],
        ["C","D","E","G","H","K","N","R","S","Y"]
    ),
    "adp-riboxanation": ModificationType(
        "adp-riboxanation", version,
        [
            ModificationDatabase(uniprot, [r"ADP[-_\s]*ribox[-\w_]*"])
        ],
        ['R']
    ), # NEW TYPE
    "amidation": ModificationType(
        "amidation", version,
        [
        #ModificationDatabase(dbptm, "Amidation"),
            ModificationDatabase(uniprot, [r"[-\w_]*(?<!de)amid[-\w_]*"])
        ],
        ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    ),
    "ampylation" : ModificationType(
        "ampylation", version,
        [
            # ModificationDatabase(dbptm, "AMPylation"),
            ModificationDatabase(uniprot, [r"(?<!S)AMP[-\w_]*", r"[-\w_]*adenylate[-\w_]*"])
        ],
        ["S","T","Y"]
    ),
    "benzoylation": ModificationType(
        "benzoylation", version, 
        [
            # ModificationDatabase(cplmv4, "Benzoylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*benzoyl[-\w_]*"])
        ],
        ["K"]
    ),
    "beta-hydroxybutyrylation": ModificationType(
        "beta-hydroxybutyrylation", version, 
        [
            # ModificationDatabase(cplmv4, "β-Hydroxybutyrylation"),
            ModificationDatabase(uniprot, [r"(?:β|beta)[-_\s]*hydroxybutyryl[-\w_]*"])
        ],
        ["K"]
    ),
    "biotinylation": ModificationType(
        "biotinylation", version,
        [
            # ModificationDatabase(cplmv4, "Biotinylation"), ModificationDatabase(dbptm, "Biotinylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*biotinyl[-\w_]*"])
        ],
        ["K"]
    ),
    "blocked_amino_end": ModificationType(
        "blocked_amino_end", version,
        [
            # ModificationDatabase(dbptm, "Blocked amino end"),
            ModificationDatabase(uniprot, [r"blocked__(?:amino|carboxyl)__end"])
        ],
        ["A","C","D","E","G","H","I","L","M","N","P","Q","R","S","T","V"]
    ),
    "bromination" : ModificationType(
        "bromination", version,
        [
            ModificationDatabase(uniprot, [r"[-\w_]*bromo[-\w_]*"])
        ],
        ['H', 'W', 'Y']
    ), # NEW TYPE
    "butyrylation": ModificationType(
        "butyrylation", version,
        [
            # ModificationDatabase(cplmv4, "Butyrylation"), ModificationDatabase(dbptm, "Butyrylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*butyryl[_\w-]*"])
        ],
        ["K"]
    ),
    # "carbamidation": ModificationType(
    #     "carbamidation", 1.4,
    #     [ModificationDatabase(dbptm, "Carbamidation")],
    #     ["C"]
    # ),
    "carboxyethylation": ModificationType(
        "carboxyethylation", version,
        [
            # ModificationDatabase(cplmv4, "Carboxyethylation"), ModificationDatabase(dbptm, "Carboxyethylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*carboxyethyl[_\w-]*"])
        ],
        ["K"]
    ),
    "carboxylation": ModificationType(
        "carboxylation", version,
        [
            # ModificationDatabase(cplmv4, "Carboxylation"), ModificationDatabase(dbptm, "Carboxylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*(?<!de)carboxy[_\w-]*"])
        ],
        ["K"]
    ),
    # "carboxymethylation": ModificationType(
    #     "carboxymethylation", 1.4, [ModificationDatabase(cplmv4, "Carboxymethylation")],
    #     ["K"]
    # ),
    "cholesterol_ester": ModificationType(
        "cholesterol_ester", version,
        [
            # ModificationDatabase(dbptm, "Cholesterol ester"),
            ModificationDatabase(uniprot, [r"[-\w_]*cholesterol(?:__glycine)?__ester"])
        ],
        ["G"]
    ),
    "citrullination": ModificationType(
        "citrullination", version,
        [
            # ModificationDatabase(dbptm, "Citrullination"),
            ModificationDatabase(uniprot, [r"[-\w_]*citrulline"])
        ],
        ["R"]
    ),
    "crotonylation": ModificationType(
        "crotonylation", version,
        [
            # ModificationDatabase(cplmv4, "Crotonylation"), ModificationDatabase(dbptm, "Crotonylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*crotonyl[-\w_]*"])
        ],
        ["K"]
    ),
    "cyclopeptide" : ModificationType(
        "cyclopeptide", version,
        [
            ModificationDatabase(uniprot, [r"[-\w_]*cyclo(?:peptide)?[-\w_]*"])
        ],
        ['L', 'R', 'K', 'W', 'M', 'S', 'A', 'E', 'I', 'N', 'C', 'F', 'P', 'D', 'G']
    ), # NEW TYPE
    "cysteinylation": ModificationType(
        "cysteinylation", version,
        [
            ModificationDatabase(uniprot, [r"[_\w-]*cysteinyl[_\w-]*"])
        ],
        ['C', 'Y', 'T', 'M', 'E', 'W', 'R', 'H', 'F', 'D']
    ), # NEW TYPE
    "c-linked_glycosylation": ModificationType(
        "c-linked_glycosylation", version,
        [
            # ModificationDatabase(dbptm, "C-linked Glycosylation"),
            ModificationDatabase(uniprot, [r"c[-_\s]*linked[-\w_]*"])    
        ],
        ["W"]
    ),
    "deamidation": ModificationType(
        "deamidation", version,
        [
            # ModificationDatabase(dbptm, "Deamidation"),
            ModificationDatabase(uniprot, [r"[-\w_]*deamidat[-\w_]*"])
        ],
        ["N","Q"]
    ),
    "deamination": ModificationType(
        "deamination", version,
        [
            # ModificationDatabase(dbptm, "Deamination"),
            ModificationDatabase(uniprot, [r"[-\w_]*allysine[-\w_]*"])
        ],
        ["K"]
    ),
    "decanoylation": ModificationType(
        "decanoylation", version,
        [
            # ModificationDatabase(dbptm, "Decanoylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*decanoyl[-\w_]*"])
        ],
        ["S","T"]
    ),
    "decarboxylation": ModificationType(
        "decarboxylation", version,
        [
            # ModificationDatabase(dbptm, "Decarboxylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*decarboxylat[-\w_]*"])
        ],
        ["D","T"]
    ),
    "dehydration": ModificationType(
        "dehydration", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*dehydro[-\w_]*"])],
        ['D', 'Y', 'S', 'C', 'T']
    ), # NEW TYPE
    "dephosphorylation": ModificationType(
        "dephosphorylation", version,
        [
            # ModificationDatabase(dbptm, "Dephosphorylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*dephospho[-\w_]*"])
        ],
        ["S","T","Y"]
    ),
    # "dietylphosphorylation": ModificationType(
    #     "dietylphosphorylation", 1.4, [ModificationDatabase(cplmv4, "Dietylphosphorylation")],
    #     ["K"]
    # ),
    "disulfide_bond": ModificationType(
        "disulfide_bond", version,
        [
            # ModificationDatabase(dbptm, "Disulfide bond"),
            ModificationDatabase(uniprot, [r"[-\w_]*disulfide[-\w_]*"])
        ],
        ["C"]
    ),
    "d-glucuronoylation": ModificationType(
        "d-glucuronoylation", version,
        [
            # ModificationDatabase(dbptm, "D-glucuronoylation"),
            ModificationDatabase(uniprot, [r"d[-_\s]*glucuronoyl[-\w_]*"])
        ],
        ["G"]
    ),
    "fadylation" : ModificationType(
        "fadylation", version,
        [ModificationDatabase(uniprot, [r"FAD[-\w_]*"])],
        ['C', 'H']
    ), # NEW TYPE
    "farnesylation": ModificationType(
        "farnesylation", version,
        [
            # ModificationDatabase(dbptm, "Farnesylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*farnesyl[-\w_]*"])
        ],
        ["C"]
    ),
    "formation_of_an_isopeptide_bond": ModificationType(
        "formation_of_an_isopeptide_bond", version,
        [
            # ModificationDatabase(dbptm, "Formation of an isopeptide bond"),
            ModificationDatabase(uniprot, [r"[_\w-]*isopeptide[_\w-]*"])
        ],
        ["E","Q"]
    ),
    "formylation": ModificationType(
        "formylation", version,
        [
            # ModificationDatabase(cplmv4, "Formylation"), ModificationDatabase(dbptm, "Formylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*formyl[-\w_]*"])
        ],
        ["G","K","M"]
    ),
    # "gamma-carboxyglutamic_acid": ModificationType(
    #     "gamma-carboxyglutamic_acid", version,
    #     [
    #         # ModificationDatabase(dbptm, "Gamma-carboxyglutamic acid"),
    #         ModificationDatabase(uniprot, [r"(?:gamma|γ)-?carboxyglutamic__acid[-\w_]*"])
    #     ],
    #     ["E"]
    # ),
    "genarylation": ModificationType(
        "genarylation", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*geranyl[_\w-]*"])],
        ['C', 'W']
    ), # NEW TYPE
    "geranylgeranylation": ModificationType(
        "geranylgeranylation", version,
        [
            # ModificationDatabase(dbptm, "Geranylgeranylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*geranylgeranyl[-\w_]*"])
        ],
        ["C"]
    ),
    "glutamylation" : ModificationType(
        "glutamylation", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*glutamyl[-\w_]*"])],
        ['E', 'G', 'K', 'C', 'Q']
    ), # NEW TYPE
    "glutarylation": ModificationType(
        "glutarylation", version,
        [
            # ModificationDatabase(cplmv4, "Glutarylation"), ModificationDatabase(dbptm, "Glutarylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*glutaryl[-\w_]*"])
        ],
        ["K"]
    ),
    "glutathionylation": ModificationType(
        "glutathionylation", version,
        [
            # ModificationDatabase(dbptm, "Glutathionylation"),
            ModificationDatabase(uniprot, [r"[-\w_]*glutathionyl(?:__cysteine)?[-\w_]*"])
        ],
        ["C"]
    ),
    "glycation": ModificationType(
        "glycation", version, 
        [
            # ModificationDatabase(cplmv4, "Glycation"),
            ModificationDatabase(uniprot, [r"[-\w_]*glycation[-\w_]*"])
        ],
        ["K"]
    ),
    "gmpylation": ModificationType(
        "gmpylation", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*GMP[-\w_]*"])],
        ['C']
    ), # NEW TYPE
    "gpi-anchor": ModificationType(
        "gpi-anchor", version, 
        [
            # ModificationDatabase(dbptm, "GPI-anchor"),
            ModificationDatabase(uniprot, [r"[-\w_]*gpi[-\s]?(?:like-)?anchor[-\w_]*"])
        ],
        ["A","C","D","G","N","S","T"]
    ),
    "histidylation": ModificationType(
        "histidylation", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*histidyl[-\w_]*"])],
        ['C', 'H', 'Y']
    ), # NEW TYPE
    # "hmgylation": ModificationType(
    #     "hmgylation", 1.4, [ModificationDatabase(cplmv4, "HMGylation")],
    #     ["K"] #OK
    # ),
    # "hydroxyceramide_ester": ModificationType(
    #     "hydroxyceramide_ester", version, 
    #     [
    #         # ModificationDatabase(dbptm, "Hydroxyceramide ester"),
    #         ModificationDatabase(uniprot, [r"[-\w_]*hydroxyceramide(?:__glutamate)?__ester[-\w_]*"])    
    #     ],
    #     ["Q"]
    # ),
    "hydroxylation": ModificationType(
        "hydroxylation", version,
        [
            # ModificationDatabase(cplmv4, "Hydroxylation"), ModificationDatabase(dbptm, "Hydroxylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*(?<!de)hydroxy[_\w-]*"])
        ],
        ["C","D","E","F","H","I","K","L","N","P","R","S","T","V","W","Y"]
    ),
    "hypusine": ModificationType(
        "hypusine", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*hypusine[-\w_]*"])],
        ['K']
    ), # NEW TYPE
    "imidazolation": ModificationType(
        "imidazolation", version,
        [ModificationDatabase(uniprot, [r"[-\w_]*imidazol[-\w_]*"])],
        ['M', 'A', 'C', 'K', 'S', 'G', 'R', 'Q']
    ), # NEW TYPE
    "iodination": ModificationType(
        "iodination", version,
        [
            # ModificationDatabase(dbptm, "Iodination"),
            ModificationDatabase(uniprot, [r"[-\w_]*iodo[_\w-]*", r"[_\w-]*thyroxine[_\w-]*"])
        ],
        ["Y"]
    ),
    "isomerization": ModificationType(
        "isomerization", version,
        [ModificationDatabase(uniprot, [r"(?:^|-)[dl]-[_\w-]*"])],
        ['R', 'K', 'G', 'I', 'M', 'N', 'C', 'F', 'V', 'A', 'W', 'L', 'T', 'S']
    ), # NEW TYPE
    "lactoylation": ModificationType(
        "lactoylation", version,
        [
            # ModificationDatabase(dbptm, "Lactoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*lactoyl[_\w-]*"])
        ],
        ["K"]
    ),
    # "lactylation": ModificationType(
    #     "lactylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Lactylation"), ModificationDatabase(dbptm, "Lactylation")],
    #     ["K"]
    # ),
    "lipoylation": ModificationType(
        "lipoylation", version,
        [
            # ModificationDatabase(cplmv4, "Lipoylation"), ModificationDatabase(dbptm, "Lipoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*lipoyl[_\w-]*"])
        ],
        ["K"]
    ),
    "malonylation": ModificationType(
        "malonylation", version,
        [
            # ModificationDatabase(cplmv4, "Malonylation"), ModificationDatabase(dbptm, "Malonylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*malonyl[_\w-]*"])
        ],
        ["K"]
    ),
    "methylation": ModificationType(
        "methylation", version,
        [
            # ModificationDatabase(cplmv4, "Methylation"), ModificationDatabase(dbptm, "Methylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*methyl[_\w-]*"])
        ],
        ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","Y"]
    ),
    # "mgcylation": ModificationType(
    #     "mgcylation", 1.4, [ModificationDatabase(cplmv4, "MGcylation")],
    #     ["K"]
    # ),
    # "mgylation": ModificationType(
    #     "mgylation", 1.4, [ModificationDatabase(cplmv4, "MGylation")],
    #     ["K"]
    # ),
    "myristoylation": ModificationType(
        "myristoylation", version, 
        [
            # ModificationDatabase(dbptm, "Myristoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*myristoyl[_\w-]*"])
        ],
        ["C","G","K"]
    ),
    "neddylation": ModificationType(
        "neddylation", version,
        [
            # ModificationDatabase(cplmv4, "Neddylation"), ModificationDatabase(dbptm, "Neddylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*NEDD8[_\w-]*"])
        ],
        ["K"]
    ),
    "nitration": ModificationType(
        "nitration", version, 
        [
            # ModificationDatabase(dbptm, "Nitration"),
            ModificationDatabase(uniprot, [r"[_\w-]*nitro[_\w-]*", r"[_\w-]*nitrated[_\w-]*"])
        ],
        ["Y"]
    ),
    "n-carbamoylation": ModificationType(
        "n-carbamoylation", version, 
        [
            # ModificationDatabase(dbptm, "N-carbamoylation"),
            ModificationDatabase(uniprot, [r"n[-_\s]*carbamoyl[_\w-]*"])
        ],
        ["A"]
    ),
    "n-linked_glycosylation": ModificationType(
        "n-linked_glycosylation", version,
        [
            # ModificationDatabase(dbptm, "N-linked Glycosylation"),
            ModificationDatabase(uniprot, [r"n[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"])
        ],
        ["D","I","K","N","R","S","T","V","W"]
    ),
    "n-palmitoylation": ModificationType(
        "n-palmitoylation", version, 
        [
            # ModificationDatabase(dbptm, "N-palmitoylation"),
            ModificationDatabase(uniprot, [r"n[-_\s]*palmitoyl[_\w-]*"])    
        ],
        ["C","G","K"]
    ),
    "octanoylation": ModificationType(
        "octanoylation", version,
        [
            # ModificationDatabase(dbptm, "Octanoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*octanoyl[_\w-]*"])    
        ],
        ["S","T"]
    ),
    "oxidation": ModificationType(
        "oxidation", version,
        [
            # ModificationDatabase(dbptm, "Oxidation"),
            ModificationDatabase(uniprot, [r"[_\w-]*oxo[_\w-]*", r"[_\w-]*sulf[ei]nic__acid[_\w-]*"])
        ],
        ["C","L","M","S","W"]
    ),
    "o-linked_glycosylation": ModificationType(
        "o-linked_glycosylation", version,
        [
            # ModificationDatabase(dbptm, "O-linked Glycosylation"),
            ModificationDatabase(uniprot, [r"o[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"])    
        ],
        ["K","P","S","T","Y"]
    ),
    "o-palmitoleoylation": ModificationType(
        "o-palmitoleoylation", version,
        [
            # ModificationDatabase(dbptm, "O-palmitoleoylation"),
            ModificationDatabase(uniprot, [r"o[-_\s]*palmitoleoyl[_\w-]*"])
        ],
        ["S"]
    ),
    "o-palmitoylation": ModificationType(
        "o-palmitoylation", version,
        [
            # ModificationDatabase(dbptm, "O-palmitoylation"),
            ModificationDatabase(uniprot, [r"o[-_\s]*palmitoyl[_\w-]*"])
        ],
        ["S","T"]
    ),
    "palmitoylation": ModificationType(
        "palmitoylation", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*palmitoyl[_\w-]*"])],
        ['L', 'K', 'S', 'G', 'T', 'A', 'C']
    ), # NEW TYPE
    "phosphatidylethanolamine_amidation": ModificationType(
        "phosphatidylethanolamine_amidation", version, 
        [
            # ModificationDatabase(dbptm, "Phosphatidylethanolamine amidation"),
            ModificationDatabase(uniprot, [r"[_\w-]*phosphatidylethanolamine__amidated[_\w-]*"])    
        ],
        ["G"]
    ),
    "phosphoglycerylation": ModificationType(
        "phosphoglycerylation", version,
        [
            # ModificationDatabase(cplmv4, "Phosphoglycerylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*glycerophospho[_\w-]*"])
        ],
        ["K"]
    ),
    "phosphorylation": ModificationType(
        "phosphorylation", version,
        [
            # ModificationDatabase(dbptm, "Phosphorylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*(?<!de)phospho[_\w-]*", r"[_\w-]*aspartylphosphate[_\w-]*"])    
        ],
        ["A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","W","Y"]
    ),
    "prenylation": ModificationType(
        "prenylation", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*prenyl[_\w-]*"])],
        ['W']
    ), # NEW TYPE
    "propionylation": ModificationType(
        "propionylation", version,
        [
            # ModificationDatabase(cplmv4, "Propionylation"), ModificationDatabase(dbptm, "Propionylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*propionyl[_\w-]*"])
        ],
        ["K"]
    ),
    "pupylation": ModificationType(
        "pupylation", version, 
        [
            # ModificationDatabase(cplmv4, "Pupylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*pup[_\w-]*"])    
        ],
        ["K"]
    ),
    "pyridoxal_phosphate_addition": ModificationType(
        "pyridoxal_phosphate_addition", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*pyridoxal__phosphate[_\w-]*"])],
        ['K', 'D']
    ), # NEW TYPE
    "pyrrolidone_carboxylic_acid": ModificationType(
        "pyrrolidone_carboxylic_acid", version, 
        [
            # ModificationDatabase(dbptm, "Pyrrolidone carboxylic acid"),
            ModificationDatabase(uniprot, [r"[_\w-]*pyrrolidone__carboxylic__acid[_\w-]*"])    
        ],
        ["E","Q"]
    ),
    "pyrrolylation": ModificationType(
        "pyrrolylation", version, 
        [
            # ModificationDatabase(dbptm, "Pyrrolylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*pyrrolyl[_\w-]*"])    
        ],
        ["C"]
    ),
    "pyruvate": ModificationType(
        "pyruvate", version, 
        [
            # ModificationDatabase(dbptm, "Pyruvate"),
            ModificationDatabase(uniprot, [r"[_\w-]*pyruvic__acid[_\w-]*"])
        ],
        ["C","S"]
    ),
    "quinone": ModificationType(
        "quinone", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*quinone[_\w-]*"])],
        ['Y', 'W', 'C', 'K']
    ), # NEW TYPE
    "serotonylation" : ModificationType(
        "serotonylation", version,
        [
            # ModificationDatabase(dbptm, "Serotonylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*seroton[_\w-]*"])
        ],
        ["Q"]
    ),
    "stearoylation" : ModificationType(
        "stearoylation", version,
        [
            # ModificationDatabase(dbptm, "Stearoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*stearoyl[_\w-]*"])
        ],
        ["C"]
    ),
    "succinylation": ModificationType(
        "succinylation", version,
        [
            # ModificationDatabase(cplmv4, "Succinylation"), ModificationDatabase(dbptm, "Succinylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*succin[iy][_\w-]*"])    
        ],
        ["C","K","W"]
    ),
    "sulfation" : ModificationType(
        "sulfation", version,
        [
            # ModificationDatabase(dbptm, "Sulfation"),
            ModificationDatabase(uniprot, [r"[_\w-]*sulfo[_\w-]*"])
        ],
        ["C","S","T","Y"]
    ),
    "sulfhydration" : ModificationType(
        "sulfhydration", version,
        [
            # ModificationDatabase(dbptm, "Sulfhydration"),
            ModificationDatabase(uniprot, [r"[_\w-]*cysteine__persulfide[_\w-]*"])
        ],
        ["C"]
    ),
    "sulfilimine_crosslink": ModificationType(
        "sulfilimine_crosslink", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*sulfilimine[_\w-]*"])],
        ['M', 'K']
    ), # NEW TYPE
    "sulfoxidation" : ModificationType(
        "sulfoxidation", version,
        [
            # ModificationDatabase(dbptm, "Sulfoxidation"),
            ModificationDatabase(uniprot, [r"[_\w-]*sulfoxide[_\w-]*"])
        ],
        ["M"]
    ),
    "sumoylation": ModificationType(
        "sumoylation", version,
        [
            # ModificationDatabase(cplmv4, "Sumoylation"), ModificationDatabase(dbptm, "Sumoylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*SUMO[_\w-]*"])    
        ],
        ["K"]
    ),
    "s-archaeol" : ModificationType(
        "s-archaeol", version,
        [
            # ModificationDatabase(dbptm, "S-archaeol"),
            ModificationDatabase(uniprot, [r"s[_\s-]?archaeol[_\w-]*"])
        ],
        ["C"]
    ),
    "s-carbamoylation" : ModificationType(
        "s-carbamoylation", version,
        [
            # ModificationDatabase(dbptm, "S-carbamoylation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*carbamoyl[_\w-]*"])
        ],
        ["C"]
    ),
    "s-cyanation" : ModificationType(
        "s-cyanation", version,
        [
            # ModificationDatabase(dbptm, "S-Cyanation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*cyano[_\w-]*"])
        ],
        ["C"]
    ),
    "s-cysteinylation" : ModificationType(
        "s-cysteinylation", version,
        [
            # ModificationDatabase(dbptm, "S-cysteinylation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*cysteinyl[_\w-]*"])    
        ],
        ["C"]
    ),
    "s-diacylglycerol" : ModificationType(
        "s-diacylglycerol", version,
        [
            # ModificationDatabase(dbptm, "S-diacylglycerol"),
            ModificationDatabase(uniprot, [r"s[-_\s]*diacylglycerol[_\w-]*"])    
        ],
        ["C"]
    ),
    "s-linked_glycosylation": ModificationType(
        "s-linked_glycosylation", version,
        [
            # ModificationDatabase(dbptm, "S-linked Glycosylation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"])
        ],
        ["C"]
    ),
    "s-nitrosylation" : ModificationType(
        "s-nitrosylation", version,
        [
            # ModificationDatabase(dbptm, "S-nitrosylation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*nitro[_\w-]*"])
        ],
        ["C"]
    ),
    "s-palmitoylation" : ModificationType(
        "s-palmitoylation", version,
        [
            # ModificationDatabase(dbptm, "S-palmitoylation"),
            ModificationDatabase(uniprot, [r"s[-_\s]*palmitoyl[_\w-]*"])
        ],
        ["C"]
    ),
    "thiocarboxylation" : ModificationType(
        "thiocarboxylation", version,
        [
            # ModificationDatabase(dbptm, "Thiocarboxylation"),
            ModificationDatabase(uniprot, [r"[_\w-]*thioglycine[_\w-]*"])
        ],
        ["G"]
    ),
    "thioester_crosslink": ModificationType(
        "thioester_crosslink", version,
        [ModificationDatabase(uniprot, [r"[_\w-]*thioester[_\w-]*"])],
        ['Q', 'C', 'G']
    ), # NEW TYPE
    "ubiquitination": ModificationType(
        "ubiquitination", version,
        [
            # ModificationDatabase(cplmv4, "Ubiquitination"), ModificationDatabase(dbptm, "Ubiquitination"),
            ModificationDatabase(uniprot, [r"[_\w-]*ubiquitin[_\w-]*"])
        ],
        ["C","K","R","S"]
    ),
    "umpylation" : ModificationType(
        "umpylation", version,
        [
            # ModificationDatabase(dbptm, "UMPylation"),
            ModificationDatabase(uniprot, [r"UMP[_\w-]*"])
        ],
        ["S","T","Y"]
    ),
    "2-hydroxyisobutyrylation": ModificationType(
        "2-hydroxyisobutyrylation", version, 
        [
            # ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation"),
            ModificationDatabase(uniprot, [r"2[-_\s]*hydroxyisobutyryl[_\w-]*"])    
        ],
        ["K"]
    ),
    }


def update_db_for_modifications(list_of_mods_to_check: List[str]):
    """
    This function updates the local BLASTDB for a given modification.

    Parameters
    ----------
    list_of_mods_to_check: List[str]
        List of modifications (which are keys to any of the ModificationType's stored in the MODIFICATIONS dictionary),
        for which the database should be updated.

    """
    for m in list_of_mods_to_check:
        _generate_blastdb_if_not_up_to_date(MODIFICATIONS[m])


def _generate_blastdb_if_not_up_to_date(modification: ModificationType):
    """
    This function generates a new local BLASTDB when a newer database version for a specific modification is available.

    Parameters
    ----------
    modification: ModificationType
        ModificationType for which a BLASTDB will be generated

    """
    data_dir = get_data_dir()

    BLASTDB_PATH = get_blastdb_name_for_modification(
        modification.type, modification.version
    )

    # If an up-to-date BLASTDB for the given modification already exists, do nothing.
    if Path(f"{data_dir}/{BLASTDB_PATH}.pdb").exists():
        logging.info(f"BLAST database for modification {modification.type} already exists.")
        return

    # If no up-to-date BLASTDB exists, check whether a FASTA file exists with the modification data
    fasta_location_mod = f"{data_dir}/{modification.type}-{modification.version}.fasta"

    if not Path(fasta_location_mod).exists():
        logging.info(f"Fasta file for modification {modification.type} doesn't exist. Proceeding with uniprot download.")
        # # If this does not exist, create a FASTA file with the modification data
        # try:
        #     _get_fasta_from_zenodo(modification, data_dir, fasta_location_mod)
        # except requests.HTTPError:
        #     logging.error(f"Could not fetch FLAMS {modification.type} Database {modification.version} from Zenodo. Please try again later. Exiting FLAMS...")
        #     logging.error("If needed, you can also download the databases from CPLM and dbPTM yourself, see the docs. (not recommended: very slow!)")
        #     sys.exit()

# FOR OWN DATABASE DOWNLOAD - on a fresh install:
# Comment out try/except above (lines 516-521), uncomment the following line of code
        # downloads uniprot database
        get_fasta_from_uniprot(data_dir)

    # Generate local BLASTDB from FASTA in fasta_location_mod
    _generate_blastdb(data_dir, modification)


def _generate_blastdb(data_dir, modification: ModificationType):
    """
    This function generates a local BLASTDB for a given modification.

    Parameters
    ----------
    data_dir: directory
        Platform-specific directory that stores app data. The local BLAST database will be stored here.
    modification: ModificationType
        ModificationType for which a local BLAST database will be generated

    """
    try:
        # We presume that the FASTA is stored in a file {modification.type}.fasta inside the data_dir.
        # We will write the local BLASTDB to out_path
        out_db_name = get_blastdb_name_for_modification(
            modification.type, modification.version
        )
        subprocess.call(
            f'cd "{data_dir}" && makeblastdb -in {modification.type}-{modification.version}.fasta '
            f'-dbtype prot -input_type fasta -parse_seqids'
            f" -out {out_db_name}",
            shell=True,
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(
            "You need to install BLAST and include it in system PATH."
        ) from e


def get_blastdb_name_for_modification(modification: str, version=None):
    """
    This function gets the name of the local BLASTDB for a given modification.

    Parameters
    ----------
    modification: str
        Description of a specific modification,
        must be the key to any of the ModificationType's stored in the MODIFICATIONS dictionary.
    version: float
        Database version

    """
    # If version was not specified, get the current
    if not version:
        version = MODIFICATIONS[modification].version

    return f"{modification}-{version}"


# def _get_fasta_from_zenodo(modification, data_dir, fasta_location_mod):
#     """
#     This function downloads the .zip file containing the fasta file with information from all databases containing information on the specified modification.
#     It will create and store a .fasta file containing all entries of all databases related to this to modification file in the fasta_location.

#     Parameters
#     ----------
#     modification: ModificationType
#         ModificationType for which a fasta file will be generated
#     data_dir: directory
#         Platform-specific directory that stores app data. The local BLAST database will be stored here.
#     fasta_location_mod:
#         Output .fasta file containing all modification entries of all databases for the specified modification

#     """
#     URL = version_urls.get(modification.version)
#     req = requests.get(URL.format(modification.type, modification.version), stream=True)
#     # Raise an exception if HTTP request failed.
#     req.raise_for_status()
#     size_in_mb = int(req.headers.get("content-length")) / 1048576
#     logging.info(f"Downloading FLAMS {modification.type} Database {modification.version}, please wait. Size: {size_in_mb:.1f} MB")

#     with ZipFile(BytesIO(req.content)) as myzip:
#         # Extract the single txt file and return as UTF-8 string
#         ptm = myzip.read(myzip.namelist()[0]).decode("UTF-8")

#     filename = Path(fasta_location_mod)
#     with filename.open("w+", encoding="UTF-8") as f:
#         f.write(ptm)


#############################################
# Optional: downloading your own databases #
#############################################

def get_fasta_from_uniprot(data_dir):
    """
    This function downloads the UniProt databases for all modification types.
    To update change the version in this .py file.

    Parameters
    ----------
    data_dir: directory
        Platform-specific directory that stores app data. The UniProt databases will be stored here.

    """
    uniprot.get_fasta(MODIFICATIONS, data_dir)
    

# def _get_fasta_for_blast(modification: ModificationType, data_dir, fasta_location_mod):
#     """
#     This function creates a fasta file combining the information from all databases containing information on the specified modification.
#     It will create and store a fasta file containing all entries of all databases related to this to modification file in the fasta_location.

#     Parameters
#     ----------
#     modification: ModificationType
#         ModificationType for which a fasta file will be generated
#     data_dir: directory
#         Platform-specific directory that stores app data. The local BLAST database will be stored here.
#     fasta_location_mod:
#         Output .fasta file containing all modification entries of all databases for the specified modification

#     """
#     # create fasta
#     for db in modification.dbs:
#         module_name = db.module.__name__.split('.')[-1]
#         fasta_location_per_database = f"{data_dir}/{modification.type}-{modification.version}-{module_name}.fasta"
#         file_before_deduplication = f"{fasta_location_mod.removesuffix('.fasta')}-before_deduplication.fasta"
#         if not Path(fasta_location_per_database).exists():
#             _get_fasta_from_dbs(modification, db, fasta_location_per_database)
#         fileMod = open(file_before_deduplication, "a")
#         fileDb = open(fasta_location_per_database, "r")
#         fileMod.write(fileDb.read())
#         fileMod.close()
#         fileDb.close()

#     # remove duplicate entries, that will create errors when creating BLAST DB
#         # should only happen in rare cases, dbPTM and CPLM try to avoid duplicates, but they do happen
#     seen = set()
#     records = []
#     for record in SeqIO.parse(file_before_deduplication, "fasta"):
#         if record.id not in seen:
#             seen.add(record.id)
#             records.append(record)
#         else:
#             logging.warning(f'Your database contained a duplicate entry for {record.id}, only the first one was added to your BLAST database.')
#     SeqIO.write(records, fasta_location_mod, "fasta")


# def _get_fasta_from_dbs(modification: ModificationType, db, fasta_location):
#     """
#     This function calls on the get_fasta function from a specific module for a modification.
#     It will create and store a fasta file containing all entries in the specified database related to this to modification file in the fasta_location.

#     Parameters
#     ----------
#     modification: ModificationType
#         ModificationType for which a fasta file will be generated
#     db:
#         PTM database for which all entries will be converted to a fasta file (either dbptm or cplmv4)
#     fasta_location:
#         Output .fasta file containing all modification entries of the specified database for the specified modification

#     """
#     db.module.get_fasta(db.descriptor, fasta_location)
