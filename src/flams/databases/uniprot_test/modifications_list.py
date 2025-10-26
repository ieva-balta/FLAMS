
"""
taking the modification list from setup.py
making a dictionary with the uniprot query terms for each modification

to run on bash cd into the uniprot_test folder and run this file
"""
from uniprot import get_fasta_rest


MODIFICATIONS = {
    "acetylation": "ft_mod_res:acetyl*",
    "adp-ribosylation": "(ft_mod_res:ADP-ribosyl* OR ft_mod_res:ADP ribosyl*)",
    "adp-riboxanation": "ft_mod_res:ADP-riboxanated", # new type
    "amidation": "ft_mod_res:amide",
    "ampylation" : "ft_mod_res:AMP",
    "benzoylation": "ft_mod_res:benzoyl*",
    "beta-hydroxybutyrylation": "ft_mod_res:beta-hydroxybutyryl",
    "biotinylation": "ft_mod_res:biotinyl*",
    "blocked_amino_end": "ft_mod_res:Blocked amino end",
    "bromination": "ft_mod_res:bromo*", # new type
    "butyrylation": "ft_mod_res:butyryl*", # overlaps with beta-hydroxybutyrylation
    # "carbamidation": "",
    "carboxyethylation": "ft_mod_res:carboxyethyl",
    "carboxylation": "ft_mod_res:carboxy*", #carboxylic acid and carboxylysine and carboxyglutamate
    # "carboxymethylation": "",
    "cholesterol_ester": "ft_lipid:Cholesterol glycine ester",
    "citrullination": "ft_mod_res:citrulline",
    "crotonylation": "ft_mod_res:crotonyllysine",
    "cyclopeptide": "ft_crosslnk:cyclopeptide",
    "cysteinylation": "(ft_mod_res:cysteinyl OR ft_crosslnk:cysteinyl)", # new type, overlap with s-cysteinylation etc
    "c-linked_glycosylation": "ft_carbohyd:c-linked",
    "deamidation": "ft_mod_res:deamidated",
    "deamination": "ft_mod_res:allysine",
    "decanoylation": "ft_lipid:*decanoyl*",
    "decarboxylation": "ft_mod_res:decarboxylated",
    "dehydration": "(ft_mod_res:dehydro* OR ft_mod_res:didehydro*)", # new type
    "dephosphorylation": "ft_mod_res:dephospho",
    # "dietylphosphorylation": "",
    "disulfide_bond": "ft_disulfid:*", # maybe ft_chain: disulfide bond as well? or smotheing liek it
    "d-enantiomers": "ft_mod_res:d-", # new type, is applicable? (reported as modified residues on uniprot)
    "d-glucuronoylation": "ft_mod_res:d-glucuronoyl",
    "farnesylation": "ft_lipid:farnesyl",
    "formation_of_an_isopeptide_bond": "ft_crosslnk:isopeptide", # lots of overlap with the sumo, ubiq, nedd etc.
    "formylation": "ft_mod_res:formyl*",
    "gamma-carboxyglutamic_acid": "ft_chain:gamma-carboxyglutamic acid",
    "geranylation": "ft_lipid:geranyl", # new type
    "geranylgeranylation": "ft_lipid:geranylgeranyl", 
    "glutamylation": "ft_mod_res:glutamyl", # new type
    "glutarylation": "ft_mod_res:glutaryl*",
    "glutathionylation": "ft_mod_res:glutathionyl cysteine",
    "glycation": "ft_carbohyd:glycation",
    "gmpylation": "ft_mod_res:GMP", # new type
    "gpi-anchor": "ft_lipid:gpi-anchor",
    "histidylation": "ft_crosslnk:histidyl", # new type
    # "hmgylation": "",
    "hydroxyceramide_ester": "ft_lipid:hydroxyceramide glutamate ester",
    "hydroxylation": "(ft_mod_res:hydroxy* OR ft_mod_res:dihydroxy*)", # doesn't catch all dihydroxy
    "hypusine": "ft_mod_res:hypusine", # new type
    "imidazolation": "ft_crosslnk:imidazol", # new type
    "iodination": "ft_mod_res:iodotyrosine",
    "lactoylation": "ft_mod_res:lactoyllysine",
    # "lactylation": "",
    "lipoylation": "ft_mod_res:lipoyllysine",
    "malonylation": "ft_mod_res:malonyl*",
    "methylation": "ft_mod_res:methyl*",
    # "mgcylation": "",
    # "mgylation": "",
    "myristoylation": "ft_lipid:myristoyl*",
    "neddylation": "ft_crosslnk:NEDD8", 
    "nitration": "ft_mod_res:nitrotyrosine",
    "n-carbamoylation": "ft_mod_res:N-carbamoylalanine",
    "n-linked_glycosylation": "(ft_carbohyd:n-linked OR ft_carbohyd:n-beta-linked OR ft_carbohyd:n-alpha-linked)",
    "n-palmitoylation": "ft_lipid:n-palmitoyl",
    "octanoylation": "ft_lipid:octanoyl",
    "oxidation": "ft_mod_res:oxo* OR ft_mod_res:suldoxide OR ft_mod_res:sulfenic acid", #Tryptophylquinone? Methionine sulfoxide? Cysteine sulfenic acid?
    "o-linked_glycosylation": "(ft_carbohyd:o-linked OR ft_carbohyd:o-alpha-linked)",
    "o-palmitoleoylation": "ft_lipid:o-palmitoleoyl",
    "o-palmitoylation": "ft_lipid:o-palmitoyl", # overlap with palmitoylation
    "palmitoylation": "ft_lipid:palmitoyl", # new type
    "phosphatidylethanolamine_amidation": "ft_lipid:phosphatidylethanolamine amidated",
    # "phosphoglycerylation": "",
    "phosphorylation": "ft_mod_res:phospho*",
    "propionylation": "ft_mod_res:propionyl*",
    "pupylation": "ft_crosslnk:Pup",
    "pyridoxal_phosphate_addition": "ft_mod_res:pyridoxal phosphate",
    "pyrrolidone_carboxylic_acid": "ft_mod_res:pyrrolidone carboxylic acid",
    "pyrrolylation": "ft_mod_res:*pyrrolyl*",
    "pyruvate": "ft_mod_res:pyruvic acid",
    "serotonylation" : "ft_mod_res:serotonin",
    "stearoylation" : "ft_lipid:stearoyl",
    "succinylation": "ft_mod_res:succiny*", # previously succinyl*
    "sulfation" : "ft_mod_res:sulfo*",
    "sulfhydration" : "ft_mod_res:cysteine persulfide",
    "sulfoxidation" : "ft_mod_res:sulfoxide",
    "sumoylation": "ft_crosslnk:SUMO", # doesn't query "Modified residue (large scale data)" : "Sumoylated lysine"
    "s-archaeol" : "ft_lipid:s-archaeol",
    "s-carbamoylation" : "ft_mod_res:s-carbamoylcysteine",
    "s-cyanation" : "ft_mod_res:s-cyanocysteine",
    "s-cysteinylation" : "(ft_mod_res:s-cysteinyl OR ft_crosslnk:s-cysteinyl)",
    "s-diacylglycerol" : "ft_lipid:s-diacylglycerol",
    "s-linked_glycosylation": "ft_carbohyd:s-linked",
    "s-nitrosylation" : "ft_mod_res:s-nitrosocysteine", # overlaps with nitro*
    "s-palmitoylation" : "ft_lipid:s-palmitoyl",
    "thiocarboxylation" : "ft_mod_res:thioglycine",
    "ubiquitination": "ft_crosslnk:ubiquitin", # same problem as sumoylation
    "umpylation" : "ft_mod_res:UMP",
    "2-hydroxyisobutyrylation": "ft_mod_res:2-hydroxyisobutyryl", # overlap with butyryl*
    }

# code to run through the list and create the databases
for m, desc in MODIFICATIONS.items():
    #get_fasta_rest(desc, f"./dbs/{m}.fasta")
    get_fasta_rest(desc, f"/scratch/leuven/368/vsc36826/IBP/run4/{m}.fasta")