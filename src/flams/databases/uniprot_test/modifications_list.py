
"""
taking the modification list from setup.py
making a dictionary with the uniprot query terms for each modification

to run on bash cd into the uniprot_test folder and run this file
"""
from uniprot import get_fasta_rest


# Here we store a dict of modifications that can be queried for.
    # sorted alphabetically
MODIFICATIONS = {
    "acetylation": "ft_mod_res:acetyl*",
    #     ModificationType(
    #     "acetylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Acetylation"), ModificationDatabase(dbptm,"Acetylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:acetyl*"], #acetyl*
    #     ["A","C","D","E","G","K","M","P","R","S","T","V","Y"]
    # ), # works both with keywrod and ft_mod_res
    "adp-ribosylation": "ft_mod_res:ADP",
    #     ModificationType(
    #     "adp-ribosylation", 1.4,
    #     #[ModificationDatabase(dbptm, "ADP-ribosylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:ADP-ribosyl"],  
    #     ["C","D","E","G","H","K","N","R","S","Y"]
    # ), # works with ft_mod_res:ADP-rybosil, doesn't work with ADP-rybosil*
    "amidation": "ft_mod_res:amide",
        #     ModificationType(
        #     "amidation", 1.4,
        #     #[ModificationDatabase(dbptm, "Amidation")],
        #     [ModificationDatabase(uniprot), "ft_mod_res:amide"],     # KW-0027, amide
        #     ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
        # ), # ft_mod_res:amide
    # "ampylation" : ModificationType(
    #     "ampylation", 1.4,
    #     #[ModificationDatabase(dbptm, "AMPylation")],
    #     [ModificationDatabase(uniprot), ""],  
    #     ["S","T","Y"]
    # ), # ????
    "benzoylation": "ft_mod_res:benzoyl*",
        #     ModificationType(
        #     "benzoylation", 1.4, 
        #     #[ModificationDatabase(cplmv4, "Benzoylation")],
        #     [ModificationDatabase(uniprot), "ft_mod_res:benzoyl*"], #benzoyl
        #     ["K"]
        # ), # benzoyl* works (finds 1)
    # example - P10412 doesnt match
    "beta-hydroxybutyrylation": "ft_mod_res:beta-hydroxybutyryl",
    #     ModificationType(
    #     "beta-hydroxybutyrylation", 1.4, 
    #     #[ModificationDatabase(cplmv4, "β-Hydroxybutyrylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:beta-hydroxybutyryl"], #beta-hydroxybutyryl
    #     ["K"]
    # ), #beta-hydroxybutyryl works
    "biotinylation": "ft_mod_res:biotinyl*",
    #     ModificationType(
    #     "biotinylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Biotinylation"), ModificationDatabase(dbptm, "Biotinylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:biotinyl*"],
    #     ["K"]
    # ), # biotinyl* works
    "blocked_amino_end": "ft_mod_res:Blocked amino end",
        # ModificationType(
        # "blocked_amino_end", 1.4,
        # #[ModificationDatabase(dbptm, "Blocked amino end")],
        # [ModificationDatabase(uniprot), "ft_mod_res:Blocked amino end"],  #blocked amino end
        # ["A","C","D","E","G","H","I","L","M","N","P","Q","R","S","T","V"]
    # ), #"Blocked amino end" works
    "butyrylation": "ft_mod_res:butyryl*", # includes the beta-hydroxybutyrylation
    #     ModificationType(
    #     "butyrylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Butyrylation"), ModificationDatabase(dbptm, "Butyrylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:butyryl*"],
    #     ["K"]
    # ), # butyryl* works, could overlap with the hydroxybutyryl..
    # "carbamidation": ModificationType(
    #     "carbamidation", 1.4,
    #     #[ModificationDatabase(dbptm, "Carbamidation")],
    #     [ModificationDatabase(uniprot), ""],
    #     ["C"]
    # ), # ??? don't think uniprot has this one
    "carboxyethylation": "ft_mod_res:carboxyethyl",
    #     ModificationType(
    #     "carboxyethylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Carboxyethylation"), ModificationDatabase(dbptm, "Carboxyethylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:carboxyethyl"],
    #     ["K"]
    # ), # works but does ity find all?
    "carboxylation": "ft_mod_res:carboxylysine",
    #     ModificationType(
    #     "carboxylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Carboxylation"), ModificationDatabase(dbptm, "Carboxylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:carboxylysine"],
    #     ["K"]
    # ),
    # "carboxymethylation": "ft_mod_res:carboxymethyl",
    # #     ModificationType(
    # #     "carboxymethylation", 1.4, 
    # #     #[ModificationDatabase(cplmv4, "Carboxymethylation")],
    # #     [ModificationDatabase(uniprot), "ft_mod_res:carboxymethyl"],
    # #     ["K"]
    # # ), #??? doesnt exist?
    "cholesterol_ester": "ft_lipid:Cholesterol glycine ester",
    #     ModificationType(
    #     "cholesterol_ester", 1.4,
    #     #[ModificationDatabase(dbptm, "Cholesterol ester")],
    #     [ModificationDatabase(uniprot), "ft_lipid:Cholesterol glycine ester"],
    #     ["G"]
    # ),
    "citrullination": "ft_mod_res:citrulline",
    #     ModificationType(
    #     "citrullination", 1.4,
    #     #[ModificationDatabase(dbptm, "Citrullination")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:citrulline"],
    #     ["R"]
    # ),
    "crotonylation": "ft_mod_res:crotonyllysine",
    #     ModificationType(
    #     "crotonylation", 1.4,
    #     #[ModificationDatabase(cplmv4, "Crotonylation"), ModificationDatabase(dbptm, "Crotonylation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:crotonyllysine"],
    #     ["K"]
    # ),
    "c-linked_glycosylation": "ft_carbohyd:c-linked",
    "deamidation": "ft_mod_res:deamidated",
    #     ModificationType(
    #     "deamidation", 1.4,
    #     #[ModificationDatabase(dbptm, "Deamidation")],
    #     [ModificationDatabase(uniprot), "ft_mod_res:deamidated"],
    #     ["N","Q"]
    # ),
    "deamination": "ft_mod_res:allysine",
        #ModificationType(
    #     "deamination", 1.4,
    #     [ModificationDatabase(dbptm, "Deamination")],
    #     ["K"]
    # ),
    "decanoylation": "ft_lipid:*decanoyl*",
    #ModificationType(
    #     "decanoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Decanoylation")],
    #     ["S","T"]
    # ),
    "decarboxylation": "ft_mod_res:decarboxylated",
        #ModificationType(
    #     "decarboxylation", 1.4,
    #     [ModificationDatabase(dbptm, "Decarboxylation")],
    #     ["D","T"]
    # ),
    "dephosphorylation": "ft_mod_res:dephospho",
        #ModificationType(
    #     "dephosphorylation", 1.4,
    #     [ModificationDatabase(dbptm, "Dephosphorylation")],
    #     ["S","T","Y"]
    # ),
    # "dietylphosphorylation": "ft_mod_res:dietylphosphorylated",
    # #ModificationType(
    # #     "dietylphosphorylation", 1.4, [ModificationDatabase(cplmv4, "Dietylphosphorylation")],
    # #     ["K"] #OK
    # # ), #??? not on uniprot
    "disulfide_bond": "ft_disulfid:*", # maybe ft_chain: disulfide bond as well? or smotheing liek it
        #ModificationType(
    #     "disulfide_bond", 1.4,
    #     [ModificationDatabase(dbptm, "Disulfide bond")],
    #     ["C"]
    # ),
    "d-glucuronoylation": "ft_mod_res:d-glucuronoyl",
    #ModificationType(
    #     "d-glucuronoylation", 1.4,
    #     [ModificationDatabase(dbptm, "D-glucuronoylation")],
    #     ["G"]
    # ),
    "farnesylation": "ft_lipid:farnesyl",
    #ModificationType(
    #     "farnesylation", 1.4,
    #     [ModificationDatabase(dbptm, "Farnesylation")],
    #     ["C"]
    # ),
    "formation_of_an_isopeptide_bond": "ft_crosslnk:isopeptide", # lots of overlap with the sumo, ubiq, nedd
    #ModificationType(
    #     "formation_of_an_isopeptide_bond", 1.4,
    #     [ModificationDatabase(dbptm, "Formation of an isopeptide bond")],
    #     ["E","Q"]
    # ),
    "formylation": "ft_mod_res:formyl*",
    #ModificationType(
    #     "formylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Formylation"), ModificationDatabase(dbptm, "Formylation")],
    #     ["G","K","M"]
    # ),
    "gamma-carboxyglutamic_acid": "ft_chain:gamma-carboxyglutamic acid",
    #ModificationType(
    #     "gamma-carboxyglutamic_acid", 1.4,
    #     [ModificationDatabase(dbptm, "Gamma-carboxyglutamic acid")],
    #     ["E"]
    # ),
    "geranylgeranylation": "ft_lipid:geranylgeranyl cysteine",
    #ModificationType(
    #     "geranylgeranylation", 1.4,
    #     [ModificationDatabase(dbptm, "Geranylgeranylation")],
    #     ["C"]
    # ),
    "glutarylation": "ft_mod_res:glutaryllysine",
    #ModificationType(
    #     "glutarylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Glutarylation"), ModificationDatabase(dbptm, "Glutarylation")],
    #     ["K"]
    # ),
    "glutathionylation": "ft_mod_res:glutathionyl cysteine",
    #ModificationType(
    #     "glutathionylation", 1.4,
    #     [ModificationDatabase(dbptm, "Glutathionylation")],
    #     ["C"]
    # ),
    "glycation": "ft_carbohyd:glycation",
    #ModificationType(
    #     "glycation", 1.4, [ModificationDatabase(cplmv4, "Glycation")],
    #     ["K"]
    # ),
    "gpi-anchor": "ft_lipid:gpi-anchor",
    #ModificationType(
    #     "gpi-anchor", 1.4, [ModificationDatabase(dbptm, "GPI-anchor")],
    #     ["A","C","D","G","N","S","T"]
    # ),
    # "hmgylation": 
    # #ModificationType(
    # #     "hmgylation", 1.4, [ModificationDatabase(cplmv4, "HMGylation")],
    # #     ["K"] #OK
    # # ),
    "hydroxyceramide_ester": "ft_lipid:hydroxyceramide glutamate ester",
    # ModificationType(
    #     "hydroxyceramide_ester", 1.4, [ModificationDatabase(dbptm, "Hydroxyceramide ester")],
    #     ["Q"]
    # ),
    "hydroxylation": "ft_mod_res:hydroxy*",
    # ModificationType(
    #     "hydroxylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Hydroxylation"), ModificationDatabase(dbptm, "Hydroxylation")],
    #     ["C","D","E","F","H","I","K","L","N","P","R","S","T","V","W","Y"]
    # ),
    "iodination": "ft_mod_res:iodotyrosine",
    # ModificationType(
    #     "iodination", 1.4,
    #     [ModificationDatabase(dbptm, "Iodination")],
    #     ["Y"]
    # ),
    "lactoylation": "ft_mod_res:lactoyllysine",
    # ModificationType(
    #     "lactoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Lactoylation")],
    #     ["K"]
    # ),
    # "lactylation": "",
    # # ModificationType(
    # #     "lactylation", 1.4,
    # #     [ModificationDatabase(cplmv4, "Lactylation"), ModificationDatabase(dbptm, "Lactylation")],
    # #     ["K"]
    # # ),
    "lipoylation": "ft_mod_res:lipoyllysine",
    # ModificationType(
    #     "lipoylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Lipoylation"), ModificationDatabase(dbptm, "Lipoylation")],
    #     ["K"]
    # ),
    "malonylation": "ft_mod_res:malonyl*",
    # ModificationType(
    #     "malonylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Malonylation"), ModificationDatabase(dbptm, "Malonylation")],
    #     ["K"]
    # ),
    "methylation": "ft_mod_res:methyl*",
    # ModificationType(
    #     "methylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Methylation"), ModificationDatabase(dbptm, "Methylation")],
    #     ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","Y"]
    # ),
    # "mgcylation": "",
    # # ModificationType(
    # #     "mgcylation", 1.4, [ModificationDatabase(cplmv4, "MGcylation")],
    # #     ["K"]
    # # ),
    # "mgylation": "",
    # # ModificationType(
    # #     "mgylation", 1.4, [ModificationDatabase(cplmv4, "MGylation")],
    # #     ["K"]
    # # ),
    "myristoylation": "ft_lipid:myristoyl*",
    # ModificationType(
    #     "myristoylation", 1.4, [ModificationDatabase(dbptm, "Myristoylation")],
    #     ["C","G","K"]
    # ),
    "neddylation": "ft_crosslnk:NEDD8", # might not be complete
    # ModificationType(
    #     "neddylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Neddylation"), ModificationDatabase(dbptm, "Neddylation")],
    #     ["K"]
    # ),
    "nitration": "ft_mod_res:nitrotyrosine",
    # ModificationType(
    #     "nitration", 1.4, [ModificationDatabase(dbptm, "Nitration")],
    #     ["Y"]
    # ),
    "n-carbamoylation": "ft_mod_res:N-carbamoylalanine",
    # ModificationType(
    #     "n-carbamoylation", 1.4, [ModificationDatabase(dbptm, "N-carbamoylation")],
    #     ["A"]
    # ),
    "n-linked_glycosylation": "ft_carbohyd:n-linked",
    # ModificationType(
    #     "n-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "N-linked Glycosylation")],
    #     ["D","I","K","N","R","S","T","V","W"]
    # ),
    "n-palmitoylation": "ft_lipid:n-palmitoyl",
    # ModificationType(
    #     "n-palmitoylation", 1.4, [ModificationDatabase(dbptm, "N-palmitoylation")],
    #     ["C","G","K"]
    # ),
    "octanoylation": "ft_lipid:octanoyl",
    # ModificationType(
    #     "octanoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Octanoylation")],
    #     ["S","T"]
    # ),
    # "oxidation": "", #Tryptophylquinone? Methionine sulfoxide? Cysteine sulfenic acid?
    # # ModificationType(
    # #     "oxidation", 1.4,
    # #     [ModificationDatabase(dbptm, "Oxidation")],
    # #     ["C","L","M","S","W"]
    # # ),
    "o-linked_glycosylation": "ft_carbohyd:o-linked",
    # ModificationType(
    #     "o-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-linked Glycosylation")],
    #     ["K","P","S","T","Y"]
    # ),
    "o-palmitoleoylation": "ft_lipid:o-palmitoleoyl",
    # ModificationType(
    #     "o-palmitoleoylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-palmitoleoylation")],
    #     ["S"]
    # ),
    "o-palmitoylation": "ft_lipid:o-palmitoyl",
    # ModificationType(
    #     "o-palmitoylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-palmitoylation")],
    #     ["S","T"]
    # ),
    "phosphatidylethanolamine_amidation": "ft_lipid:phosphatidylethanolamine amidated",
    # ModificationType(
    #     "phosphatidylethanolamine_amidation", 1.4, [ModificationDatabase(dbptm, "Phosphatidylethanolamine amidation")],
    #     ["G"]
    # ),
    # "phosphoglycerylation": "";
    # # ModificationType(
    # #     "phosphoglycerylation", 1.4,
    # #     [ModificationDatabase(cplmv4, "Phosphoglycerylation")],
    # #     ["K"]
    # # ),
    "phosphorylation": "ft_mod_res:phospho*",
    # ModificationType(
    #     "phosphorylation", 1.4,
    #     [ModificationDatabase(dbptm, "Phosphorylation")],
    #     ["A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","W","Y"]
    # ),
    "propionylation": "ft_mod_res:propionyl*",
    # ModificationType(
    #     "propionylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Propionylation"), ModificationDatabase(dbptm, "Propionylation")],
    #     ["K"]
    # ),
    "pupylation": "ft_crosslnk:Pup",
    # ModificationType(
    #     "pupylation", 1.4, [ModificationDatabase(cplmv4, "Pupylation")],
    #     ["K"]
    # ),
    "pyrrolidone_carboxylic_acid": "ft_mod_res:pyrrolidone carboxylic acid",
    # ModificationType(
    #     "pyrrolidone_carboxylic_acid", 1.4, [ModificationDatabase(dbptm, "Pyrrolidone carboxylic acid")],
    #     ["E","Q"]
    # ),
    "pyrrolylation": "ft_mod_res:*pyrrolyl*",
    # ModificationType(
    #     "pyrrolylation", 1.4, [ModificationDatabase(dbptm, "Pyrrolylation")],
    #     ["C"]
    # ),
    "pyruvate": "ft_mod_res:pyruvic acid",
    # ModificationType(
    #     "pyruvate", 1.4, [ModificationDatabase(dbptm, "Pyruvate")],
    #     ["C","S"]
    # ),
    "serotonylation" : "ft_mod_res:serotonin",
    # ModificationType(
    #     "serotonylation", 1.4,
    #     [ModificationDatabase(dbptm, "Serotonylation")],
    #     ["Q"]
    # ),
    "stearoylation" : "ft_lipid:stearoyl",
    # ModificationType(
    #     "stearoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Stearoylation")],
    #     ["C"]
    # ),
    "succinylation": "succinyl*",
    # ModificationType(
    #     "succinylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Succinylation"), ModificationDatabase(dbptm, "Succinylation")],
    #     ["C","K","W"]
    # ),
    "sulfation" : "ft_mod_res:sulfo*",
    # ModificationType(
    #     "sulfation", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfation")],
    #     ["C","S","T","Y"]
    # ),
    "sulfhydration" : "ft_mod_res:cysteine persulfide",
    # ModificationType(
    #     "sulfhydration", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfhydration")],
    #     ["C"]
    # ),
    "sulfoxidation" : "ft_mod_res:sulfoxide",
    # ModificationType(
    #     "sulfoxidation", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfoxidation")],
    #     ["M"]
    # ),
    "sumoylation": "ft_crosslnk:SUMO", # how to query the sumoylated lysine???
    # ModificationType(
    #     "sumoylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Sumoylation"), ModificationDatabase(dbptm, "Sumoylation")],
    #     ["K"]
    # ),
    "s-archaeol" : "ft_lipid:s-archaeol",
    # ModificationType(
    #     "s-archaeol", 1.4,
    #     [ModificationDatabase(dbptm, "S-archaeol")],
    #     ["C"]
    # ),
    "s-carbamoylation" : "ft_mod_res:s-carbamoylcysteine",
    # ModificationType(
    #     "s-carbamoylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-carbamoylation")],
    #     ["C"]
    # ),
    "s-cyanation" : "ft_mod_res:s-cyanocysteine",
    # ModificationType(
    #     "s-cyanation", 1.4,
    #     [ModificationDatabase(dbptm, "S-Cyanation")],
    #     ["C"]
    # ),
    "s-cysteinylation" : "ft_mod_res:s-cysteinyl",
    # ModificationType(
    #     "s-cysteinylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-cysteinylation")],
    #     ["C"]
    # ),
    "s-diacylglycerol" : "ft_lipid:s-diacylglycerol",
    # ModificationType(
    #     "s-diacylglycerol", 1.4,
    #     [ModificationDatabase(dbptm, "S-diacylglycerol")],
    #     ["C"]
    # ),
    "s-linked_glycosylation": "ft_carbohyd:s-linked",
    # ModificationType(
    #     "s-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-linked Glycosylation")],
    #     ["C"]
    # ),
    "s-nitrosylation" : "ft_mod_res:s-nitrosocysteine", # overlaps with nitro*?
    # ModificationType(
    #     "s-nitrosylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-nitrosylation")],
    #     ["C"]
    # ),
    "s-palmitoylation" : "ft_lipid:s-palmitoyl",
    # ModificationType(
    #     "s-palmitoylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-palmitoylation")],
    #     ["C"]
    # ),
    "thiocarboxylation" : "ft_mod_res:thioglycine",
    # ModificationType(
    #     "thiocarboxylation", 1.4,
    #     [ModificationDatabase(dbptm, "Thiocarboxylation")],
    #     ["G"]
    # ),
    "ubiquitination": "ft_crosslnk:ubiquitin", # same problem as sumoylation
    # # ModificationType(
    # #     "ubiquitination", 1.4,
    # #     [ModificationDatabase(cplmv4, "Ubiquitination"), ModificationDatabase(dbptm, "Ubiquitination")],
    # #     ["C","K","R","S"]
    # # ),
    "umpylation" : "ft_mod_res:UMP",
    # ModificationType(
    #     "umpylation", 1.4,
    #     [ModificationDatabase(dbptm, "UMPylation")],
    #     ["S","T","Y"]
    # ),
    "2-hydroxyisobutyrylation": "ft_mod_res:2-hydroxyisobutyryl",
    # ModificationType(
    #     "2-hydroxyisobutyrylation", 1.4, [ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation")],
    #     ["K"]
    # ),
    }

# code to run through the list and create the databases
for m, desc in MODIFICATIONS.items():
    #get_fasta_rest(desc, f"./dbs/{m}.fasta")
    get_fasta_rest(desc, f"/scratch/leuven/368/vsc36826/IBP/run2/{m}.fasta")