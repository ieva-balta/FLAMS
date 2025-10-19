
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
    "adp-ribosylation": "ft_mod_res:ADP-ribosyl",
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
    "butyrylation": "ft_mod_res:butyryl*",
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
    "formation_of_an_isopeptide_bond": "ft_crosslnk:isopeptide",
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
    # "glycation": ModificationType(
    #     "glycation", 1.4, [ModificationDatabase(cplmv4, "Glycation")],
    #     ["K"]
    # ),
    # "gpi-anchor": ModificationType(
    #     "gpi-anchor", 1.4, [ModificationDatabase(dbptm, "GPI-anchor")],
    #     ["A","C","D","G","N","S","T"]
    # ),
    # "hmgylation": ModificationType(
    #     "hmgylation", 1.4, [ModificationDatabase(cplmv4, "HMGylation")],
    #     ["K"] #OK
    # ),
    # "hydroxyceramide_ester": ModificationType(
    #     "hydroxyceramide_ester", 1.4, [ModificationDatabase(dbptm, "Hydroxyceramide ester")],
    #     ["Q"]
    # ),
    # "hydroxylation": ModificationType(
    #     "hydroxylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Hydroxylation"), ModificationDatabase(dbptm, "Hydroxylation")],
    #     ["C","D","E","F","H","I","K","L","N","P","R","S","T","V","W","Y"]
    # ),
    # "iodination": ModificationType(
    #     "iodination", 1.4,
    #     [ModificationDatabase(dbptm, "Iodination")],
    #     ["Y"]
    # ),
    # "lactoylation": ModificationType(
    #     "lactoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Lactoylation")],
    #     ["K"]
    # ),
    # "lactylation": ModificationType(
    #     "lactylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Lactylation"), ModificationDatabase(dbptm, "Lactylation")],
    #     ["K"]
    # ),
    # "lipoylation": ModificationType(
    #     "lipoylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Lipoylation"), ModificationDatabase(dbptm, "Lipoylation")],
    #     ["K"]
    # ),
    # "malonylation": ModificationType(
    #     "malonylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Malonylation"), ModificationDatabase(dbptm, "Malonylation")],
    #     ["K"]
    # ),
    # "methylation": ModificationType(
    #     "methylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Methylation"), ModificationDatabase(dbptm, "Methylation")],
    #     ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","Y"]
    # ),
    # "mgcylation": ModificationType(
    #     "mgcylation", 1.4, [ModificationDatabase(cplmv4, "MGcylation")],
    #     ["K"]
    # ),
    # "mgylation": ModificationType(
    #     "mgylation", 1.4, [ModificationDatabase(cplmv4, "MGylation")],
    #     ["K"]
    # ),
    # "myristoylation": ModificationType(
    #     "myristoylation", 1.4, [ModificationDatabase(dbptm, "Myristoylation")],
    #     ["C","G","K"]
    # ),
    # "neddylation": ModificationType(
    #     "neddylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Neddylation"), ModificationDatabase(dbptm, "Neddylation")],
    #     ["K"]
    # ),
    # "nitration": ModificationType(
    #     "nitration", 1.4, [ModificationDatabase(dbptm, "Nitration")],
    #     ["Y"]
    # ),
    # "n-carbamoylation": ModificationType(
    #     "n-carbamoylation", 1.4, [ModificationDatabase(dbptm, "N-carbamoylation")],
    #     ["A"]
    # ),
    # "n-linked_glycosylation": ModificationType(
    #     "n-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "N-linked Glycosylation")],
    #     ["D","I","K","N","R","S","T","V","W"]
    # ),
    # "n-palmitoylation": ModificationType(
    #     "n-palmitoylation", 1.4, [ModificationDatabase(dbptm, "N-palmitoylation")],
    #     ["C","G","K"]
    # ),
    # "octanoylation": ModificationType(
    #     "octanoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Octanoylation")],
    #     ["S","T"]
    # ),
    # "oxidation": ModificationType(
    #     "oxidation", 1.4,
    #     [ModificationDatabase(dbptm, "Oxidation")],
    #     ["C","L","M","S","W"]
    # ),
    # "o-linked_glycosylation": ModificationType(
    #     "o-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-linked Glycosylation")],
    #     ["K","P","S","T","Y"]
    # ),
    # "o-palmitoleoylation": ModificationType(
    #     "o-palmitoleoylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-palmitoleoylation")],
    #     ["S"]
    # ),
    # "o-palmitoylation": ModificationType(
    #     "o-palmitoylation", 1.4,
    #     [ModificationDatabase(dbptm, "O-palmitoylation")],
    #     ["S","T"]
    # ),
    # "phosphatidylethanolamine_amidation": ModificationType(
    #     "phosphatidylethanolamine_amidation", 1.4, [ModificationDatabase(dbptm, "Phosphatidylethanolamine amidation")],
    #     ["G"]
    # ),
    # "phosphoglycerylation": ModificationType(
    #     "phosphoglycerylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Phosphoglycerylation")],
    #     ["K"]
    # ),
    # "phosphorylation": ModificationType(
    #     "phosphorylation", 1.4,
    #     [ModificationDatabase(dbptm, "Phosphorylation")],
    #     ["A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","W","Y"]
    # ),
    # "propionylation": ModificationType(
    #     "propionylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Propionylation"), ModificationDatabase(dbptm, "Propionylation")],
    #     ["K"]
    # ),
    # "pupylation": ModificationType(
    #     "pupylation", 1.4, [ModificationDatabase(cplmv4, "Pupylation")],
    #     ["K"]
    # ),
    # "pyrrolidone_carboxylic_acid": ModificationType(
    #     "pyrrolidone_carboxylic_acid", 1.4, [ModificationDatabase(dbptm, "Pyrrolidone carboxylic acid")],
    #     ["E","Q"]
    # ),
    # "pyrrolylation": ModificationType(
    #     "pyrrolylation", 1.4, [ModificationDatabase(dbptm, "Pyrrolylation")],
    #     ["C"]
    # ),
    # "pyruvate": ModificationType(
    #     "pyruvate", 1.4, [ModificationDatabase(dbptm, "Pyruvate")],
    #     ["C","S"]
    # ),
    # "serotonylation" : ModificationType(
    #     "serotonylation", 1.4,
    #     [ModificationDatabase(dbptm, "Serotonylation")],
    #     ["Q"]
    # ),
    # "stearoylation" : ModificationType(
    #     "stearoylation", 1.4,
    #     [ModificationDatabase(dbptm, "Stearoylation")],
    #     ["C"]
    # ),
    # "succinylation": ModificationType(
    #     "succinylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Succinylation"), ModificationDatabase(dbptm, "Succinylation")],
    #     ["C","K","W"]
    # ),
    # "sulfation" : ModificationType(
    #     "sulfation", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfation")],
    #     ["C","S","T","Y"]
    # ),
    # "sulfhydration" : ModificationType(
    #     "sulfhydration", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfhydration")],
    #     ["C"]
    # ),
    # "sulfoxidation" : ModificationType(
    #     "sulfoxidation", 1.4,
    #     [ModificationDatabase(dbptm, "Sulfoxidation")],
    #     ["M"]
    # ),
    # "sumoylation": ModificationType(
    #     "sumoylation", 1.4,
    #     [ModificationDatabase(cplmv4, "Sumoylation"), ModificationDatabase(dbptm, "Sumoylation")],
    #     ["K"]
    # ),
    # "s-archaeol" : ModificationType(
    #     "s-archaeol", 1.4,
    #     [ModificationDatabase(dbptm, "S-archaeol")],
    #     ["C"]
    # ),
    # "s-carbamoylation" : ModificationType(
    #     "s-carbamoylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-carbamoylation")],
    #     ["C"]
    # ),
    # "s-cyanation" : ModificationType(
    #     "s-cyanation", 1.4,
    #     [ModificationDatabase(dbptm, "S-Cyanation")],
    #     ["C"]
    # ),
    # "s-cysteinylation" : ModificationType(
    #     "s-cysteinylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-cysteinylation")],
    #     ["C"]
    # ),
    # "s-diacylglycerol" : ModificationType(
    #     "s-diacylglycerol", 1.4,
    #     [ModificationDatabase(dbptm, "S-diacylglycerol")],
    #     ["C"]
    # ),
    # "s-linked_glycosylation": ModificationType(
    #     "s-linked_glycosylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-linked Glycosylation")],
    #     ["C"]
    # ),
    # "s-nitrosylation" : ModificationType(
    #     "s-nitrosylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-nitrosylation")],
    #     ["C"]
    # ),
    # "s-palmitoylation" : ModificationType(
    #     "s-palmitoylation", 1.4,
    #     [ModificationDatabase(dbptm, "S-palmitoylation")],
    #     ["C"]
    # ),
    # "thiocarboxylation" : ModificationType(
    #     "thiocarboxylation", 1.4,
    #     [ModificationDatabase(dbptm, "Thiocarboxylation")],
    #     ["G"]
    # ),
    # "ubiquitination": ModificationType(
    #     "ubiquitination", 1.4,
    #     [ModificationDatabase(cplmv4, "Ubiquitination"), ModificationDatabase(dbptm, "Ubiquitination")],
    #     ["C","K","R","S"]
    # ),
    # "umpylation" : ModificationType(
    #     "umpylation", 1.4,
    #     [ModificationDatabase(dbptm, "UMPylation")],
    #     ["S","T","Y"]
    # ),
    # "2-hydroxyisobutyrylation": ModificationType(
    #     "2-hydroxyisobutyrylation", 1.4, [ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation")],
    #     ["K"]
    # ),
    }

# code to run through the list and create the databases
for m, desc in MODIFICATIONS.items():
    get_fasta_rest(desc, f"./dbs/{m}.fasta")