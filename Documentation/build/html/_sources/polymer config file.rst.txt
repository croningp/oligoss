Polymer Config File
###################


Section 1: General MS1 and Chemistry Rules
==========================================

- "MONOMERS" :

    - **Description:** dictionary of monomer one letter codes, along with associated neutral monoisotopic masses and reactivity classes (functional groups). This dictionary should specify the neutral monoisotopic mass for each monomer, its reactive functional groups, and the number of each reactive functional group. This information is essential to ensure that all theoretical sequences generated are chemically feasible.

    - **Format:**  {'X': [mass, [[rxn_class1, n], [rxn_class2, y]]} 
                   where X = monomer one letter code; mass = neutral monoisotopic mass; rxn_class1 and rxn_class2 = reactivity classes (functional groups - e.g. amine, aldehyde); n = number of rxn_class1 functional groups; y = number of rxn_class2 functional groups.

- "MASS_DIFF" :

    - **Description:** the mass difference when adding an additional monomer on to a polymer chain.

    - **Format:** float.

    - **Options:** For condensation polymer = H2O.

- "ELONGATION_UNIT" :

    - **Description:** the number of additional monomer units typically added when elongating a polymer.

    - **Format:** integer.

    - **Options:** This will typically be 1 for most polymers.

- "REACTIVITY_CLASSES" :

    - **Description:** dictionary of reactivity classes with associated compatible classes and monomers.

    - **Format:** {'classA' : [['classX', 'classY'], ['A', 'B']} 
                  where classA = reaction class, classX and classY = classes that are cross-reactive with classA, and A and B are monomers within reaction class classA.

- "SYMMETRY" :

    - **Description:** bool to define whether polymer is identical at both ends.

    - **Options:** true or false.

    - **Recommended:** Set to false for polymers with different termini (e.g. N- and C- termini for peptides), true for polymers with identical functional groups at both termini in linear chains.

- "CHAIN_TERMINATORS" :

    - **Description:** list of monomers that terminate chain elongation.

    - **Options:** null unless you are certain a certain monomer will stop elongation of the chain in your reactions.

- "LOSS_PRODUCTS" :

    -  **Description:** dictionary of monomer one letter codes and associated side chain neutral loss products, i.e. masses that can be lost from the monomer side chain.

- "IONIZABLE_SIDECHAINS" :

    - **Description:** dictionary of monomers that can be ionised with extra adducts at the SIDE CHAIN, with associated adducts, minimum and maximum absolute charge states in both positive and negative mode.

    - **Format:** {"X": {"pos": (adduct, a, b), "neg": (adduct, a, b)}}
                  where "X" = monomer one letter code, adduct = adduct string (must be found in either CATIONS OR ANIONS in GlobalChemicalConstantss), a = min side chain charge, b = max side chain charge for ionized form.

- "INTRINSICALLY_CHARGED_MONOMERS" :

    - **Description:** dictionary of monomers that have an intrinsic charge (i.e. charged without addition of adducts), with associated lists of permissible adducts.

Section 2: MS2 Fragmentation
============================

- **Description:** This section should contain all the information required to construct a basic MS2 fragment series for linear polymers. Each fragment type is defined as a key in the FRAG_SERIES dict. Typically, fragment keys are one letter codes used to denote fragments - e.g. the standard 'b' and 'y' fragment series for peptides.

- **Output:** When building fragment series, the fragment generator will use the following convention: f'{frag}{n}' where frag = fragment one letter code, n = number of monomers in fragment (e.g. b1, b2, b3, y1, y2, y3 etc...).

- "FRAG_SERIES" :

    - **Description:** dictionary of fragment series one letter codes and associated properties.

    - **Format:** Each fragment series key = fragment one letter code, value = subdictionary containing all relevant fragment properties.

    - "terminus" :

        - **Description:** end in which fragmentation starts.

        - **Options:** -1 (for end of the sequence) or 0 (for start of the sequence).

    - "mass_diff" :

         between the fragment and corresponding subsequence. This may vary depending on whether cations or anions are being fragmented, therefore separate mass_diffs are specified for positive mode ('pos') and negative mode ('neg') mass spec. Example: the 'y3' fragment of peptide sequence 'AGVS' = the mass of (GVS+H) in positive mode, and (GVS-H) in negative mode.

    - "fragmentation_unit" : 

        - **Description:** the minimum number of monomer units typically added and / or removed at a time when building a fragment series.

        - **Recommended:** Default is ELONGATION_UNIT.

- "MS2_SIGNATURE_IONS" :

    - **Description:** MS2 fragments which can be used as markers for monomers and / or small subsequences.

    - **Format:** nested dictionary of format {"type": {{"X"}: [signature_mass]}} 
                  where "type" is the signature ion type, "X" is the monomer and signature_mass is a float which corresponds to the mass of the signature ion for the monomer ("X").