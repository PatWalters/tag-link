#!/usr/bin/env python

"""Usage: tag.py --in INFILE_NAME --out OUTFILE_NAME --smarts SMARTS --tag TAG

--in INFILE_NAME input SMILES file
--out OUTFILE_NAME output SMILES file
--smarts SMARTS SMARTS expression for atom to be tagged
--tag TAG tag to add
"""

import sys
from rdkit import Chem
from docopt import docopt
from tqdm import tqdm


def tag_molecule(mol, query, replacement_type):
    result_smiles_list = []
    match_list = mol.GetSubstructMatches(query)
    for idx in [x[0] for x in match_list]:
        new_mol = Chem.RWMol(mol)
        new_idx = new_mol.AddAtom(Chem.Atom(replacement_type))
        new_mol.AddBond(idx, new_idx, order=Chem.rdchem.BondType.SINGLE)
        result_smiles_list.append(Chem.MolToSmiles(new_mol, True))
    return result_smiles_list


def main(arg_string):
    cmd_input = docopt(__doc__)
    input_file = cmd_input.get("--in")
    output_file = cmd_input.get("--out")
    smarts = cmd_input.get("--smarts")
    tag = cmd_input.get("--tag")

    tag_mol = Chem.MolFromSmiles(tag)
    if tag_mol is None:
        print(f"Could not parse tag {tag}", file=sys.stderr)
        sys.exit(1)
    tag_atomic_num = -1
    for tag_atm in tag_mol.GetAtoms():
        tag_atomic_num = tag_atm.GetAtomicNum()
        break
    if tag_atomic_num == -1:
        print(f"Could not parse tag atom {tag}")
        sys.exit(1)

    suppl = Chem.SmilesMolSupplier(input_file)
    query_mol = Chem.MolFromSmarts(smarts)
    if query_mol is None:
        print(f"Could not parse SMARTS {smarts}")
        sys.exit(0)

    ofs = open(output_file, "w")
    tagged_smiles_set = set()
    for in_mol in tqdm(suppl):
        tagged_smiles_list = tag_molecule(in_mol, query_mol, 79)
        for smiles in tagged_smiles_list:
            if smiles not in tagged_smiles_set:
                tagged_smiles_set.add(smiles)
                print(smiles, file=ofs)


main(__doc__)