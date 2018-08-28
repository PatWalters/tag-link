#!/usr/bin/env python

"""Usage:
tag.py --in INFILE_NAME --out OUTFILE_NAME --smarts SMARTS --tag TAG [--replace]

--in INFILE_NAME input SMILES file
--out OUTFILE_NAME output SMILES file
--smarts SMARTS SMARTS expression for atom to be tagged
--tag TAG tag to add
"""

import sys
from rdkit import Chem
from docopt import docopt
from tqdm import tqdm


def tag_molecule(mol, query, atomic_num, replace):
    if replace:
        return replace_fragment(mol, query, atomic_num)
    else:
        return attach_atom(mol, query, atomic_num)


def attach_atom(mol, query, atomic_num):
    result_smiles_list = []
    match_list = mol.GetSubstructMatches(query)
    for idx in [x[0] for x in match_list]:
        new_mol = Chem.RWMol(mol)
        new_idx = new_mol.AddAtom(Chem.Atom(atomic_num))
        new_mol.AddBond(idx, new_idx, order=Chem.rdchem.BondType.SINGLE)
        result_smiles_list.append(Chem.MolToSmiles(new_mol, True))
    return result_smiles_list


def get_neighbor_atoms(atm):
    return [x.GetOtherAtom(atm).GetIdx() for x in atm.GetBonds()]


def replace_fragment(mol,query,atomic_num):
    new_mol = Chem.RWMol(mol)
    match_list = new_mol.GetSubstructMatches(query)
    if len(match_list) == 0:
        return None
    for idx in match_list[0]:
        atm = new_mol.GetAtomWithIdx(idx)
        nbr_list = get_neighbor_atoms(atm)
        good_neighbors = [x for x in nbr_list if x not in match_list[0]]
        if len(good_neighbors):
            nbr_atm = good_neighbors[0]
            new_idx = new_mol.AddAtom(Chem.Atom(atomic_num))
            new_mol.AddBond(nbr_atm, new_idx, order=Chem.rdchem.BondType.SINGLE)
            new_mol = Chem.DeleteSubstructs(new_mol,query)
            return [Chem.MolToSmiles(new_mol)]
    return None


def main(arg_string):
    cmd_input = docopt(__doc__)
    input_file = cmd_input.get("--in")
    output_file = cmd_input.get("--out")
    smarts = cmd_input.get("--smarts")
    tag = cmd_input.get("--tag")
    do_replace = cmd_input.get("--replace")

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
        tagged_smiles_list = tag_molecule(in_mol, query_mol, 79, do_replace)
        for smiles in tagged_smiles_list:
            if smiles not in tagged_smiles_set:
                tagged_smiles_set.add(smiles)
                print(smiles, file=ofs)


main(__doc__)