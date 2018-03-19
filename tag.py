#!/usr/bin/env python

"""Usage: tag.py --in INFILE_NAME --out OUTFILE_NAME --smarts SMARTS --tag TAG

--in INFILE_NAME input SMILES file
--out OUTFILE_NAME output SMILES file
--smarts SMARTS SMARTS expression for atom to be tagged
--tag TAG tag to add
"""

import sys
from openeye.oechem import *
from docopt import docopt
from tqdm import tqdm


def adjust_hydrogens(mol):
    for atm in mol.GetAtoms():
        h_count = OEDefaultMDLHCount(atm)
        if h_count < 0:
            h_count = 0
        atm.SetImplicitHCount(h_count)


def match_smarts(mol, pat):
    match_list = []
    for p in pat.Match(mol):
        for pair in p.GetAtoms():
            match_list.append(pair.target)
    return match_list


def has_open_valence(atm):
    return atm.GetImplicitHCount() >= 1


def add_tag_atom(mol,tag):
    new_mol = OEGraphMol(mol)
    for atm in mol.GetAtoms():
        atm.SetMapIdx(0)
    base_atm = [atm for atm in new_mol.GetAtoms(OEHasMapIdx(1))][0]
    tag_mol = OEGraphMol()
    OEParseSmiles(tag_mol,tag)
    tag_atoms = [atm for atm in tag_mol.GetAtoms()]
    link_atm = tag_atoms[0]
    link_atm.SetMapIdx(1)
    OEAddMols(new_mol,tag_mol)
    atm1, atm2 = [atm for atm in new_mol.GetAtoms(OEHasMapIdx(1))]
    new_mol.NewBond(atm1, atm2, 1)
    for atm in new_mol.GetAtoms():
        atm.SetMapIdx(0)
    adjust_hydrogens(new_mol)
    return OECreateIsoSmiString(new_mol)


def tag_molecule(mol, smarts, tag):
    used = set()
    pat = OESubSearch(smarts)
    matches = match_smarts(mol, pat)
    for atm in matches:
        if has_open_valence(atm):
            atm.SetMapIdx(1)
            new_smi = add_tag_atom(mol, tag)
            if new_smi not in used:
                used.add(new_smi)
    return list(used)


if __name__ == "__main__":
    input = docopt(__doc__)
    input_file = input.get("--in")
    output_file = input.get("--out")
    smarts = input.get("--smarts")
    tag = input.get("--tag")

    ifs = oemolistream(input_file)
    ofs = open(output_file,"a")

    tagged_smiles = set()
    for in_mol in tqdm(ifs.GetOEGraphMols()):
        [tagged_smiles.add(x) for x in tag_molecule(in_mol, smarts, tag)]
    for smi in tagged_smiles:
        print(smi,file=ofs)
