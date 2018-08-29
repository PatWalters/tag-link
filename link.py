#!/usr/bin/env python

"""Usage: link.py --file1 FILE1 --link1 LINK1 --file2 FILE2 --link2 LINK2 --out OUTFILE

--file1 FILE1 first file to be linked
--link1 LINK1 atom type in file1
--file2 FILE2 second file to be linked
--link2 LINK2 atom type in file2
--out OUTFILE output file
"""

import sys
from openeye.oechem import *
from docopt import docopt
from tqdm import tqdm


def adjust_hydrogens(mol):
    for atom in mol.GetAtoms():
        h_count = OEDefaultMDLHCount(atom) - atom.GetExplicitValence()
        if h_count < 0:
            h_count = 0
        atom.SetImplicitHCount(h_count)
    OEAssignMDLHydrogens(mol)


def link_molecules(mol1, type1, mol2, type2):
    num1 = OEGetAtomicNum(type1)
    num2 = OEGetAtomicNum(type2)
    match_atom_list1 = [x for x in mol1.GetAtoms(OEHasAtomicNum(num1))]
    if len(match_atom_list1) == 0:
        print("Error: molecule1 does not contain atom type %s" % type1)
        sys.exit(0)
    atm1 = match_atom_list1[0]
    nbr1 = [x for x in atm1.GetAtoms()][0]

    match_atom_list2 = [x for x in mol2.GetAtoms(OEHasAtomicNum(num2))]
    if len(match_atom_list2) == 0:
        print("Error: molecule2 does not contain atom type %s" % type2)
        sys.exit(0)
    atm2 = match_atom_list2[0]
    nbr2 = [x for x in atm2.GetAtoms()][0]

    nbr1.SetMapIdx(1)
    nbr2.SetMapIdx(1)
    mol1.DeleteAtom(atm1)
    mol2.DeleteAtom(atm2)
    OEAddMols(mol1, mol2)

    link_atms = [x for x in mol1.GetAtoms() if x.GetMapIdx() == 1]
    l1, l2 = link_atms
    mol1.NewBond(l1, l2, 1)
    adjust_hydrogens(mol1)

    for atm in mol1.GetAtoms():
        atm.SetMapIdx(0)
    return mol1

def main(doc_string):
    input = docopt(__doc__)
    file_name1 = input.get("--file1")
    ifs1 = oemolistream(file_name1)
    file_name2 = input.get("--file2")
    ifs2 = oemolistream(file_name2)
    link1 = input.get("--link1")
    link2 = input.get("--link2")
    ofs = oemolostream(input.get("--out"))

    mol_list1 = [OEGraphMol(x) for x in ifs1.GetOEGraphMols()]
    mol_list2 = [OEGraphMol(x) for x in ifs2.GetOEGraphMols()]

    for m1 in mol_list1:
        for m2 in tqdm(mol_list2):
            n1 = OEGraphMol(m1)
            n2 = OEGraphMol(m2)
            new_mol = link_molecules(n1, link1, n2, link2)
            new_title = m2.GetTitle()
            if len(m1.GetTitle()):
                new_title = m1.GetTitle() + "_" + m2.GetTitle()
            new_mol.SetTitle(new_title)
            OEWriteMolecule(ofs, new_mol)


if __name__ == "__main__":
    main(__doc__)