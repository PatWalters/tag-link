#!/usr/bin/env python

"""Usage: link.py --file1 FILE1 --link1 LINK1 --file2 FILE2 --link2 LINK2 --out OUTFILE

--file1 FILE1 first file to be linked
--link1 LINK1 linking SMARTS in file1
--file2 FILE2 second file to be linked
--link2 LINK2 linking SMARTS in file2
--out OUTFILE output file
"""

import sys
import os
from rdkit import Chem
from docopt import docopt
from tqdm import tqdm


def map_attached_atom(in_mol, qmol, map_idx):
    """
    Read an input molecule, identify atoms mapped by a query molecule (SMARTS), identify the atom attaching
    that group to the remainder of the molecule, tag the attached atom with the AtomMapNum specified by
    map_idx, delete the atoms matched by the query
    :param in_mol: input mmolecule
    :param qmol: query molecule
    :param map_idx: atom map to be assigned
    :return: molecule with atom tagged and group matching qmol removed
    """
    mol = Chem.Mol(in_mol)
    match_list = mol.GetSubstructMatches(qmol)
    if len(match_list) == 0:
        return None
    for m in match_list[0]:
        atm = mol.GetAtomWithIdx(m)
        nbr_list = [x.GetOtherAtom(atm).GetIdx() for x in atm.GetBonds()]
        not_list = [x for x in nbr_list if x not in match_list[0]]
        if len(not_list):
            attach_idx = not_list[0]
            attach_atm = mol.GetAtomWithIdx(attach_idx)
            attach_atm.SetAtomMapNum(map_idx)
    return Chem.DeleteSubstructs(mol, qmol)


def link_molecules(in_mol_1, qmol_1, in_mol_2, qmol_2):
    """
    Link two molecules with a single bond
    E.g. [Au]CCC + [Au]c1ccccc1 -> CCC-c1ccccc1
    :param in_mol_1: the first molecule
    :param qmol_1: query to match in the first molecule
    :param in_mol_2: the second molecule
    :param qmol_2: query to match in the second molecule
    :return: modified molecule
    """
    mol_1 = map_attached_atom(in_mol_1, qmol_1, 1)
    mol_2 = map_attached_atom(in_mol_2, qmol_2, 2)
    new_mol = Chem.RWMol(Chem.CombineMols(mol_1, mol_2))
    m1 = -1
    m2 = -1
    for atm in new_mol.GetAtoms():
        if atm.GetAtomMapNum() == 1:
            m1 = atm.GetIdx()
            atm.SetAtomMapNum(0)
        if atm.GetAtomMapNum() == 2:
            m2 = atm.GetIdx()
            atm.SetAtomMapNum(0)
    if m1 >= 0 and m2 >= 0:
        new_mol.AddBond(m1, m2, order=Chem.rdchem.BondType.SINGLE)
    return new_mol


def get_supplier(file_name):
    """
    Get an RDKit supplier based on a file name
    :param file_name: input file name
    :return: supplier
    """
    ext = os.path.splitext(file_name)[-1]
    try:
        if ext == ".smi":
            return Chem.SmilesMolSupplier(file_name)
        elif ext == ".sdf":
            return Chem.SDMolSupplier(file_name)
        else:
            print(f"{ext} is not a support input file type")
    except FileNotFoundError:
        print(f"Could not open {file_name}",file=sys.stderr)
        sys.exit(1)


def get_writer(file_name):
    """
    Get an RDKit writer based on a file name
    :param file_name: output file name
    :return: writer
    """
    ext = os.path.splitext(file_name)[-1]
    try:
        if ext == ".smi":
            return Chem.SmilesWriter(file_name)
        elif ext == ".sdf":
            return Chem.SDWriter(file_name)
        else:
            print(f"{ext} is not a support input file type")
    except FileNotFoundError:
        print(f"Could not open {file_name}", file=sys.stderr)
        sys.exit(1)


def main(cmd_string):
    cmd_input = docopt(cmd_string)
    file_name1 = cmd_input.get("--file1")
    file_name2 = cmd_input.get("--file2")
    link1 = cmd_input.get("--link1")
    link2 = cmd_input.get("--link2")
    writer = get_writer(cmd_input.get("--out"))

    suppl_1 = get_supplier(file_name1)
    suppl_2 = get_supplier(file_name2)

    mol_list1 = [x for x in suppl_1]
    mol_list2 = [x for x in suppl_2]

    query1 = Chem.MolFromSmarts(link1)
    query2 = Chem.MolFromSmarts(link2)

    for m1 in mol_list1:
        name1 = m1.GetProp("_Name")
        for m2 in tqdm(mol_list2):
            name2 = m2.GetProp("_Name")
            linked_mol = link_molecules(m1, query1, m2, query2)
            Chem.SanitizeMol(linked_mol)
            linked_mol.SetProp("_Name",name1+"_"+name2)
            writer.write(linked_mol)


if __name__ == "__main__":
    main(__doc__)
