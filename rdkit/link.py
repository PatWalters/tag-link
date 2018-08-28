#!/usr/bin/env python

"""Usage: link.py --file1 FILE1 --link1 LINK1 --file2 FILE2 --link2 LINK2 --out OUTFILE

--file1 FILE1 first file to be linked
--link1 LINK1 atom type in file1
--file2 FILE2 second file to be linked
--link2 LINK2 atom type in file2
--out OUTFILE output file
"""

import sys
from rdkit import Chem
from docopt import docopt


def tag_attached_atom(mol, qmol, map_idx):
    match_list = mol.GetSubstructMatches(qmol)
    if len(match_list) == 0:
        return None
    for m in match_list[0]:
        atm = mol.GetAtomWithIdx(m)
        nbr_list = [x.GetOtherAtom(atm).GetIdx() for x in atm.GetBonds()]
        not_list = [x for  x in nbr_list if x not in match_list[0]]
        if len(not_list):

            attach_atm = mol.GetAtomWithIdx()
    sys.exit(0)



def map_attached_atom(mol, qmol, map_idx):
    """
    Assign an AtomMapNum to the atom attached to an atom adjacent to an atom with a specified atomic number
    Also deletes the atom with the specified atomic number
    E.g. with [Au]CCC as input generated [C:1]CC
    :param mol: input molecule
    :param atomic_num: atomic number of the atom to be removed
    :param map_idx: map index to assign to the the adjacent atom
    :return: modified molecule
    """
    match_list = mol.GetSubstructureMatches(qmol)
    if len(match_list) == 0:
        return None
    for atm in rw_mol.GetAtoms():
        if atm.GetAtomicNum() == atomic_num:
            to_remove = atm.GetIdx()
            nbr = [x.GetOtherAtom(atm) for x in atm.GetBonds()][0]
            nbr.SetAtomMapNum(map_idx)
    if to_remove >= 0:
        rw_mol.RemoveAtom(to_remove)
    return rw_mol


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
    print(m1,m2)
    if m1 >= 0 and m2 >= 0:
        new_mol.AddBond(m1, m2, order=Chem.rdchem.BondType.SINGLE)
    return new_mol


def get_atomic_num(smiles):
    mol = Chem.MolFromSmiles(smiles)
    atm = [x for x in mol.GetAtoms()][0]
    return atm.GetAtomicNum()


mol = Chem.MolFromSmiles("CCCC(=O)O")
query = Chem.MolFromSmarts("C(=O)O")
tag_attached_atom(mol,query,1)


#
# if __name__ == "__main__":
#     cmd_input = docopt(__doc__)
#     file_name1 = cmd_input.get("--file1")
#     file_name2 = cmd_input.get("--file2")
#     link1 = cmd_input.get("--link1")
#     link2 = cmd_input.get("--link2")
#     ofs = open(cmd_input.get("--out"), "w")
#
#     suppl_1 = Chem.SmilesMolSupplier(file_name1)
#     suppl_2 = Chem.SmilesMolSupplier(file_name2)
#
#     mol_list1 = [x for x in suppl_1]
#     mol_list2 = [x for x in suppl_2]
#
#     type_1 = get_atomic_num(link1)
#     type_2 = get_atomic_num(link2)
#
#     for m1 in mol_list1:
#         for m2 in tqdm(mol_list2):
#             print(link_molecules(m1, link1, m2, link2))
