import parmed
import unittest
import sys

class PDBPrepTestCase(unittest.TestCase):

    def test_res_equal(self, pdb_file_one, pdb_file_two):
        # Get the structures
        structure_one = parmed.load_file(pdb_file_one)
        structure_two = parmed.load_file(pdb_file_two)
        # Do they have the same length
        print(pdb_file_one + ": "+str(len(structure_one.residues)))
        print(pdb_file_two + ": "+str(len(structure_two.residues)))
        # Do they have the same member
        print(list(structure_one.residues))
        print(list(structure_two.residues))

if __name__ == "__main__":
    pdb_prep = PDBPrepTestCase()
    pdb_prep.test_res_equal(pdb_file_one=sys.argv[1], pdb_file_two=sys.argv[2])