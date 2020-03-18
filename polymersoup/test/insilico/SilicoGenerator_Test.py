from polymersoup.insilico.SilicoGenerator import *
import unittest

class test_sequence_and_mass_generators(unittest.TestCase):
    """
    Class for testing basic functions from helpers.helpers.py that are used
    to generate lists of sequences from input monomers and calculate sequence
    neutral masses

    Args:
        unittest ([type]): [description]
    """
    def test_sequence_generator(self):
        """
        tests generate_all_sequence function in helpers
        """
        monomers = ['A', 'V']
        sequences = sorted(generate_all_sequences(monomers, 3, 3))
        compositions = generate_all_sequences(monomers, 3, 3, False)

        self.assertEqual(
            sequences,
            sorted(generate_all_sequences(['V', 'A'], 3, 3)
            ))
        self.assertNotEqual(compositions, sequences)

    def test_sequence_mass(self):
        """
        tests find_sequence_mass function in helpers
        """
        sequences = ['AGS', 'GSA', 'SGA']
        calc_masses = [find_sequence_mass(sequence) for sequence in sequences]
        real_masses = [233.1012]

        unique_calc_masses = list(set(calc_masses))

        self.assertEqual(real_masses, unique_calc_masses)

        A_monomer, A_dimer, A_trimer = 'A', 'AA', 'AAA'

        real_monomer_mass = round(MONOMERS["A"][0],4)
        calc_monomer_mass = find_sequence_mass(A_monomer)

        dimer_mass = find_sequence_mass(A_dimer)
        trimer_mass = find_sequence_mass(A_trimer)

        self.assertEqual(real_monomer_mass, calc_monomer_mass)
        self.assertNotEqual(real_monomer_mass, dimer_mass)
        self.assertNotEqual(real_monomer_mass, find_sequence_mass('S'))
        self.assertNotEqual(calc_monomer_mass, dimer_mass)
        self.assertNotEqual(calc_monomer_mass, find_sequence_mass('G'))

        real_massdiff = round(-MASS_DIFF, 4)
        calc_massdiff = round(trimer_mass - dimer_mass - real_monomer_mass,4)
        calc_massdiff - round(calc_massdiff, 4)

        self.assertEqual(real_massdiff, calc_massdiff)

class test_simple_adduct_addition(unittest.TestCase):
    """
    Class for testing addition of single adducts to sequence masses - does
    NOT test for functions that build adducts for multiple metal centre
    complexes, as those functions have not yet been written (as of 18.07.2019)

    Args:
        unittest ([type]): [description]
    """
    def test_add_adducts_sequence_mass(self):

        cations = ['H', 'Fe', 'Mn']
        anions = ['-H', 'Cl', 'ClO4']

        sequence = 'AAVG'
        neutral_mass = find_sequence_mass(sequence)
        min_z, max_z = 1, None

        pos_adduct_masses = add_adducts_sequence_mass(
            neutral_mass,
            cations,
            min_z,
            max_z,
            'pos')

        neg_adduct_masses = add_adducts_sequence_mass(
            neutral_mass,
            anions,
            min_z,
            max_z,
            'neg')

        self.assertIs(list, type(pos_adduct_masses))
        self.assertNotEqual(list(set(pos_adduct_masses)), pos_adduct_masses)


if __name__ == '__main__':
    unittest.main()
