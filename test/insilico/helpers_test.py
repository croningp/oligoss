from polymersoup.insilico.helpers.helpers import *
import unittest

class test_generate_dict_isobaric_sequences(unittest.TestCase):
    """
    Class for testing basic functions from helpers.helpers.py that are used
    to generate dictionary of isobaric sets of sequences from
    sequence lists

    Args:
        unittest ([type]): [description]
    """
    def test_isobaric(self):
        """
        tests generate_all_sequence function in helpers
        """
        seq_group = ['AAV', 'AVA', "VAA", "VVA", "AVV", "VAV"]
        isobaric_dict = generate_dict_isobaric_sequences(seq_group)

        self.assertIn("AAV", isobaric_dict.keys())
        self.assertNotIn("VAA", isobaric_dict.keys())

        self.assertIn("AVV", isobaric_dict.keys())
        self.assertNotIn("VVA", isobaric_dict.keys())

        self.assertEqual(isobaric_dict["AAV"], ["AAV", "AVA", "VAA"])
        self.assertEqual(isobaric_dict["AVV"], ["AVV", "VAV", "VVA"])

class test_generate_reading_frames_sequence(unittest.TestCase):
    """
    Testing the generate_reading_frames_sequence(..) function
    """
    def test_reading_frame_shift(self):

        test_sequence = "ABCDEFG"
        reading_frame_shifts = ["ABCDEFG", "GABCDEF", "FGABCDE", "EFGABCD", "DEFGABC", "CDEFGAB", "BCDEFGA"]
        self.assertEqual(generate_reading_frames_sequence(test_sequence), reading_frame_shifts)

if __name__ == '__main__':
    unittest.main()
