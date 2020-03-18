from polymersoup.insilico.MS1_silico import *
import unittest

class test_generate_mass_dictionary(unittest.TestCase):
	"""
	Class for testing the generate_mass_dictionary(..) function
	"""
	def test_glycine_3(self):
		glycine_3_dict = generate_mass_dictionary(["G"],
													3)
		actual_masses = {"G": 75.0320, "GG":132.0535, "GGG": 189.0750}
		self.assertEqual(glycine_3_dict, actual_masses)

	def test_glycine_tyrosine_3(self):
		glycine_tyrosine_3 = generate_mass_dictionary(["G", "Y"],
															3)
		actual_masses = {"G": 75.0320, "GG":132.0535, "GGG": 189.0750,
							"Y":181.0739, "YY":344.1372, "YYY":507.2005,
							"GY":238.0954, "YG":238.0954,
							"GGY":295.1168, "GYG":295.1168, "YGG":295.1168,
							"YYG":401.1587, "YGY":401.1587, "GYY":401.1587}
		self.assertEqual(actual_masses.keys(), glycine_tyrosine_3.keys())
		for k in glycine_tyrosine_3.keys():
			self.assertAlmostEqual(glycine_tyrosine_3[k], actual_masses[k], places = 3)

if __name__ == '__main__':
    unittest.main()
