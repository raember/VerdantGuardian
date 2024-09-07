import unittest
from time import sleep

from verdantguardian.services.open_elevation import OpenElevationService

class OpenElevationServiceTest(unittest.TestCase):
    def test_simple(self):
        zurich_coords = (47.35953600, 8.63564520)
        self.assertEqual(582, OpenElevationService.get_elevation(zurich_coords[0], zurich_coords[1]))

if __name__ == '__main__':
    unittest.main()
