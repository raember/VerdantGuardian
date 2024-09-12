import unittest
from time import sleep

from common import ZURICH_COORDS, ZURICH_ALTITUDE
from verdantguardian.services.open_elevation import OpenElevationService

class OpenElevationServiceTest(unittest.TestCase):
    def test_simple(self):
        lat, lon = ZURICH_COORDS
        self.assertEqual(ZURICH_ALTITUDE, OpenElevationService.get_elevation(lat, lon))

if __name__ == '__main__':
    unittest.main()
