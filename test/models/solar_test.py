import unittest
from datetime import datetime, timezone, timedelta

from astropy import units as u
from astropy.coordinates import Latitude, Longitude, get_sun

from common import ZURICH_COORDS, ZURICH_TIMESTAMP, ZURICH_SUN_ALT_AZ
from verdantguardian.models.solar import get_loc, get_sun_pos, calc, observe, get_times
from astropy.time import Time


class SolarTestCase(unittest.TestCase):
    def test_get_loc(self):
        lat, lon = ZURICH_COORDS
        loc = get_loc(lat, lon)
        self.assertAlmostEqual(lat, float(loc.lat / u.deg), places=6)
        self.assertAlmostEqual(lon, float(loc.lon / u.deg), places=6)

    def test_get_sun_pos(self):
        # https://www.suncalc.org/#/47.3662,8.5415,18/2024.08.01/12:00/404/3
        lat, lon = ZURICH_COORDS
        loc = get_loc(lat, lon)
        timestamp = Time(ZURICH_TIMESTAMP)
        sun_pos = get_sun_pos(loc, timestamp)
        alt, az = ZURICH_SUN_ALT_AZ
        self.assertAlmostEqual(alt, float(sun_pos.alt / u.deg), places=1)
        self.assertAlmostEqual(az, float(sun_pos.az / u.deg), places=2)

    def test_sun_pos_2(self):
        date = datetime(2003, 10, 17, 12, 30, 30, tzinfo=timezone(timedelta(hours=-7)))
        loc = get_loc(39.742476, -105.1786)
        timestamp = Time(date)
        sun_pos = get_sun_pos(loc, timestamp)

        c = sun_pos.icrs
        ra_rad = c.ra.wrap_at(180 * u.deg).radian
        dec_rad = c.dec.radian
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 4.2))
        plt.subplot(111, projection="aitoff")
        plt.title("Aitoff projection of sun position")
        plt.grid(True)
        plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
        plt.subplots_adjust(top=0.95, bottom=0.0)
        plt.show()
        self.assertAlmostEqual(0.0, float(sun_pos.alt / u.deg), places=1)
        self.assertAlmostEqual(0.0, float(sun_pos.az / u.deg), places=2)

    def test_brightest_sun(self):
        # https://en.wikipedia.org/wiki/Lux#Illuminance
        DIRECT_SUNLIGHT_LUX = (32000, 100000) * u.lx
        loc = get_loc(Latitude(0.0*u.deg), Longitude(0.0*u.deg))
        date = datetime(2024, 9, 1, 22, tzinfo=timezone.utc)
        time = Time(date)
        observer = observe(loc, date)
        sun = get_sun(time)
        times = get_times(observer, sun, time)

    def test_calc(self):
        date = datetime(2003, 10, 17, 12, 30, 30, tzinfo=timezone(timedelta(hours=-7)))
        calc(date)


if __name__ == '__main__':
    unittest.main()
