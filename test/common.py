from datetime import datetime

from astropy.coordinates import Latitude
from astropy.time import TimezoneInfo
from astropy import units as u

ZURICH_COORDS = (Latitude(47.36616*u.deg), Latitude(8.54150*u.deg))
ZURICH_ALTITUDE = 404.0
ZURICH_TZ = TimezoneInfo(utc_offset=2*u.hour)
ZURICH_TIMESTAMP = datetime(2024, 8, 1, 12, 0, 0, 0, tzinfo=ZURICH_TZ)
ZURICH_SUN_ALT_AZ = (54.96, 139.55)  # Altitude(deg), Azimuth(deg; N=0°, E=90°, S=180°, W=260°)
