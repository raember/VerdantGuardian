import json
from datetime import datetime, timezone, timedelta

import pytz
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder

from common import ZURICH_COORDS
from verdantguardian.models.solar import calc, get_sun_pos, get_loc, observe

if __name__ == '__main__':
    tzf = TimezoneFinder()
    lat, lng = ZURICH_COORDS
    # lng = 150.0
    # Erkert, Tübingen:
    lat, lng = 48.32, 9.035
    tz = pytz.timezone(tzf.timezone_at(lat=lat, lng=lng))
    now = datetime.now(tz)
    now = datetime(1969, 4, 2, 12, tzinfo=tz)
    loc = get_loc(lat*u.deg, lng*u.deg, height_over_ground=6*u.m)
    observe(loc, now)
    calc(datetime.now())
