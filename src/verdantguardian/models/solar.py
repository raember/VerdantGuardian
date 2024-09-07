# import astropy.coordinates as coord
from astropy.coordinates import EarthLocation, AltAz, get_sun as sun2altaz, SkyCoord
from astropy.time import Time
import astropy.units as u

from verdantguardian.services.open_elevation import OpenElevationService


def get_loc(lat: float, lon: float) -> EarthLocation:
    height = OpenElevationService.get_elevation(lat, lon)
    return EarthLocation.from_geodetic(lon, lat, height)

def get_sun_pos(loc: EarthLocation, time: Time = Time.now()) -> SkyCoord:
    local_altaz = AltAz(location=loc, obstime=time)
    sun: SkyCoord = sun2altaz(time)
    return sun.transform_to(local_altaz)
