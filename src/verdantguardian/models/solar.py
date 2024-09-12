# import astropy.coordinates as coord
from datetime import datetime, timezone, timedelta, tzinfo

import numpy as np
from ambiance import Atmosphere
from astroplan import Observer, FixedTarget
from astropy.coordinates import EarthLocation, AltAz, get_sun as sun2altaz, SkyCoord, Latitude, Longitude
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity, temperature
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter

from common import ZURICH_COORDS
from verdantguardian.services.open_elevation import OpenElevationService
from verdantguardian.utils.constants import TSI_MIN, SUN_LUX_CONVERSION_FACTOR


def get_loc(lat: Latitude, lon: Longitude, height_over_ground: Quantity = 0.0*u.m) -> EarthLocation:
    height = OpenElevationService.get_elevation(lat.value, lon.value) * u.m
    return EarthLocation.from_geodetic(lon, lat, height + height_over_ground)

def get_sun_pos(loc: EarthLocation, time: Time = Time.now()) -> SkyCoord:
    local_altaz = AltAz(location=loc, obstime=time)
    sun: SkyCoord = sun2altaz(time)
    return sun.transform_to(local_altaz)

MAR_EQUINOX = Time(datetime(2024, 3, 20, 3, 7, tzinfo=timezone(timedelta(hours=0), 'UTC')))
JUN_SOLSTICE = Time(datetime(2024, 6, 20, 20, 51, tzinfo=timezone(timedelta(hours=0), 'UTC')))
SEP_EQUINOX = Time(datetime(2024, 9, 22, 12, 44, tzinfo=timezone(timedelta(hours=0), 'UTC')))
DEC_SOLSTICE = Time(datetime(2024, 12, 21, 9, 20, tzinfo=timezone(timedelta(hours=0), 'UTC')))

ATMOSPHERE = Atmosphere(0)


def observe(loc: EarthLocation, dt: datetime) -> Observer:
    _ = timezone.utc if dt.tzinfo == timezone.utc else dt.tzinfo._utcoffset
    atmosphere = Atmosphere(loc.height)
    pressure = atmosphere.pressure[0]
    temp_c = atmosphere.temperature_in_celsius[0]
    observer = Observer(location=loc, timezone=dt.tzinfo, pressure=pressure * u.Pa, temperature=temp_c * u.deg_C)
    return observer
    time = Time(datetime(dt.year, dt.month, dt.day, 12, tzinfo=dt.tzinfo))
    sun = sun2altaz(time)

    civil_deg = 6*u.deg
    sun_radius = 0.5*u.deg
    dusk_time = observer.target_rise_time(time, sun, which='previous', horizon=-civil_deg)
    sunrise_time = observer.target_rise_time(time, sun, which='previous', horizon=-sun_radius)
    full_morning_time = observer.target_rise_time(time, sun, which='previous', horizon=sun_radius)
    full_evening_time = observer.target_set_time(time, sun, which='next', horizon=sun_radius)
    sunset_time = observer.target_set_time(time, sun, which='next', horizon=-sun_radius)
    dawn_time = observer.target_set_time(time, sun, which='next', horizon=-civil_deg)
    zenith_time = sunrise_time + (sunset_time - sunrise_time) / 2
    zenith_pos = observer.sun_altaz(zenith_time)

    sun_pos = observer.sun_altaz(time)

    print(f"Time:         {time.to_datetime(dt.tzinfo)}")
    print(f"Dusk:         {dusk_time.to_datetime(dt.tzinfo)}")
    print(f"Sunrise:      {sunrise_time.to_datetime(dt.tzinfo)}")
    print(f"Full morning: {full_morning_time.to_datetime(dt.tzinfo)}")
    print(f"Zenith:       {zenith_time.to_datetime(dt.tzinfo)}")
    print(f"Full evening: {full_evening_time.to_datetime(dt.tzinfo)}")
    print(f"Sunset:       {sunset_time.to_datetime(dt.tzinfo)}")
    print(f"Dawn:         {dawn_time.to_datetime(dt.tzinfo)}")
    print(f"Highest altitude: {zenith_pos.alt}")

    time_from = time - TimeDelta(timedelta(hours=12))
    time_to = time + TimeDelta(timedelta(hours=12))
    times = np.linspace(time_from, time_to, 100)
    # Calculate airmass
    poss = observer.altaz(times, sun)
    airmass = poss.secz
    # Mask out nonsense airmasses
    am = np.ma.array(airmass, mask=airmass < 1)
    dni = 1361 * np.exp(-0.14 * am) * u.W / u.m ** 2
    dni = np.power(0.7 * np.ones_like(am), am * 0.678) * TSI_MIN
    dhi = np.zeros_like(am)
    ghi = dhi + dni * np.cos(poss.alt)
    lux = ghi * SUN_LUX_CONVERSION_FACTOR * 1.0
    ax = plt.gca()
    ax.plot(times.plot_date, lux, linestyle='-', linewidth=1.5)
    date_formatter = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    xlo, xhi = (times[0]), (times[-1])
    ax.set_xlim([xlo.plot_date, xhi.plot_date])
    ax.set_yscale('log')
    plt.show()

    from astroplan.plots import plot_airmass
    ax = plot_airmass(sun, observer, time)
    plt.show()
    return

def get_times(observer: Observer, sun: SkyCoord, time: Time) -> dict[str, Time]:
    civil_deg = 6 * u.deg
    sun_radius = 0.5 * u.deg
    day = time.to_datetime().day
    dusk_time = observer.target_rise_time(time, sun, horizon=-civil_deg)
    sunrise_time = observer.target_rise_time(time, sun, horizon=-sun_radius)
    full_morning_time = observer.target_rise_time(time, sun, horizon=sun_radius)
    full_evening_time = observer.target_set_time(time, sun, horizon=sun_radius)
    sunset_time = observer.target_set_time(time, sun, horizon=-sun_radius)
    dawn_time = observer.target_set_time(time, sun, horizon=-civil_deg)
    zenith_time = sunrise_time + (sunset_time - sunrise_time) / 2
    zenith_pos = observer.sun_altaz(zenith_time)
    return {
        'dusk_time': dusk_time,
        'sunrise_time': sunrise_time,
        'full_morning_time': full_morning_time,
        'full_evening_time': full_evening_time,
        'sunset_time': sunset_time,
        'dawn_time': dawn_time,
        'zenith_time': zenith_time,
        'zenith_pos': zenith_pos,
    }

def calc(dt: datetime):
    # Solar Position Algorithm - Michalsky, Joseph J. 1988
    # https://www.nrel.gov/docs/fy08osti/34302.pdf
    date = Time(dt)
    # dt = date.to_datetime()
    # jd_2 = int(365.25*(dt.year + 4716)) + int(30.6001*(dt.month + 1)) + dt.day + (2 - int(dt.year/100) + int(int(dt.year/100)/4)) - 1524.5
    # jde_2 = jd_2 + (dt.hour*60*60 + dt.minute*60 + dt.second)/86400
    date.format = 'jd'

    # 3.1: Calculate the Julian and Julian Ephemeris Day, Century, and Millennium
    jde = date.value
    jd = np.uint32(jde)
    # jc = (jd - 2451545)/36525
    jce = (jde - 2451545)/36525
    jc = np.uint32(jce)
    jme = jce/10

    # 3.2: Calculate the Earth heliocentric longitude, latitude, and radius vector (L, B, and R):
    L0 = np.array([
        [175347046, 0, 0],
        [3341656, 4.6692568, 6283.07585],
        [34894, 4.6261, 12566.1517],
        [3497, 2.7441, 5753.3849],
        [3418, 2.8289, 3.5231],
        [3136, 3.6277, 77713.7715],
        [2676, 4.4181, 7860.4194],
        [2343, 6.1352, 3930.2097],
        [1324, 0.7425, 11506.7698],
        [1273, 2.0371, 529.691],
        [1199, 1.1096, 1577.3435],
        [990, 5.233, 5884.927],
        [902, 2.045, 26.298],
        [857, 3.508, 398.149],
        [780, 1.179, 5223.694],
        [753, 2.533, 5507.553],
        [505, 4.583, 18849.228],
        [492, 4.205, 775.523],
        [357, 2.92, 0.067],
        [317, 5.849, 11790.629],
        [284, 1.899, 796.298],
        [271, 0.315, 10977.079],
        [243, 0.345, 5486.778],
        [206, 4.806, 2544.314],
        [205, 1.869, 5573.143],
        [202, 2.458, 6069.777],
        [156, 0.833, 213.299],
        [132, 3.411, 2942.463],
        [126, 1.083, 20.775],
        [115, 0.645, 0.98],
        [103, 0.636, 4694.003],
        [102, 0.976, 15720.839],
        [102, 4.267, 7.114],
        [99, 6.21, 2146.17],
        [98, 0.68, 155.42],
        [86, 5.98, 161000.69],
        [85, 1.3, 6275.96],
        [85, 3.67, 71430.7],
        [80, 1.81, 17260.15],
        [79, 3.04, 12036.46],
        [75, 1.76, 5088.63],
        [74, 3.5, 3154.69],
        [74, 4.68, 801.82],
        [70, 0.83, 9437.76],
        [62, 3.98, 8827.39],
        [61, 1.82, 7084.9],
        [57, 2.78, 6286.6],
        [56, 4.39, 14143.5],
        [56, 3.47, 6279.55],
        [52, 0.19, 12139.55],
        [52, 1.33, 1748.02],
        [51, 0.28, 5856.48],
        [49, 0.49, 1194.45],
        [41, 5.37, 8429.24],
        [41, 2.4, 19651.05],
        [39, 6.17, 10447.39],
        [37, 6.04, 10213.29],
        [37, 2.57, 1059.38],
        [36, 1.71, 2352.87],
        [36, 1.78, 6812.77],
        [33, 0.59, 17789.85],
        [30, 0.44, 83996.85],
        [30, 2.74, 1349.87],
        [25, 3.16, 4690.48],
    ], dtype=np.float32)
    L1 = np.array([
        [628331966747, 0, 0],
        [206059, 2.678235, 6283.07585],
        [4303, 2.6351, 12566.1517],
        [425, 1.59, 3.523],
        [119, 5.796, 26.298],
        [109, 2.966, 1577.344],
        [93, 2.59, 18849.23],
        [72, 1.14, 529.69],
        [68, 1.87, 398.15],
        [67, 4.41, 5507.55],
        [59, 2.89, 5223.69],
        [56, 2.17, 155.42],
        [45, 0.4, 796.3],
        [36, 0.47, 775.52],
        [29, 2.65, 7.11],
        [21, 5.34, 0.98],
        [19, 1.85, 5486.78],
        [19, 4.97, 213.3],
        [17, 2.99, 6275.96],
        [16, 0.03, 2544.31],
        [16, 1.43, 2146.17],
        [15, 1.21, 10977.08],
        [12, 2.83, 1748.02],
        [12, 3.26, 5088.63],
        [12, 5.27, 1194.45],
        [12, 2.08, 4694],
        [11, 0.77, 553.57],
        [10, 1.3, 6286.6],
        [10, 4.24, 1349.87],
        [9, 2.7, 242.73],
        [9, 5.64, 951.72],
        [8, 5.3, 2352.87],
        [6, 2.65, 9437.76],
        [6, 4.67, 4690.48],
    ], dtype=np.float32)
    L2 = np.array([
        [52919, 0, 0],
        [8720, 1.0721, 6283.0758],
        [309, 0.867, 12566.152],
        [27, 0.05, 3.52],
        [16, 5.19, 26.3],
        [16, 3.68, 155.42],
        [10, 0.76, 18849.23],
        [9, 2.06, 77713.77],
        [7, 0.83, 775.52],
        [5, 4.66, 1577.34],
        [4, 1.03, 7.11],
        [4, 3.44, 5573.14],
        [3, 5.14, 796.3],
        [3, 6.05, 5507.55],
        [3, 1.19, 242.73],
        [3, 6.12, 529.69],
        [3, 0.31, 398.15],
        [3, 2.28, 553.57],
        [2, 4.38, 5223.69],
        [2, 3.75, 0.98],
    ], dtype=np.float32)
    L3 = np.array([
        [289, 5.844, 6283.076],
        [35, 0, 0],
        [17, 5.49, 12566.15],
        [3, 5.2, 155.42],
        [1, 4.72, 3.52],
        [1, 5.3, 18849.23],
        [1, 5.97, 242.73],
    ], dtype=np.float32)
    L4 = np.array([
        [114, 3.142, 0],
        [8, 4.13, 6283.08],
        [1, 3.84, 12566.15],
    ], dtype=np.float32)
    L5 = np.array([
        [1, 3.14, 0],
    ], dtype=np.float32)
    B0 = np.array([
        [280, 3.199, 84334.662],
        [102, 5.422, 5507.553],
        [80, 3.88, 5223.69],
        [44, 3.7, 2352.87],
        [32, 4, 1577.34],
    ], dtype=np.float32)
    B1 = np.array([
        [9, 3.9, 5507.55],
        [6, 1.73, 5223.69],
    ], dtype=np.float32)
    R0 = np.array([
        [100013989, 0, 0],
        [1670700, 3.0984635, 6283.07585],
        [13956, 3.05525, 12566.1517],
        [3084, 5.1985, 77713.7715],
        [1628, 1.1739, 5753.3849],
        [1576, 2.8469, 7860.4194],
        [925, 5.453, 11506.77],
        [542, 4.564, 3930.21],
        [472, 3.661, 5884.927],
        [346, 0.964, 5507.553],
        [329, 5.9, 5223.694],
        [307, 0.299, 5573.143],
        [243, 4.273, 11790.629],
        [212, 5.847, 1577.344],
        [186, 5.022, 10977.079],
        [175, 3.012, 18849.228],
        [110, 5.055, 5486.778],
        [98, 0.89, 6069.78],
        [86, 5.69, 15720.84],
        [86, 1.27, 161000.69],
        [65, 0.27, 17260.15],
        [63, 0.92, 529.69],
        [57, 2.01, 83996.85],
        [56, 5.24, 71430.7],
        [49, 3.25, 2544.31],
        [47, 2.58, 775.52],
        [45, 5.54, 9437.76],
        [43, 6.01, 6275.96],
        [39, 5.36, 4694],
        [38, 2.39, 8827.39],
        [37, 0.83, 19651.05],
        [37, 4.9, 12139.55],
        [36, 1.67, 12036.46],
        [35, 1.84, 2942.46],
        [33, 0.24, 7084.9],
        [32, 0.18, 5088.63],
        [32, 1.78, 398.15],
        [28, 1.21, 6286.6],
        [28, 1.9, 6279.55],
        [26, 4.59, 10447.39],
    ], dtype=np.float32)
    R1 = np.array([
        [103019, 1.10749, 6283.07585],
        [1721, 1.0644, 12566.1517],
        [702, 3.142, 0],
        [32, 1.02, 18849.23],
        [31, 2.84, 5507.55],
        [25, 1.32, 5223.69],
        [18, 1.42, 1577.34],
        [10, 5.91, 10977.08],
        [9, 1.42, 6275.96],
        [9, 0.27, 5486.78],
    ], dtype=np.float32)
    R2 = np.array([
        [4359, 5.7846, 6283.0758],
        [124, 5.579, 12566.152],
        [12, 3.14, 0],
        [9, 3.63, 77713.77],
        [6, 1.87, 5573.14],
        [3, 5.47, 18849.23],
    ], dtype=np.float32)
    R3 = np.array([
        [145, 4.273, 6283.076],
        [7, 3.92, 12566.15],
    ], dtype=np.float32)
    R4 = np.array([
        [4, 2.56, 6283.08],
    ], dtype=np.float32)
    def calc_l(arr: np.ndarray) -> np.ndarray:
        vals = arr[:, 0] * np.cos(arr[:, 1] + (arr[:, 2] * jme))
        return vals.sum()
    def calc_arr(arrs: list[np.ndarray]):
        val = 0
        for i, a in enumerate(arrs):
            val += calc_l(a) * jme ** i
        return val / 1e8
    def limit_deg(deg: np.float32) -> np.float32:
        f = (360 / deg) % 1
        return f * 360 if deg >= 0.0 else 360 - 360 * f
    e_lon = limit_deg(calc_arr([L0, L1, L2, L3, L4, L5]) * 180 / np.pi)
    e_lat = limit_deg(calc_arr([B0, B1]) * 180 / np.pi)
    e_r = calc_arr([R0, R1, R2, R3, R4])

    # 3.3: Calculate the geocentric longitude and latitude (theta and beta)
    theta = limit_deg(e_lon + 180)
    beta = -e_lat

    # 3.4: Calculate the nutation in longitude and obliquity (del_psi and del_eps)
    x0 = 297.85036 + 445267.111480 * jce - 0.0019142 * jce**2 + jce**3 / 189474
    x1 = 357.52772 + 35999.050340 * jce - 0.0001603 * jce**2 - jce**3 / 300000
    x2 = 134.96298 + 477198.867398 * jce + 0.0086972 * jce**2 + jce**3 / 56250
    x3 = 93.27191 + 483202.017538 * jce - 0.0036825 * jce**2 + jce**3 / 327270
    x4 = 125.04452 - 1934.136261 * jce + 0.0020708 * jce**2 + jce**3 / 450000
    sin_coeffs = np.array([
        [0, 0, 0, 0, 1],
        [-2, 0, 0, 2, 2],
        [0, 0, 0, 2, 2],
        [0, 0, 0, 0, 2],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [-2, 1, 0, 2, 2],
        [0, 0, 0, 2, 1],
        [0, 0, 1, 2, 2],
        [-2, -1, 0, 2, 2],
        [-2, 0, 1, 0, 0],
        [-2, 0, 0, 2, 1],
        [0, 0, -1, 2, 2],
        [2, 0, 0, 0, 0],
        [0, 0, 1, 0, 1],
        [2, 0, -1, 2, 2],
        [0, 0, -1, 0, 1],
        [0, 0, 1, 2, 1],
        [-2, 0, 2, 0, 0],
        [0, 0, -2, 2, 1],
        [2, 0, 0, 2, 2],
        [0, 0, 2, 2, 2],
        [0, 0, 2, 0, 0],
        [-2, 0, 1, 2, 2],
        [0, 0, 0, 2, 0],
        [-2, 0, 0, 2, 0],
        [0, 0, -1, 2, 1],
        [0, 2, 0, 0, 0],
        [2, 0, -1, 0, 1],
        [-2, 2, 0, 2, 2],
        [0, 1, 0, 0, 1],
        [-2, 0, 1, 0, 1],
        [0, -1, 0, 0, 1],
        [0, 0, 2, -2, 0],
        [2, 0, -1, 2, 1],
        [2, 0, 1, 2, 2],
        [0, 1, 0, 2, 2],
        [-2, 1, 1, 0, 0],
        [0, -1, 0, 2, 2],
        [2, 0, 0, 2, 1],
        [2, 0, 1, 0, 0],
        [-2, 0, 2, 2, 2],
        [-2, 0, 1, 2, 1],
        [2, 0, -2, 0, 1],
        [2, 0, 0, 0, 1],
        [0, -1, 1, 0, 0],
        [-2, -1, 0, 2, 1],
        [-2, 0, 0, 0, 1],
        [0, 0, 2, 2, 1],
        [-2, 0, 2, 0, 1],
        [-2, 1, 0, 2, 1],
        [0, 0, 1, -2, 0],
        [-1, 0, 1, 0, 0],
        [-2, 1, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 1, 2, 0],
        [0, 0, -2, 2, 2],
        [-1, -1, 1, 0, 0],
        [0, 1, 1, 0, 0],
        [0, -1, 1, 2, 2],
        [2, -1, -1, 2, 2],
        [0, 0, 3, 2, 2],
        [2, -1, 0, 2, 2],
    ],dtype=np.float32)
    del_psi_coeffs = np.array([
        [-171996, -174.2],
        [-13187, -1.6],
        [-2274, -0.2],
        [2062, 0.2],
        [1426, -3.4],
        [712, 0.1],
        [-517, 1.2],
        [-386, -0.4],
        [-301, 0],
        [217, -0.5],
        [-158, 0],
        [129, 0.1],
        [123, 0],
        [63, 0],
        [63, 0.1],
        [-59, 0],
        [-58, -0.1],
        [-51, 0],
        [48, 0],
        [46, 0],
        [-38, 0],
        [-31, 0],
        [29, 0],
        [29, 0],
        [26, 0],
        [-22, 0],
        [21, 0],
        [17, -0.1],
        [16, 0],
        [-16, 0.1],
        [-15, 0],
        [-13, 0],
        [-12, 0],
        [11, 0],
        [-10, 0],
        [-8, 0],
        [7, 0],
        [-7, 0],
        [-7, 0],
        [-7, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [-6, 0],
        [-6, 0],
        [5, 0],
        [-5, 0],
        [-5, 0],
        [-5, 0],
        [4, 0],
        [4, 0],
        [4, 0],
        [-4, 0],
        [-4, 0],
        [-4, 0],
        [3, 0],
        [-3, 0],
        [-3, 0],
        [-3, 0],
        [-3, 0],
        [-3, 0],
        [-3, 0],
        [-3, 0],
    ], dtype=np.float32)
    del_eps_coeffs = np.array([
        [92025, 8.9],
        [5736, -3.1],
        [977, -0.5],
        [-895, 0.5],
        [54, -0.1],
        [-7, 0],
        [224, -0.6],
        [200, 0],
        [129, -0.1],
        [-95, 0.3],
        [0, 0],
        [-70, 0],
        [-53, 0],
        [0, 0],
        [-33, 0],
        [26, 0],
        [32, 0],
        [27, 0],
        [0, 0],
        [-24, 0],
        [16, 0],
        [13, 0],
        [0, 0],
        [-12, 0],
        [0, 0],
        [0, 0],
        [-10, 0],
        [0, 0],
        [-8, 0],
        [7, 0],
        [9, 0],
        [7, 0],
        [6, 0],
        [0, 0],
        [5, 0],
        [3, 0],
        [-3, 0],
        [0, 0],
        [3, 0],
        [3, 0],
        [0, 0],
        [-3, 0],
        [-3, 0],
        [3, 0],
        [3, 0],
        [0, 0],
        [3, 0],
        [3, 0],
        [3, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    ], dtype=np.float32)
    print(jd)
