from verdantguardian.services.session import SESS


class OpenElevationService:
    @staticmethod
    def get_elevation(lat: float, lon: float) -> float:
        assert isinstance(lat, float) & isinstance(lon, float)
        resp = SESS.get(f"https://api.open-elevation.com/api/v1/lookup?locations={lat},{lon}")
        resp.raise_for_status()
        return resp.json()['results'][0]['elevation']
