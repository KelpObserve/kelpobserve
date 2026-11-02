"""Sun-glint geometry for Sentinel-2 tiles, before anything is downloaded.

Sun glint is the specular reflection of the sun off the sea surface. Its
strength falls off with the angle between the sensor's view and the direction
of perfect specular reflection: the glint angle. A glint angle below about 20°
puts the scene inside the glint pattern of a moderately rough sea, and no
atmospheric correction fully removes what that adds over water.

The angles are read from the tile metadata (``MTD_TL.xml``), which Planetary
Computer serves as the ``granule-metadata`` asset of each L2A item: a few
kilobytes, not the 800 MB bundle.

Azimuth convention
------------------
Both azimuths are clockwise from north. The sun azimuth points from the target
toward the sun. For the viewing azimuth the Sentinel-2 metadata is not explicit,
and the two readings put a scene near or far from the specular direction; see
:data:`VIEW_AZIMUTH_POINTS_TO_SENSOR` for how that was settled with measured
glint rather than assumed.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any

#: Whether the Sentinel-2 viewing azimuth points from the target to the sensor.
#: Settled on 23 Sesimbra scenes from 2026 whose L2A B11 glint was measured for
#: the Sesimbra benchmark: with this reading the glint angle ranks the measured
#: glint with Spearman -0.65, and the orbit measured at 0.11-0.34 in June and July
#: (R080) sits at 13.5-15° while the one measured at 0.01-0.04 (R037) sits at
#: 28-30°. The opposite reading gives -0.08 and swaps the two orbits.
VIEW_AZIMUTH_POINTS_TO_SENSOR = True


def glint_angle(
    sun_zenith_deg: float,
    sun_azimuth_deg: float,
    view_zenith_deg: float,
    view_azimuth_deg: float,
    *,
    view_azimuth_points_to_sensor: bool = VIEW_AZIMUTH_POINTS_TO_SENSOR,
) -> float:
    """Angle in degrees between the view direction and the specular direction.

    Zero means the sensor looks straight at the sun's mirror image. With the sun
    azimuth pointing to the sun and the view azimuth pointing to the sensor, the
    specular view lies at the sun's zenith and the opposite azimuth.
    """
    view_az = view_azimuth_deg if view_azimuth_points_to_sensor else view_azimuth_deg + 180.0
    ts, tv = math.radians(sun_zenith_deg), math.radians(view_zenith_deg)
    relative = math.radians(sun_azimuth_deg - view_az)
    cosine = math.cos(ts) * math.cos(tv) - math.sin(ts) * math.sin(tv) * math.cos(relative)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


@dataclass(frozen=True)
class TileAngles:
    """Tile-mean solar and viewing angles, in degrees."""

    sun_zenith_deg: float
    sun_azimuth_deg: float
    #: Band id (0 = B1 … 12 = B12, Sentinel-2 numbering) to (zenith, azimuth).
    view_by_band: dict[int, tuple[float, float]]

    def view(self, band_id: int | None = None) -> tuple[float, float]:
        """Viewing angles for one band, or the circular mean over all bands."""
        if band_id is not None:
            return self.view_by_band[band_id]
        zeniths = [z for z, _ in self.view_by_band.values()]
        sin_sum = sum(math.sin(math.radians(a)) for _, a in self.view_by_band.values())
        cos_sum = sum(math.cos(math.radians(a)) for _, a in self.view_by_band.values())
        return sum(zeniths) / len(zeniths), math.degrees(math.atan2(sin_sum, cos_sum)) % 360.0

    def glint_angle(self, band_id: int | None = None) -> float:
        zenith, azimuth = self.view(band_id)
        return glint_angle(self.sun_zenith_deg, self.sun_azimuth_deg, zenith, azimuth)


def parse_tile_angles(xml_text: str) -> TileAngles:
    """Mean sun and per-band viewing angles from a Sentinel-2 ``MTD_TL.xml``."""
    root = ET.fromstring(xml_text)

    def local(tag: str) -> str:
        return tag.rsplit("}", 1)[-1]

    sun = next((e for e in root.iter() if local(e.tag) == "Mean_Sun_Angle"), None)
    if sun is None:
        raise ValueError("No Mean_Sun_Angle in the tile metadata.")
    sun_values = {local(c.tag): float(c.text or "nan") for c in sun}
    views: dict[int, tuple[float, float]] = {}
    for element in root.iter():
        if local(element.tag) != "Mean_Viewing_Incidence_Angle":
            continue
        values = {local(c.tag): float(c.text or "nan") for c in element}
        views[int(element.attrib["bandId"])] = (values["ZENITH_ANGLE"], values["AZIMUTH_ANGLE"])
    if not views:
        raise ValueError("No Mean_Viewing_Incidence_Angle in the tile metadata.")
    return TileAngles(sun_values["ZENITH_ANGLE"], sun_values["AZIMUTH_ANGLE"], views)


def item_tile_angles(
    item: Any, asset: str = "granule-metadata", timeout_s: float = 60.0
) -> TileAngles:
    """Tile angles for a signed Planetary Computer Sentinel-2 L2A item."""
    import requests

    response = requests.get(item.assets[asset].href, timeout=timeout_s)
    response.raise_for_status()
    return parse_tile_angles(response.text)
