"""Which Sentinel-2 overpasses are worth processing at a site: image quality only.

Four checks, each read from the L2A product on Planetary Computer before any
L1C bundle is downloaded:

1. **Tile cloud**, from metadata. A cost control, not a judgement: a cloudy
   tile can be clear over a small site.
2. **Glint geometry**: the angle between the view and the sun's specular
   reflection, from the tile metadata (:mod:`.geometry`). At Sesimbra one orbit
   sits inside 15° every summer and every one of its scenes was glinted.
3. **Measured glint**: median 1610 nm (B11) reflectance over water. Water
   reflects almost nothing at 1610 nm, so what remains is glint, haze or
   adjacency. The Sesimbra benchmark registered 0.02 as its screen.
4. **Haze**: median over-water red (:func:`.planetary.measure_over_water_red`).

Every overpass keeps its measured values and all its failure reasons, whatever
the verdict, so a threshold can be revisited without re-reading the imagery.
Scenes are never chosen by how well a retrieval scores on them.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from datetime import date, datetime
from typing import Any

import numpy as np

from kelpobserve.coastal.acquire.geometry import item_tile_angles
from kelpobserve.coastal.acquire.planetary import (
    HAZY_MIN,
    BBox,
    _coverage_fraction,
    measure_over_water,
    mgrs_tile,
    open_catalog,
)

logger = logging.getLogger(__name__)

#: B11 is band id 11 in the Sentinel-2 tile metadata.
SWIR_BAND_ID = 11


@dataclass(frozen=True)
class UsabilityThresholds:
    #: Glint angles below this put the view inside the glint pattern. 20° is
    #: conservative against Sesimbra, where glinted scenes sat at 13.5–15°.
    glint_angle_min_deg: float = 20.0
    #: Median B11 over water above this is glint or haze. The Sesimbra
    #: benchmark's registered screen.
    swir_glint_max: float = 0.02
    #: Median over-water red at or above this is haze, suspended sediment or a
    #: bright shallow bed; see :mod:`.planetary` before reusing at turbid sites.
    red_hazy_min: float = HAZY_MIN
    #: Recorded as ``tile_cloudy`` without reading the imagery.
    max_tile_cloud_pct: float = 60.0
    #: NIR cut for water when no water geometry is supplied. Looser than the
    #: haze screen's, because glint raises NIR too and a tight cut would drop
    #: exactly the glinted pixels this is meant to measure.
    glint_water_nir_max: float = 0.25


@dataclass(frozen=True)
class SceneUsability:
    sensing_utc: datetime
    platform: str | None
    tile: str | None
    relative_orbit: int | None
    tile_cloud_pct: float | None
    glint_angle_deg: float | None
    swir_glint: float | None
    over_water_red: float | None
    n_water_pixels: int
    #: ``usable`` or the first failing check, in the order of the module docstring.
    verdict: str
    reasons: list[str] = field(default_factory=list)

    @property
    def date(self) -> date:
        return self.sensing_utc.date()

    @property
    def usable(self) -> bool:
        return self.verdict == "usable"

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["sensing_utc"] = self.sensing_utc.isoformat()
        return payload


def assess(
    item: Any,
    bbox: BBox,
    thresholds: UsabilityThresholds = UsabilityThresholds(),
    *,
    water_geometry: dict[str, Any] | None = None,
) -> SceneUsability:
    """All four checks on one L2A item."""
    props = item.properties
    tile_cloud = props.get("eo:cloud_cover")
    base = {
        "sensing_utc": item.datetime,
        "platform": props.get("platform"),
        "tile": mgrs_tile(item),
        "relative_orbit": props.get("sat:relative_orbit"),
        "tile_cloud_pct": None if tile_cloud is None else round(float(tile_cloud), 1),
    }
    if tile_cloud is not None and float(tile_cloud) > thresholds.max_tile_cloud_pct:
        return SceneUsability(
            **base,
            glint_angle_deg=None,
            swir_glint=None,
            over_water_red=None,
            n_water_pixels=0,
            verdict="tile_cloudy",
            reasons=["tile_cloudy"],
        )
    try:
        angle = item_tile_angles(item).glint_angle(SWIR_BAND_ID)
        swir, n_water = measure_over_water(
            item,
            bbox,
            "B11",
            water_nir_max=thresholds.glint_water_nir_max,
            water_geometry=water_geometry,
        )
        red, _ = measure_over_water(item, bbox, "B04", water_geometry=water_geometry)
    except Exception as exc:  # noqa: BLE001 - one unreadable scene is not fatal
        logger.warning("%s unreadable: %s", item.id, exc)
        return SceneUsability(
            **base,
            glint_angle_deg=None,
            swir_glint=None,
            over_water_red=None,
            n_water_pixels=0,
            verdict="unreadable",
            reasons=["unreadable"],
        )
    return classify(base, angle, swir, red, n_water, thresholds)


def classify(
    base: dict[str, Any],
    angle: float,
    swir: float,
    red: float,
    n_water: int,
    thresholds: UsabilityThresholds,
) -> SceneUsability:
    reasons = []
    if angle < thresholds.glint_angle_min_deg:
        reasons.append("glint_geometry")
    if not np.isfinite(swir):
        reasons.append("no_water")
    elif swir > thresholds.swir_glint_max:
        reasons.append("glint_measured")
    if np.isfinite(red) and red >= thresholds.red_hazy_min:
        reasons.append("hazy")
    return SceneUsability(
        **base,
        glint_angle_deg=round(angle, 2),
        swir_glint=round(swir, 4) if np.isfinite(swir) else None,
        over_water_red=round(red, 4) if np.isfinite(red) else None,
        n_water_pixels=n_water,
        verdict=reasons[0] if reasons else "usable",
        reasons=reasons,
    )


def overpasses(items: list[Any], bbox: BBox) -> list[Any]:
    """One item per overpass: the tile covering most of ``bbox``, then least cloud.

    Keyed by platform and sensing minute, not date: two orbits can pass the same
    site on the same day, as R037 and R080 do at Sesimbra.
    """
    groups: dict[tuple[str, str], list[Any]] = {}
    for item in items:
        if item.datetime is None:
            continue
        key = (str(item.properties.get("platform", "")), item.datetime.strftime("%Y-%m-%dT%H:%M"))
        groups.setdefault(key, []).append(item)
    chosen = [
        max(
            group,
            key=lambda it: (
                round(_coverage_fraction(it, bbox), 3),
                -float(it.properties.get("eo:cloud_cover", 100.0)),
            ),
        )
        for group in groups.values()
    ]
    return sorted(chosen, key=lambda it: it.datetime)


def screen_usability(
    bbox: BBox,
    start: date | str,
    end: date | str,
    thresholds: UsabilityThresholds = UsabilityThresholds(),
    *,
    water_geometry: dict[str, Any] | None = None,
    catalog: Any = None,
    collection: str = "sentinel-2-l2a",
) -> list[SceneUsability]:
    """Assess every Sentinel-2 overpass of ``bbox`` between two dates.

    ``bbox`` should be a tight window over the site's water, not the processing
    extent. ``water_geometry`` (GeoJSON, longitude and latitude), for instance
    from the site's lidar, makes the water mask independent of the image.
    """
    catalog = catalog or open_catalog()
    items = list(
        catalog.search(
            collections=[collection],
            bbox=list(bbox),
            datetime=f"{date.fromisoformat(str(start)).isoformat()}/"
            f"{date.fromisoformat(str(end)).isoformat()}",
        ).items()
    )
    logger.info("assessing %d STAC item(s) over %s", len(items), bbox)
    return [
        assess(item, bbox, thresholds, water_geometry=water_geometry)
        for item in overpasses(items, bbox)
    ]

