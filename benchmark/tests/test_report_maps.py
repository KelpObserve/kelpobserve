"""The report's map helpers: geometry written for the story map must be well formed."""

import numpy as np
import pytest

from benchmark import report_maps as maps

rasterio = pytest.importorskip("rasterio")


def test_rounded_keeps_nesting():
    assert maps.rounded([[1.23456789, 2.0], [3.000000049, 4]], 6) == [[1.234568, 2.0], [3.0, 4.0]]


def test_polygon_utm_is_a_closed_lonlat_ring_near_the_site():
    ring = maps.polygon_utm(482000, 4251000, 482010, 4251010)["coordinates"][0]
    assert len(ring) == 5 and ring[0] == ring[-1]
    longitudes, latitudes = zip(*ring)
    assert -9.21 < min(longitudes) < max(longitudes) < -9.20
    assert 38.40 < min(latitudes) < max(latitudes) < 38.42


def test_overlay_corners_run_clockwise_from_the_top_left():
    # MapLibre's image source expects top-left, top-right, bottom-right, bottom-left.
    top_left, top_right, bottom_right, bottom_left = maps.overlay_corners()
    assert top_left[0] < top_right[0] and bottom_left[0] < bottom_right[0]
    assert top_left[1] > bottom_left[1] and top_right[1] > bottom_right[1]


def depth_raster(path):
    """Depth rising from 0 m in the north to 30 m in the south, on a lon/lat grid."""
    from rasterio.transform import from_origin

    depth = np.repeat(np.linspace(0, 30, 60)[:, None], 80, axis=1).astype("float32")
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=60,
        width=80,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(-9.3, 38.5, 0.001, 0.001),
    ) as target:
        target.write(depth, 1)


def test_contours_follow_the_depth_and_respect_the_clip(tmp_path):
    path = tmp_path / "depth.tif"
    depth_raster(path)
    lines = maps.contour_lines(path, (10.0, 20.0), "test")
    assert sorted(line["properties"]["depth_m"] for line in lines) == [10.0, 20.0]
    ten = next(line for line in lines if line["properties"]["depth_m"] == 10.0)
    latitudes = [point[1] for point in ten["geometry"]["coordinates"]]
    assert max(latitudes) - min(latitudes) < 1e-4  # the 10 m contour runs east-west
    clipped = maps.contour_lines(path, (10.0,), "test", clip=(-9.28, 38.4, -9.26, 38.5))
    longitudes = [point[0] for line in clipped for point in line["geometry"]["coordinates"]]
    assert -9.28 - 1e-6 <= min(longitudes) and max(longitudes) <= -9.26 + 1e-6


def test_report_links_preserve_published_evidence_and_external_sources(tmp_path):
    source = tmp_path / "report.md"
    body = ('See <a\nhref="results.md#july">the results</a>, '
            '<a href="missing.md">unpublished notes</a>, '
            '<a href="https://x.org/method.md">the method</a> and <a href="#limits">limits</a>.')
    published = {(tmp_path / "results.md").resolve(): "supporting/results.html"}
    cleaned = maps.document_links(body, source, published)
    assert cleaned == ('See <a href="supporting/results.html#july">the results</a>, '
                       'unpublished notes, <a href="https://x.org/method.md">the method</a> '
                       'and <a href="#limits">limits</a>.')
