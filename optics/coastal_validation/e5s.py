"""E5s: kelp spectral contrast and detectability from published libraries.

Registration: ``registration/e5s-kelp-spectra.json``. Characterisation, with decision
rules fixed before any contrast was computed. Stages: ``freeze``, ``run``.

Rowan & Phinn (UQ eSpace, doi 10.48610/c5f745e) may be re-used with acknowledgement
but not redistributed: it is read from the git-ignored cache, and only derived
statistics leave this module.
"""

from __future__ import annotations

import csv
import datetime as dt
import glob
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from coastal_validation import common

ROOT = common.ROOT
REGISTRATION = ROOT / "registration" / "e5s-kelp-spectra.json"
LOCK = ROOT / "registration" / "e5s-kelp-spectra-lock.json"
OUTPUT = ROOT / "reports" / "e5s-kelp-spectra"
SPECTRA = ROOT / ".cache" / "spectra"
ACOLITE_RSR = Path.home() / "acolite" / "data" / "RSR"
E2C = ROOT / "reports" / "e2c-sesimbra"
E2C_DATES = ("2025-12-24", "2026-02-22")
S2_BANDS = ("B1", "B2", "B3", "B4", "B5")
S2_10M = ("B2", "B3", "B4")
NEO_BANDS = ("B1", "B2", "B3", "B4", "B5")
IOP_NM = (442.0, 492.0, 559.0, 665.0, 704.0)
MIN_SPECTRA = 5
N_DRAWS, SEED = 5000, 42
SOLAR_ZENITH_DEG = 35.0
DEPTHS_M = (1.0, 3.0, 5.0, 8.0, 12.0, 16.0)
SIGMAS = (0.0015, 0.003, 0.006)
PRIMARY_SIGMA = 0.003
T_AW = 0.52
AUC_SEPARABLE = 0.8
FLAT_TEST_DEPTH_M = 5.0
UNDERSTORY = ("red_foliose", "turf")
BARE = ("rock", "sand")

SOURCE_FILES = {
    "rowan_phinn": SPECTRA / "rowan_phinn" / "spectra.csv",
    "vahtmae": SPECTRA / "pangaea_sav_spectra.tab",
    "fuller_table": SPECTRA / "hawaii" / "RS_Algae_SpeciesData.csv",
    "glovers": SPECTRA / "glovers_824861.tab",
}


def load_registration() -> dict[str, Any]:
    return json.loads(REGISTRATION.read_text())


# -- Libraries and classes ---------------------------------------------------------


def _norm(text: str) -> str:
    return " ".join(str(text).split()).casefold()


def _tag(library: Any, source: str, taxon: Any, geometry: Any) -> Any:
    """Add source, taxon and geometry to every spectrum's metadata."""
    from kelpobserve.coastal.spectra import SpectralLibrary

    metadata = tuple(
        {**m, "source": source, "taxon": taxon(m), "geometry": geometry(m)}
        for m in library.metadata
    )
    return SpectralLibrary(
        library.wavelengths_nm, library.values, metadata, library.quantity, source
    )


def _fuller_files() -> list[Path]:
    return sorted((SPECTRA / "hawaii" / "spectra" / "1_Spectra_RS_Algae_Paper").glob("*.asd"))


def _vigo_files() -> list[Path]:
    return sorted(Path(p) for p in glob.glob(str(SPECTRA / "vigo" / "CSV_Files" / "*.csv")))


def _wallis_files() -> list[tuple[Path, Path]]:
    pairs = []
    for raw in sorted(SPECTRA.glob("wallis/*/*/*_rawdata.txt")):
        if "representative" in raw.name:  # a subset of WallisLake_200805, not new spectra
            continue
        meta = raw.with_name(raw.name.replace("_rawdata.txt", "_metadata.txt"))
        if meta.exists():
            pairs.append((raw, meta))
    return pairs


def load_libraries() -> dict[str, Any]:
    """Every registered library as reflectance, with source, taxon and geometry."""
    from kelpobserve.coastal import spectra

    out: dict[str, Any] = {}
    rp = spectra.read_rows_csv(SOURCE_FILES["rowan_phinn"])
    out["rowan_phinn"] = _tag(rp, "rowan_phinn", lambda m: m["target"], lambda m: "air")

    va = spectra.read_pangaea_wide(SOURCE_FILES["vahtmae"]).as_reflectance()

    def va_taxon(m: dict[str, str]) -> str:
        latin = m.get("Scientific name (Target in Latin)", "").strip()
        if latin and not latin.lower().startswith("mostly"):
            return latin.split()[0]
        return " ".join(m.get("Description (Target description)", "").split())

    out["vahtmae"] = _tag(va, "vahtmae", va_taxon, lambda m: "air")

    rows = list(csv.DictReader(SOURCE_FILES["fuller_table"].open(encoding="cp1252")))
    by_id = {r["Spectral ID"]: r for r in rows}
    # Rows without a division repeat the sample above them (three scans of one sample).
    current: dict[str, str] = {}
    for r in rows:
        if r["Division"].strip():
            current = r
        by_id[r["Spectral ID"]] = {**current, "Spectral ID": r["Spectral ID"]}
    files = [f for f in _fuller_files() if f.stem in by_id]
    fu = spectra.read_asd_library(
        files, [{k: v or "" for k, v in by_id[f.stem].items() if k} for f in files]
    )
    out["fuller"] = _tag(
        fu,
        "fuller",
        lambda m: m["Genus"].strip(),
        lambda m: "air" if m["Method"].strip() == "ex-situ" else "in_situ_unstated",
    )

    gl = spectra.read_pangaea_long(
        SOURCE_FILES["glovers"],
        wavelength_column="Lambda [nm]",
        value_column="Refl tot",
        spectrum_columns=("Event", "Image no/name"),
    )
    out["glovers"] = _tag(gl, "glovers", lambda m: m["Comment"].strip(), lambda m: "underwater")

    vigo = [spectra.read_columns_csv(f) for f in _vigo_files()]
    out["vigo"] = [
        _tag(
            lib,
            "vigo",
            lambda m: m["file"]
            .replace("ReflectanceSpectralSignatureLibrary_", "")
            .replace(".csv", ""),
            lambda m: "underwater",
        )
        for lib in vigo
    ]

    wallis = []
    for raw, meta in _wallis_files():
        lib = spectra.read_columns_csv(
            raw, wavelength_column="Wavelength [nm]", delimiter="\t", encoding="latin-1"
        )
        info = {}
        for r in csv.DictReader(meta.open(encoding="latin-1"), delimiter="\t"):
            row = {_norm(k): (v or "") for k, v in r.items() if k}
            info[row.get("spectrum name", "").split("xls_")[-1]] = row

        def substrate(m: dict[str, str], info: dict[str, Any] = info) -> str:
            return info.get(m["name"], {}).get("substrate specification", "").strip()

        wallis.append(_tag(lib, "wallis", substrate, lambda m: "air"))
    out["wallis"] = wallis
    return out


def class_members(libraries: dict[str, Any]) -> dict[str, list[tuple[Any, list[int]]]]:
    """Per class, the (library, row indices) holding its registered labels."""
    classes = load_registration()["classes"]
    members: dict[str, list[tuple[Any, list[int]]]] = {}
    for name, sources in classes.items():
        members[name] = []
        for source, labels in sources.items():
            if source == "note":
                continue
            wanted = {_norm(label) for label in labels}
            libs = libraries[source] if isinstance(libraries[source], list) else [libraries[source]]
            for lib in libs:
                rows = [
                    i
                    for i, m in enumerate(lib.metadata)
                    if _norm(m["taxon"]) in wanted
                    or (source == "vigo" and any(_norm(w) in _norm(m["taxon"]) for w in wanted))
                ]
                if rows:
                    members[name].append((lib, rows))
    return members


def label_check() -> dict[str, Any]:
    """Spectra and taxa per class and source: a label check, before any contrast."""
    members = class_members(load_libraries())
    report: dict[str, Any] = {}
    for name, parts in members.items():
        counts: dict[str, int] = {}
        taxa: set[tuple[str, str]] = set()
        for lib, rows in parts:
            for i in rows:
                m = lib.metadata[i]
                counts[m["source"]] = counts.get(m["source"], 0) + 1
                taxa.add((m["source"], m["taxon"]))
        report[name] = {"n": sum(counts.values()), "by_source": counts, "n_taxa": len(taxa)}
    return report


# -- Band values and taxon-balanced statistics ---------------------------------------


def band_responses() -> dict[str, list[Any]]:
    from kelpobserve.coastal.spectra import read_rsr

    s2 = read_rsr(ACOLITE_RSR / "S2A_MSI_V4.txt")
    neo = read_rsr(ACOLITE_RSR / "PNEO3.txt")
    return {"sentinel2": [s2[b] for b in S2_BANDS], "pleiades_neo": [neo[b] for b in NEO_BANDS]}


def class_values(
    members: dict[str, list[tuple[Any, list[int]]]], bands: list[Any], keep: Any = None
) -> dict[str, dict[str, np.ndarray]]:
    """Per class: band values ``(n, n_bands)`` and the taxon, source and geometry of each row."""
    from kelpobserve.coastal.spectra import resample

    out: dict[str, dict[str, np.ndarray]] = {}
    for name, parts in members.items():
        values, taxa, sources, geometry = [], [], [], []
        for lib, rows in parts:
            chosen = [i for i in rows if keep is None or keep(lib.metadata[i])]
            if not chosen:
                continue
            sub = lib.select(lambda m, ids=set(id(lib.metadata[i]) for i in chosen): id(m) in ids)
            values.append(resample(sub, bands))
            taxa += [f"{m['source']}:{m['taxon']}" for m in sub.metadata]
            sources += [m["source"] for m in sub.metadata]
            geometry += [m["geometry"] for m in sub.metadata]
        if values:
            v = np.vstack(values)
            # Blue, green and red must be present; an edge band (B1 or B5 here, the
            # Pléiades Neo coastal blue) may be missing where a library stops short
            # of it, and is then skipped band by band.
            ok = np.all(np.isfinite(v[:, 1:4]), axis=1)
            out[name] = {
                "values": v[ok],
                "taxa": np.array(taxa)[ok],
                "sources": np.array(sources)[ok],
                "geometry": np.array(geometry)[ok],
            }
    if "rock" in out and "sand" in out:
        out["bare"] = {k: np.concatenate([out["rock"][k], out["sand"][k]]) for k in out["rock"]}
    return out


def taxon_weights(taxa: np.ndarray) -> np.ndarray:
    names, inverse, counts = np.unique(taxa, return_inverse=True, return_counts=True)
    return 1.0 / (counts[inverse] * names.size)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    finite = np.isfinite(values)
    if not finite.any():
        return float("nan")
    values, weights = values[finite], weights[finite]
    order = np.argsort(values)
    cum = np.cumsum(weights[order]) - 0.5 * weights[order]
    return float(np.interp(q, cum / weights.sum(), values[order]))


def class_table(data: dict[str, dict[str, np.ndarray]], band_names: tuple[str, ...]) -> dict:
    table: dict[str, Any] = {}
    for name, d in data.items():
        w = taxon_weights(d["taxa"])
        table[name] = {
            "n": int(d["values"].shape[0]),
            "n_taxa": int(np.unique(d["taxa"]).size),
            "scored": bool(d["values"].shape[0] >= MIN_SPECTRA),
            "bands": {
                b: {
                    "n": int(np.isfinite(d["values"][:, j]).sum()),
                    "median": round(weighted_quantile(d["values"][:, j], w, 0.5), 5),
                    "p25": round(weighted_quantile(d["values"][:, j], w, 0.25), 5),
                    "p75": round(weighted_quantile(d["values"][:, j], w, 0.75), 5),
                }
                for j, b in enumerate(band_names)
            },
        }
    return table


def draw(d: dict[str, np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """Taxon-balanced draws, ``(n_bands, N_DRAWS)``: a taxon uniformly, then a spectrum."""
    names = np.unique(d["taxa"])
    picks = rng.integers(0, names.size, N_DRAWS)
    rows = np.empty(N_DRAWS, dtype=int)
    for k, name in enumerate(names):
        members = np.flatnonzero(d["taxa"] == name)
        chosen = picks == k
        rows[chosen] = members[rng.integers(0, members.size, int(chosen.sum()))]
    return d["values"][rows].T


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    p, q = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if p.size == 0 or q.size == 0:
        return float("nan")
    ranks = np.argsort(np.argsort(np.concatenate([p, q]))) + 1
    return float((ranks[: p.size].sum() - p.size * (p.size + 1) / 2) / (p.size * q.size))


def nd(x: np.ndarray, i: int, j: int) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return (x[i] - x[j]) / (x[i] + x[j])


def features(x: np.ndarray, band_names: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Every single band and every normalised difference of two bands."""
    out = {b: x[j] for j, b in enumerate(band_names)}
    for i, j in itertools.combinations(range(len(band_names)), 2):
        out[f"ND({band_names[i]},{band_names[j]})"] = nd(x, i, j)
    return out


def orientation(kelp: np.ndarray, other: np.ndarray, band_names: tuple[str, ...]) -> dict:
    """+1 where kelp has the higher median in air, else -1, per feature."""
    fk, fo = features(kelp, band_names), features(other, band_names)
    return {k: 1.0 if np.nanmedian(fk[k]) >= np.nanmedian(fo[k]) else -1.0 for k in fk}


def separability(
    kelp: np.ndarray, other: np.ndarray, sign: dict[str, float], band_names: tuple[str, ...]
) -> dict[str, float]:
    fk, fo = features(kelp, band_names), features(other, band_names)
    return {k: round(auc(sign[k] * fk[k], sign[k] * fo[k]), 4) for k in fk}


# -- Water column --------------------------------------------------------------------


def load_iops() -> dict[str, Any]:
    """Forced-560 QAA at Sesimbra per date at the five bands, and their median."""
    from kelpobserve.coastal.optics.qaa import SceneIOPs

    per_date = {}
    for date in E2C_DATES:
        raw = json.loads((E2C / date / "forced_560" / "iops.json").read_text())["iops"]
        wl = list(raw["wavelengths_nm"])
        idx = [wl.index(w) for w in IOP_NM]
        per_date[date] = (np.array(raw["a"])[idx], np.array(raw["bb"])[idx])
    sets = dict(per_date)
    sets["median"] = (
        np.median([v[0] for v in per_date.values()], axis=0),
        np.median([v[1] for v in per_date.values()], axis=0),
    )
    return {
        name: SceneIOPs(
            wavelengths_nm=np.array(IOP_NM),
            a=a,
            bb=bb,
            kd=np.full(len(IOP_NM), np.nan),
            a_cdm_443=float("nan"),
            bbp_555=float("nan"),
            y=float("nan"),
            n_deep_pixels=0,
            ref_wavelength_nm=560.0,
        )
        for name, (a, bb) in sets.items()
    }


def _attenuation(iops: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from kelpobserve.coastal.inversion import lee

    u = iops.bb / (iops.a + iops.bb)
    inv_cos = 1 / np.cos(np.deg2rad(lee._subsurface_solar_zenith_deg(SOLAR_ZENITH_DEG)))
    k_c = (inv_cos + lee._du_column(u)) * (iops.a + iops.bb)
    k_b = (inv_cos + lee._du_bottom(u)) * (iops.a + iops.bb)
    return u, k_c, k_b


def forward(rho_b: np.ndarray, depth: float, iops: Any) -> np.ndarray:
    """Above-water Rrs ``(n_bands, n)`` from bottom reflectance, as lee.forward_model."""
    from kelpobserve.coastal.optics import water

    u, k_c, k_b = _attenuation(iops)
    deep = (0.089 * u + 0.125 * u * u)[:, None]
    sub = deep * (1 - np.exp(-k_c[:, None] * depth)) + rho_b / np.pi * np.exp(-k_b[:, None] * depth)
    return np.asarray(water.rrs_subsurface_to_above(sub))


def retrieve(rrs: np.ndarray, depth: float, iops: Any) -> np.ndarray:
    from kelpobserve.coastal.inversion.lee import invert_scene

    n = rrs.shape[1]
    return invert_scene(
        rrs[:, :, None],
        iops.wavelengths_nm,
        np.full((n, 1), depth),
        iops,
        np.ones((n, 1), bool),
        SOLAR_ZENITH_DEG,
    )[:, :, 0]


def through_water(rho_b: np.ndarray, depth: float, iops: Any, sigma: float, rng: Any) -> np.ndarray:
    rrs = forward(rho_b, depth, iops) + rng.normal(0, sigma / np.pi, rho_b.shape)
    return retrieve(rrs, depth, iops)


# -- Analyses ------------------------------------------------------------------------


def _scored(data: dict[str, Any]) -> list[str]:
    return [c for c in data if c != "kelp" and data[c]["values"].shape[0] >= MIN_SPECTRA]


def analyse(data: dict[str, Any], iop_sets: dict[str, Any], *, full: bool) -> dict[str, Any]:
    """Air and water separability of kelp against every scored class.

    ``full`` runs every IOP set and sigma; otherwise only the primary ones.
    """
    rng = np.random.default_rng(SEED)
    classes = [c for c in _scored(data) if c not in ("sediment", "sediment_film")]
    draws = {c: draw(data[c], rng) for c in ["kelp", *classes]}
    signs = {c: orientation(draws["kelp"], draws[c], S2_BANDS) for c in classes}
    result: dict[str, Any] = {"classes": classes, "air": {}, "water": {}}
    for sigma in SIGMAS if full else (PRIMARY_SIGMA,):
        noisy = {c: v + rng.normal(0, sigma, v.shape) for c, v in draws.items()}
        result["air"][f"{sigma:g}"] = {
            c: separability(noisy["kelp"], noisy[c], signs[c], S2_BANDS) for c in classes
        }
    for iop_name in iop_sets if full else ("median",):
        iops = iop_sets[iop_name]
        for sigma in SIGMAS if full else (PRIMARY_SIGMA,):
            key = f"{iop_name}|{sigma:g}"
            result["water"][key] = {}
            for depth in DEPTHS_M:
                seen = {c: through_water(v, depth, iops, sigma, rng) for c, v in draws.items()}
                result["water"][key][f"{depth:g}"] = {
                    c: separability(seen["kelp"], seen[c], signs[c], S2_BANDS) for c in classes
                }
    return result


def neo_air(data: dict[str, Any]) -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    classes = [c for c in _scored(data) if c not in ("sediment", "sediment_film")]
    draws = {c: draw(data[c], rng) for c in ["kelp", *classes]}
    noisy = {c: v + rng.normal(0, PRIMARY_SIGMA, v.shape) for c, v in draws.items()}
    return {
        c: separability(
            noisy["kelp"], noisy[c], orientation(draws["kelp"], draws[c], NEO_BANDS), NEO_BANDS
        )
        for c in classes
    }


def z_max_table(table: dict[str, Any], iop_sets: dict[str, Any]) -> dict[str, Any]:
    from kelpobserve.coastal.detectability import z_max_from_k

    out: dict[str, Any] = {}
    kelp = table["kelp"]["bands"]
    for iop_name, iops in iop_sets.items():
        _, _, k_b = _attenuation(iops)
        out[iop_name] = {
            "k_bottom_per_m": dict(zip(S2_BANDS, np.round(k_b, 4).tolist(), strict=True))
        }
        for c, row in table.items():
            if c in ("kelp", "sediment", "sediment_film") or not row["scored"]:
                continue
            out[iop_name][c] = {
                f"{sigma:g}": {
                    b: round(
                        float(
                            z_max_from_k(
                                k_b[j],
                                sigma,
                                abs(kelp[b]["median"] - row["bands"][b]["median"]),
                                T_AW,
                            )
                        ),
                        2,
                    )
                    for j, b in enumerate(S2_BANDS)
                }
                for sigma in SIGMAS
            }
    return out


def geometry_table(members: dict[str, Any], bands: list[Any]) -> dict[str, Any]:
    """Fuller genera with ex-situ and in-situ spectra: in-situ over ex-situ medians."""
    from kelpobserve.coastal.spectra import resample

    lib = next(lb for parts in members.values() for lb, _ in parts if lb.source == "fuller")
    values = resample(lib, bands)
    out: dict[str, Any] = {}
    for genus in sorted({m["taxon"] for m in lib.metadata}):
        rows = {
            g: [i for i, m in enumerate(lib.metadata) if m["taxon"] == genus and m["geometry"] == g]
            for g in ("air", "in_situ_unstated")
        }
        if min(len(r) for r in rows.values()) < 2:
            continue
        air = np.nanmedian(values[rows["air"]], axis=0)
        wet = np.nanmedian(values[rows["in_situ_unstated"]], axis=0)
        out[genus] = {
            "n_ex_situ": len(rows["air"]),
            "n_in_situ": len(rows["in_situ_unstated"]),
            "ratio_in_over_ex": dict(zip(S2_BANDS, np.round(wet / air, 3).tolist(), strict=True)),
        }
    return out


def deepest(water: dict[str, Any], key: str, cls: str, feature: str) -> float:
    reached = [
        float(d) for d, v in water[key].items() if v.get(cls, {}).get(feature, 0) >= AUC_SEPARABLE
    ]
    return max(reached) if reached else 0.0


def decisions(pooled: dict, air_only: dict, table: dict, cc_table: dict, zmax: dict) -> dict:
    from kelpobserve.coastal.detectability import SUBSTRATE_ALBEDO

    key = f"median|{PRIMARY_SIGMA:g}"
    # D1: E4t-B index among ND of the 10 m bands, kelp against bare.
    candidates = {}
    for i, j in itertools.combinations(S2_10M, 2):
        f = f"ND({i},{j})"
        a = pooled["air"][f"{PRIMARY_SIGMA:g}"]["bare"][f]
        w = pooled["water"][key][f"{FLAT_TEST_DEPTH_M:g}"]["bare"][f]
        candidates[f] = {"auc_air": a, "auc_5m": w, "score": min(a, w)}
    chosen = max(candidates, key=lambda f: candidates[f]["score"])
    reported = {
        f: {
            "auc_air": pooled["air"][f"{PRIMARY_SIGMA:g}"]["bare"][f],
            "auc_5m": pooled["water"][key][f"{FLAT_TEST_DEPTH_M:g}"]["bare"][f],
        }
        for f in features(np.zeros((len(S2_BANDS), 1)), S2_BANDS)
        if f.startswith("ND(") and f not in candidates
    }
    # D2: brown against red.
    d2: dict[str, Any] = {}
    for cls in UNDERSTORY:
        best = {}
        for f in pooled["water"][key]["1"][cls]:
            dp = deepest(pooled["water"], key, cls, f)
            da = deepest(air_only["water"], key, cls, f) if cls in air_only["classes"] else 0.0
            best[f] = {"pooled_m": dp, "air_only_m": da, "both_m": min(dp, da)}
        top = max(best, key=lambda f: best[f]["both_m"])
        d2[cls] = {"best_feature": top, **best[top], "all": best}
    revise = any(d2[c]["both_m"] >= 1.0 for c in UNDERSTORY)
    # D3: CC-BY class medians against the library defaults.
    groups = {
        "bare": ("sand", "rock"),
        "brown": ("kelp", "brown_other"),
        "red": ("red_foliose", "turf"),
    }
    d3 = {}
    for g, cls in groups.items():
        rows = [cc_table[c] for c in cls if c in cc_table]
        d3[g] = {}
        for j, b in enumerate(S2_BANDS[:4]):
            lib_value = list(SUBSTRATE_ALBEDO[g].values())[j]
            ours = float(np.mean([r["bands"][b]["median"] for r in rows]))
            d3[g][b] = {
                "library": lib_value,
                "cc_by_median": round(ours, 5),
                "ratio": round(ours / lib_value, 3),
            }
    return {
        "D1_flat_test_index": {
            "chosen": chosen,
            "candidates": candidates,
            "reported_not_chosen": reported,
        },
        "D2_brown_against_red": {"revise_statement": revise, "by_class": d2},
        "D3_library_defaults": d3,
        "D4_detectability": {
            "table": "z_max",
            "baltic_statement": "brown against bare separable to about 20 m at 444 nm",
        },
    }


# -- Stages --------------------------------------------------------------------------


def _files() -> list[Path]:
    import kelpobserve.coastal.detectability as detectability
    import kelpobserve.coastal.inversion.lee as lee
    import kelpobserve.coastal.optics.qaa as qaa
    import kelpobserve.coastal.optics.water as water
    import kelpobserve.coastal.spectra as spectra

    files = {REGISTRATION, Path(__file__), *SOURCE_FILES.values()}
    files |= set(_fuller_files()) | set(_vigo_files())
    files |= {p for pair in _wallis_files() for p in pair}
    files |= {ACOLITE_RSR / "S2A_MSI_V4.txt", ACOLITE_RSR / "PNEO3.txt"}
    files |= {E2C / d / "forced_560" / "iops.json" for d in E2C_DATES}
    files |= {Path(m.__file__) for m in (spectra, detectability, lee, qaa, water)}
    return sorted(p.resolve() for p in files)


def freeze() -> Path:
    if LOCK.exists():
        raise FileExistsError("E5s inputs are already frozen.")
    common._write_json(
        LOCK,
        {
            "schema_version": "1.0",
            "experiment_id": load_registration()["experiment_id"],
            "frozen_at": dt.datetime.now(dt.UTC).isoformat(),
            "files": {str(p): common.digest(p) for p in _files()},
        },
    )
    return LOCK


def verify_lock() -> dict[str, Any]:
    lock = json.loads(LOCK.read_text())
    for name, expected in lock["files"].items():
        if common.digest(Path(name)) != expected:
            raise ValueError(f"Locked input changed: {name}")
    return {"path": str(LOCK), "sha256": common.digest(LOCK), "n_files": len(lock["files"])}


def run() -> Path:
    lock = verify_lock()
    libraries = load_libraries()
    members = class_members(libraries)
    bands = band_responses()
    iop_sets = load_iops()
    is_air = lambda m: m["geometry"] == "air"  # noqa: E731
    not_rp = lambda m: m["source"] != "rowan_phinn"  # noqa: E731
    s2 = class_values(members, bands["sentinel2"])
    s2_air = class_values(members, bands["sentinel2"], keep=is_air)
    neo = class_values(members, bands["pleiades_neo"])
    table = class_table(s2, S2_BANDS)
    table_air = class_table(s2_air, S2_BANDS)
    cc_table = class_table(class_values(members, bands["sentinel2"], keep=not_rp), S2_BANDS)
    pooled = analyse(s2, iop_sets, full=True)
    air_only = analyse(s2_air, iop_sets, full=False)
    zmax = z_max_table(table, iop_sets)
    result = {
        "experiment_id": load_registration()["experiment_id"],
        "ran_at": dt.datetime.now(dt.UTC).isoformat(),
        "input_lock": lock,
        "citation_rowan_phinn": load_registration()["libraries"]["rowan_phinn"]["citation"],
        "labels": label_check(),
        "band_centres_nm": {
            s: {r.name: round(r.centre_nm, 1) for r in rs} for s, rs in bands.items()
        },
        "class_bands_sentinel2": {"pooled": table, "air_only": table_air, "cc_by_only": cc_table},
        "class_bands_pleiades_neo": class_table(neo, NEO_BANDS),
        "iops": {
            n: {"a": np.round(v.a, 5).tolist(), "bb": np.round(v.bb, 5).tolist()}
            for n, v in iop_sets.items()
        },
        "separability_pooled": pooled,
        "separability_air_only": air_only,
        "separability_pleiades_neo_air": neo_air(neo),
        "z_max_m": zmax,
        "geometry_fuller": geometry_table(members, bands["sentinel2"]),
        "decisions": decisions(pooled, air_only, table, cc_table, zmax),
    }
    common._write_json(OUTPUT / "results.json", result)
    verify_lock()
    return OUTPUT / "results.json"
