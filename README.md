# KelpObserve Sesimbra 2026: reproduction archive of the satellite kelp-cover benchmark, its revisions and the optics prerequisites

Code, registrations, input locks and archived runs behind the Sesimbra analyses of KelpObserve, an AQUASERV Transnational Access project (PID 42082) hosted by CCMAR (Centro de Ciências do Mar do Algarve, Universidade do Algarve, Faro, Portugal). The benchmark at the centre of these analyses asks how well Sentinel-2 and Pléiades Neo imagery predict the kelp cover recorded by divers at one rocky site on the Arrábida coast near Sesimbra (38.409° N, 9.202° W), and whether a correction for the water column improves that prediction. The command `kelpobserve replay` reruns an experiment in a fresh workspace and compares every output, file by file, with the archived run.

The technical report published on 3 October 2026 is at [kelp-observe.oceanstream.io/sesimbra/technical-report.html](https://kelp-observe.oceanstream.io/sesimbra/technical-report.html). Its source, `benchmark/docs/technical-report.md`, is a working draft updated on 10 October 2026, whose sections 6.6 to 6.10 report exploratory experiments of 9 and 10 October; their code and outputs are not part of this repository. The field observations are published as the dataset *KelpObserve Sesimbra 2026: diver quadrat surveys of* Saccorhiza polyschides *kelp and water-column sonde profiles, Arrábida coast, Portugal* ([doi:10.5281/zenodo.23282075](https://doi.org/10.5281/zenodo.23282075)); the transect and mapping videos are in a companion record ([doi:10.5281/zenodo.23350632](https://doi.org/10.5281/zenodo.23350632)), and the software that built the dataset is [KelpObserve/arrabida-dataset](https://github.com/KelpObserve/arrabida-dataset) ([doi:10.5281/zenodo.23454705](https://doi.org/10.5281/zenodo.23454705)).

## The analyses

The registered benchmark, v1, was frozen on 20 September 2026. Its reference is the cover of *Saccorhiza polyschides* that CCMAR divers estimated in 0.5 × 0.5 m quadrats: the models were developed on the July campaign by leave-one-transect-out cross-validation, sealed, and applied once to the September campaign, which had been held out. Five exploratory revisions followed on 2 October 2026: r1, r3, r4, r5 and r6. The draft protocol r2 was never executed; r3 ran under the name r2, and its run records carry that name (`benchmark/registration/revision-r3.json` records the rename).

The optics prerequisite is E5s, which carried kelp spectral libraries through Sesimbra's water, with the water optics that the registered analysis E2c fitted on two winter Sentinel-2 scenes using the quasi-analytical algorithm (QAA) with a forced 560 nm reference. E2c's results are inputs here; E2c itself is not replayed. A draft registration, r7, would give the water-column correction an attenuation taken from Copernicus Marine ocean colour.

## Replays

Each experiment was replayed on 2 November 2026 in a fresh workspace, with Python 3.12.12 and the package versions pinned in `requirements-replay.txt`. Every output was compared file by file with the archived run.

| Experiment | Report section | Stages replayed | Result |
|---|---|---|---|
| `v1` | 5, Registered results | `fit`, `evaluate` | Every result file byte-identical, including the two sealed model files |
| `v1-readiness` | 4.1 and 4.2 | `scene-qc`, `optics-reference`, `cross-sensor` | Byte-identical: `scene_qc.json`, `reference_iops.json` (with the Kd(486) of 1.32 and 1.26 m⁻¹) and `cross_sensor.json` |
| `r1` | 6.1 | `qc`, `models`, `models --variant V1`, `radiometry`, `injection` | All result files identical; `injection.json` equal (it differs only in the recorded digest of the E5s results) |
| `r3` | 6.2 | `mask`, `models`, `injection`, `windows` | All result files identical or equal |
| `r4` | 6.3 | `run --draws 100` | Identical: `draws.csv` of 26,000 rows; `envelope.json` equal (it differs only in the recorded digest of the E5s results) |
| `r5` | 6.4 | `run` | Identical: `diagnosis.json` |
| `r6` | 6.5 | `bias`, `windows`, `pneo`, `s2`, `offset` | All result files identical, including the three `draws.csv` tables of 22,400, 22,800 and 41,600 rows |
| `e5s` | [E5s report](optics/reports/e5s-kelp-spectra.md) | `run` | Equal: `results.json` differs only in run time and lock record |

*Identical* means the same bytes. *Equal* means the same values once recorded paths are mapped to the same roots (see *Recorded paths*) and provenance fields are set aside: creation times, git state, runner and registration checksums, and the recorded digests of inputs this repository stores with root tokens, such as the E5s results; for r3, the r2 names in the archived records are mapped to r3. The provenance differences are listed in each replay's report, `out/replay/<experiment>-<time>/replay.json`.

Two r1 stages, `screen` and `windows`, read Sentinel-2 assets from the Planetary Computer and run only when named: `kelpobserve replay r1 --stages "revision screen" "revision windows"`. Their replay of 18 October 2026 against the live service was identical: `screen.json`, the 374 rows of `overpasses.csv` and `windows.json`. The service's holdings can change, so a later replay of these two stages may differ for reasons outside the archive.

### Recorded paths

The registrations, input locks and archived runs write each path they record relative to a named root:

| Root | Resolves to |
|---|---|
| `$DATA` | The `data` path of `kelpobserve.toml`; `$DATA/external`, the two EMODnet rasters, resolves to `benchmark/inputs/emodnet/` |
| `$CACHE` | The `cache` path, or for E5s the `optics_cache` path |
| `$ACOLITE` | The `acolite` path |
| `$REG` | The experiment's registration: `benchmark/registration/` or `optics/registration/` |
| `$OUT` | The run outputs: `benchmark/archive/out/` for the archived runs, the workspace's `out/` in a replay |
| `$ROOT` | The experiment's folder: `benchmark/` or `optics/`, or the replay's workspace |
| `$LIB` | The package holding the coastal library, `src/kelpobserve/` |
| `$E5S` | The E5s results, `optics/reports/e5s-kelp-spectra/results.json` |

A replay resolves each root before any stage runs (`src/kelpobserve/relocation.py` and `src/kelpobserve/replay.py`), and maps the paths a replayed run records back to the same roots before comparing it with the archived run.

### Checks before the stages run

A replay checks the inputs it can check before any stage runs, and stops if they do not verify. The input locks record the SHA-256 of every input of v1 and E5s: the external inputs as they were frozen before each experiment ran, and the files of this repository, the coastal library among them, as published. For v1, all 91 files of the input lock verify. For v1-readiness, the 250 prepared image products match the hashes in the prepare manifest. For e5s, all 648 files of the E5s lock verify.

### Provenance gaps

The runner that produced the v1 runs of 20 September 2026 was not preserved; `benchmark/benchmark/run_benchmark.py` stands in for it, and with it every v1 result file replays byte-identical. Some revision modules and registrations were edited after their runs; the replays use the current files, their results reproduce, and the checksums that differ from those in the run manifests are listed as provenance differences.

All experiments run on one version of the coastal library, `src/kelpobserve/coastal`. It tabulates pure-water absorption to 2300 nm (Pope and Fry 1997 to 700 nm, Röttgers et al., WOPP v3, beyond); the v1 runs used an earlier version that held it flat at 0.624 m⁻¹ beyond 700 nm, so the v1 replays select that tabulation, and with it the near-infrared bands of `reference_iops.json` reproduce.

### Not replayed

| Stages | Note |
|---|---|
| `r1 acolite`, `r3 acolite`, `r6 check` | They write into the shared imagery cache. |
| `r3 fetch` | It downloads from the Copernicus Data Space. |
| v1 `inventory`, `glint-screen`, `water-levels`, `fetch`, `prepare`, `reference-region` | They lie upstream of the input lock; their archived outputs are verified by hash. |
| E2c | Its results (`optics/reports/e2c-sesimbra/`) are the water optics E5s reads. It is not replayed: the 56 ACOLITE products of the two winter scenes it fitted are no longer in the optics cache. |

## Running a replay

### Inputs the repository does not carry

The replays read inputs from outside the repository, at paths set in `kelpobserve.toml`, a git-ignored copy of `kelpobserve.example.toml`. Relative paths resolve against the repository root, and the environment variable `KELPOBSERVE_CONFIG` can name another file. A replay that needs a path that is unset or missing stops and names it.

| Key | Content | Needed by |
|---|---|---|
| `data` | The CCMAR deliveries as delivered and the Pléiades Neo archives, in the layout `ccmar-diving/`, `pneo/`. No replayed stage reads them; v1's input lock pins them | v1 |
| `cache` | The imagery cache, about 57 GB: ACOLITE products made with ACOLITE at revision 64a02ff, Sentinel-2 Level-1C archives from the Copernicus Data Space and the extracted Pléiades Neo images | v1, v1-readiness, r1, r3, r4, r5, r6 |
| `optics_cache` | The optics cache, with the spectral libraries | e5s |
| `acolite` | A checkout of [ACOLITE](https://github.com/acolite/acolite) at revision 64a02ff. The benchmark stages read its sensor defaults (the water-mask band and threshold), and E5s reads its sensor response functions | v1-readiness, r1, r3, r4, r5, r6, e5s |
| `workspace` | Where replays write their workspaces; `out/replay` by default | optional |

The 193 Sentinel-2 ACOLITE products that the replayed stages read are published byte-identical as *KelpObserve Sesimbra 2026: ACOLITE-corrected Sentinel-2 scenes read by the benchmark replays* ([doi:10.5281/zenodo.23844565](https://doi.org/10.5281/zenodo.23844565)). Unzipping its 17 zips into one folder recreates the cache's `acolite/` and `acolite-r1/V1/` folders, and its `manifest.json` gives the SHA-256 of each file. Set `cache` to that folder; with the ACOLITE checkout, r3 then replays in full. The other experiments also read Pléiades Neo products. On macOS, unzipping onto an exFAT, FAT or network volume writes a `._` file beside each GeoTIFF, which the stages would read as a raster; delete them with `find <folder> -name '._*' -delete`.

The v1 input lock pins the CCMAR deliveries as they were delivered. The Zenodo record carries sanitised release copies whose bytes differ, and a replay from those copies stops at the lock check. The spectral libraries, among them the kelp spectra of Rowan et al. (2025, [doi:10.48610/c5f745e](https://doi.org/10.48610/c5f745e)), are used under their terms and are not redistributed; `THIRD_PARTY_NOTICES.md` names each of them. Pléiades Neo imagery © Airbus DS was supplied through ESA Third Party Missions project PP0110832 and cannot be redistributed; the repository holds no Pléiades Neo imagery.

### Installation

```bash
git clone https://github.com/KelpObserve/kelpobserve
cd kelpobserve
uv venv --python 3.12.12
source .venv/bin/activate
uv pip install -r requirements-replay.txt -e .
cp kelpobserve.example.toml kelpobserve.toml   # then set the paths
```

`requirements-replay.txt` pins the package versions installed in the environment that ran the benchmark and its revisions. `.env.example` names the credentials of the commands that download; a copy saved as `.env` is git-ignored. `kelpobserve kd490 fetch` needs the Copernicus Marine pair, from the environment or from the file given with `--env-file`. The benchmark's own download stages, which the replays do not run, also take `--env-file` and need the Copernicus Marine and Copernicus Data Space pairs.

### Commands

| Command | Action |
|---|---|
| `kelpobserve experiments` | Lists the eight experiments with their stages and the stages left out, each with its reason. |
| `kelpobserve replay <experiment> [--stages …]` | Replays `v1`, `v1-readiness`, `r1`, `r3`, `r4`, `r5`, `r6` or `e5s` in a new workspace under `out/replay/` and compares every output with the archived run. `--stages` limits the replay to the stages named, by label (as `experiments` prints them) or by output folder. |
| `kelpobserve recompare <workspace>` | Compares a finished replay again without rerunning it. |
| `kelpobserve kd490 fetch [--env-file .env]` | Downloads the CMEMS KD490 subset and records its version and SHA-256. |
| `kelpobserve kd490 table` | Writes the attenuation table of the two Pléiades Neo scenes to `out/kd490/<time>/attenuation.json`. |

```bash
kelpobserve experiments
kelpobserve replay v1
kelpobserve replay r1 --stages "revision models"
```

`replay` prints a summary of its checks and of each stage, and exits with status 0 when every stage reproduced. It builds the workspace from the experiment's registration and the archived runs it reads, runs each stage in a separate interpreter and writes its runs only into that workspace; the archive is never modified. The experiment's modules run unchanged. A start-up hook, `src/kelpobserve/replay_site/sitecustomize.py`, points the paths they are anchored to into the workspace and answers their own lock verification with the check the replay has already made against this checkout.

### Running time

On the machine that made the replays, v1 and v1-readiness each took under a minute, e5s a few seconds, r4 about 2 minutes, r3 about 3 minutes, r1 about 6 minutes (its two network stages another 6, when replayed on 18 October 2026), and r6 about 6 minutes.

## CMEMS KD490 as an external attenuation (exploratory)

The water-column-corrected branch of v1 could not run: the in-scene QAA fit over the Pléiades Neo deep-water reference raised plausibility flags on both dates, with Kd(486) of 1.32 and 1.26 m⁻¹. The module `src/kelpobserve/attenuation` takes an attenuation from outside the image, from the Copernicus Marine Atlantic ocean-colour product OCEANCOLOUR_ATL_BGC_L3_MY_009_113 ([doi:10.48670/moi-00286](https://doi.org/10.48670/moi-00286)), dataset `cmems_obs-oc_atl_bgc-transp_my_l3-multi-1km_P1D`: daily KD490 at 1 km, with a per-pixel uncertainty.

`kelpobserve kd490 fetch` freezes a subset of it (the Sesimbra region, 1 December 2025 to 30 September 2026, dataset version 202603, 2.8 MB) in `.cache/cmems-kd490/` and records its version and SHA-256 in `benchmark/registration/evidence/cmems-kd490-subset.json`. The subset itself is not in the repository. The record also holds a digest of the subset's values and the dates the catalogue gives for the dataset part: released on 31 March 2026, to be retired on 24 November 2026. A fetch never replaces the record. When the subset is present it is only checked; a new download is accepted when its values or its bytes match the record, and is otherwise recorded beside it and reported. `kelpobserve kd490 table` opens a subset only when it matches the record.

`kelpobserve kd490 table` takes, for each Pléiades Neo scene, the valid pixel nearest the quadrat cluster within 5 km and 3 days, and extends its Kd(490) to the scene's bands with Sesimbra's own shape of the inherent optical properties: the median non-water absorption and particulate backscatter of E5s, from the forced-560 QAA fits to the winter Sentinel-2 scenes of 24 December 2025 and 22 February 2026, scaled by the one factor at which their Lee et al. (2005) Kd(490) equals the CMEMS value. For both scenes the nearest valid pixel lies 0.29 km away, on the day itself.

| Scene | CMEMS KD490 (m⁻¹) | Uncertainty | Implied Kd at 486, 563 and 655 nm (m⁻¹) | In-scene QAA Kd(486) (m⁻¹) |
|---|---|---|---|---|
| Pléiades Neo, 8 July 2026 | 0.086 | 39 % | 0.087, 0.112, 0.439 | 1.32 |
| Pléiades Neo, 3 September 2026 | 0.153 | 84 % | 0.154, 0.153, 0.482 | 1.26 |

The rule that decides when CMEMS stands in for the in-scene fit, `choose()` in `src/kelpobserve/attenuation/fallback.py`, is not yet written, and the table records no choice for either scene. The experiment that would use it, r7, is a draft registration (`benchmark/registration/revision-r7.json`) that lists the decisions still open; no corrected-branch result has been computed. Among the limits the draft records, the 1 km pixel lies near cliffs and partly over optically shallow water, and a correct attenuation does not remove the radiometric residual of 0.01 to 0.04 that r1 to r6 found as large as the seabed contrast.

## Layout

| Path | Content |
|---|---|
| `benchmark/benchmark/` | The analysis modules: `run_benchmark.py` with the v1 stages, `revision.py` and `revision_r3.py` to `revision_r6.py` with the revisions, and `report_maps.py` and `habitat_map.py`, which draw the report's maps and charts. |
| `benchmark/registration/` | Registrations and the input lock. `experiment.json` is the v1 registration with its amendments, and `input-lock.json` pins every input of v1; `revision-r1.json` to `revision-r6.json` register the revisions; `evidence/` holds the record of the frozen CMEMS subset; `revision-r7.json` is a draft that has not been registered. |
| `benchmark/inputs/` | The registered habitat area and processing extent, and in `emodnet/` the two EMODnet depth rasters the stages read, cut to the Sesimbra area (see `THIRD_PARTY_NOTICES.md`). |
| `benchmark/archive/` | `out/` holds the archived run outputs behind sections 2 to 6.5 of the published report. `index.json` lists each of their files with its SHA-256; the six it marks as not copied, the two sealed model files and four `draws.csv` tables, are compared by hash. |
| `benchmark/docs/` | The technical report and its figures, the registered protocol (`plan.md`), the supporting reports on the September inputs, the scene selection, the optics readiness, the registered results and the revisions r1 and r3 to r6. |
| `benchmark/tests/` | The benchmark's tests. |
| `optics/coastal_validation/` | The E5s module, `e5s.py`, and `common.py` with the folder it reads and writes. |
| `optics/registration/` | The registrations of E2c and E5s, and the E5s input lock. |
| `optics/reports/` | The reports of E2c and E5s, with their archived results. |
| `src/kelpobserve/coastal/` | The coastal library the experiments run on: water optics, QAA and the Lee et al. inversion, masks, detectability, bathymetry, scene and raster IO, the ACOLITE driver, scene acquisition and run provenance; `src/kelpobserve/storage/` holds its object-storage helpers. |
| `src/kelpobserve/` | The replay tool and the attenuation module. |
| `tests/` | Tests of the replay tool and the attenuation module. |

## Tests

```bash
pytest
```

runs the nine tests of the replay tool and the attenuation module, in `tests/`; they need no external inputs. The benchmark's 92 tests run from its folder:

```bash
cd benchmark
PYTHONPATH=. python -m pytest tests
```

## Licence and citation

The repository's own code and documents are released under the Apache License 2.0 (see `LICENSE`). The data the analyses used and the map backgrounds of the report figures are listed with their terms in `THIRD_PARTY_NOTICES.md`. The field dataset and the Sentinel-2 products are licensed CC BY 4.0 under their Zenodo records.

The archive is deposited on Zenodo as [doi:10.5281/zenodo.23478585](https://doi.org/10.5281/zenodo.23478585), and `CITATION.cff` gives its citation. Cite the field data as:

> Rusu, A., Botezat, O., Paulo, D., Encarnação, J., Rosado, G. (2026). KelpObserve Sesimbra 2026: diver quadrat surveys of Saccorhiza polyschides kelp and water-column sonde profiles, Arrábida coast, Portugal (Version 1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23282075

## Credits

The diver surveys, quadrat photographs and videos were made by the CCMAR Scientific Dive Centre (Diogo Paulo, João Encarnação, Guilherme Rosado).

## Funding

The project leading to these results has received funding from the European Union's Horizon Europe research and innovation programme under grant agreement No 101131121 (AQUASERV), through Transnational Access project PID 42082 (KelpObserve) hosted by CCMAR. This study received support from the Portuguese node of EMBRC-ERIC, specifically EMBRC.PT ALG-01-0145-FEDER-022121. This output reflects only the authors' view and the European Union cannot be held responsible for any use that may be made of the information contained therein.

## Contact

Questions about the archive can be sent through [oceanstream.io/contact](https://oceanstream.io/contact).
