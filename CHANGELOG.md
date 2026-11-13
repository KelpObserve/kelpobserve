# Changelog

## 1.1.0

- The r7 development experiment, executed on 2 November 2026: its frozen registration and evidence, the `local-cover` command and the result bundles in `benchmark/reports/`. No primary Sentinel-2 or Pléiades Neo retrieval was spectrally compatible with the observed reflectance.
- The CMEMS KD490 matchup is repaired. The selector requires flags=0, positive in-range KD490 and finite, uncensored uncertainty, and stops on unfamiliar flag metadata or units. The pixel chosen before, 0.29 km from the reef, was flagged as land on all seven dates checked; the repaired candidate lies 0.891 km away. `kd490 qualify` writes the matchup audit.
- New commands `optics-inputs fetch`, `optics-inputs audit` and `pneo-optical-review`, with the bundles of the optical diagnosis, the quadrat inspection, the Pléiades Neo glint and water review and the July atmospheric comparison.
- Technical report and results: the September water item cites Copernicus Marine KD490 for the change in attenuation (revised on 13 November 2026).
- The E5s report names the spectral-library loader by its path in this repository.
- The registrations, archive and replays of v1 and r1 to r6 are unchanged.

## 1.0.0

- The registered benchmark v1 with its readiness check, the exploratory revisions r1 and r3 to r6, and the E5s optics analysis with the E2c results it reads: registrations, input locks, archived runs and the `kelpobserve replay` command.
- The coastal library the analyses ran with, as `kelpobserve.coastal`, and the two EMODnet depth rasters the stages read.
- An exploratory module for attenuation from CMEMS KD490, and r7 as a draft registration.
