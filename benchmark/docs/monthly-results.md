# Monthly leave-one-quadrat-out results

All values come from [metrics.csv](../results/monthly/metrics.csv). MAE and RMSE are cover percentage points; lower is better. R² is reported in the CSV. Predictions are not clipped. Models are trained independently within each month; settings were not tuned inside this evaluation. These are exploratory results from repeatedly inspected data.

## Identical samples for comparison with depth

July n=30; September n=25. The README figures use the median, depth and 5×5 statistics rows.

| Model | July MAE | July RMSE | September MAE | September RMSE |
|---|---:|---:|---:|---:|
| Mean baseline | 35.31 | 37.91 | 9.08 | 15.93 |
| Median baseline | 39.00 | 41.91 | 6.80 | 15.62 |
| Depth-only ridge | 24.91 | 32.43 | 9.08 | 15.64 |
| Spectral-statistics ridge, 5×5 | 35.24 | 40.08 | 12.40 | 17.56 |
| Single-pixel six-band ridge | 42.51 | 48.38 | 13.06 | 21.77 |

## All measured covers

July n=30; September n=29. These September scores cannot be compared directly with depth on n=25. Additional regressors use 5×5 statistics only.

| Method | July MAE | July RMSE | September MAE | September RMSE |
|---|---:|---:|---:|---:|
| Mean baseline | 35.31 | 37.91 | 12.97 | 18.58 |
| Median baseline | 39.00 | 41.91 | 9.66 | 18.98 |
| Band means, 1×1, ridge | 42.51 | 48.38 | 15.84 | 23.58 |
| Lyzenga, 1×1, ridge | 39.45 | 43.75 | 12.99 | 19.14 |
| Band means + Lyzenga, 1×1, ridge | 40.79 | 45.96 | 17.29 | 25.52 |
| Band means, 3×3, ridge | 35.78 | 39.11 | 15.24 | 23.90 |
| Lyzenga, 3×3, ridge | 38.00 | 42.44 | 15.34 | 22.99 |
| Band means + Lyzenga, 3×3, ridge | 37.18 | 41.76 | 21.46 | 32.59 |
| Spectral statistics, 3×3, ridge | 31.20 | 37.88 | 18.28 | 27.38 |
| Statistics + Lyzenga, 3×3, ridge | 32.71 | 38.95 | 21.24 | 32.12 |
| Band means, 5×5, ridge | 32.62 | 37.07 | 13.61 | 22.34 |
| Lyzenga, 5×5, ridge | 32.52 | 38.66 | 14.62 | 20.24 |
| Band means + Lyzenga, 5×5, ridge | 34.07 | 39.95 | 14.75 | 23.28 |
| Spectral statistics, 5×5, ridge | 35.24 | 40.08 | 12.11 | 17.84 |
| Statistics + Lyzenga, 5×5, ridge | 35.85 | 42.51 | 12.57 | 18.60 |
| Spectral statistics, 5×5, random forest | 34.24 | 37.37 | 13.37 | 18.98 |
| Spectral statistics, 5×5, Extra Trees | 33.76 | 36.36 | 11.93 | 17.72 |
| Spectral statistics, 5×5, RBF SVR | 37.49 | 41.04 | 10.61 | 17.58 |
| Spectral statistics, 5×5, PLS (2 components) | 35.52 | 39.63 | 12.99 | 21.33 |
| Spectral statistics, 5×5, ridge, global empty subtraction | 35.24 | 40.08 | 12.11 | 17.84 |
| Spectral statistics, 5×5, Extra Trees, global empty subtraction | 33.76 | 36.36 | 11.93 | 17.72 |
| Spectral statistics, 5×5, ridge, local empty subtraction | 37.15 | 42.67 | 14.37 | 22.80 |
| Spectral statistics, 5×5, Extra Trees, local empty subtraction | 35.68 | 38.29 | 13.05 | 18.63 |

## September quality-filtered sensitivity

The same n=21 quadrats are used for each method below. These retain the input inventory quality flags, including position/depth concerns. No additional regressors or background-subtraction variants were tested on this subset.

| Method | MAE | RMSE | R² |
|---|---:|---:|---:|
| Mean baseline | 10.57 | 17.46 | -0.102 |
| Median baseline | 7.86 | 17.01 | -0.046 |
| Depth-only ridge | 10.28 | 17.16 | -0.065 |
| Band means, 1×1, ridge | 13.82 | 23.59 | -1.012 |
| Lyzenga, 1×1, ridge | 13.65 | 21.95 | -0.743 |
| Band means + Lyzenga, 1×1, ridge | 17.10 | 29.12 | -2.066 |
| Band means, 3×3, ridge | 14.62 | 24.24 | -1.126 |
| Lyzenga, 3×3, ridge | 14.01 | 22.50 | -0.830 |
| Band means + Lyzenga, 3×3, ridge | 16.08 | 24.07 | -1.095 |
| Spectral statistics, 3×3, ridge | 10.93 | 18.38 | -0.221 |
| Statistics + Lyzenga, 3×3, ridge | 11.49 | 18.93 | -0.296 |
| Band means, 5×5, ridge | 11.64 | 18.53 | -0.241 |
| Lyzenga, 5×5, ridge | 11.30 | 17.97 | -0.168 |
| Band means + Lyzenga, 5×5, ridge | 11.92 | 18.67 | -0.261 |
| Spectral statistics, 5×5, ridge | 14.33 | 19.52 | -0.378 |
| Statistics + Lyzenga, 5×5, ridge | 14.91 | 20.93 | -0.584 |

For methods and study context, see the [technical report](monthly-evaluation.md).

## Interpretation

Depth has lower MAE than the median baseline for July, but not September. For September, median and depth RMSE are nearly equal. Global empty subtraction did not improve results; local subtraction increased MAE. Window and model rankings vary, and none of these comparisons isolates a kelp-specific optical signal from depth, substrate or other scene effects. More independent observations and validation are needed. See the [monthly README](monthly/README.md) for limitations and reproduction commands.
