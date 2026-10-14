# QUEST gravity: verified run results

Run: 2026-10-14T19:09:11.541324+00:00. All figures use actual supplier data.

- Parsed observations: 1,049,547; line sections: 340.
- Per-line block means: 79,791; planned traverse/infill lines: 246.
- Nonfinite required rows excluded: 0.
- Duplicate flight/line/time records: 0.
- Published Bouguer range: -183.926 to -72.725 mGal.
- Geographic/projected coordinate discrepancy: median 1.303 m; max 1.312 m.

## Held-line comparison

| Method | RMSE (mGal) | MAE (mGal) | Bias (mGal) | Supported blocks |
|---|---:|---:|---:|---:|
| linear | 1.912 | 1.455 | -0.011 | 98.8% |
| idw | 1.934 | 1.476 | -0.007 | 98.8% |
| nearest | 2.981 | 2.273 | -0.002 | 98.8% |

Prediction minus observation defines error. The methods share the same geometry mask.
Coverage is the fraction of held-out block means with a supported prediction.
Scores are weighted by block count, not independent geophysical observations.
See `tables/cv_summary.csv` for median per-line RMSE and `cv_by_fold.csv` for fold variation.
Linear interpolation was selected in advance as a transparent baseline, not optimized against these folds.

## What these results establish

The workflow reads a million-record release, resolves practical metadata conflicts,
checks source grid/channel consistency, exports properly oriented GeoTIFFs, and
evaluates interpolation with whole planned lines withheld. The 250 m supplier
grid spacing does not imply 250 m geophysical resolving power. The source is
already filtered and levelled. GRVBGL100 and the final 3.5 km products have
different filtering histories and must not be treated as interchangeable.

## Limits on interpretation

Validation is conditional on the supplier's prior corrections and levelling,
which used survey-wide information. It is not independent external accuracy,
and close infill lines can make prediction easier. The small grid-versus-line
residuals are consistency checks of related products, not accuracy estimates.
The Gaussian residual depends on the chosen scale; edge cells are excluded.
High/low gravity regions describe relative mass-distribution signals but do not
uniquely identify lithology, ore bodies or depth. The gradient is only an
exploratory map-texture diagnostic, not a navigation performance result.
