"""Run from the repository root: python -m quest_gravity.pipeline."""
import argparse
import hashlib
import json
import logging
import platform
from datetime import datetime, timezone
from pathlib import Path
import importlib.metadata

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

from .core import (read_quest_grid, metrics, LocalInterpolator,
                   regional_residual, gradient_mgal_km, write_geotiff)
from .download import sha256

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'data/raw/release/GBCReport2008-8_Gravity_Data'
OUT = ROOT/'outputs'
CHANNELS = ('SECOND YEAR DAY FLT LINE LATWGS84 LONGWGS84 UTMXNAD83 UTMYNAD83 '
            'GPSZNAD83 MSLZ RALT FZ AZ GRVRAW GRVLAT GRVEOT GRVFRA GRVFRA100 '
            'GRVFRAL100 BTERRAIN SBGCOR GRVBGL100 ISOCORR ISOBGL3.5 '
            'GRVFRAL3.5 GRVBGL3.5').split()


def load_lines():
    """Preserve headers/line markers as an audit, parse exactly 27 columns."""
    source = SOURCE/'Grav.xyz'
    numeric = ROOT/'data/raw/numeric.tmp'
    markers = []
    count = 0
    with source.open(encoding='ascii') as f, numeric.open('w') as out:
        for n, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith('/'):
                continue
            fields = s.split()
            if fields[0] in ('Line', 'Tie_line'):
                markers.append({'type': fields[0], 'line': fields[1], 'source_row': n})
                continue
            if len(fields) != len(CHANNELS):
                raise ValueError(f'Malformed source row {n}: {len(fields)} columns')
            out.write(s+'\n')
            count += 1
    try:
        df = pd.read_csv(numeric, sep=r'\s+', names=CHANNELS, na_values=['*', '**', '***'])
    finally:
        numeric.unlink(missing_ok=True)
    if len(df) != count:
        raise ValueError('Parser row-count mismatch')
    pd.DataFrame(markers).to_csv(OUT/'tables/line_markers.csv', index=False)
    return df


def block_by_line(df, size):
    """Means within projected cells AND exact line section, before CV split.

    Averaging reduces oversampling along tracks. No across-line mixing here.
    Coordinates are means of observations, not arbitrary block corners.
    """
    work = df[['LINE', 'UTMXNAD83', 'UTMYNAD83', 'GRVBGL100', 'GRVFRAL100']].copy()
    work['bx'] = np.floor(work.UTMXNAD83/size).astype(int)
    work['by'] = np.floor(work.UTMYNAD83/size).astype(int)
    b = work.groupby(['LINE', 'bx', 'by'], as_index=False).agg(
        easting_m=('UTMXNAD83', 'mean'), northing_m=('UTMYNAD83', 'mean'),
        bouguer_mgal=('GRVBGL100', 'mean'), free_air_mgal=('GRVFRAL100', 'mean'),
        count=('GRVBGL100', 'size'))
    b['planned_line'] = np.floor(b.LINE).astype(int)
    return b


def fit_interpolator(df, cfg):
    # Shared endpoints can occur across line sections; aggregate exact XY only.
    unique = df.groupby(['easting_m', 'northing_m'], as_index=False).bouguer_mgal.mean()
    return LocalInterpolator(unique[['easting_m', 'northing_m']].values,
                             unique.bouguer_mgal.values,
                             cfg['max_nearest_distance_m'], cfg['max_triangle_edge_m'])


def validate_lines(blocks, cfg):
    traverse = blocks[blocks.planned_line >= 1000].copy()
    ids = np.sort(traverse.planned_line.unique())
    mapping = {int(line): int(i % cfg['cv_folds']) for i, line in enumerate(ids)}
    traverse['fold'] = traverse.planned_line.map(mapping)
    traverse[['LINE', 'planned_line', 'fold']].drop_duplicates().to_csv(
        OUT/'tables/cv_line_assignment.csv', index=False)
    records, all_predictions = [], []
    for fold in range(cfg['cv_folds']):
        train, test = traverse[traverse.fold != fold], traverse[traverse.fold == fold]
        assert set(train.planned_line).isdisjoint(set(test.planned_line))
        model = fit_interpolator(train, cfg)
        xy = test[['easting_m', 'northing_m']].values
        observed = test.bouguer_mgal.values
        for method in ['linear', 'idw', 'nearest']:
            pred, dist = model.predict(xy, method)
            error = pred-observed
            records.append({'fold': fold, 'method': method, 'test_blocks': len(test),
                            'coverage_fraction': float(np.isfinite(pred).mean()), **metrics(error)})
            all_predictions.append(pd.DataFrame({
                'fold': fold, 'method': method, 'planned_line': test.planned_line.values,
                'easting_m': xy[:, 0], 'northing_m': xy[:, 1],
                'observed_mgal': observed, 'predicted_mgal': pred,
                'error_mgal': error, 'nearest_train_m': dist}))
        logging.info('Completed held-line fold %d/%d', fold+1, cfg['cv_folds'])
    predictions = pd.concat(all_predictions, ignore_index=True)
    scores = pd.DataFrame(records)
    scores.to_csv(OUT/'tables/cv_by_fold.csv', index=False)
    predictions.to_csv(OUT/'tables/cv_predictions.csv.gz', index=False, compression='gzip')
    summary = []
    for method, group in predictions.groupby('method'):
        valid = group.dropna(subset=['error_mgal'])
        line_rmse = valid.groupby('planned_line').error_mgal.apply(lambda a: np.sqrt(np.mean(a*a)))
        summary.append({'method': method, **metrics(valid.error_mgal),
                        'coverage_fraction': len(valid)/len(group),
                        'median_line_rmse_mgal': float(line_rmse.median()),
                        'lines_scored': len(line_rmse)})
    summary = pd.DataFrame(summary).sort_values('rmse_mgal')
    summary.to_csv(OUT/'tables/cv_summary.csv', index=False)
    return summary, predictions


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 12, 'axes.titleweight': 'bold',
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white'})


def map_panel(fig, ax, x, y, z, title, units='mGal', diverge=False, cmap=None):
    finite = z[np.isfinite(z)]
    if not len(finite):
        ax.text(.5, .5, 'No supported cells', ha='center', transform=ax.transAxes)
        return
    options = {}
    if diverge:
        limit = float(np.max(np.abs(finite)))
        options['norm'] = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    dx, dy = x[1]-x[0], y[1]-y[0]
    im = ax.imshow(z, origin='lower', extent=[(x[0]-dx/2)/1000, (x[-1]+dx/2)/1000,
                  (y[0]-dy/2)/1000, (y[-1]+dy/2)/1000],
                  cmap=cmap or ('RdBu_r' if diverge else 'viridis'),
                  interpolation='nearest', **options)
    ax.set(title=title, xlabel='Easting (km)', ylabel='Northing (km)', facecolor='#edf1f4')
    ax.set_aspect('equal')
    ax.tick_params(labelsize=9)
    fig.colorbar(im, ax=ax, shrink=.72, pad=.025, label=units)


def finish(fig, filename, title, subtitle):
    fig.suptitle(title, x=.06, ha='left', fontsize=20, fontweight='bold', color='#153449')
    fig.text(.06, .90, subtitle, fontsize=10, color='#48606c')
    fig.text(.06, .025, 'Source: Geoscience BC Report 2008-08 / Sander Geophysics  |  NAD83 / UTM 10N (EPSG:26910)',
             fontsize=9, color='#48606c')
    fig.subplots_adjust(top=.82, bottom=.12, left=.06, right=.98, wspace=.32)
    fig.savefig(OUT/'figures'/filename, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default=str(ROOT/'config.json'))
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    if cfg['crs'] != 'EPSG:26910' or cfg['cv_channel'] != 'GRVBGL100':
        raise ValueError('This audited workflow supports only EPSG:26910 / GRVBGL100')
    for name in ['figures', 'grids', 'tables']:
        (OUT/name).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(OUT/'processing.log', mode='w'), logging.StreamHandler()])
    style()
    expected = json.loads((ROOT/'docs/source_checksums.json').read_text())
    for name, checksum in expected.items():
        if sha256(ROOT/'data/raw'/name) != checksum:
            raise ValueError(f'Archive checksum mismatch: {name}')
    logging.info('Parsing all released line observations')
    df = load_lines()
    qc = {'rows_parsed': len(df), 'line_sections': int(df.LINE.nunique()),
          'missing_by_channel': {k: int(v) for k, v in df.isna().sum().items()},
          'duplicate_time_line_records': int(df.duplicated(['YEAR','DAY','FLT','LINE','SECOND']).sum())}
    required = ['UTMXNAD83','UTMYNAD83','LATWGS84','LONGWGS84','GRVBGL100','GRVFRAL100','LINE']
    good = np.isfinite(df[required]).all(axis=1)
    qc['rows_excluded_nonfinite_required'] = int((~good).sum())
    df = df[good].copy()
    # No automatic sigma clipping: geological extrema are not presumed errors.
    sample = df.iloc[::cfg['consistency_stride']]
    trans = Transformer.from_crs('EPSG:4326', cfg['crs'], always_xy=True)
    xx, yy = trans.transform(sample.LONGWGS84.values, sample.LATWGS84.values)
    position_error = np.hypot(xx-sample.UTMXNAD83, yy-sample.UTMYNAD83)
    qc['crs_check'] = {'n': len(sample), 'median_difference_m': float(np.median(position_error)),
                       'max_difference_m': float(np.max(position_error)),
                       'note': 'Modern EPSG transformation versus supplied WGS84/NAD83 realization; not survey positional accuracy.'}
    if np.max(position_error) > 20:
        raise ValueError('Geographic/projected coordinate mismatch requires review')
    df[CHANNELS].describe().to_csv(OUT/'tables/channel_statistics.csv')
    qc['isostatic_identity'] = metrics(df['ISOBGL3.5']-(df['GRVBGL3.5']-df.ISOCORR))
    qc['bouguer100_minus_freeair100_minus_slab'] = metrics(df.GRVBGL100-(df.GRVFRAL100-df.SBGCOR))
    qc['utc_seconds_above_86400'] = int((df.SECOND >= 86400).sum())
    qc['ranges'] = {c: [float(df[c].min()),float(df[c].max())] for c in
                    ['LATWGS84','LONGWGS84','UTMXNAD83','UTMYNAD83','GRVBGL100']}
    grids, consistency = {}, []
    for file, channel, label in [
        ('GRVBG3.5KM.grd', 'GRVBGL3.5', 'Supplier Bouguer anomaly'),
        ('GRVFA3.5KM.grd', 'GRVFRAL3.5', 'Supplier free-air gravity'),
        ('ISOBG3.5KM.grd', 'ISOBGL3.5', 'Supplier isostatic Bouguer residual')]:
        x, y, z = read_quest_grid(SOURCE/file)
        grids[channel] = (x, y, z)
        write_geotiff(OUT/'grids'/file.replace('.grd','.tif'), x, y, z, label)
        sampled = RegularGridInterpolator((y,x), z, bounds_error=False)(
            sample[['UTMYNAD83','UTMXNAD83']].values)
        score = {'file': file, 'line_channel': channel, **metrics(sampled-sample[channel].values)}
        consistency.append(score)
        if score['n'] < .9*len(sample) or score['rmse_mgal'] > .1:
            raise ValueError(f'Grid/channel consistency failed: {file}')
    pd.DataFrame(consistency).to_csv(OUT/'tables/supplier_grid_consistency.csv', index=False)
    logging.info('Supplier grid orientation and channel consistency verified')
    x, y, bouguer = grids['GRVBGL3.5']
    qc['supplier_bouguer_range_mgal'] = [float(np.nanmin(bouguer)), float(np.nanmax(bouguer))]
    qc['supplier_valid_grid_nodes'] = int(np.isfinite(bouguer).sum())
    blocks = block_by_line(df, cfg['block_size_m'])
    blocks.to_csv(OUT/'tables/line_blocks.csv.gz', index=False, compression='gzip')
    qc['blocks'] = len(blocks)
    qc['planned_traverse_lines'] = int(blocks.loc[blocks.planned_line>=1000,'planned_line'].nunique())
    summary, predictions = validate_lines(blocks, cfg)
    # Fix the method in advance; report alternatives rather than tuning to test scores.
    model = fit_interpolator(blocks[blocks.planned_line >= 1000], cfg)
    step = cfg['output_grid_m']
    gx = np.arange(np.floor(df.UTMXNAD83.min()/step)*step, df.UTMXNAD83.max()+step, step)
    gy = np.arange(np.floor(df.UTMYNAD83.min()/step)*step, df.UTMYNAD83.max()+step, step)
    xx, yy = np.meshgrid(gx, gy)
    prediction, distance = model.predict(np.c_[xx.ravel(), yy.ravel()])
    own_grid = prediction.reshape(xx.shape)
    distance = distance.reshape(xx.shape)
    distance[~np.isfinite(own_grid)] = np.nan
    write_geotiff(OUT/'grids/bouguer100_linear_1000m.tif', gx, gy, own_grid,
                  'Linear interpolation of supplier GRVBGL100; no new physical reductions')
    write_geotiff(OUT/'grids/nearest_support_distance_m.tif', gx, gy, distance,
                  'Distance to nearest block mean, not uncertainty', 'm')
    # Work at 1 km node spacing on the already low-pass-filtered supplier grid.
    sx, sy, sz = x[::4], y[::4], bouguer[::4,::4]
    regional, residual, support = regional_residual(sz, 1000, cfg['regional_sigma_m'], cfg['regional_min_support'])
    gradient = gradient_mgal_km(sz, 1000)
    write_geotiff(OUT/'grids/regional_gaussian.tif', sx, sy, regional, 'Exploratory Gaussian regional component; sigma in metadata')
    write_geotiff(OUT/'grids/residual_gaussian.tif', sx, sy, residual, 'Exploratory Gaussian residual, not isostatic anomaly; sigma in metadata')
    write_geotiff(OUT/'grids/horizontal_gradient_mgal_km.tif', sx, sy, gradient, 'Horizontal gradient magnitude', 'mGal/km')
    # 30 km has no 95%-supported cells in this narrow survey. Do not relax
    # the edge criterion to force a result; compare supported 10/15/20 km.
    residuals = {s: regional_residual(sz,1000,s*1000,cfg['regional_min_support'])[1] for s in [10,15,20]}
    common = np.logical_and.reduce([np.isfinite(r) for r in residuals.values()])
    if not common.any():
        raise ValueError('No common cells satisfy regional-filter support criteria')
    sensitivity = []
    for sigma, r in residuals.items():
        sensitivity.append({'sigma_km':sigma, 'common_cells':int(common.sum()),
                            'std_mgal':float(np.std(r[common])),
                            'min_mgal':float(np.min(r[common])), 'max_mgal':float(np.max(r[common]))})
    pd.DataFrame(sensitivity).to_csv(OUT/'tables/regional_scale_sensitivity.csv',index=False)

    fig, axes = plt.subplots(1,3,figsize=(17,8))
    for ax, channel, title, div in zip(axes, grids,
            ['Bouguer anomaly','Free-air gravity','Isostatic Bouguer residual'],[False,True,True]):
        a,b,c = grids[channel]
        map_panel(fig,ax,a,b,c,title,diverge=div)
    finish(fig,'01_gravity_atlas.png','QUEST | Regional airborne gravity',
           'Published fields converted in Python  |  250 m grid nodes; supplier 3.5 km half-wavelength filtering  |  Values in mGal')

    fig, axes = plt.subplots(1,3,figsize=(17,8))
    map_panel(fig,axes[0],sx,sy,sz,'Published Bouguer field')
    map_panel(fig,axes[1],sx,sy,residual,f"Exploratory residual\nGaussian sigma = {cfg['regional_sigma_m']/1000:g} km",diverge=True)
    map_panel(fig,axes[2],sx,sy,gradient,'Horizontal gradient magnitude','mGal/km',cmap='magma')
    finish(fig,'02_spatial_analysis.png','QUEST | Wavelengths and spatial structure',
           'Residual requires 95% kernel support  |  Gradient is a map-texture indicator, not a positioning-accuracy estimate')

    fig, axes = plt.subplots(1,3,figsize=(17,8))
    for kind,color,label in [(False,'#267e91','Traverse / infill'),(True,'#da7b25','Control')]:
        sub = blocks[(blocks.planned_line<1000)==kind]
        axes[0].scatter(sub.easting_m/1000,sub.northing_m/1000,s=.4,c=color,label=label,rasterized=True)
    axes[0].set(title='Survey coverage',xlabel='Easting (km)',ylabel='Northing (km)',aspect='equal')
    axes[0].legend(markerscale=8,fontsize=9,loc='upper right')
    map_panel(fig,axes[1],gx,gy,own_grid,'Our line-data interpolation\nGRVBGL100 / linear / 1 km nodes')
    map_panel(fig,axes[2],gx,gy,distance/1000,'Distance to nearest data block','km',cmap='cividis')
    finish(fig,'03_line_processing.png','QUEST | From line observations to GIS rasters',
           '500 m per-line block means  |  Outside-hull, long-triangle and distance masking  |  No extrapolation')

    fig, axes = plt.subplots(1,3,figsize=(17,6))
    axes[0].bar(summary.method, summary.rmse_mgal,color=['#247c8c','#559fa7','#b5d4d6'])
    axes[0].set(ylabel='Held-line RMSE (mGal)',title='Interpolation comparison')
    p = predictions[(predictions.method=='linear') & predictions.error_mgal.notna()]
    h = axes[1].hexbin(p.observed_mgal,p.predicted_mgal,gridsize=55,mincnt=1,cmap='viridis',bins='log')
    lo,hi=p.observed_mgal.min(),p.observed_mgal.max()
    axes[1].plot([lo,hi],[lo,hi],color='#d35c50',lw=1)
    axes[1].set(xlabel='Held-out observation (mGal)',ylabel='Prediction (mGal)',title='Linear: held-line predictions')
    fig.colorbar(h,ax=axes[1],label='Block count',shrink=.8)
    axes[2].hist(p.error_mgal,bins=75,color='#247c8c')
    axes[2].axvline(0,color='#d35c50',lw=1)
    axes[2].set(xlabel='Prediction - observation (mGal)',ylabel='Block count',title='Linear: error distribution')
    finish(fig,'04_validation.png','QUEST | Whole-line interpolation checks',
           'Five folds; all sections/reflights of a planned line stay together  |  Control lines excluded  |  Conditional on supplier processing')

    # Representative geographically central traverse profile, deterministic selection.
    ids = np.sort(blocks.loc[(blocks.planned_line>=1000)&(blocks.planned_line<2000),'planned_line'].unique())
    line = int(ids[len(ids)//2])
    profile = blocks[blocks.planned_line==line].sort_values('easting_m')
    fig,ax=plt.subplots(figsize=(12,5))
    ax.plot(profile.easting_m/1000,profile.bouguer_mgal,color='#247c8c',label='GRVBGL100 block means')
    ax.set(xlabel='Easting (km), NAD83 / UTM 10N',ylabel='Bouguer channel (mGal)',title=f'Planned traverse {line}')
    ax.grid(alpha=.2);ax.legend()
    finish(fig,'05_profile.png','QUEST | A representative 2D profile',
           'Supplier line-filtered Bouguer channel  |  Horizontal axis is easting, not cumulative flight distance')
    profile.to_csv(OUT/'tables/representative_profile.csv',index=False)

    metadata = {'created_utc':datetime.now(timezone.utc).isoformat(),
                'project':'Original QUEST airborne gravity, Geoscience BC Report 2008-08',
                'source_url':'https://www.geosciencebc.com/major-projects/quest/',
                'source_checksums':expected,'config':cfg,'qc':qc,
                'software':{'python':platform.python_version(), **{p:importlib.metadata.version(p) for p in
                             ['numpy','pandas','scipy','matplotlib','rasterio','pyproj']}},
                'validation_scope':'Interpolation consistency of supplier-corrected/levelled data; not independent instrument accuracy.',
                'not_performed':['new terrain correction','normal gravity reduction','remove-compute-restore',
                                 'independent terrestrial validation','full navigation simulation','magnetic processing']}
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2))
    summary_lines = '\n'.join(f"| {r.method} | {r.rmse_mgal:.3f} | {r.mae_mgal:.3f} | {r.bias_mgal:.3f} | {r.coverage_fraction:.1%} |" for r in summary.itertuples())
    (OUT/'RESULTS.md').write_text(f'''# QUEST gravity: verified run results

Run: {metadata['created_utc']}. All figures use actual supplier data.

- Parsed observations: {qc['rows_parsed']:,}; line sections: {qc['line_sections']}.
- Per-line block means: {len(blocks):,}; planned traverse/infill lines: {qc['planned_traverse_lines']}.
- Nonfinite required rows excluded: {qc['rows_excluded_nonfinite_required']}.
- Duplicate flight/line/time records: {qc['duplicate_time_line_records']}.
- Published Bouguer range: {qc['supplier_bouguer_range_mgal'][0]:.3f} to {qc['supplier_bouguer_range_mgal'][1]:.3f} mGal.
- Geographic/projected coordinate discrepancy: median {qc['crs_check']['median_difference_m']:.3f} m; max {qc['crs_check']['max_difference_m']:.3f} m.

## Held-line comparison

| Method | RMSE (mGal) | MAE (mGal) | Bias (mGal) | Supported blocks |
|---|---:|---:|---:|---:|
{summary_lines}

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

Read `docs/METHODS.md` and `docs/METADATA_ISSUES.md` before presenting the work.
''')
    logging.info('Completed pipeline; results written to %s', OUT)


if __name__ == '__main__':
    main()
