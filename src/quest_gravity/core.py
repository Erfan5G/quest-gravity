"""Small, explicit numerical functions used by the analysis and its tests."""
from pathlib import Path
import struct
import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import gaussian_filter, binary_erosion


def read_quest_grid(path):
    """Read ONLY uncompressed little-endian float32, unrotated QUEST GRD files.

    Header layout checked against GMT.jl readgeosoft implementation; not a
    general Geosoft reader. Returns ascending node-centre x, y and z[y,x].
    Unsupported encodings fail explicitly, rather than silently misplotting.
    """
    with Path(path).open('rb') as f:
        h = f.read(512)
        if len(h) != 512:
            raise ValueError('Truncated Geosoft header')
        es, sf, nx, ny, kx = struct.unpack('<5i', h[:20])
        dx, dy, x0, y0, rotation, base, mult = struct.unpack('<7d', h[20:76])
        if (es, sf, kx) != (4, 2, 1) or rotation != 0 or base != 0 or mult != 1:
            raise ValueError('Unsupported Geosoft encoding, orientation or scaling')
        if not (0 < nx < 10000 and 0 < ny < 10000 and dx > 0 and dy > 0):
            raise ValueError('Invalid grid geometry')
        data = np.fromfile(f, dtype='<f4', count=nx*ny)
        if data.size != nx*ny:
            raise ValueError('Truncated grid payload')
    z = data.reshape(ny, nx).astype(float)
    z[z < -1e30] = np.nan
    return x0 + dx*np.arange(nx), y0 + dy*np.arange(ny), z


def metrics(error):
    e = np.asarray(error, dtype=float)
    e = e[np.isfinite(e)]
    if not len(e):
        return {'n': 0, 'bias_mgal': None, 'mae_mgal': None,
                'rmse_mgal': None, 'p95_abs_mgal': None}
    return {'n': len(e), 'bias_mgal': float(e.mean()),
            'mae_mgal': float(np.abs(e).mean()),
            'rmse_mgal': float(np.sqrt(np.mean(e**2))),
            'p95_abs_mgal': float(np.percentile(np.abs(e), 95))}


class LocalInterpolator:
    """Linear triangulation, IDW and nearest baseline, with common support.

    Mask outside convex hull, where nearest data exceed max_distance, or
    where the longest enclosing triangle edge exceeds max_edge. Distance
    and triangle checks limit bridging of internal coverage gaps.
    """
    def __init__(self, xy, z, max_distance=3000., max_edge=7000.):
        self.xy = np.asarray(xy)
        self.z = np.asarray(z)
        self.tree = cKDTree(self.xy)
        self.tri = Delaunay(self.xy)
        self.linear = LinearNDInterpolator(self.tri, self.z)
        self.max_distance = max_distance
        self.max_edge = max_edge
        corners = self.xy[self.tri.simplices]
        self.edge = np.maximum.reduce([
            np.linalg.norm(corners[:, i]-corners[:, j], axis=1)
            for i, j in [(0, 1), (1, 2), (0, 2)]])

    def predict(self, targets, method='linear'):
        targets = np.asarray(targets)
        d, ix = self.tree.query(targets, k=min(12, len(self.xy)))
        if d.ndim == 1:
            d, ix = d[:, None], ix[:, None]
        simplex = self.tri.find_simplex(targets)
        supported = ((simplex >= 0) & (d[:, 0] <= self.max_distance)
                     & (self.edge[np.maximum(simplex, 0)] <= self.max_edge))
        if method == 'linear':
            pred = np.asarray(self.linear(targets))
        elif method == 'nearest':
            pred = self.z[ix[:, 0]].copy()
        elif method == 'idw':
            w = np.where(d <= self.max_distance, 1/np.maximum(d, 1e-9)**2, 0.)
            total = w.sum(axis=1)
            pred = np.divide((w*self.z[ix]).sum(axis=1), total,
                             out=np.full(len(targets), np.nan), where=total > 0)
            exact = d[:, 0] < 1e-8
            pred[exact] = self.z[ix[exact, 0]]
        else:
            raise ValueError(method)
        pred[~supported] = np.nan
        return pred, d[:, 0]


def regional_residual(z, spacing, sigma_m=20000., min_support=.95):
    """Exploratory Gaussian regional trend, NOT an isostatic correction.

    Normalize by valid-data support; omit residuals with less than 95%
    kernel support. sigma is a standard deviation, not a cutoff wavelength.
    """
    valid = np.isfinite(z)
    sig = sigma_m/spacing
    support = gaussian_filter(valid.astype(float), sig, mode='constant', cval=0)
    numerator = gaussian_filter(np.where(valid, z, 0.), sig, mode='constant', cval=0)
    regional = np.divide(numerator, support, out=np.full_like(z, np.nan), where=support > 0)
    regional[~valid | (support < min_support)] = np.nan
    return regional, z-regional, support


def gradient_mgal_km(z, spacing):
    """Finite-difference horizontal gradient; mask one-cell boundary halo."""
    gy, gx = np.gradient(z, spacing/1000., spacing/1000.)
    g = np.hypot(gx, gy)
    g[~binary_erosion(np.isfinite(z), structure=np.ones((3, 3)))] = np.nan
    return g


def write_geotiff(path, x, y, z, description, units='mGal'):
    import rasterio
    from rasterio.transform import from_origin
    dx, dy = float(x[1]-x[0]), float(y[1]-y[0])
    transform = from_origin(x[0]-dx/2, y[-1]+dy/2, dx, dy)
    # Array has ascending northing. GeoTIFF requires northern row first.
    array = np.where(np.isfinite(z[::-1]), z[::-1], -9999.).astype('float32')
    with rasterio.open(path, 'w', driver='GTiff', height=len(y), width=len(x),
                       count=1, dtype='float32', crs='EPSG:26910',
                       transform=transform, nodata=-9999., compress='deflate') as ds:
        ds.write(array, 1)
        ds.set_band_description(1, description)
        ds.update_tags(AREA_OR_POINT='Point', units=units,
                       source='Geoscience BC Report 2008-08; Sander Geophysics',
                       processing='See README and outputs/metadata.json')
