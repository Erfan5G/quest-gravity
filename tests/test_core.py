"""Tests target consequential scientific/GIS mistakes, not plotting details."""
import unittest
import struct
import tempfile
from pathlib import Path
import numpy as np
import rasterio
from quest_gravity.core import (read_quest_grid, LocalInterpolator, metrics,
                                regional_residual, gradient_mgal_km, write_geotiff)


class ScientificChecks(unittest.TestCase):
    def test_plane_recovered_and_no_extrapolation(self):
        xy=np.array([[0.,0.],[1000,0],[0,1000],[1000,1000]])
        z=2*xy[:,0]/1000-3*xy[:,1]/1000+10
        model=LocalInterpolator(xy,z)
        p,_=model.predict(np.array([[250.,750.],[2000.,2000.]]))
        self.assertAlmostEqual(p[0],8.25)
        self.assertTrue(np.isnan(p[1]))

    def test_large_triangle_mask(self):
        model=LocalInterpolator(np.array([[0.,0.],[10000,0],[0,10000]]),np.array([1.,2.,3.]))
        p,_=model.predict([[100,100]])
        self.assertTrue(np.isnan(p[0]))

    def test_idw_exact_observation(self):
        xy=np.array([[0.,0.],[1000,0],[0,1000]])
        model=LocalInterpolator(xy,[5.,6.,7.])
        p,_=model.predict([[0.,0.]],'idw')
        self.assertAlmostEqual(p[0],5.)

    def test_geotiff_orientation_crs_and_pixel_centres(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'grid.tif'
            x=np.array([400000.,401000.]); y=np.array([5900000.,5901000.,5902000.])
            z=np.array([[1.,2.],[3.,4.],[5.,np.nan]])
            write_geotiff(path,x,y,z,'test')
            with rasterio.open(path) as ds:
                self.assertEqual(ds.crs.to_epsg(),26910)
                np.testing.assert_allclose(ds.xy(0,0),(x[0],y[-1]))
                a=ds.read(1,masked=True)
                self.assertEqual(a[0,0],5.)
                self.assertTrue(a.mask[0,1])
                self.assertEqual(a[-1,0],1.)

    def test_geosoft_binary_and_nodata(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'test.grd'
            h=struct.pack('<5i7d',4,2,2,2,1,250.,250.,400000.,5900000.,0.,0.,1.)
            p.write_bytes(h+bytes(512-len(h))+np.array([1.,2.,3.,-1e32],dtype='<f4').tobytes())
            x,y,z=read_quest_grid(p)
            self.assertEqual(x[1],400250.)
            self.assertEqual(y[1],5900250.)
            self.assertEqual(z[0,1],2.)
            self.assertTrue(np.isnan(z[1,1]))
            p.write_bytes(p.read_bytes()[:520])
            with self.assertRaises(ValueError): read_quest_grid(p)

    def test_gaussian_constant_no_false_residual(self):
        z=np.full((101,101),37.)
        regional,residual,support=regional_residual(z,1000.,5000.)
        self.assertTrue(np.isnan(residual[0,0]))
        self.assertLess(np.nanmax(np.abs(residual)),1e-10)

    def test_gradient_units_and_mask(self):
        y,x=np.mgrid[:10,:10]
        z=2.*x+3.*y
        g=gradient_mgal_km(z,1000.)
        self.assertAlmostEqual(g[4,4],np.sqrt(13.))
        self.assertTrue(np.isnan(g[0,0]))

    def test_metric_sign_and_rmse(self):
        m=metrics([1.,-3.,np.nan])
        self.assertEqual(m['n'],2)
        self.assertEqual(m['bias_mgal'],-1.)
        self.assertAlmostEqual(m['rmse_mgal'],np.sqrt(5.))


if __name__=='__main__': unittest.main()
