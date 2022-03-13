# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:44:23 2021

@author: Boorn
"""

import xarray as xray
import os
import atmos.data as dat
import cmaps
import numpy as np
from atmos.constants import const as constants
from atmos.data import get_coord
import atmos.xrhelper as xr
import pandas as pd

#%%
def divergence_spherical_2d(Fx, Fy, lat=None, lon=None):
    """Return the 2-D spherical divergence.
    Parameters
    ----------
    Fx, Fy : ndarrays or xray.DataArrays
        Longitude, latitude components of a vector function in
        spherical coordinates.  Latitude and longitude should be the
        second-last and last dimensions, respectively, of Fx and Fy.
        Maximum size 5-D.
    lat, lon : ndarrays, optional
        Longitude and latitude in degrees.  If these are omitted, then
        Fx and Fy must be xray.DataArrays with latitude and longitude
        in degrees within the coordinates.
    Returns
    -------
    d, d1, d2 : ndarrays or xray.DataArrays
        d1 = dFx/dx, d2 = dFy/dy, and d = d1 + d2.
        d: Vertically integrated moisture flux divergence
        d1: Vertically integrated moisture flux divergence - x
        d2: Vertically integrated moisture flux divergence - y
    Reference
    ---------
    Atmospheric and Oceanic Fluid Dynamics: Fundamentals and
    Large-Scale Circulation, by Geoffrey K. Vallis, Cambridge
    University Press, 2006 -- Equation 2.30.
    """
    
    dims0=Fx.dims[-2:]
    nmax = 5
    ndim = Fx.ndim
    if ndim > nmax:
        raise ValueError('Input data has too many dimensions. Max 5-D.')

    if isinstance(Fx, xray.DataArray):
        i_DataArray = True
        name, attrs, coords, _ = xr.meta(Fx)
        if lat is None:
            lat = get_coord(Fx, 'lat')
        if lon is None:
            lon = get_coord(Fx, 'lon')
    else:
        i_DataArray = False
        if lat is None or lon is None:
            raise ValueError('Lat/lon inputs must be provided when input '
                'data is an ndarray.')

    R = constants.radius_earth.values
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)

    # Add singleton dimensions for looping, if necessary
    for i in range(ndim, nmax):
        Fx = np.expand_dims(Fx, axis=0)
        Fy = np.expand_dims(Fy, axis=0)

    dims = Fx.shape
    nlon = dims[-1]
    nlat = dims[-2]

    d1 = np.zeros(dims, dtype=float)
    d2 = np.zeros(dims, dtype=float)

    for i in range(nlat):
        dx = np.cumsum(np.gradient(lon_rad))
        coslat = np.cos(lat_rad[i])
        for k1 in range(dims[0]):
            for k2 in range(dims[1]):
                for k3 in range(dims[2]):
                    sub = Fx[k1,k2,k3,i,:]
                    d1[k1,k2,k3,i,:] = np.gradient(sub, dx) / (R*coslat)
    for j in range(nlon):
        dy = np.cumsum(np.gradient(lat_rad))
        coslat = np.cos(lat_rad)
        # Set to NaN at poles to keep from blowing up
        coslat[abs(lat) > 89] = np.nan
        for k1 in range(dims[0]):
            for k2 in range(dims[1]):
                for k3 in range(dims[2]):
                    sub = Fy[k1,k2,k3,:,j] * coslat
                    d2[k1,k2,k3,:,j] = np.gradient(sub, dy) / (R*coslat)

    # Collapse any additional dimensions that were added
    for i in range(ndim, d1.ndim):
        d1, d2 = d1[0], d2[0]

    d = d1 + d2

    if i_DataArray:
        d = xray.DataArray(d, coords=coords, dims=dims0)
        d1 = xray.DataArray(d1, coords=coords, dims=dims0)
        d2 = xray.DataArray(d2, coords=coords, dims=dims0)
        d.attrs['units'] = 'kg/m2/s'
        d1.attrs['units'] = 'kg/m2/s'
        d2.attrs['units'] = 'kg/m2/s'
        
    return d, d1, d2

def moisture_flux_conv(uq, vq, lat=None, lon=None, plev=None, pdim=-3,
                       pmin=0, pmax=1e6, return_comp=False, already_int=False):
    """Return the vertically integrated moisture flux convergence.
    Parameters
    ----------
    uq : ndarray or xray.DataArray
        Zonal moisture flux, with latitude as the second-last dimension,
        longitude as the last dimension.
        i.e. zonal wind (m/s) * specific humidity (kg/kg)
        If already_int is True, then uq is already vertically integrated.
        Otherwise, uq is on vertical pressure levels.
    vq : ndarray or xray.DataArray
        Meridional moisture flux, with latitude as the second-last dimension,
        longitude as the last dimension.
        i.e. meridional wind (m/s) * specific humidity (kg/kg)
        If already_int is True, then vq is already vertically integrated.
        Otherwise, vq is on vertical pressure levels.
    lat, lon : ndarray, optional
        Latitudes and longitudes in degrees.  If omitted, then uq and
        vq must be xray.DataArrays and the coordinates are extracted
        from them.
    plev : ndarray, optional
        Pressure levels in Pascals.  If omitted, then extracted from
        DataArray inputs.
    pdim : int, optional
        Dimension of pressure levels in uq and vq.
    pmin, pmax : float, optional
        Lower and upper bounds (inclusive) of pressure levels (Pa)
        to include in integration.
    return_comp : bool, optional
        If True, return additional components, otherwise just total
        moisture flux convergence.
    already_int : bool, optional
        If True, then uq and vq inputs have already been vertically
        integrated.  Otherwise, the vertical integration is
        calculated here.
    Returns
    -------
    If return_comp is False:
    mfc : ndarray or xray.DataArray
        Vertically integrated moisture flux convergence in mm/day.
    If return_comp is True:
    mfc, mfc_x, mfc_y, uq_int, vq_int : ndarrays or xray.DataArrays
        Vertically integrated moisture flux convergence in mm/day
        (total, x- and y- components) and vertically integrated
        moisture fluxes.
    """

    if already_int:
        uq_int, vq_int = uq, vq
    else:
        uq_int = dat.int_pres(uq, plev, pdim=pdim, pmin=pmin, pmax=pmax)
        vq_int = dat.int_pres(vq, plev, pdim=pdim, pmin=pmin, pmax=pmax)

    mfc, mfc_x, mfc_y = divergence_spherical_2d(uq_int, vq_int, lat, lon)

    # Convert from divergence to convergence, and to mm/day
    mfc = -mfc * 1e5
    mfc_x = -mfc_x * 1e5
    mfc_y = -mfc_y * 1e5

    if isinstance(mfc, xray.DataArray):
        mfc.name = 'Vertically integrated moisture flux convergence'
        mfc.attrs['units'] = '10-5 kg/m2/s'
        
    print(mfc.name)
    
    if return_comp:
        return mfc, mfc_x, mfc_y, uq_int, vq_int
    else:
        return mfc

def calc_vimfc(hus_name, ua_name, va_name, out_name):
    with xray.open_dataset(hus_name) as ds:
        q = ds.hus
    with xray.open_dataset(ua_name) as ds:
        u = ds.ua
    with xray.open_dataset(va_name) as ds:
        v = ds.va
    
    mfc = xray.concat([
        moisture_flux_conv(u.sel(time=time0)*q.sel(time=time0),
                           v.sel(time=time0)*q.sel(time=time0)) 
        for time0 in q.time.values], 
        pd.Index(q.time.values, name="time"))
    mfc = xr.Dataset({'mfc':mfc})
    mfc.to_netcdf(out_name)
    
# calculation
hus_name = 'hus_Amon_historical_1955-2004.nc'
ua_name = os.path.join(dir2, 'ua_Amon_historical_1955-2004.nc')
va_name = os.path.join(dir2, 'va_Amon_historical_1955-2004.nc')
out_name = va_name.replace('va','vimfc')