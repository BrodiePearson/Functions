#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:22:28 2018

@author: JennaPalmer
"""

# Documentation: http://unidata.github.io/netcdf4-python/http://unidata.github.io/netcdf4-python/
#
# INPUTS:
#   path = path to data with quotes and end slash - Ex: '/users/user/Desktop/'
#   filename = name of dataset with quotes and .nc - Ex: 'testdata.nc'
#   datadict = dictionary containing all key-value pairs

def data2nc(path,filename,datadict):
    
    import xarray as xr
    
    outputdata = path + filename
 
    xrdict = xr.Dataset(dict)
    xrdict.to_netcdf(outputdata)
    