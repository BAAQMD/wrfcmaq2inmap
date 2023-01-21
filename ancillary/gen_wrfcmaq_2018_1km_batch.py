#!/usr/bin/env python
# Preprocesses WRF/MCIP/CMAQ output for use in InMap
# <beidler.james@epa.gov>

from __future__ import division
from __future__ import print_function
from builtins import range
from builtins import object
import numpy as np
import netCDF4 as ncf
import sys
from fauxioapi import Grid
import pandas as pd

########

class GridBounds(object):
    '''
    Set the boundaries of the new grid and the old grid based on the new grid
    '''
    def __init__(self, in_grid, out_grid):
        self._set_bounds(in_grid, out_grid)

    def _set_bounds(self, in_grid, out_grid):
        # Calculate the raw difference for the distance between the origination points on the new and old grids
        #print(in_grid.XORIG)
        x_dist = (out_grid.XORIG - in_grid.XORIG)
        y_dist = (out_grid.YORIG - in_grid.YORIG)
        # Calculate the starting points for copying on both grids
        # Keep in mind that the grid starts at the SW corner
        if x_dist < 0: # If the old x origination is inside the new grid
            self.icol_o = 0 # Start at the first old grid cell 
            self.ocol_o = int(abs(x_dist/out_grid.XCELL)) # Offset some new out_grid cells
        else: # If the old x origination is outside of equal to the orig on new grid
            self.icol_o = int(abs(x_dist/in_grid.XCELL))
            self.ocol_o = 0
        if y_dist < 0:
            self.irow_o = 0
            self.out_row_orig = int(abs(y_dist/out_grid.YCELL))
        else:
            self.irow_o = int(abs(y_dist/in_grid.XCELL))
            self.out_row_orig = 0
        # Calculate the end points to form boundary fields
        in_x_end = in_grid.XORIG + (in_grid.NCOLS * in_grid.XCELL)
        in_y_end = in_grid.YORIG + (in_grid.NROWS * in_grid.XCELL)
        out_x_end = out_grid.XORIG + (out_grid.NCOLS * out_grid.XCELL)
        out_y_end = out_grid.YORIG + (out_grid.NROWS * out_grid.YCELL)
        x_dist = out_x_end - in_x_end
        y_dist = out_y_end - in_y_end
        if x_dist > 0: # If the endpoint of the old is inside the new
            self.icol_e = in_grid.NCOLS
            self.ocol_e = int(out_grid.NCOLS - abs(x_dist/out_grid.cell)) # Calculate the last out_grid cell to aggregate
            if self.ocol_e > int(self.ocol_e): # Check to see if there is anything after the decimal ie. a partial column calculated.
                self.ocol_e = int(self.ocol_e) + 1  # If there is a partial column then add a column to be processed to hold the extra input data
        else:
            self.icol_e = int(in_grid.NCOLS - abs(x_dist/in_grid.XCELL))
            self.ocol_e = out_grid.NCOLS
        if y_dist > 0:
            self.irow_e = in_grid.NROWS
            self.orow_e = int(out_grid.NROWS - abs(y_dist/out_grid.YCELL))
            if self.orow_e > int(self.orow_e):
                self.orow_e = int(self.orow_e) + 1
        else:
            self.irow_e = int(in_grid.NROWS - abs(y_dist/in_grid.XCELL))
            self.orow_e = out_grid.NROWS
        self.chunk_dim = float(out_grid.XCELL)/float(in_grid.XCELL)  # How many times larger is the output cell versus the input

# Variables to read from the WRF outputs
metvars = {
'U': ('Time','bottom_top','south_north','west_east_stag'),
'V': ('Time','bottom_top','south_north_stag','west_east'),
'W': ('Time','bottom_top_stag','south_north','west_east'),
'PH': ('Time','bottom_top_stag','south_north','west_east'),
'PHB': ('Time','bottom_top_stag','south_north','west_east'),
'T': ('Time','bottom_top','south_north','west_east'),
'P': ('Time','bottom_top','south_north','west_east'),
'PB': ('Time','bottom_top','south_north','west_east'),
'QRAIN': ('Time','bottom_top','south_north','west_east'),
'QCLOUD': ('Time','bottom_top','south_north','west_east'),
'CLDFRA': ('Time','bottom_top','south_north','west_east'),
'GLW': ('Time','south_north','west_east'),
'SWDOWN': ('Time','south_north','west_east'),
'HFX': ('Time','south_north','west_east'),
'UST': ('Time','south_north','west_east'),
'PBLH': ('Time','south_north','west_east'),
'LU_INDEX': ('Time','south_north','west_east')}

class NCF(ncf.Dataset):
    """
    NCF subclass
    """
    def __init__(self, file_name, mode='r'):
        #print('Opening %s' %file_name)
        ncf.Dataset.__init__(self, file_name, mode, format='NETCDF3_64BIT')

    def _set_dims(self, in_ncf, out_grid, layers):
        '''
        Set up the outfile dimensions and attributes based on the new grid
        and the input file
        '''
        out_dims = {'Time': 24, 'bottom_top': layers, 'bottom_top_stag': layers+1,
          'south_north': out_grid.NROWS, 'south_north_stag': out_grid.NROWS+1, 
          'west_east': out_grid.NCOLS, 'west_east_stag': out_grid.NCOLS+1} 
        for dim, value in out_dims.items(): 
            self.createDimension(dim, int(value))
        for att_name in in_ncf.ncattrs():
            att_val = getattr(in_ncf, att_name)
            setattr(self, att_name, att_val) 
        self.LAYERS = layers

    def find_date(self, times, rundate):
        dates = []
        for dt in times:
            dt = ''.join([c.decode() for c in dt])
            dates.append(dt[:4]+dt[5:7]+dt[8:10]+dt[11:13])
        try:
            start_time = dates.index(str(rundate)+'00')
        except ValueError:
            raise ValueError('Could not find %s in WRF' %rundate)
        else:
            return slice(start_time, start_time+24)

    def regrid(self, in_ncf, in_grid, out_grid, rundate, layers_fn, layers=28):
        """
        the main regridding section
        sets loop over variables and decides how to regrid
        """
        time_slice = self.find_date(in_ncf.variables['Times'], rundate) 
        self._set_dims(in_ncf, out_grid, layers)
        bounds = GridBounds(in_grid, out_grid)
        print('wrf nc')
        if in_ncf.dimensions['bottom_top'] == 50:
            print("50")
            layer_idx = self.layer_map(layers_fn)
        else:
            print("not50")
            layer_idx = [x for x in range(layers)]
        # Loop through and subset each species variable
        for varname, dims in metvars.items():
            print(varname)
            var = in_ncf.variables[varname]
            if 'bottom_top' in dims:
                #print(layer_idx)
                idx=[x-1 for x in layer_idx]
                arr = var[time_slice,idx,:,:]
            elif 'bottom_top_stag' in dims:
                stag_idx = [0,]+[x for x in layer_idx]
                arr = var[time_slice,stag_idx,:,:]          
            else:
                #print('t1')             
                arr = var[time_slice,:,:]
            var_out = self.createVariable(varname, np.float32, dims)
            for att in ['FieldType','MemoryOrder','description','units','stagger','coordinates']:
                setattr(var_out, att, getattr(var, att))
            if 'south_north_stag' in dims:
                #print('t2')
                arr = self._subset(arr, bounds, row_stag=True)
            elif 'west_east_stag' in dims:
                #print('t3')
                arr = self._subset(arr, bounds, col_stag=True)
            else:
                #print('t4')
                arr = self._subset(arr, bounds)
            #print('hhhhh')
            #print(len(arr))
            #print(len(var_out[:]))          
            var_out[:] = arr
            self.sync() 

    def layer_map(self, fn='vert_layers_50_to_28.csv'):
        df = pd.read_csv(fn, usecols=['lay1','lay2'])
        return [x for x in df['lay1'].values]

    def _subset(self, data_in, bounds, col_stag=False, row_stag=False):
        """
        Window WRF emissions, assume row and col are last two dimensions
        """
        col_end = bounds.ocol_e
        if col_stag:
            col_end = col_end + 1
        col_slice = slice(bounds.icol_o, bounds.icol_o + col_end)
        row_end = bounds.orow_e
        if row_stag:
            row_end = row_end + 1
        row_slice = slice(bounds.irow_o, bounds.irow_o + row_end)
        if len(data_in.shape) == 3:
            data_out = data_in[:,row_slice,col_slice]
        elif len(data_in.shape) == 4:
            data_out = data_in[:,:,row_slice,col_slice]
        return data_out

    def append_alt(self, dens):
        '''
        Append the inverse density from the MCIP
        '''
        dims = ('Time','bottom_top','south_north','west_east')
        var_out = self.createVariable('ALT', np.float32, dims)
        var_out.description = 'Inverse MCIP DENS'
        var_out.units = 'm**3/kg'
        arr = dens[:24,:]
        if arr.shape == var_out.shape:
            var_out[:] = 1/arr
            self.sync()
        else:
            raise ValueError('Input shape of DENS does not match output dimensions')

    def append_cmaq(self, cmaq, dens, layers_fn):
        '''
        Append the CMAQ concentrations
        '''
        layer_idx = self.layer_map(layers_fn)
        dims = ('Time','bottom_top','south_north','west_east')
        # MW
        mw = {'PAR': 72.1, 'ETH': 28, 'ETHY': 26, 'MEOH': 32, 'ETOH': 46.1, 'OLE': 42.1, 'TOL': 92.1, 
          'XYLMN': 106.2, 'FORM': 30, 'ALD2': 44, 'ETHA': 30.1, 'IOLE': 56.1, 'ALDX': 58.1, 'NAPH': 128.2, 
          'PRPA': 44.1, 'KET': 72.1, 'ISOP': 68.1, 'TERP': 136, 'SESQ': 204, 'NH3': 17.031, 'NO': 30.01,
          'NO2': 46, 'SO2': 64, 'SULF': 98, 'NO3': 62, 'N2O5': 108, 'HONO': 47, 'HNO3': 63, 'PNA': 79,
          'CRON': 153, 'CLNO2': 81.5, 'PAN': 121, 'PANX': 121, 'OPAN': 161, 'NTR1': 119.1, 'NTR2': 135.1, 
          'INTR': 147.1, 'OH': 17.1, 'HO2H': 34}
        # Output variable definitions for gasses and aerosols
        cmaq_vars = {'aVOC': {'type': 'gas', 'units': 'ug/m3', 'species': ['PAR','ETH','ETHY','MEOH',
            'ETOH','OLE','TOL','XLYMN','FORM','ALD2','ETHA','IOLE','ALDX','NAPH','PRPA','KET']},
          'bVOC': {'type': 'gas', 'units': 'ug/m3', 'species': ['ISOP','TERP','SESQ']},
          'aSOA': {'type': 'aero', 'units': 'ug/m3', 'species': ['AXYL1J','AXYL2J','AXYL3J','ATOL1J',
             'ATOL2J','ATOL3J','ABNZ1J','ABNZ2J','ABNZ3J','AALK1J','AALK2J','AOLGAJ','APAH1J','APAH2J','APAH3J']},
          'bSOA': {'type': 'aero', 'units': 'ug/m3', 'species': ['AISO1J','AISO2J','AISO3J','ATRP1J',
             'ATRP2J','ASQTJ','AOLGBJ']},
          'TotalPM25': {'type': 'aero', 'units': 'ug/m3', 'species': ['ASO4I','ANO3I','ANH4I','ANAI',
             'ACLI','AECI','AOTHRI','ASO4J','ANO3J','ANH4J','ANAJ','ACLJ','AECJ','AOTHRJ','AFEJ',
             'ASIJ','ATIJ','ACAJ','AMGJ','AMNJ','AALJ','AKJ','ALVPO1I','ASVPO1I','ASVPO2I','ALVPO1J',
             'ASVPO1J','ASVPO2J','ASVPO3J','AIVPO1J','ALVOO1I','ALVOO2I','ASVOO1I','ASVOO2I',
             'AXYL1J','AXYL2J','AXYL3J','ATOL1J','ATOL2J','ATOL3J','ABNZ1J','ABNZ2J','ABNZ3J',
             'AISO1J','AISO2J','AISO3J','ATRP1J','ATRP2J','ASQTJ','AALK1J','AALK2J','APAH1J',
             'APAH2J','APAH3J','AORGCJ','AOLGBJ','AOLGAJ','ALVOO1J','ALVOO2J','ASVOO1J','ASVOO2J',
             'ASVOO3J','APCSOJ']},
          'gNH': {'type': 'gas', 'units': 'ug/m3', 'species': ['NH3',]},
          'oh': {'type': 'gas', 'units': 'ug/m3', 'species': ['OH',]},
# fake variable, change in the future
          'h2o2': {'type': 'gas', 'units': 'ug/m3', 'species': ['HO2H',]},

          'gNO': {'type': 'gas', 'units': 'ug/m3', 'species': ['NO','NO2']},
          'gS': {'type': 'gas', 'units': 'ug/m3', 'species': ['SULF','SO2']},
          'pNH': {'type': 'aero', 'units': 'ug/m3', 'species': ['ANH4I','ANH4J']},
          'pNO': {'type': 'aero', 'units': 'ug/m3', 'species': ['ANO3I','ANO3J']},
          'pS': {'type': 'aero', 'units': 'ug/m3', 'species': ['ASO4I','ASO4J']},
          'gN': {'type': 'gas', 'units': 'ug/m3', 'species': ['NO3','N2O5','N2O5','HONO','HNO3','PNA',
             'CRON','CLNO2','PAN','PANX','OPAN','NTR1','NTR2','INTR']}}
        for varname, desc in cmaq_vars.items():
            print(varname)
#            print(dims)
 #           arr_out = np.zeros(dims)
            var_out = self.createVariable(varname, np.float32, dims)
            var_out.description = '+'.join(desc['species'])
            var_out.units = desc['units']
  #added
            arr_out = None

            for spec in desc['species']:
                if spec in cmaq.variables:
                    if desc['type'] == 'gas':
                        mult = mw[spec] * 1000 * dens[:24,:] / 28.9647
                    else:
                        mult = 1
                    var = cmaq.variables[spec]
                    # This is setup to work with 35/25 layer files
                    if len(cmaq.dimensions['LAY']) == 50:
                        if arr_out is None:
                            arr_out = var[:24,layer_idx,:] * mult
                        else:
                            arr_out += var[:24,layer_idx,:] * mult
                    else:
                        if arr_out is None:
                            arr_out = var[:24,:] * mult
                        else: 
                            arr_out += var[:24,:] * mult
            var_out[:] = arr_out[:]
            self.sync()
        # Calc NO/NO2 partition
        var_out = self.createVariable('NO_NO2partitioning', np.float32, dims)
        var_out.description = 'NO/(NO+NO2)'
        var_out.units = 'fraction'
        if len(cmaq.dimensions['LAY']) == 50:
            var_out[:] = cmaq.variables['NO'][:24,layer_idx,:] / (cmaq.variables['NO'][:24,layer_idx,:] + cmaq.variables['NO2'][:24,layer_idx,:])
        else:
            var_out[:] = cmaq.variables['NO'][:24,:] / (cmaq.variables['NO'][:24,:] + cmaq.variables['NO2'][:24,:])
        self.sync()
        # Calc partitions from previously calculated vars
        partitions = {'bOrgPartitioning': {'num': 'bSOA', 'den': ['bSOA','bVOC']},
          'aOrgPartitioning': {'num': 'aSOA', 'den': ['aSOA','aVOC']},
          'NHPartitioning': {'num': 'pNH', 'den': ['gNH','pNH']},
          'NOPartitioning': {'num': 'pNO', 'den': ['gNO','pNO','gN']},
          'SPartitioning': {'num': 'pS', 'den': ['gS','pS']}}
        for varname, desc in partitions.items():
            print(varname)
            var_out = self.createVariable(varname, np.float32, dims)
            var_out.description = '%s/(%s)' %(desc['num'], '+'.join(desc['den']))
            var_out.units = 'fraction'
            #arr_out = np.zeros(dims)
            arr_out = None
            for poll in desc['den']:
                print(poll)
                if arr_out is None:
                    arr_out = self.variables[poll][:]
                else:
                    arr_out += self.variables[poll][:]
            var_out[:] = self.variables[desc['num']][:] / arr_out
            self.sync()

def main():
#    if len(sys.argv) != 7:
#        raise ValueError('./gen_wrfcmaq.py wrfoutput metcro3d cmaq_conc rundate outputfilename outputgrid')
    # Command line options
    wrf = sys.argv[1]       # Path to WRF output file
    mcip = sys.argv[2]      # Path to MCIP METCRO3D file
    cmaq = sys.argv[3]      # Path to CMAQ output concentration file
    rundate = sys.argv[4]   # Current date to process
    inmap_out = sys.argv[5] # Output file name
    # Grid description file (https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html)
    griddesc = sys.argv[7]
    # WRF grid name as it appears in the grid description file
    wrf_grid = 'WRF_1_CA'
    out_grid = 'CMAQ_1_CA'
    # Mapping of 50 layer files to 28 layer files
    layers_fn = sys.argv[6]
    in_grid = Grid(wrf_grid, griddesc)
    out_grid = Grid(out_grid, griddesc)
    with NCF(wrf) as in_ncf, NCF(inmap_out, 'w') as out_ncf, NCF(mcip) as mcip, NCF(cmaq) as cmaq:
        # Regrid the WRF input to the CMAQ grid and domain
        out_ncf.regrid(in_ncf, in_grid, out_grid, rundate, layers_fn)
        # Insert the ALT variable from the MCIP DENS
        out_ncf.append_alt(mcip.variables['DENS'])
        # Append the CMAQ concentrations
        out_ncf.append_cmaq(cmaq,mcip.variables['DENS'], layers_fn)

if __name__ == '__main__':
	main()
