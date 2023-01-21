import numpy as np
import netCDF4 as ncf
import pandas as pd
from wrfcmaq2inmap.vardefs import * 
from wrfcmaq2inmap.gridtools import * 

vardefs = VarDefs()

# Variables to read from the WRF outputs
class InMAP(ncf.Dataset):
    """
    NCF subclass with functions for processing to the InMAP file
    """
    def __init__(self, file_name, mode='r'):
        print('Opening %s' %file_name, flush=True)
        ncf.Dataset.__init__(self, file_name, mode, format='NETCDF4_CLASSIC') #format='NETCDF3_64BIT')
        self.LAYERS = 0

    def regrid(self, in_ncf, bounds, rundate, layers_fn):
        '''
        the main regridding section
        sets loop over variables and decides how to regrid
        '''
        time_slice = self.find_date(in_ncf.variables['Times'], rundate) 
        # Define the layer mapping if the WRF layers > MCIP layers
        if in_ncf.dimensions['bottom_top'].size != self.LAYERS:
            layer_idx = self.layer_map(layers_fn)
        else:
            layer_idx = [x for x in range(layers)]
        # Staggered layer index
        stag_idx = [0,]+[x+1 for x in layer_idx]
        # Loop through and subset each species variable
        for varname, dims in vardefs.metvars.items():
            print(varname, flush=True)
            var = in_ncf.variables[varname]
            var_out = self.createVariable(varname, np.float32, dims)
            for att in ['description','units','stagger','coordinates']:
                setattr(var_out, att, getattr(var, att))
            # Differentiate staggered and unstaggered col/rows
            if 'south_north_stag' in dims:
                row_slice, col_slice = self._cell_slice(bounds, row_stag=True)
            elif 'west_east_stag' in dims:
                row_slice, col_slice = self._cell_slice(bounds, col_stag=True)
            else:
                row_slice, col_slice = self._cell_slice(bounds)
            # Differentiate staggered and unstaggered layers
            if 'bottom_top' in dims:
                lay_slice = layer_idx
            elif 'bottom_top_stag' in dims:
                lay_slice = stag_idx
            if len(var.shape) == 3:
                var_out[:] = var[time_slice,row_slice,col_slice]
            else:
                var_out[:] = var[time_slice,lay_slice,row_slice,col_slice]
            self.sync() 

    def layer_map(self, fn):
        '''
        Map the layers from one definition to another
        Where lay1 is the layer set you map from and lay2 is the target layer set
        '''
        df = pd.read_csv(fn, usecols=['lay1','lay2'])
        return [x-1 for x in df['lay1'].values]

    def _cell_slice(self, bounds, col_stag=False, row_stag=False):
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
        return (row_slice, col_slice)

    def find_date(self, times, rundate):
        '''
        Find the index for the first and last hour of the run date in the wrf output file
        '''
        dates = []
        for dt in times:
            dt = ''.join([c.decode() for c in dt])
            dates.append(dt[:4]+dt[5:7]+dt[8:10]+dt[11:13])
        try:
            start_time = dates.index(str(rundate)+'00')
        except ValueError:
            raise ValueError('Could not find %s in WRF' %rundate)
        try:
            end_time = dates.index(str(rundate)+'23')
        except ValueError:
            raise ValueError('Could not find %s in WRF' %rundate)
        else:
            return slice(start_time, end_time+1)

    def set_dims(self, in_ncf, out_grid):
        '''
        Set up the outfile dimensions and attributes based on the new grid
        and the input file
        '''
        out_dims = {'Time': 24, 'bottom_top': self.LAYERS, 'bottom_top_stag': self.LAYERS+1,
          'south_north': out_grid.NROWS, 'south_north_stag': out_grid.NROWS+1, 
          'west_east': out_grid.NCOLS, 'west_east_stag': out_grid.NCOLS+1} 
        for dim, value in out_dims.items(): 
            self.createDimension(dim, int(value))
        inatts_ = in_ncf.__dict__
        for att_name in in_ncf.ncattrs():
            setattr(self, att_name, inatts_.get(att_name)) 

    def append_alt(self, dens):
        '''
        Append the inverse density from the MCIP
        '''
        print('ALT', flush=True)
        dims = ['Time','bottom_top','south_north','west_east']
        var_out = self.createVariable('ALT', np.float32, dims)
        var_out.description = 'Inverse MCIP DENS'
        var_out.units = 'm**3/kg'
        arr = dens[:24,:]
        if arr.shape == var_out.shape:
            var_out[:] = 1/arr
            self.sync()
        else:
            print(arr.shape, var_out.shape)
            raise ValueError('Input shape of DENS does not match output dimensions')

    def append_cmaq(self, cmaq):
        '''
        Append the CMAQ concentrations if they are already in the defined CMAQ output file
        '''
        dims = ['Time','bottom_top','south_north','west_east']
        for varname in vardefs.cmaq_vars:
            print(varname)
            var = cmaq.variables[varname]
            var_out = self.createVariable(varname, np.float32, dims)
            var_out.description = var.var_desc
            var_out.units = var.units
            var_out[:] = var[:24,:]
            self.sync()
     
    def append_calc_cmaq(self, cmaq, dens, mech):
        '''
        Calculate the partitioning variables from the CMAQ concentrations and append to the netCDF
        '''
        vardefs.set_mech(mech)
        dims = ['Time','bottom_top','south_north','west_east']
        for varname, desc in vardefs.cmaq_map.items():
            print(varname, flush=True)
            var_out = self.createVariable(varname, np.float32, dims)
            arr_out = np.zeros(var_out.shape)
            var_out.description = '+'.join(desc['species'])
            var_out.units = desc['units']
            for spec in desc['species']:
                if desc['type'] == 'gas':
                    if desc['units'] == 'ppbC':
                        coeff = vardefs.mw[spec] * 1000
                    else:
                        coeff = vardefs.mw[spec] * 1000 * dens[:24,:] / 28.9647
                else:
                    coeff = 1
                try:
                    var = cmaq.variables[spec]
                except KeyError:
                    print('WARNING: Missing %s in CMAQ conc' %spec)
                else:
                    arr_out += var[:24,:] * coeff
            var_out[:] = arr_out[:]
            self.sync()
        self.calc_no_part(cmaq, dims)
        self.calc_other_part(cmaq, dims)

    def calc_no_part(self, cmaq, dims):
        '''
        Calculate and fill the NO/NO2 partition
        '''
        var_out = self.createVariable('NO_NO2partitioning', np.float32, dims)
        var_out.description = 'NO/(NO+NO2)'
        var_out.units = 'fraction'
        var_out[:] = cmaq.variables['NO'][:24,:] / (cmaq.variables['NO'][:24,:] +\
          cmaq.variables['NO2'][:24,:])
        self.sync()

    def calc_other_part(self, cmaq, dims):
        '''
        Calc partitions from previously calculated vars
        '''
        partitions = {'bOrgPartitioning': {'num': 'bSOA', 'den': ['bSOA','bVOC']},
          'aOrgPartitioning': {'num': 'aSOA', 'den': ['aSOA','aVOC']},
          'NHPartitioning': {'num': 'pNH', 'den': ['gNH','pNH']},
          'NOPartitioning': {'num': 'pNO', 'den': ['gNO','pNO','gN']},
          'SPartitioning': {'num': 'pS', 'den': ['gS','pS']}}
        for varname, desc in partitions.items():
            print(varname, flush=True)
            var_out = self.createVariable(varname, np.float32, dims)
            var_out.description = '%s/(%s)' %(desc['num'], '+'.join(desc['den']))
            var_out.units = 'fraction'
            arr_out = np.zeros(var_out.shape)
            for poll in desc['den']:
                arr_out += self.variables[poll][:]
            var_out[:] = self.variables[desc['num']][:] / arr_out
            self.sync()

if __name__ == '__main__':
	main()
