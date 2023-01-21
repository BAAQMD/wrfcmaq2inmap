# Library containing routines to define gridded modeling domains and calculate boundaries/offsets

from pyproj import Proj

class GridDef:
    """
    Define the gridded modeling domains from a wrf or ioapi file
    """
    def __init__(self):
        self.grid_atts = ['GDNAM','GDTYP','P_ALP','P_BET','P_GAM','XCENT','YCENT','XORIG',
          'YORIG','XCELL','YCELL','NCOLS','NROWS','NTHIK']

    def wrf_grid(self, ncf):
        '''
        Find the WRF grid definition
        '''
        att_pairs = {'TRUELAT1': 'P_ALP', 'TRUELAT2': 'P_BET', 'STAND_LON': 'XCENT',
          'MOAD_CEN_LAT': 'YCENT', 'DX': 'XCELL', 'DY': 'YCELL'}
        for wrfatt, gridatt in att_pairs.items():
            setattr(self, gridatt, round(float(getattr(ncf, wrfatt)), 4))
        self.P_GAM = self.XCENT
        self.NTHIK = 1
        self.GDNAM = 'WRF'
        self.NCOLS = int(ncf.dimensions['west_east'].size)
        self.NROWS = int(ncf.dimensions['south_north'].size)
        if ncf.MAP_PROJ == 1:
            self.GDTYP = 2
        else:
            self.GDTYP = int(ncf.MAP_PROJ)
        proj = Proj(self.proj4())
        # XLONG and XLAT are at the centroids
        lon = float(ncf.variables['XLONG'][0,0,0])
        lat = float(ncf.variables['XLAT'][0,0,0])
        x, y = proj(lon, lat)
        # Round the ORIG values to the closest projection value that is a multiple
        #  of a grid cell
        self.XORIG = round((x - (0.5 * self.XCELL))/self.XCELL, 0) * self.XCELL
        self.YORIG = round((y - (0.5 * self.YCELL))/self.YCELL, 0) * self.YCELL

    def io_grid(self, ncf):
        '''
        Define a grid based on an IOAPI file
        '''
        for att in self.grid_atts:
            setattr(self, att, getattr(ncf, att))

    def check_grids(self, grid2):
        '''
        Compare the projection and cell size for this grid to another
        Raises a fatal value error for mismatches
        '''
        err = ''
        for projatt in ['GDTYP','P_ALP','P_BET','P_GAM','XCENT','YCENT']:
            if getattr(self, projatt) != getattr(grid2, projatt):
                err = 'Projection mismatch between the WRF input domain and the output domain'
                break
        if self.XCELL != grid2.XCELL or self.YCELL != grid2.YCELL:
            err = 'Grid cell size mismatch'
        if err:
            raise ValueError(err)

    def proj4(self):
        '''
        Return the proj4 string for the projection used with this gridding domain
        '''
        if self.GDTYP == 1:
            proj = '+proj=latlon'
        elif self.GDTYP == 2:
            proj_vars = (self.P_ALP, self.P_BET, self.XCENT, self.YCENT)
            proj = '+proj=lcc +lat_1=%s +lat_2=%s +lon_0=%s +lat_0=%s +a=6370000 +b=6370000 +units=m +no_defs' %proj_vars[:]
        elif self.GDTYP == 6:
            proj_vars = (self.P_BET, self.YCENT, self.XCENT)
            proj = '+proj=stere +lat_ts=%s +lat_0=%s +lon_0=%s' %proj_vars[:]
        return proj

class GridBounds:
    '''
    Set the boundaries of the new grid and the old grid based on the new grid
    '''
    def __init__(self, in_grid, out_grid):
        in_grid.check_grids(out_grid)
        self._set_bounds(in_grid, out_grid)

    def _set_bounds(self, in_grid, out_grid):
        '''
        Define the boundaries and the nesting offsets of the two grids
        '''
        # Calculate the raw difference for the distance between the origination points on the new and old grids
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

