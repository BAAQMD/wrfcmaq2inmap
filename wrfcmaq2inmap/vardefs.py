# Variable definitions related to the CMAQ and WRF files for InMAP processing

class VarDefs:
    def __init__(self):
        self._init_vars()

    def set_mech(self, mech):
        if mech.lower() == 'cb6':
            self.mw = self.cb6_coeff.copy()
            self.cmaq_map = self.cb6_map.copy()
        elif mech.lower() == 'saprc':
            self.mw = self.saprc_coeff.copy()
            self.cmaq_map = self.saprc_map.copy()
        self.mw.update(self.nonvoc_coeff)
        self.cmaq_map.update(self.nonvoc_map)

    def _init_vars(self):
        # Set the input WRF variables
        self.metvars = {
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
        # CMAQ variables to write to InMAP file
        self.cmaq_vars = ['TotalPM25','gS','pS','aVOC','bVOC','aSOA','bSOA','oh','h2o2','pNO','gNO','pNH','gNH']
        # Mapping for non-VOC clumps
        self.nonvoc_coeff = {'NH3': 17.031, 'NO': 30.01, 'INTR': 147.1,
          'NO2': 46, 'SO2': 64, 'SULF': 98, 'NO3': 62, 'N2O5': 108, 'HONO': 47, 'HNO3': 63, 'PNA': 79,
          'CRON': 153, 'CLNO2': 81.5, 'PAN': 121, 'PANX': 121, 'OPAN': 161, 'NTR1': 119.1, 'NTR2': 135.1}
        self.nonvoc_map = {'aSOA': {'type': 'aero', 'units': 'ug/m3', 'species': ['AXYL1J','AXYL2J','AXYL3J','ATOL1J',
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
          'gNO': {'type': 'gas', 'units': 'ug/m3', 'species': ['NO','NO2']},
          'gS': {'type': 'gas', 'units': 'ug/m3', 'species': ['SULF','SO2']},
          'pNH': {'type': 'aero', 'units': 'ug/m3', 'species': ['ANH4I','ANH4J']},
          'pNO': {'type': 'aero', 'units': 'ug/m3', 'species': ['ANO3I','ANO3J']},
          'pS': {'type': 'aero', 'units': 'ug/m3', 'species': ['ASO4I','ASO4J']},
          'gN': {'type': 'gas', 'units': 'ug/m3', 'species': ['NO3','N2O5','N2O5','HONO','HNO3','PNA',
             'CRON','CLNO2','PAN','PANX','OPAN','NTR1','NTR2','INTR']}}
        # Coefficients for CB6 calculations (MW)
        self.cb6_coeff = {'PAR': 72.1, 'ETH': 28, 'ETHY': 26, 'MEOH': 32, 'ETOH': 46.1, 'OLE': 42.1, 'TOL': 92.1, 
          'XYLMN': 106.2, 'FORM': 30, 'ALD2': 44, 'ETHA': 30.1, 'IOLE': 56.1, 'ALDX': 58.1, 'NAPH': 128.2, 
          'PRPA': 44.1, 'KET': 72.1, 'ISOP': 68.1, 'TERP': 136, 'SESQ': 204}
        # CB6 CMAQ output variable definitions
        self.cb6_map = {'aVOC': {'type': 'gas', 'units': 'ug/m3', 'species': ['PAR','ETH','ETHY','MEOH',
            'ETOH','OLE','TOL','XLYMN','FORM','ALD2','ETHA','IOLE','ALDX','NAPH','PRPA','KET']},
          'bVOC': {'type': 'gas', 'units': 'ug/m3', 'species': ['ISOP','TERP','SESQ']}}
        # Coefficients for CB6 calculations (MW)
        self.saprc_coeff = {'IPRD': 5, 'MACR': 4, 'ISOPRENE': 5, 'APIN': 10, 'TERP': 10, 'SESQ': 15,
          'NPHE': 6, 'CRES': 7, 'BALD': 7, 'BENZENE': 6, 'TOLUENE': 7, 'MXYL': 8, 'OXYL': 8,
          'PXYL': 8, 'ARO1': 7, 'ARO2MN': 8, 'NAPHTHAL': 10, 'ALK1': 2, 'ALK2': 3, 'ALK3': 4, 
          'ALK4': 5, 'ALK5': 8, 'RCOOH': 3, 'CCOOH': 2, 'HCOOH': 1, 'ACETONE': 3, 'ACETYLENE': 2,
          'ACROLEIN': 3, 'BACL': 4, 'BUTADIENE13': 4, 'CCHO': 2, 'CRES': 7, 'ETHENE': 2, 'ETOH': 2,
          'GLY': 2, 'HCHO': 1, 'MEK': 4, 'MOEH': 1, 'MGLY': 3, 'MVK': 4, 'OLE1': 5, 'OLE2': 5,
          'PRD2': 6, 'PROPENE': 3, 'RCHO': 3, 'RNO3': 6}
        self.saprc_map = {'aVOC': {'type': 'gas', 'units': 'ppbC', 'species': ['NPHE','CRES',
            'BALD','BENZENE','TOLUENE','MXYL','OXYL','PXYL','ARO1','ARO2MN','NAPHTHAL','ALK1','ALK2',
            'ALK3','ALK4','ALK5']},
          'bVOC': {'type': 'gas', 'units': 'ppbC', 'species': ['IPRD','MACR','ISOPRENE','APIN',
            'TERP','SESQ']}}

