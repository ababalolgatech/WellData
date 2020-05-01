# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:42:59 2019

@author: Ayodeji Babalola

"""

import numpy as np
import matlab_funcs as mfun
import well as wd

well = getattr(wd,"well") # extracting well class
class func:
#----------------------------------------------------
    def read(self):
       if (self.well_flag is None):
           self.init()
      
       self.Inline  = np.zeros(self.nwells)
       self.Crossline = np.zeros(self.nwells)
       self.Xloc = np.zeros(self.nwells)
       self.Yloc = np.zeros(self.nwells)   
       
       for i in range(self.nwells):
           self.wells[i] = mfun.load_obj(self.wellnames[i]) # m is welldata class
           self.wells[i].wellname = self.wellnames[i]
           if(self.wells[i].Inline is not None):
               self.Inline[i] = self.wells[i].Inline
               
           if(self.wells[i].Crossline is not None):
               self.Crossline[i] = self.wells[i].Crossline
           
           if(self.wells[i].Xloc is not None):
               self.Xloc[i] = self.wells[i].Xloc
               
           if(self.wells[i].Yloc is not None):
               self.Yloc[i] = self.wells[i].Yloc        
#----------------------------------------------------
    def init(self):
        if (self.wellnames is None):
            raise Exception ("Please add wellnames")        
        self.nwells = self.wellnames.size
        self.wells  = mfun.create_obj(well,self.nwells)  # create nwells welldata 
        self.flag = True
#----------------------------------------------------   
    def extract_logs(self,wlog):
        LOGS = mfun.cell(self.nwells,2)  # pseudo cell_arrays        
        for i in range(self.nwells):
            LOGS[i][0] = self.wells[i].time
            LOGS[i][1] = getattr(self.wells[i],wlog)        
        return LOGS
#----------------------------------------------------
    def calc_AI(self):
        if (self.well_flag is None):
           self.init()        
        for i in range(self.nwells):
           self.wells[i].calc_AI()

#----------------------------------------------------
    def calc_SI(self):
        if (self.well_flag is None):
           self.init()        
        for i in range(self.nwells):
           self.wells[i].calc_SI()
#----------------------------------------------------
    def calc_Vp(self):
        if (self.well_flag is None):
           self.init()        
        for i in range(self.nwells):
           self.wells[i].calc_Vp()           
#----------------------------------------------------
    def calc_Vs(self):
        if (self.well_flag is None):
           self.init()        
        for i in range(self.nwells):
           self.wells[i].calc_Vs()           