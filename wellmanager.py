# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:34:21 2019

@author: Ayodeji Babalola
"""
import well as wd
import numpy as np
import matlab_funcs as mfun
#from matlab_funcs import save_obj,load_obj

well = getattr(wd,"well") # extracting well class
class wellmanager():
    def __init__(self,wellnames):
        self.wells = well()
        self.nwells = None
        self.wellnames = wellnames
        self.Inline = None
        self.Crossline = None
        self.Xloc = None
        self.Yloc = None
        self.well_flag = None
        self.domain = None  
        self.init()
        
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
            print ("To do :Please add wellnames")        
        self.nwells = self.wellnames.size
        self.wells  = mfun.create_obj(well,self.nwells)  # create nwells welldata 
        self.flag = True
#----------------------------------------------------   
    def extract_logs(self,wlog):
        LOGS = mfun.cell(self.nwells,2)  # pseudo cell_arrays        
        for i in range(self.nwells):
            LOGS[i,0] = self.wells[i].time
            LOGS[i,1] = getattr(self.wells[i],wlog)        
        return LOGS
#----------------------------------------------------
    def calc_AI(self):
        for i in range(self.nwells):
           self.wells[i].calc_AI()

#----------------------------------------------------
    def calc_SI(self):       
        for i in range(self.nwells):
           self.wells[i].calc_SI()
#----------------------------------------------------
    def calc_Vp(self):     
        for i in range(self.nwells):
           self.wells[i].calc_Vp()           
#----------------------------------------------------
    def calc_Vs(self):       
        for i in range(self.nwells):
           self.wells[i].calc_Vs()         
#-----------------------------------------------------------------------------
#                                   TESTS        
#*****************************************************************************                 
 
#----------------------------------------------- 
if __name__ == '__main__':   
    """
    fn = ('depth','gr','por','rho','sw','vp','vs')
    fn_arr = np.array(['depth','gr','por','rho','sw','vp','vs'])
    mat = mfun.load_ascii('x1001.las',12)
    xx1 = mfun.load_ascii('x1001.las',12,fn)
    xx2 = mfun.load_ascii('x1001.las',12,fn_arr)
    
    x1001 = well()
    
    x1001.depth = mat[:,0]
    x1001.time = mat[:,0]
    x1001.GR = mat[:,1]
    x1001.Por = mat[:,2]
    x1001.Rhob = mat[:,3]
    x1001.Sw = mat[:,4]
    x1001.Vp = mat[:,5]
    x1001.Vs = mat[:,6]
    
    x1002 = well()
    x1002.time = mat[:,0]*1000
    x1002.depth = mat[:,0]*1000
    x1002.Sw = mat[:,4]*1000
    x1002.Vp = mat[:,5]*1000
    x1002.Vs = mat[:,6]*1000
    
    # loading and saving objects
    mfun.save_obj("x1001_saved.obj",x1001)
    mfun.save_obj("x1002_saved.obj",x1002)
    x1001_loaded = mfun.load_obj("x1001_saved.obj")
    """
    # wellmanager
    
    wellnames = np.array(['wells_files\\w1002','wells_files\\w1003'])
    wm = wellmanager(wellnames)
    wm.read()
    wm.calc_AI()
    Vp = wm.extract_logs("Vp")

            