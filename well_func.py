# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 12:56:24 2019

@author: Chris
"""
import matlab_funcs as mfun
import matplotlib.pyplot as plt
import numpy as np
class func:
#---------------------------------
    def calc_Vp(self,units  = None):
        """
        if (self.init_flag is None):
            self.init()
        """

        if (self.DT is None):
            raise Exception ('Sonic log is not loaded')

        if(units is not None):
            self.units = units
        elif(units is not None and self.units is None):
            self.units = 'm'
            print('Vp is calculated in m/s')
            print('sonic is assumed to be in microft/s')

        if (self.units == 'ft'):
            self.Vp = 1e6/self.DT  ;
        else:
            self.Vp = 304800/self.DT  ;
#---------------------------------
    def calc_Vs(self,units = None):
        """
        if (self.init_flag is None):
            self.init()
        """
        if (self.DTS is None):
            raise Exception ('Sonic log is not loaded')

        if(units is not None):
            self.units = units
        elif(units is not None and self.units is None):
            self.units = 'm'
            print('Vp is calculated in m/s')
            print('sonic is assumed to be in microft/s')

        if (self.units == 'ft'):
            self.Vs = 1e6/self.DTS 
        else:
            self.Vs = 304800/self.DTS  
#--------------------------------- 
    def calc_AI(self,units  = None):
        """
        if (self.init_flag is None):
            self.init()
        """
        if (self.Vp is not None and self.Rhob is not None):
            self.AI = self.Vp* self.Rhob
        else:
            raise Exception ('Log unavailable')

#---------------------------------
    def calc_SI(self,units = None):
        """
        if (self.init_flag is None):
            self.init()
        """

        if (self.Vs is not None and self.Rhob is not None):
            self.SI = self.Vp* self.Rhob
        else:
            raise Exception ('Log unavailable')
#---------------------------------             
    def read(self,filepath,nheaders,nan = None):
        if (nan is None):
            nan =  -999.25
        
        mat = mfun.load_ascii(filepath,nheaders)
        mat[mat == nan] = np.nan        
        return mat
        
#---------------------------------            
    def save_obj(self,fname,obj):
        mfun.save_obj(fname,obj)
#---------------------------------            
    def QC_plots(self,Type = None):
        if (Type is None):
            Type = 'Elastic'
        if (Type == 'Elastic'):    
            #plt.figure()
            fig, (ax1, ax2,ax3) = plt.subplots(3,1)
            ax1.plot(self.time,self.Vp)
            ax1.set_title('Vp')
            
            ax2.plot(self.time,self.Vs)
            ax2.set_title('Vs') 
            
            ax3.plot(self.time,self.Rhob)
            ax3.set_title('Rhob') 
        elif(Type == 'Sonic'):
                        #plt.figure()
            fig, (ax1, ax2,ax3) = plt.subplots(3,1)
            ax1.plot(self.time,self.DT)
            ax1.set_title('DT')
            
            ax2.plot(self.time,self.DTS,)
            ax2.set_title('DTS') 
            
            ax3.plot(self.time,self.Rhob)
            ax3.set_title('Rhob') 
        
        
        
               

                
        

