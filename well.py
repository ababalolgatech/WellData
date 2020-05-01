# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 12:53:22 2019

@author: Chris
"""
import numpy as np
import math
import matlab_funcs as mfun
import Filter as ft
from scipy import interpolate
import matplotlib.pyplot as plt
import crewes_resample as crs
from scipy.interpolate import InterpolatedUnivariateSpline
# add plotting routine for Vp, Vs and Rho
class well() :
    # http://www.qtrac.eu/pyclassmulti.html#definitive
    def __init__(self):
        self.wellname  = None
        self.Xloc = None
        self.Yloc= None
        self.KB= None
        self.Inline= None
        self.Crossline= None
        self.X = None
        self.Y = None        
        # deviation surveys
        self. MD_dev= None
        self.Xdev= None
        self.Ydev= None
        self.TVD_dev= None        
        self.dip= None
        self.azm = None
        self.dx = None
        self.dy = None
        self.dz = None
        # Timd- depth curves
        self.TD_time = None
        self.TD_MD= None
        self.TD_TVD= None
        self.TD_X= None
        self.TD_Y= None
        self.dtold= None
        self.dtnew= None
        # well-logs
        self.MD = None
        self.TVD= None
        self.DT= None
        self.DTS= None
        self.Vp = None
        self.Vs = None
        self.AI = None
        self.SI= None
        self.EI = None 
        self.EEI = None 
        self.Vsh = None
        self.Por = None 
        self.Sw = None
        self.GR = None
        self.Rhob = None
        self.SW = None
        self.DT= None
        self.DTS= None
        self.Vp = None
        self.Vs = None        
        self.units = None
        self.depth = None
        self.time = None
        self.rockpar = None        
        # flags        
        self.zone_logs = None 
        self.prist_logs = None 
        self.zone_data = None
        # self.wellpath_flag = None
        self.init_flag = None
        self.tdomain_flag = None    # True False
        self.applydev_flag = None  # True , False
        self.wellpath = None   # deviated, vertical
        self.TD_flag = None   # MD,TVD
        self.lognames = np.array(['time','MD','TVD','X','Y','DT','DTS',\
                                  'GR','Rhob','depth','Vp','Vs','AI','SI',\
                                  'EI','EEI','Vsh','Por','Sw']) 
        #self.init()
        
#----------------------------         
    def __repr__(self):
            return repr("Well Class")  
    
#----------------------------
    def init(self):
        nlogs = len(self.lognames)
        if(self.MD is not None or self.tdomain_flag == False):
            z = self.MD
            self.tdomain_flag = False
        elif(self.time is not None or self.tdomain_flag == True):
            z = self.time
            self.tdomain_flag = True            

        self.prist_logs = np.zeros((z.size,nlogs))
        self.zone_data = np.zeros((z.size,nlogs))
        
        if (self.zone_logs is None):
            self.zone_logs = np.zeros((z.size,nlogs))
            for i in range(nlogs):  
                tmp = getattr(self,self.lognames[i])    
                if(tmp is not None):
                    self.zone_logs[:,i] = tmp
                    self.prist_logs[:,i] = tmp
                else:
                    continue                    
        # whenever i load time logs
        if(self.depth is None and self.time is not None):
            self.tdomain_flag = True
        self.init_flag = True
        
        if(self.MD_dev is not None and self.dx is not None and self.dy is not None):
            self.wellpath = 'dxy'  
        elif(self.MD_dev is not None and self.dip is not None and self.azm is not None):
            self.wellpath = 'inc_azm'          
        elif(self.MD_dev is not None and self.Xdev is not None and self.Ydev is not None):
            self.wellpath = 'XY'   
        elif(self.MD_dev is None and self.Xdev is None and self.Ydev is None and \
             self.dip is None and self.azm is None and self.dx is None and self.dy is  None):
            self.wellpath = 'vertical'           
 #----------------------------
    def apply_TD(self,dtnew = None):
        if (self.init_flag is None or self.init_flag is False):
            self.init()
            
        if(self.tdomain_flag is None and  self.tdomain_flag == True):
            self.depth_domain()
        # check for imported checkshot
        if(self.TD_time is None):
            raise Exception ('import checkshot')        
        if(self.TD_time is not None and self.TD_TVD is not None):
            self.TD_flag = 'TVD'  # set depth = TVD for time conversion
        if(self.TD_time is not None and self.TD_MD is not None): # Prefered because depth is in MD - no changes
            self.TD_flag = 'MD'  # set depth = MD for TD
        """    
        if(self.TD_MD is not None and self.TD_time is not None and self.TD_TVD is not None \
           and self.TD_X is not None and self.TD_Y is not None):
            self.TD_flag = 'HRS'
        """
            
        self.dt = 1  # Prefered
        if (dtnew is None and self.dtnew is None):
            dtnew=1
            self.dtnew =1
        elif(self.dtnew is None):
            self.dtnew = dtnew  
            
        # check if deviation exists
        if(self.wellpath == 'vertical'):
            self.depth = self.MD
            self.TVD = self.MD
        else:
            self.apply_dev()
            if(self.TD_flag == 'MD'):
                self.depth = self.MD
                TD_depth = self.TD_MD
            elif(self.TD_flag == 'TVD'):
                self.depth = self.TVD
                TD_depth = self.TD_TVD
                #self.recoordinate_logs_MD2TVD()
        """
        # assume checkshot is from KB
        if(self.KB is None):
            self.KB = 0 
            print ('Beware : KB is forced to zero')   
        else:
            self.TVD = self.TVD - self.KB 
            print( 'debug : self.TVD - self.KB : is this right ?')
         """  
         
        nlogs = len(self.lognames)
        for i in range(nlogs):
            tmp = getattr(self,self.lognames[i])    
            if(tmp is None or self.lognames[i] == 'depth' or self.lognames[i] == 'time' \
               or self.lognames[i] == 'MD' or self.lognames[i] == 'TVD'  or \
               self.lognames[i] == 'TD_MD' or self.lognames[i] == 'TD_TVD' or self.lognames[i] == 'TD_time'):
                continue
            else:
                if (self.lognames[i] == 'X' or self.lognames[i] == 'Y'):
                    new_time,wlog = apply_td(self.TD_time,TD_depth,self.depth,tmp,dtnew,'track') 
                    #new_time,wlog = apply_td(self.TD_time,TD_depth,self.depth,tmp,dtnew) 
                    wlog = mfun.remove_nan_interp(wlog) # temporary BUG fix
                    setattr(self,self.lognames[i],wlog)
                else:
                     new_time,wlog = apply_td(self.TD_time,TD_depth,self.depth,tmp,dtnew)
                     setattr(self,self.lognames[i],wlog)  
                     
                     
        self.time = new_time
        self.tdomain_flag = True                  
        # adding depth time 
        tmp = self.depth
        new_time,wlog = apply_td(self.TD_time,TD_depth,self.depth,tmp,dtnew,'track') 
        wlog = mfun.remove_nan_interp(wlog)
        self.depth = wlog
#---------------------------------
    def recoordinate_logs_MD2TVD(self,units  = None): 
        
        nlogs = len(self.lognames)
        for i in range(nlogs):
            tmp = getattr(self,self.lognames[i])    
            if(tmp is None or self.lognames[i] == 'depth' or self.lognames[i] == 'time' \
               or self.lognames[i] == 'MD' or self.lognames[i] == 'TVD'  or \
               self.lognames[i] == 'TD_MD' or self.lognames[i] == 'TD_TVD' or self.lognames[i] == 'TD_time'):
                continue
            else:
                f_tvd = interpolate.interp1d(self.MD_dev,tmp, kind ='cubic',fill_value="extrapolate")
                wlog = f_tvd(self.MD)
                setattr(self,self.lognames[i],wlog)       
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
            self.SI = self.Vs* self.Rhob
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
            
            if (self.Vs is not None):
                ax2.plot(self.time,self.Vs)
                ax2.set_title('Vs') 
                
            if (self.Rhob is not None) :
                ax3.plot(self.time,self.Rhob)
                ax3.set_title('Rhob') 
        elif(Type == 'Sonic'):
                        
            fig, (ax1, ax2,ax3) = plt.subplots(3,1)
            ax1.plot(self.time,self.DT)
            ax1.set_title('DT')
            
            if (self.DTS is not None):            
                ax2.plot(self.time,self.DTS)
                ax2.set_title('DTS') 
            
            if(self.Rhob is not None):
                ax3.plot(self.time,self.Rhob)
                ax3.set_title('Rhob')         
        
        
#---------------------------------            
    def apply_dev(self):
        """
        Apply dev will:
        - compute the MD,X,Y,Z from MD,azm,inc
        - convert the depth values and angle information to current well depth range
        """
        #1
        if (self.wellpath == 'inc_azm'):
            ns = self.azm.size
            dx,dy,self.TVD_dev = tvd_minc(self.MD_dev,self.azm,self.dip)
            self.Xdev = self.Xloc + dx
            self.Ydev = self.Yloc + dy
            #self.TVD_dev = self.MD_dev + dz            
        if(self.wellpath == 'vertical'):
            ns = self.MD.size
            self.TVD_dev = self.MD
            self.MD_dev = self.MD
            self.X = np.ones(ns)*self.Xloc
            self.Y = np.ones(ns)*self.Yloc
            self.TVD = self.MD
        else:
        #2  
            if (self.TVD_dev is None and self.MD_dev is not None and self.Xdev is not None ):
                self.calc_tvd_dev()
            f_tvd = interpolate.interp1d(self.MD_dev,self.TVD_dev, kind ='cubic',fill_value="extrapolate") 
            f_x = interpolate.interp1d(self.MD_dev,self.Xdev, kind ='cubic',fill_value="extrapolate")
            f_y = interpolate.interp1d(self.MD_dev,self.Ydev, kind ='cubic',fill_value="extrapolate")
            self.TVD = f_tvd(self.MD)
            self.X = f_x(self.MD)
            self.Y = f_y(self.MD)
        
        #self.depth = self.TVD  # Please note 
        # it seemed not be correct because using TVD checkshot.. it does not conform         
        self.applydev_flag = True

#---------------------------------
    def calc_tvd_dev(self):
        npts = self.Xdev.size
        tvd = np.zeros(npts)
        tvd[0] = self.MD_dev[0]
        for i in range(npts-1):
            dist = mfun.calc_dist(self.Xdev[i],self.Ydev[i],self.Xdev[i+1],self.Ydev[i+1])
            hyp = self.MD_dev[i+1] - self.MD_dev[i]
            if(math.isnan(dist)):
                dist = 0.00001
            if(math.isnan(hyp)):
                hyp = 0.00001
            dtvd = np.sqrt(hyp**2 - dist**2)            
            if (math.isnan(dtvd)):
                dtvd = 0
            if(i==0):
                tvd[i+1] = self.MD_dev[i]+dtvd
            else:
                tvd[i+1] = tvd[i]+dtvd
            if(i==263):
                pass
        self.TVD_dev = tvd
   
#---------------------------------
def calc_td(td_time,td_depth,depth,wlog,dt = None):
    if (dt is None):
        dt = 1
    f_t2d = interpolate.interp1d(td_time,td_depth, kind = 'linear',fill_value="extrapolate")
    f_d2t = interpolate.interp1d(td_depth,td_time,kind = 'linear', fill_value="extrapolate")
    
    tmp_time = f_d2t(depth)
    new_time = np.arange(np.int(tmp_time[0]),np.int(tmp_time[-1]),dt,dtype=int)
    znew = f_t2d(new_time)    
    #nwlog = ft.interp3(wlog,depth,znew) 
    
    nwlog = crs.sincnan(wlog,depth,znew)
    #nwlog = ft.sinc_interp(wlog,depth,znew)
    return nwlog,new_time,f_d2t

#---------------------------------
def calc_td_track(td_time,td_depth,depth,wlog,dt = None):
    """
    modified for continous logs like X, Y, and TD
    ver1 did not work well
    """
    if (dt is None):
        dt = 1
    f_t2d = InterpolatedUnivariateSpline(td_time,td_depth)
    f_d2t = InterpolatedUnivariateSpline(td_depth,td_time)
    
    tmp_time = f_d2t(depth)
    new_time = np.arange(np.int(tmp_time[0]),np.int(tmp_time[-1]),dt,dtype=int)
    znew = f_t2d(new_time) 
    npts = znew.size
    nwlog = ft.resamp_by_newlength(wlog,npts)   
    return nwlog,new_time,f_d2t

#---------------------------------
def calc_td_track1(td_time,td_depth,depth,wlog,dt = None):
    """
    modified for continous logs like X, Y, and TD
    ver1 and ver2 did not work well
    """
    if (dt is None):
        dt = 1
    f_t2d = interpolate.interp1d(td_time,td_depth, kind = 'linear',fill_value="extrapolate")
    f_d2t = interpolate.interp1d(td_depth,td_time,kind = 'linear', fill_value="extrapolate")
    
    tmp_time = f_d2t(depth)
    new_time = np.arange(np.int(tmp_time[0]),np.int(tmp_time[-1]),dt,dtype=int)
    znew = f_t2d(new_time) 
    nwlog = crs.sincnan(wlog,depth,znew)
    return nwlog,new_time,f_d2t
#---------------------------------  
def apply_td(td_time,td_depth,depth,wlog,dt = None,flag = None): 
    if (dt is None):
        dt=1
    if (flag is None):
        flag = 'wlog'
    """
    Apply TD will:
        - remove nan from top and base
        - interpolate nan in-between
        - apply td to clean data
        - recreate the origin time and depth and
        - add the nan values back
    """
    ndepth = depth
    depth,wlog = mfun.remove_nan_topbase(depth,wlog)   
    wlog = mfun.remove_nan_interp(wlog)
    if (flag == 'wlog'):
        twlog,tmp_time,func_d2t  = calc_td(td_time,td_depth,depth,wlog,dt)
    else:
        twlog,tmp_time,func_d2t  = calc_td_track1(td_time,td_depth,depth,wlog,dt)
    tmin = func_d2t(ndepth[0])
    tmax = func_d2t(ndepth[-1])
    time = np.arange(tmin,tmax,dt,dtype=int)
    if (time[0]<tmp_time[0]):
        tmp = np.ones(tmp_time[0]-time[0])*np.nan
        twlog = np.append(tmp,twlog)
    if(time[-1]>tmp_time[-1]):
        tmp = np.ones(time[-1]-tmp_time[-1])*np.nan
        twlog = np.append(twlog,tmp)
    return time, twlog 
#---------------------------------
def tvd_minc(zmd,azmi,incl):
    """
     TVD_MINC uses the minimum curvature method to compute true 3-D coordinates
     of a wellbore path given information from a deviation survey. The method is
     documented in 'Directional Survey Calculation', by J.T. Craig and B.V. Randal
     found in 'Petroleum Engineer',March, 1976.
    
     zmd ... vector of measured depths
     azmi ... vector of azimuth angles in degrees
     incl ... vector of inclination angles in degrees
    """
    pad = 0
    if (zmd[0] !=0):
        zmd = np.append(0,zmd)
        azmi = np.append(0,azmi)
        incl = np.append(0,incl)
        pad = 1
    torad= np.pi/180     
    phi=torad*incl
    theta=torad*azmi
    npts = zmd.size
    k = np.arange(1,npts)
    kminus = np.arange(0,npts-1)
    cosd = np.cos(phi[k] -phi[kminus]) - np.sin(phi[kminus])*np.sin(phi[k])*\
		(1- np.cos(theta[k]-theta[k-1]))
    tand = np.sqrt(cosd**(-2) -1)
    dl = np.arctan(tand)
    ind = dl>0.00001
    not_ind = dl<0.00001  # returns bool
    fc = np.zeros(dl.size)
    #fc[not_ind] = np.ones(np.sum[not_ind])
    fc[not_ind] = 1
    fc[ind] = 2*np.tan(dl[ind]/2)/dl[ind]
    sinphi = fc*(np.sin(phi[k]) * np.sin(theta[k]) + np.sin(phi[kminus]) * np.sin(theta[kminus]))/2 
    cosphi = fc*(np.sin(phi[k]) * np.cos(theta[k]) + np.sin(phi[kminus]) * np.cos(theta[kminus]))/2 
    x=np.zeros(zmd.size)
    y=np.zeros(zmd.size)
    z=np.zeros(zmd.size)
    x[k]= np.cumsum( (zmd[k]-zmd[kminus]) * cosphi )
    y[k]= np.cumsum( (zmd[k]-zmd[kminus]) * sinphi )
    z[k]= np.cumsum( fc * (zmd[k]-zmd[kminus]) * (np.cos(phi[k])+np.cos(phi[kminus]))/2)
    
    if(pad == 1):
    	x=x[k]
    	y=y[k]
    	z=z[k]

    tmp = x
    x = y
    y = tmp
    return x,y,z   


#******************************************************************************            
if __name__ == '__main__':
    #--------------------------------------------------------
    plot = getattr(mfun,"plot")
    plotd = getattr(mfun,"plotd")
    imshow = getattr(mfun,"imshow")
    mfun.plt.close('all')
    """
    wF11A = well()
    mat = wF11A.read('well_files//volve//logs//15_9-F-11A_logs.las',35)
    wF11A.Xloc = 435049.09
    wF11A.Yloc = 6478568.18
    wF11A.Inline = 10151
    wF11A.Crossline = 2277
    wF11A.KB = 54
    wF11A.MD = mat[:,0]
    wF11A.DT = mat[:,5]
    wF11A.DTS = mat[:,3]
    wF11A.Rhob = mat[:,1]
    wF11A.Vsh = mat[:,4]
    wF11A.Por = mat[:,2]
    wF11A.SW = mat[:,6]
    wF11A.calc_Vp()
    wF11A.calc_Vs()
    
    # TD
    mat = wF11A.read('well_files//volve//TD//HRS//15_9_F_11A_well_tie.txt',16)
    wF11A.TD_time = mat[:,5]
    wF11A.TD_TVD = mat[:,3]
    wF11A.TD_MD = mat[:,0]
    wF11A.TD_X = mat[:,8]
    wF11A.TD_Y = mat[:,9]
    
    # Deviation file
    mat = wF11A.read('well_files//volve//TD//HRS//15_9_F_11A_well_tie.txt',16)
    wF11A.MD_dev = mat[:,0]
    wF11A.Xdev = mat[:,8]
    wF11A.Ydev = mat[:,9]
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(wF11A.MD,wF11A.Vp)
    
    wF11A.apply_TD()
    # Todo : Appy deviation file to to get TVD and compare the values of TVD to the checkshot
    # wF11A_time = mfun.load_obj('wells_files//volve//time//wF11A')
    wF11A_time = well()
    mat = mfun.load_ascii('well_files//volve//time//15_9-F-11A_logs_time.las',35)
    wF11A_time.time = mat[:,0]
    wF11A_time.DT = mat[:,6]
    wF11A_time.calc_Vp()
    # debug
    zone = np.zeros(2)
    if (wF11A.time [0] >= wF11A_time.time[0] ):
        zone[0] = wF11A.time[0]
    else:
        zone[0] = wF11A_time.time[0]
        
    if (wF11A.time [-1]<=wF11A_time.time[-1] ):
        zone[1] = wF11A_time.time[-1]
    else:
        zone[1] = wF11A.time[-1]
    zone[0] = 2350
    zone[1] = 2515
    
    wF11A.calc_Vp()
    ax2.plot(wF11A.time,wF11A.Vp)
    

    vp_time_zn,new_time_zn = mfun.segment_logs(zone,wF11A.time,wF11A.Vp)   
    org_vp_time_zn,org_time_zn = mfun.segment_logs(zone,wF11A_time.time,wF11A_time.Vp)
    
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(new_time_zn,vp_time_zn)
    ax2.plot(org_time_zn,org_vp_time_zn) 
    wF11A.save_obj('wF11A_jan2020', wF11A)
    
    """    
    w19_19A = well()
    mat = w19_19A.read('well_files//volve//logs//15_9-19A_logs.las',32)
    w19_19A.Xloc = 437506
    w19_19A.Yloc = 6477887.5
    w19_19A.Inline = 10151
    w19_19A.Crossline = 2277
    w19_19A.KB = 54
    w19_19A.MD = mat[:,0]
    w19_19A.DT = mat[:,2]
    w19_19A.Rhob = mat[:,1]
    w19_19A.calc_Vp()
    
    # TD
    mat = w19_19A.read('well_files//volve//TD//HRS//15_9_19A_well_tie.txt',16)
    w19_19A.TD_time = mat[:,5]
    w19_19A.TD_TVD = mat[:,3]
    w19_19A.TD_MD = mat[:,0]
    w19_19A.TD_X = mat[:,8]
    w19_19A.TD_Y = mat[:,9]

    # Deviation file
    mat = w19_19A.read('well_files//volve//TD//HRS//15_9_19A_well_tie.txt',16)
    w19_19A.MD_dev = mat[:,0]
    w19_19A.Xdev = mat[:,8]
    w19_19A.Ydev = mat[:,9]

    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(w19_19A.MD,w19_19A.Vp)
    
    w19_19A.apply_TD()
    w19_19A_time = well()
    w19_19A_time = mfun.load_obj('well_files//volve//time//w19A')

    # debug
    zone = np.zeros(2)
    if (w19_19A.time [0] >= w19_19A_time.time[0] ):
        zone[0] = w19_19A.time[0]
    else:
        zone[0] = w19_19A_time.time[0]
        
    if (w19_19A.time [-1]<=w19_19A_time.time[-1] ):
        zone[1] = w19_19A_time.time[-1]
    else:
        zone[1] = w19_19A.time[-1]
    zone[0] = 2350
    zone[1] = 2515
    
    w19_19A.calc_Vp()
    ax2.plot(w19_19A.time,w19_19A.Vp)
    

    vp_time_zn,new_time_zn = mfun.segment_logs(zone,w19_19A.time,w19_19A.Vp)   
    org_vp_time_zn,org_time_zn = mfun.segment_logs(zone,w19_19A_time.time,w19_19A_time.Vp)
    
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(new_time_zn,vp_time_zn)
    ax2.plot(org_time_zn,org_vp_time_zn)  
    
    w19_19A.save_obj('w19_19A_jan2020',w19_19A)
        