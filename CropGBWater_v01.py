
######################################################
## !/usr/bin/python                                 ##
## name : CropGBWater_v01.py                        ##
## Developed by: Mesfin M. Mekonnen                 ##
## Reads Aquacrop input  (PR and ETo)               ##
######################################################

"""
Things to do:
        - put KC, crop growing stages and soil moisture into fuction
        - put main model into function
        - write a module to be imported for writing monthly and daily values
        - error handling
"""
# import system modules
import sys
import os, os.path
import datetime, time
from datetime import date, datetime, timedelta
import pandas as pd
import re
import numpy as np

# import the cropwat modules
# import file_monthly2


path = os.getcwd()
# outpath = path + "\\perGrid" 
fpath = path + "\\Inputs2"
ofpath = path + "\\Outputs2"


# version of the model 
__version__ = '$Revision: 2023-09-30$'

# the time definition 
def localtimestrf(frmtopt="datetime"):
    ''' Return the local date-time as a formatted string. 
        frmtopt="datetime" (default)   Return a complete date-time string.
        frmtopt="timeonly"             Return Hour:Minute:Seconds only.
        
    '''
    if frmtopt == "datetime":
        return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    elif frmtopt == "timeonly":
        return "Time: " + time.strftime("%H:%M:%S", time.localtime())


print ("\n***** Global water footprint assessment model *****")  # no name yet given!!
print ("            " + __version__)
print ("            " + localtimestrf())
print ("      ==========================================")



def read_crop_parameters(df_crop):

    # latitude, longitude and ID of the grid

    X_Grid = df_crop['X_CORD']     # X-grid (latitude) in degree
    Y_Grid = df_crop['Y_CORD']     # Y-grid (longitude)
    ID = int(float(df_crop['ID'])) # gird ID
    

    # crop parameters
    KCini = float(df_crop['Kc_ini'])      # initial crop coefficient Kc
    KCmid = float(df_crop['Kc_mid'])      # mid season crop coefficient Kc
    KCend = float(df_crop['Kc_end'])      # final crop coefficient Kc
    Lini = int(df_crop['L_ini'])         # length of crop initial stage
    Ldev = int(df_crop['L_dev'])         # length of crop development stage
    Lmid = int(df_crop['L_mid'])         # length of crop mid stage
    Llate = int(df_crop['L_late'])        # length of crop late stage

    pmonth = int(df_crop['Month'])       # planting month 
    pday = int(df_crop['Day'])       # planting date 

      
    # for rainfed field
    
    Zrmin_rf = float(df_crop['Rdmin'])
    Zrmax_rf = float(df_crop['Rdmax_rf'])

    # for irrigated field
    
    Zrmin_irr = float(df_crop['Rdmin'])    # initial rooting depth, for annual crop Zrmin = 0.2 and for perinnial Zrmin = Zrmax
    Zrmax_irr = float(df_crop['Rdmax_irr'])    # float(cfields[15])   or 1.25 (average    # maximum rooting depth - attained at plant development stage
    
    # for irrigated fieled but rainfed
    Zrmin_irf = float(df_crop['Rdmin'])       # initial rooting depth, for annual crop Zrmin = 0.2 and for perinnial Zrmin = Zrmax
    Zrmax_irf = float(df_crop['Rdmax_irr'])      # or use 1.6 as average -- maximum rooting depth - attained at plant development stage

    Pstd = float(df_crop['P'])        # standard water depletion fraction for no stress for ET = 5 mm/day
    Ky = float(df_crop['Ky'])          # yield response factor - to water stress
    Ym = float(df_crop['Ym'])                # maximum yield
    
    N_ferti = float(df_crop['Nitr_kg_ha'])         # fertilizer application rate 

 
    return X_Grid, Y_Grid, ID, KCini, KCmid, KCend, Lini, Ldev, \
           Lmid, Llate, pmonth, pday, Zrmin_rf, Zrmax_rf, Zrmin_irr, \
           Zrmax_irr, Zrmin_irf, Zrmax_irf, Pstd, Ky, Ym, N_ferti




def rooting_depth(J, Jmid, Zrmin, Zrmax, KC, KCini, KCmid):
    # calculation of the root depth. following pp. 279 (Annex 8) of FAO-56 
    
    if J <= Jmid:                                                       
            Zri = Zrmin + (Zrmax - Zrmin)*((KC - KCini)/(KCmid - KCini))
    else: 
            Zri = Zrmax

    return Zri


def runoff(PR, CN):
    
    

    # Runoff routing replaced with USDA SCS curve number - adopted from Aquacrop manual v
    # see Chapter 3, page 3-48 eqn (3.7e)
    """
    RO - amount of water lost by surface runoff [mm]
    PR - rainfall amount [mm]
    CN - curve number
    S - potential maximum soil water retention [mm]
    Ia = initial abstraction [mm] or the amount of water that can infiltrate before runoff occurs
    """
    S = 254*(100/CN-1)      
    Ia = 0.05*S 
    
    RO = ((PR - Ia)**2)/(PR + S - Ia)
           
    
    
    return RO

def KS_value (Dri, RAW, St, Pi, TAW):

    #Soil water depletion coefficient as function of soil water depletion (Dri), soil moisture(St),
    # total available water (TAW)
    
    if Dri > RAW:
        KS = max(0,((St)/((1-Pi)*TAW)))
    else:
        KS = 1

    return KS

def deep_percolation(Dri, PR, RO, ETc, IRR):

    # deep percolation (DP) calculation
    # when the soil is at field capacity Dri = 0 so DP follows this equation.
    # The Dr in the equation is Dr,i-1

    if Dri > 0:                                                     
        DP = 0
    else:                                                           
        DP = max(((PR - RO) + IRR - ETc - Dri),0)  

    return DP

def greyWU (N_ferti, Cmax, LR):
    
    # the grey water use (m3/ha)
    # LR -> leaching rate assumed 10%
    # Cmax -> maximum permissiable limit NO3-N [kg/m3] -
    # this corresponds to 10 mg/l (NO3-N) emission standard
    # L -> nitrogen leaching to the ground [kg/ha]
    
    LR = 0.1            
    Cmax = 0.01         
   
    L = LR * N_ferti     
 
    CWU_gy= (L/Cmax)

    return CWU_gy

def actual_yield(Ym, Ky, ETa, ETc):

    # yield reduction and actual yield calculation
    # the lower limit is set at 20% of Ym to avoid zero values - if crop earea
    # exists Y should be different from zero
    Ya = max((Ym*(1-Ky*(1-(ETa/ETc)))),0.2*Ym)

    return Ya



croplist = ['Rice', 'Wheat']

firstY = 2000
lastY = 2014
pYears = range(firstY,lastY)


#df header
wf_header=['CWR', 'CWU_gn_rf', 'CWU_gn_irr', 'CWU_gb', 'CWU_bl', 'CWU_gy', 'IRReq', 'Ya_rf', 'Ya_irr']

#header for grey WF (mm) and yield (t/ha) to be added in the monthly output
y_header = ['CWU_gy', 'Ya_rf', 'Ya_irr']




for cropi in croplist:
    
    
    # cropFile1 = cropfile
    # (ShortName, Extension) = os.path.splitext(cropFile1)
    gridIDs=[]


    for filename in os.listdir(fpath):
        
        # check files and collect the grid id
        if filename.startswith(cropi):
            #print (filename)
            #get the id and add it to the list
            gridIDcr = (os.path.splitext(os.path.basename(filename))[0])
            gridID = re.split("[_]", gridIDcr)[-1]
            gridIDs.append(gridID)
    
    
    
    
    
    # to collect the result per crop
    a = (len(gridIDs), 9)
    wf_result = np.zeros(a)
    
    c = (len(gridIDs), 3)
    
    y_result = np.zeros(c)
    
       
    # #creat daily list - this should have been in the year loop and use pYear and hYear
    #but b/c difficulty to create header consistent with the array, we assumed a leap year
    
    firstDate = date(2000,1,1)
    lastDate = date(2000,12,31)
    d_header = pd.date_range(firstDate, lastDate).strftime("%Y-%m-%d").tolist()
    
    b = (len(gridIDs), 366)  #accounting for leap year
    
    cwr_d = np.zeros(b); cwu_gn_rf_d = np.zeros(b); cwu_gn_ir_d = np.zeros(b);
    cwu_gb_d= np.zeros(b); cwu_bl_d = np.zeros(b); IRReq_d = np.zeros(b);
    
    i=-1
    
    #for each listed grid, do the green/blue separation        
    for gridID in gridIDs:
        
        i+=1
        
        
        # the input file with daily data
        
        crop_file = os.path.join(fpath, cropi + '_' + gridID + '.csv')
        soil_file = os.path.join(fpath, 'SOIL_' + gridID + '.csv')
        prcp_file = os.path.join(fpath, 'PRE_' + gridID + '.PLU')
        eto_file = os.path.join(fpath, 'PET_' + gridID + '.ETO')
  
      
        # crop and soil file
        df_crop = pd.read_csv(crop_file, delimiter='\s*,\s*')
        df_soil = pd.read_csv(soil_file, delimiter='\s*,\s*', usecols=['TAWC', 'CN'])
        # for precp and eto, Skip 8 rows from top and read the data to dataframe
        df_prcp = pd.read_csv(prcp_file, delimiter='\s+|\t+', skiprows=[i for i in range(0,8)], usecols=[0])
        df_eto = pd.read_csv(eto_file, delimiter='\s+|\t+', skiprows=[i for i in range(0,8)], usecols=[0])


        
        X_Grid, Y_Grid, ID, KCini, KCmid, KCend, \
                Lini, Ldev, Lmid, Llate, pmonth, pday, \
                Zrmin_rf, Zrmax_rf, Zrmin_irr, \
                Zrmax_irr, Zrmin_irf, Zrmax_irf, \
                Pstd, Ky, Ym, N_ferti = read_crop_parameters(df_crop)
    
        
    
        #soil parameters (field 3- for FAO and field 4 for ISRIC data)
        TAWC = df_soil['TAWC'][0]          # maximum soil moisture capacity in mm/m - soil water capacity 1000*(FC -WP)?
        CN = df_soil['CN'][0]
    
    
    
        ########
        
    
        j = i

        
        for pYear in pYears:
            
            
            
            # parameters initialization
        
            KS_rf = 1
            KS_irr = 1
            KS_irf = 1
        
            # for rainfed
            TAW_rf = TAWC * Zrmin_rf
            Dri_rf = TAWC * Zrmin_rf * Pstd                # ??? initial root zone depletion => Zro is the initial root depth = 0.2 (FAO-56, pp279);
                                                           #and RAW raidly available soil moisture at the begining (preseason - 5 day averag)
            # for irrigated
            TAW_irr = TAWC * Zrmin_irr
            Dri_irr = TAWC * Zrmin_irr * Pstd
        
            IRR = Dri_irr
        
            # for irrigated field but only rain-fed
            TAW_irf = TAWC * Zrmin_irr
            Dri_irf = TAWC * Zrmin_irr * Pstd
        
            
        
            St_rf = TAW_rf - Dri_rf
            St_irf = TAW_irf - Dri_irf
            St_irr = TAW_irr - Dri_irr
                
        
            Total = 0
            ETa_rf_total = 0       # total actual evapotranspiration
            ETa_irr_total = 0
            ETa_irf_total = 0
            ETc_total = 0       # total maximum evapotranspiration
            CWR = 0             # crop water requirement in m3/ha CWR = ETc_total * 10 
            CWU_gn_rf = 0            # the green component of crop water use in m3/ha CWUg = ETa_total * 10
            CWU_gn_irr = 0            # the green component of crop water use in m3/ha CWUg = ETa_total * 10
            CWU_gb = 0            # the green component of crop water use in m3/ha CWUg = ETa_total * 10
                                # under the irrigation case this value represent both blue and green which later have to be separated
            IRReq = 0           # irrigation water requirement in m3/ha IRReq = IRR * 10 
                                
            RO = 0              # run-off
            
            IRR = 0             # irrigation at each time step
            IRR_total = 0       # total irrigation over the period
            IRRnet = 0
        
             
            #identifying the row number of the climate data, based on the planting and start of the climate data
            #planting date
            #pYear = 2004
            
            plantingDate = date(pYear,pmonth,pday)
            
            
            climStartDate = date(1989,1,1)
            
            date_to_planting = (plantingDate - climStartDate)
            Jplant = date_to_planting.days 
            
            #total lenght of growing dates
            gDates = Lini + Ldev + Lmid + Llate
            #delta = datetime.timedelta(days=gDates)
            delta = timedelta(days=gDates)
            
            #harvest year - we use this for printingout the result
            hYear = (plantingDate + delta).year
            
            # annual ouput 
            # Writing out the annual results to the files
            # all file could be put in one csv file by adding the harvest year in the heading - no time for that now
            
            # ofile = os.path.join(ofpath, cropi + '%s.csv'%hYear)
            
            # output file to hold monthly and annual results
            ofile_all = os.path.join(ofpath, cropi + '_all%s.csv'%hYear)
        
            #computed dates for the stages
            Jdev = Jplant + Lini
            Jmid = Jdev + Ldev
            Jlate = Jmid + Lmid
            Jharv = Jlate + Llate
            
            
                      
            doy = plantingDate.timetuple().tm_yday - 1        
            #print(doy)
            
            t = -1
            #for loop over the grawing period - from planting to harvest
            #for J in range (Jplant, Jharv + 1, 1):  # Jharv + 1 - so the last value of Jharv can be read!!
            for J in range (Jplant, Jharv + 1, 1):  # Jharv + 1 - so the last value of Jharv can be read!!
                #print ('test', J, Jplant, Jharv)
                
                
                               
                #for using day of year as index
                t += 1
                
                
                doy += 1
                
                #print (doy)
                 
        
                # daily values of KC, based on eqn (66) pp. 132, Fig 34 pp. 126 of FAO-56
                if J <= Jdev:
                    KC = KCini
                elif Jdev < J < Jmid:
                    KC = KCini + (KCmid - KCini)*(J - Jdev)/Ldev
                elif Jmid <= J <= Jlate:
                    KC = KCmid
                elif Jlate < J < Jharv:
                    KC = KCmid + (KCend - KCmid)*(J - Jlate)/Llate
                else:
                    KC = KCini
        
                
                # calculation of the root depth. following pp. 279 (Annex 8) of FAO-56 
        
                Zri_rf = rooting_depth(J, Jmid, Zrmin_rf, Zrmax_rf, KC, KCini, KCmid)
                Zri_irr = rooting_depth(J, Jmid, Zrmin_irr, Zrmax_irr, KC, KCini, KCmid)
                
                # Daily calculation of different parameters
                PR = float(df_prcp.iloc[J])                # the precipitatin (mm/day) data from CRU has a factor of 10 so check if that factor was accounted for
                ETo = float(df_eto.iloc[J])               # reference evapotranspiration (mm/day) - got from FAO
                
                ETc = ETo * KC                        # crop evapo-transpiration (mm/day)
                Pi = min((max((Pstd + 0.04*(5 - ETc)),0.1)),0.8)       # depletion fraction 
        
                ''' soil water content calculation'''
                # rain-fed
                TAW_rf = TAWC * Zri_rf      # total available soil water in the root zone (mm); TAWC (mm/m) - total soil water capacity
                RAW_rf = TAW_rf * Pi        # raidly available soil water in the root zone (mm)
                St_rf = TAW_rf - Dri_rf	    # the actual available soil water content (mm). 
        
                # irrigated
                TAW_irr = TAWC * Zri_irr                                                
                RAW_irr = TAW_irr * Pi                                                  
                St_irr = TAW_irr - Dri_irr						     
        
                # irrigated field but rain-fed (no irrigation water)
                TAW_irf = TAWC * Zri_irr                                                
                RAW_irf = TAW_irf * Pi     
                St_irf = TAW_irf - Dri_irf						     
        
        
                # following Allen et al(1998): (Doll and Siebert (2007)- also have similar approach)
                # the model will run twice once with no irrigation and the other with irrigation  
                # final calculation of green and blue water component would be based on area coverage
                # currently no twice running - all the three options are run at once
        
                if Dri_irr >= RAW_irr:          # irrigation is applied when Dri >= RAW
                    IRR = TAW_irr - Dri_irr     #and in order to avoid deep percolation IRR <= Dri        
                                                        
                else:
                    IRR = 0
                
          
                # runoff (mm/day)
                RO = runoff(PR, CN)
        
                '''~~~~~~~~~~~~~~ calculation of Ks and actual evapotranspiration ~~~~~~~~~~~~~'''
                # calculation of Ks, actual evapotranspiration, soil water depletion on a daily basis
                # KS - the water stress coefficient calculated 
                # ETa - the the adjusted crop ET (mm/day)
                # Dri - soil water depletion (mm/day)- at end of the current period
                #       it can be 0 <= Dri <= TAW; pp 170 and Annex 8, pp 277
                # DP - deep percolation (mm/day)
        
                #~~~~~~ rain-fed ~~~~~~~
         
                KS_rf = KS_value (Dri_rf, RAW_rf, St_rf, Pi, TAW_rf)
                
                ETa_rf = ETc * KS_rf    
        
        
                # deep percolation (mm/day) 
                DP_rf = deep_percolation(Dri_rf, PR, RO, ETc, IRR=0)
                
                # soil water depletion (mm/day)- at end of the current period
                Dri_rf = min((max((Dri_rf - PR + RO + ETa_rf + DP_rf),0)),TAW_rf)        
                                                                                
        
                #~~~~~~ irrigated ~~~~~~
         
                KS_irr = KS_value (Dri_irr, RAW_irr, St_irr, Pi, TAW_irr)            
                
                ETa_irr = ETc * KS_irr  
        
                # deep percolation (mm/day)
                DP_irr = deep_percolation(Dri_irr, PR, RO, IRR, ETc)
                
                # soil water depletion (mm/day) - at end of the current period
                Dri_irr = min((max((Dri_irr - PR + RO - IRR + ETa_irr + DP_irr),0)),TAW_irr)    
               
        
                #~~~~~~ irrigated field  but rain-fed ~~~~~~
          
                KS_irf = KS_value (Dri_irf, RAW_irf, St_irf, Pi, TAW_irf)
                
                ETa_irf = ETc * KS_irf  
        
        
                # deep percolation
                DP_irf = deep_percolation(Dri_irf, PR, RO, ETc, IRR=0)
                
                # soil water depletion (mm/day) - at end of the current period
                Dri_irf = min((max((Dri_irf - PR + RO + ETa_irf + DP_irf),0)),TAW_irf)  
                
                '''~~~~~~~~~~~~~~ collecting daily values in array for aggregating at monthly level ~~~~~~~~~~~~~'''
                #daily blue WF
                ETa_bl = max((ETa_irr - ETa_irf),0)
                
                #j = i-1
                
                cwr_d[j][doy] = ETc; cwu_gn_rf_d[j][doy] = ETa_rf; cwu_gn_ir_d[j][doy] = ETa_irf;
                cwu_gb_d[j][doy] = ETa_irr; cwu_bl_d[j][doy] = ETa_bl; IRReq_d[j][doy] = IRR;
                
                
                """ ***** check j and id index"""
        
                # aggregation of the evapotranspiration over the growing period
                # values are in mm per growing period
                ETc_total += ETc
                ETa_rf_total += ETa_rf
                ETa_irr_total += ETa_irr
                ETa_irf_total += ETa_irf
                IRR_total += IRR
             
                #print (j, t, doy, round(ETc,2))
            
            #print(cwr_d[j,doy])
                        
            
            '''~~~~~~~~~~~~~~ aggregating daily values at monthly level ~~~~~~~~~~~~~'''
            #convert the array to panda dataframe and add the date as hearder
            df_cwr = pd.DataFrame(cwr_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID'); 
            df_cwu_gn_rf = pd.DataFrame(cwu_gn_rf_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID');
            df_cwu_gn_ir = pd.DataFrame(cwu_gn_ir_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID'); 
            df_cwu_gb = pd.DataFrame(cwu_gb_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID');
            df_cwu_bl = pd.DataFrame(cwu_bl_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID'); 
            df_IRReq = pd.DataFrame(IRReq_d.round(2), columns=d_header, index = gridIDs).rename_axis('ID');
            
            #print (df_cwr)           
            #Transpose the df so dates are on the rows as index and the values on columns
            #this steps could be heavy on the memory of large number of grids are run over longer periods- will find simple approach
            df_cwr = df_cwr.transpose(); df_cwu_gn_rf = df_cwu_gn_rf.transpose();
            df_cwu_gn_ir = df_cwu_gn_ir.transpose(); df_cwu_gb = df_cwu_gb.transpose();
            df_cwu_bl = df_cwu_bl.transpose(); df_IRReq = df_IRReq.transpose();
            
            #convert the date strings to datetime strings
            df_cwr.index = pd.to_datetime(df_cwr.index); df_cwu_gn_rf.index = pd.to_datetime(df_cwu_gn_rf.index);
            df_cwu_gn_ir.index = pd.to_datetime(df_cwu_gn_ir.index); df_cwu_gb.index = pd.to_datetime(df_cwu_gb.index);
            df_cwu_bl.index = pd.to_datetime(df_cwu_bl.index); df_IRReq.index = pd.to_datetime(df_IRReq.index);
            
            #aggregate at monthly time step
            df_cwr_m = df_cwr.resample('M').sum(); df_cwu_gn_rf_m = df_cwu_gn_rf.resample('M').sum();
            df_cwu_gn_ir_m = df_cwu_gn_ir.resample('M').sum(); df_cwu_gb_m = df_cwu_gb.resample('M').sum();
            df_cwu_bl_m = df_cwu_bl.resample('M').sum(); df_IRReq_m = df_IRReq.resample('M').sum();
            
            
            #aggregate annually
            df_cwr_y = df_cwr.resample('Y').sum(); df_cwu_gn_rf_y = df_cwu_gn_rf.resample('Y').sum();
            df_cwu_gn_ir_y = df_cwu_gn_ir.resample('Y').sum(); df_cwu_gb_y = df_cwu_gb.resample('Y').sum();
            df_cwu_bl_y = df_cwu_bl.resample('Y').sum(); df_IRReq_y = df_IRReq.resample('Y').sum();
            
           
            # transpose back the df to have months on the column
            df_cwr_m = df_cwr_m.transpose(); df_cwu_gn_rf_m = df_cwu_gn_rf_m.transpose();
            df_cwu_gn_ir_m = df_cwu_gn_ir_m.transpose(); df_cwu_gb_m = df_cwu_gb_m.transpose();
            df_cwu_bl_m = df_cwu_bl_m.transpose(); df_IRReq_m = df_IRReq_m.transpose();
            
            # transpose the annual value
            df_cwr_y = df_cwr_y.transpose(); df_cwu_gn_rf_y = df_cwu_gn_rf_y.transpose();
            df_cwu_gn_ir_y = df_cwu_gn_ir_y.transpose(); df_cwu_gb_y = df_cwu_gb_y.transpose();
            df_cwu_bl_y = df_cwu_bl_y.transpose(); df_IRReq_y = df_IRReq_y.transpose();
                        
                    
            
            # #creat monthly data header replace the headers, which in y-m-d format 
            #with the proper column header
            ml = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'];
            
            cwr_hd = ['cwr' + s for s in ml]; cwu_gn_rf_hd = ['cwu_gn_rf' + s for s in ml];
            cwu_gn_ir_hd = ['cwu_gn_ir' + s for s in ml]; cwu_gb_hd = ['cwu_gb' + s for s in ml];
            cwu_bl_hd = ['cwu_bl' + s for s in ml]; IRReq_hd = ['IRReq' + s for s in ml];
            
            cwr_y_hd = ['cwr_ann']; cwu_gn_rf_y_hd = ['cwu_gn_rf_ann'];
            cwu_gn_ir_y_hd = ['cwu_gn_ir_ann']; cwu_gb_y_hd = ['cwu_gb_ann'];
            cwu_bl_y_hd = ['cwu_bl_ann']; IRReq_y_hd = ['IRReq_ann'];
            
            #df_cwr_m.columns = cwr_hd
            
            df_cwr_m.columns = cwr_hd; df_cwu_gn_rf_m.columns = cwu_gn_rf_hd; df_cwu_gn_ir_m.columns = cwu_gn_ir_hd;
            df_cwu_gb_m.columns = cwu_gb_hd; df_cwu_bl_m.columns = cwu_bl_hd; df_IRReq_m.columns = IRReq_hd;
            
            df_cwr_y.columns = cwr_y_hd; df_cwu_gn_rf_y.columns = cwu_gn_rf_y_hd; df_cwu_gn_ir_y.columns = cwu_gn_ir_y_hd;
            df_cwu_gb_y.columns = cwu_gb_y_hd; df_cwu_bl_y.columns = cwu_bl_y_hd; df_IRReq_y.columns = IRReq_y_hd;
            
            
            
            '''~~~~~~~~~~~~~~ monthly ~~~~~~~~~~~~~'''
            # the grey water use in m3/ha
            
            CWU_gy = greyWU (N_ferti, Cmax = 0.01, LR = 0.1)
            
            # to be consisten with the other monthly ouput, which are in mm, 
            #we need to convert the annual grey WF into mm
            
            cwu_gy_mm = CWU_gy / 10.0
            
            # yield reduction and actual yield calculation
            # values in ton/ha
        
            Ya_rf = actual_yield(Ym, Ky, ETa_rf_total, ETc_total)
            Ya_irr = actual_yield(Ym, Ky, ETa_irr_total, ETc_total)
            
            # for the monthly files
            y_result[j][0]=cwu_gy_mm; y_result[j][1]=Ya_rf; y_result[j][2]=Ya_irr ;
            
            df_y = pd.DataFrame(y_result.round(2), columns=['cwu_gy', 'Y_rf', 'Y_irr'], index = gridIDs).rename_axis('ID');
            
            # merge all monthly value dataframes
            df_result_full = pd.concat([df_cwr_m, df_cwr_y, df_cwu_gn_rf_m, df_cwu_gn_rf_y, df_cwu_gn_ir_m, 
                                        df_cwu_gn_ir_y, df_cwu_gb_m, df_cwu_gb_y, df_cwu_bl_m, df_cwu_bl_y, df_IRReq_m, df_y], axis=1);
            
            # #creat a header text to explain content and unit
            # detail_header = ['The file contains monthly CWR, green and blue CWU in mm,',
            #                    ' and annual grey CWU in mm, irrigated and rainfed yield in t/ha ']
            # #df_result_full.columns = pd.MultiIndex.from_tuples(zip(detail_header, df_result_full.columns))
            
            #write to file
            
            df_result_full.to_csv(ofile_all)
            
            
            """ annual values are already in the monthly or the full result"""
            """
            # calculation of the crop water requirement and water use over the growing period in m3/ha
            # the factor 10 is used to convert the evapotranspiration (mm) into m3/ha
            
            CWR = 10 * ETc_total
            IRReq = 10 * IRR_total
            CWU_gn_rf = 10 * ETa_rf_total    
            CWU_gn_irr = 10 * ETa_irf_total
            CWU_gb = 10 * ETa_irr_total
        
            if CWU_gb > CWU_gn_irr:
                CWU_bl= CWU_gb - CWU_gn_irr
            else:
                CWU_bl = 0   
           
        
            #j = i-1
            #for d, k in zip(t, ind):
            wf_result[j][0] = CWR; wf_result[j][1]=CWU_gn_rf;
            wf_result[j][2]=CWU_gn_irr; wf_result[j][3]=CWU_gb; wf_result[j][4]=CWU_bl;
            wf_result[j][5]=CWU_gy; wf_result[j][6]=IRReq; wf_result[j][7]=Ya_rf;
            wf_result[j][8]=Ya_irr ;
    
            
            #convert the list to dataframe
            #df_wf = pd.DataFrame(wf_reslt.round(2), columns = wf_header, index=gridIDs)
            df_wf = pd.DataFrame(wf_result.round(2), columns = wf_header, index=gridIDs)
            
            df_wf.index.name = "ID"
  
            df_wf.to_csv(ofile)
            """
            
            
    
    print ('         completed running crop: %s file' % (cropi) )

#print ' The file = %s' % (cropFile1), CropPeriod, Zrmin_rf
print ("      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

print ("            %s\n" % localtimestrf())

