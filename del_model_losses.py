import sys
from neicio.shake import ShakeGrid
from xml.dom import minidom
import copy
import numpy as np
from neicio.gmt import GMTGrid
from neicpager.nationpop import readNationPop
from neicpager.exposure import Exposure
from neicpager.fatality import lognormal
from neicio.esri import EsriGrid
from neicutil.text import commify
from xmlmodel import readXMLModel
from econstruct import readCountryGDP


dir = '/Users/sverros/Documents/Modules/Correlation_losses/loss_data/'

def initialize_loss(shakemap, voi, dist, mag):
    popfile =    dir+'lspop2012.flt'
    isofile =    dir+'isogrid.bil'
    growthfile = dir+'WPP2012_POP_F02_POPULATION_GROWTH_RATE.XLS'

    popgrid = EsriGrid(popfile)
    isogrid = EsriGrid(isofile)
    ratedict = readNationPop(growthfile)

    outgrid = GMTGrid()
    outgrid.geodict = shakemap.geodict.copy()
    outgrid.Attributes = shakemap.getAttributes()
    outgrid.griddata = MMI_conv2(shakemap.griddata.copy(), voi, dist, mag)
    
    expobj = Exposure(outgrid,popgrid,isogrid,ratedict,2012) 

    return expobj, outgrid

def model_econ(dat, shakemap, expob):
    shakemap_mod = copy.deepcopy(shakemap)
    expobj = copy.deepcopy(expob)
    data = copy.copy(dat)

    modelfile =  dir+'economy.xml'
    nationfile = dir+'isogrid.bil'
    GDPfile =    dir+'countrygdp.xml'
    CountryModels, ModelVersion = readXMLModel(modelfile, 'economic')
    EconGDP = readCountryGDP(GDPfile)

    event = shakemap.getAttributes()

    mmimin = list(np.arange(0.5,10.5))
    mmimax = list(np.arange(1.5,11.5))

    shakemap_mod.griddata = data
    pop = expobj.popGrid.getGeoDict()
    shakemap_mod.interpolateToGrid(pop, method = 'linear')
    expobj.shakeGrid = shakemap_mod
    results = expobj.getCountryExposure(np.array(zip(mmimin, mmimax)))

    DollarsByCountry = {}
    ExposureByCountry = {}

    Results = {}
    ccodes = results.keys()
    foundCountry = False
    dollars = 0
    eventyear = event['event']['event_timestamp'].year
    for ccode in ccodes:
        # deal with weird country codes that we may have created
        if ccode:
            model = CountryModels[ccode]
            if EconGDP.has_key(ccode):
                countrygdp = EconGDP[ccode]
                minyear = min(countrygdp['years'])
                maxyear = max(countrygdp['years'])
                if eventyear < minyear:
                    raise Exception,'Inappropriate event year %i: No GDP data available for %s before %i' % (eventyear,ccode,minyear)
                if eventyear < maxyear:
                    yearidx = countrygdp['years'].index(eventyear)
                    modelgdp = countrygdp['gdp'][yearidx]
                else:
                    modelgdp = model['pcgdp']
            else:
                modelgdp = model['pcgdp']
            countrydict = EconGDP[ccode]
            foundCountry = True
        else:
            continue
        myexp = copy.deepcopy(results[ccode])
        myexp[8]['exposure'] = myexp[8]['exposure'] + myexp[9]['exposure']
        myexp[9]['exposure'] = 0
        for i in range(0,len(myexp)):
            myexp[i]['exposure'] = myexp[i]['exposure'] * modelgdp  * model['alpha']
        DollarsByCountry[ccode] = lognormal(myexp,model['theta'],model['beta'])[0]
        ExposureByCountry[ccode] = copy.deepcopy(myexp)
        dollars = dollars + DollarsByCountry[ccode]

    if not foundCountry:
        Results['BestEstimate'] = 0
        Results['LowerBound'] = 0
        Results['UpperBound'] = 0
        dollars = 0
        return Results
        
    # someday we'll add in error analysis based on (hopefully just) the CDF, norm, and G value.
    # print 'Estimated dollar damage: %f' % (dollars)
    Results['BestEstimate'] = int(dollars)
    Results['LowerBound'] = int(dollars/10)
    Results['UpperBound'] = int(dollars*10)
        
    return Results

def model_loss(data, shakemap, expob):
    # code from mhearne
    shakemap_mod = copy.deepcopy(shakemap)
    expobj = copy.deepcopy(expob)
    mmimin = list(np.arange(0.5,10.5))
    mmimax = list(np.arange(1.5,11.5))

    shakemap_mod.griddata = data
    pop = expobj.popGrid.getGeoDict()
    shakemap_mod.interpolateToGrid(pop, method = 'linear')

    expobj.shakeGrid = shakemap_mod

    results = expobj.getCountryExposure(np.array(zip(mmimin,mmimax)))

    fatfile = dir+'fatality.xml'
    ccode_out = []
    fats_out = []

    for ccode,cexp in results.iteritems():
        theta,beta = getModel(fatfile,ccode)
        if theta is None:
            if ccode == None:
                ccode = ''
            fats = 0
        else:
            cexp[8]['exposure'] = cexp[8]['exposure'] + cexp[9]['exposure']
            cexp[9]['exposure'] = 0
            fats,fatrates = lognormal(cexp,theta,beta)
        ccode_out.append(ccode)
        fats_out.append(int(fats))
    return {'ccodes':ccode_out, 'fats':fats_out}

def getModel(fatfile,ccode):
    # Code from mhearne
    
    root = minidom.parse(fatfile)
    models = root.getElementsByTagName('model')
    theta = None
    beta = None
    for model in models:
        mcode = model.getAttribute('ccode')
        if mcode == ccode:
            theta = float(model.getAttribute('theta'))
            beta = float(model.getAttribute('beta'))
            break
    root.unlink()

    return (theta,beta)


def MMI_conv2(data, key, dists = None, mag = None):

    if (key == 'PGA')or(key=='pga'):
        c = { 'C1' :  1.78, 'C2' :  1.55, 'C3' : -1.60, 'C4' : 3.70, 'C5'   : -0.91,
              'C6' :  1.02, 'C7' : -0.17, 'T1' :  1.57, 'T2' : 4.22, 'SMMI' : 0.66,
              'SPGM' : 0.35 }
        c2 = { 'C1' :  1.71, 'C2' :  2.08, 'T1' :  0.14, 'T2' : 2.0 }
    elif (key=='PGV')or(key=='pgv'):
        c =  { 'C1' :  3.78, 'C2' :  1.47, 'C3' :  2.89, 'C4' : 3.16, 'C5'   :  0.90, 
               'C6' :  0.00, 'C7' : -0.18, 'T1' :  0.53, 'T2' : 4.56, 'SMMI' : 0.63,  
               'SPGM' : 0.38 }
        c2 =  { 'C1' :  4.62, 'C2' :  2.17, 'T1' : -1.21, 'T2' : 2.0 }
    elif (key=='psa03'):
        c = { 'C1' :  1.26, 'C2' :  1.69, 'C3' : -4.15, 'C4' : 4.14, 'C5'   : -1.05,
            'C6' :  0.60, 'C7' :  0.00, 'T1' :  2.21, 'T2' : 4.99, 'SMMI' : 0.82, 
            'SPGM' : 0.44 }
        c2 = { 'C1' :  1.15, 'C2' :  1.92, 'T1' :  0.44, 'T2' : 2.0 }
    elif (key=='psa10'):
        c = { 'C1' :  2.50, 'C2' :  1.51, 'C3' :  0.20, 'C4' : 2.90, 'C5'   :  2.27,
              'C6' : -0.49, 'C7' : -0.29, 'T1' :  1.65, 'T2' : 4.98, 'SMMI' : 0.75,  
              'SPGM' : 0.47 }
        c2 =  { 'C1' :  2.71, 'C2' :  2.17, 'T1' : -0.33, 'T2' : 2.0 }
    elif (key=='psa30'):
        c =  { 'C1' :  3.81, 'C2' :  1.17, 'C3' :  1.99, 'C4' : 3.01, 'C5'   :  1.91,
               'C6' : -0.57, 'C7' : -0.21, 'T1' :  0.99, 'T2' : 4.96, 'SMMI' : 0.89, 
               'SPGM' : 0.64 } 
        c2 =  { 'C1' :  7.35, 'C2' :  3.45, 'T1' : -1.55, 'T2' : 2.0 } 

    if 'PGV' not in key:
        units = 9.81
    else:
        units = 1.0

    if dists is not None and mag is not None:
        doresid = True
        ldd = np.log10(np.clip(dists, 10, 300))
        if mag < 3.0:
            mag = 3.0
        elif mag > 7.3:
            mag = 7.3
    else:
        doresid = False
        
    # Convert (for accelerations) from %g to cm/s^2
    # then take the log10
    lamps = np.log10(data * units)
    mmi = np.zeros_like(data)
    #
    # This is the MMI 1 to 2 range that is discussed in the paper but not 
    # specifically implemented
    #
    idx = lamps < c2['T1']
    mmi[idx] = c2['C1'] + c2['C2'] * lamps[idx]
    #
    # This is the lower segment of the bi-linear fit
    #
    idx = (lamps >= c2['T1']) & (lamps < c['T1'])
    mmi[idx] = c['C1'] + c['C2'] * lamps[idx]
    #
    # This is the upper segment of the bi-linear fit
    #
    idx = lamps >= c['T1']
    mmi[idx] = c['C3'] + c['C4'] * lamps[idx]

    if doresid:
        mmi += c['C5'] + c['C6'] * ldd + c['C7'] * mag

    mmi = np.clip(mmi, 1.0, 10.0)
    return mmi


#if __name__ == '__main__':
#    smfile = sys.argv[1]
#    shakemap = ShakeGrid(smfile, variable = 'PGA')
#    data = shakemap.griddata.copy()
#
#    expobj, outgrid = initialize_loss(shakemap, 'PGA')
#
#    mmi_data = MMI_conv(data, 'PGA')
#    loss = model_loss(mmi_data, shakemap, expobj)
#    econ_est = model_econ(mmi_data, shakemap, expobj)
#    ccodes, fats = [], []
#    keep = [[int(loss['fats'][i]), loss['ccodes'][i]] for i in range(0, np.size(loss['fats'])) if loss['fats'][i] != '0']
#    for xx in range(0, np.size(keep)/2):
#        ccodes.append(keep[xx][1])
#        fats.append(keep[xx][0])
#    ubounds = econ_est['UpperBound']
#    lbounds = econ_est['LowerBound']
#    bestest = econ_est['BestEstimate']

