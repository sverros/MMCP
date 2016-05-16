#!/usr/bin/env python

#stdlib imports
from xml.dom.minidom import parse
import sys

#third party imports
import numpy

class Attribute:
    pass

class SemiData:
    inventory = Attribute()
    index = Attribute()
    workforce = Attribute()
    collapse = Attribute()
    casualty = Attribute()
    version = Attribute()
    pass

def readSemiXML(xmlfile):
    semidata = SemiData()
    semidata.index.ccodes = []
    semidata.index.btypes = []
    semidata.index.bdescs = []
    semidata.inventory.zlabels = []
    semidata.casualty.zlabels = []
    semidata.collapse.zlabels = []
    semidata.workforce.xlabels = []
        
    dom = parse(xmlfile)
    root = dom.getElementsByTagName('semi_empirical')[0]

    countries = root.getElementsByTagName('country')
    for country in countries:
        ccode = str(country.getAttribute('code'))
        semidata.index.ccodes.append(ccode)
    
    buildings = root.getElementsByTagName('global_building')
    for building in buildings:
        semidata.index.btypes.append(str(building.getAttribute('code')))
        semidata.index.bdescs.append(str(building.getAttribute('desc')))

    regions = root.getElementsByTagName('density_region')
    for region in regions:
        rname = str(region.getAttribute('region'))
        semidata.inventory.zlabels.append(rname)

    ctimes = root.getElementsByTagName('casualty_time')
    for ctime in ctimes:
        timeofday = str(ctime.getAttribute('time'))
        semidata.casualty.zlabels.append(timeofday)

    mmivalues = root.getElementsByTagName('collapse_mmi')
    for mmivalue in mmivalues:
        mmi = float(mmivalue.getAttribute('mmi'))
        semidata.collapse.zlabels.append(mmi)

    wsectors = root.getElementsByTagName('workforce_sector')
    for wsector in wsectors:
        sector = str(wsector.getAttribute('sector'))
        semidata.workforce.xlabels.append(sector)

    #All of the supplementary lists above (labels, building codes, etc.) need to be numpy arrays
    semidata.index.ccodes = numpy.array(semidata.index.ccodes)
    semidata.index.btypes = numpy.array(semidata.index.btypes)
    semidata.index.bdescs = numpy.array(semidata.index.bdescs)
    semidata.inventory.zlabels = numpy.array(semidata.inventory.zlabels)
    semidata.casualty.zlabels = numpy.array(semidata.casualty.zlabels)
    semidata.collapse.zlabels = numpy.array(semidata.collapse.zlabels)
    semidata.workforce.xlabels = numpy.array(semidata.workforce.xlabels)

    ncountries = len(semidata.index.ccodes)
    nbtypes = len(semidata.index.btypes)
    nregions = len(semidata.inventory.zlabels)
    ntimes = len(semidata.casualty.zlabels)
    nmmi = len(semidata.collapse.zlabels)
    nsectors = len(semidata.workforce.xlabels)

    semidata.inventory.data = numpy.zeros((ncountries,nbtypes,nregions))
    semidata.casualty.data = numpy.zeros((ncountries,nbtypes,ntimes))
    semidata.collapse.data = numpy.zeros((ncountries,nbtypes,nmmi))
    semidata.workforce.data = numpy.zeros((ncountries,nsectors))
        
    countries = root.getElementsByTagName('country')
    for country in countries:
        ccode = country.getAttribute('code')
        cidx = numpy.where(semidata.index.ccodes == ccode)
        inventory = country.getElementsByTagName('inventory')[0]
        buildings = inventory.getElementsByTagName('building')
        for building in buildings:
            bcode = building.getAttribute('code')
            bidx = numpy.where(semidata.index.btypes == bcode)
            for density_label in semidata.inventory.zlabels:
                density = float(building.getAttribute(density_label))
                didx = numpy.where(semidata.inventory.zlabels == density_label)
                semidata.inventory.data[cidx,bidx,didx] = density
            collapse = building.getAttribute('collapse').split(',')
            collapse = [float(c) for c in collapse]
            semidata.collapse.data[cidx,bidx,:] = numpy.array(collapse)
            for ctime in semidata.casualty.zlabels:
                tidx = numpy.where(semidata.casualty.zlabels == ctime)
                value = float(building.getAttribute(ctime))
                semidata.casualty.data[cidx,bidx,tidx] = value

        #Now read in the workforce data
        workforce = country.getElementsByTagName('workforce')[0]
        for sector in semidata.workforce.xlabels:
            sidx = numpy.where(semidata.workforce.xlabels == sector)
            value = float(workforce.getAttribute(sector))
            semidata.workforce.data[cidx,sidx] = value

    return semidata

def readXMLModel(xmlfile,modeltype):
    dom = parse(xmlfile)
    root = dom.getElementsByTagName('models')[0]
    vstr = root.getAttribute('vstr')
    mtype = root.getAttribute('type')
    if mtype != modeltype:
        raise Exception,'Wrong model type in %s (%s requested, %s found)' % (xmlfile,modeltype,mtype)
    models = root.getElementsByTagName('model')
    countrymodels = {}
    for model in models:
        cdict = {}
        atts = model.attributes
        natts = len(atts)
        for i in range(0,natts):
            node = atts.item(i)
            key = node.name
            value = node.nodeValue
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
            cdict[key] = value
        if not cdict.has_key('ccode'):
            raise Exception,'One of the models in %s is missing "ccode" key.' % xmlfile
        countrymodels[cdict['ccode']] = cdict
    dom.unlink()
    return (countrymodels,vstr)

if __name__ == '__main__':
    semidata = readSemiXML(sys.argv[1])
