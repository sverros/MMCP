#!/usr/bin/env python

from xml.dom.minidom import parse

def readCountryGDP(gdpfile):
    dom = parse(gdpfile)
    root = dom.firstChild
    countrylist = root.getElementsByTagName('country')
    countrydict = {}
    for country in countrylist:
        ccode = country.getAttribute('ccode')
        yearstrlist = country.getAttribute('years').split(',')
        years = [int(yearstr) for yearstr in yearstrlist]
        gdpstrlist = country.getAttribute('pcgdp').split(',')
        gdp = [float(gdpstr) for gdpstr in gdpstrlist]
        cdict = {}
        cdict['years'] = years
        cdict['gdp'] = gdp
        countrydict[ccode] = cdict
    dom.unlink()
    return countrydict

def readEconData(xmlfile):
    dom = parse(xmlfile)
    root = dom.firstChild
    modellist = root.getElementsByTagName('model')
    modeldict = {}
    for model in modellist:
        tdict = {}
        ccode = model.getAttribute('ccode')
        tdict['theta'] = float(model.getAttribute('theta'))
        tdict['beta'] = float(model.getAttribute('beta'))
        tdict['G'] = float(model.getAttribute('g'))
        tdict['alpha'] = float(model.getAttribute('alpha'))
        try:
            tdict['population'] = float(model.getAttribute('population'))
        except:
            pass
            #print ccode
        try:
            tdict['gdp'] = float(model.getAttribute('pcgdp'))
        except:
            print 'Failed to get gdp from %s' % ccode
        modeldict[ccode] = tdict
        
    dom.unlink()
    return modeldict
