import os.path
import urllib2
import urlparse
from xml.dom.minidom import parseString
import sys
from StringIO import StringIO
import argparse
import zipfile
import json
import datetime
from time import strptime


def getEventInfo(gridurl):
    gridfh = urllib2.urlopen(gridurl)
    gdata = gridfh.read()
    gridfh.close()

    gdata = gdata[0:gdata.find('<grid_data>')] + '</shakemap_grid>'
    xdom = parseString(gdata)
    root = xdom.getElementsByTagName('shakemap_grid')[0]
    infodict = {}
    infodict['id'] = root.getAttribute('event_id')
    event = root.getElementsByTagName('event')[0]
    infodict['lat'] = float(event.getAttribute('lat'))
    infodict['lon'] = float(event.getAttribute('lon'))
    infodict['depth'] = float(event.getAttribute('depth'))
    infodict['mag'] = float(event.getAttribute('magnitude'))
    timestr = event.getAttribute('event_timestamp')
    timestr = timestr[0:19]
    time = datetime.datetime(*strptime(timestr,"%Y-%m-%dT%H:%M:%S")[0:6])
    infodict['locstring'] = event.getAttribute('event_description')
    infodict['year'] = time.year
    infodict['month'] = time.month
    infodict['day'] = time.day
    infodict['hour'] = time.hour
    infodict['minute'] = time.minute
    infodict['second'] = time.second
    ctimestr = root.getAttribute('process_timestamp')
    ctimestr = ctimestr[0:19]
    ctime = datetime.datetime(*strptime(ctimestr,"%Y-%m-%dT%H:%M:%S")[0:6])
    infodict['created'] = ctime.strftime('%s')
    root.unlink()
    return infodict

def writeEvent(gridurl, stationurl, uncertaintyurl, eventdict, shakehome):
     datadir = os.path.join(shakehome,'data')
     inputdir = os.path.join(datadir,'input')
     if not os.path.isdir(inputdir):
         os.makedirs(inputdir)
    
     try:
         fh = urllib2.urlopen(gridurl)
         data = fh.read()
         fh.close()
         datafile = os.path.join(inputdir,'grid.xml')
         f = open(datafile,'wt')
         f.write(data)
         f.close()
     except:
         print 'No grid file found.'
         
     try:
         fh = urllib2.urlopen(stationurl)
         data = fh.read()
         fh.close()
         datafile = os.path.join(inputdir,'stationlist.xml')
         f = open(datafile,'wt')
         f.write(data)
         f.close()
     except:
         print 'No stationlist file found.'

     try:
         fh = urllib2.urlopen(uncertaintyurl)
         zf = zipfile.ZipFile(StringIO(fh.read()))
         data = zf.read('uncertainty.xml')
         fh.close()
         datafile = os.path.join(inputdir,'uncertainty.xml')
         f = open(datafile,'wt')
         f.write(data)
         f.close()
     except:
         print 'No uncertainty file found.'


def getShakeURLs(shakeurl):
    urlt = 'http://earthquake.usgs.gov/fdsnws/event/1/query?eventid=[EVENTID]&format=geojson'
    eventid = urlparse.urlparse(shakeurl).path.strip('/').split('/')[-1]
    url = urlt.replace('[EVENTID]',eventid)
    fh = urllib2.urlopen(url)
    data = fh.read()
    jdict = json.loads(data)
    fh.close()
    contentlist = jdict['properties']['products']['shakemap'][0]['contents'].keys()
    gridurl = None
    stationurl = None
    uncertaintyurlzp = None
    for content in contentlist:
        if content.endswith('grid.xml'):
            gridurl = jdict['properties']['products']['shakemap'][0]['contents'][content]['url']
        if content.find('stationlist.xml') > -1:
            stationurl = jdict['properties']['products']['shakemap'][0]['contents'][content]['url']
        if content.endswith('uncertainty.xml.zip'):
            uncertaintyurlzp = jdict['properties']['products']['shakemap'][0]['contents'][content]['url']
    return (gridurl,stationurl, uncertaintyurlzp)

def getGridfromArg(url):

    desc = ''' Get ShakeMap grid.xml, stationlist.xml, and uncertainty.xml
Example:
    %(prog)s http://earthquake.usgs.gov/earthquakes/shakemap/global/shake/capstone2014_nmsw_m7.7_se/
    from NEIC web site.
    '''
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('url', help='the URL of the desired ShakeMap.')
    args = parser.parse_args()
    if not args.url:
        print parser.print_help()
        sys.exit(1)
    shakeurl = args.url

    if shakeurl.find('#') == -1:
        endidx = len(shakeurl)
    else:
        endidx = shakeurl.find('#')
    shakeurl = shakeurl[0:endidx]
    if not shakeurl.endswith('/'):
        shakeurl += '/'

    shakehome = '/Users/sverros/Documents/Inputs'

    if shakeurl.find('_se') > -1:
        gridurl = urlparse.urljoin(shakeurl,'download/grid.xml')
        stationlisturl = urlparse.urljoin(shakeurl,'download/stationlist.xml')
        uncertaintyurlzp = urlparse.urljoin(shakeurl,'download/uncertainty.xml.zip')
    else:
        gridurl,stationlisturl, uncertaintyurlzp = getShakeURLs(shakeurl)

    eventdict = getEventInfo(gridurl)
    writeEvent(gridurl, stationlisturl, uncertaintyurlzp, eventdict, shakehome)
    print 'Files retrieved'

def main(shakeurl):
    #remove any stuff after a # sign in the url
    if shakeurl.find('#') == -1:
        endidx = len(shakeurl)
    else:
        endidx = shakeurl.find('#')
    shakeurl = shakeurl[0:endidx]
    if not shakeurl.endswith('/'):
        shakeurl += '/'

    shakehome = '/Users/sverros/Documents/Inputs/'

    if shakeurl.find('_se') > -1:
        gridurl = urlparse.urljoin(shakeurl,'download/grid.xml')
        stationlisturl = urlparse.urljoin(shakeurl,'download/stationlist.xml')
        uncertaintyurlzp = urlparse.urljoin(shakeurl,'download/uncertainty.xml.zip')
    else:
        gridurl,stationlisturl, uncertaintyurlzp = getShakeURLs(shakeurl)

    eventdict = getEventInfo(gridurl)
    
    writeEvent(gridurl, stationlisturl, uncertaintyurlzp, eventdict, shakehome)
    print 'Files retrieved'



if __name__ == '__main__':
    desc = ''' Get ShakeMap grid.xml, stationlist.xml, and uncertainty.xml

Example:
    %(prog)s http://earthquake.usgs.gov/earthquakes/shakemap/global/shake/capstone2014_nmsw_m7.7_se/ 

    from NEIC web site.
    '''
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('url', help='the URL of the desired ShakeMap.')
    args = parser.parse_args()
    if not args.url:
        print parser.print_help()
        sys.exit(1)
    main(args.url)
