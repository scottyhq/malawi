#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Web Map Tile Service time dimension demonstration
-------------------------------------------------

This example further demonstrates WMTS support within cartopy. Optional
keyword arguments can be supplied to the OGC WMTS 'gettile' method. This
allows for the specification of the 'time' dimension for a WMTS layer
which supports it.

There are 10000+ WMS services out there. Here are some compiled lists:
    http://www.skylab-mobilesystems.com/en/wms_serverlist.html 
    http://directory.spatineo.com
    http://directory.spatineo.com/service/42691/
    
# See using open street map:
#http://scitools.org.uk/cartopy/docs/v0.15/examples/tube_stations.html

# Planet Labs
    https://www.planet.com/docs/reference/tile-services/

# Google tiles:
https://ocefpaf.github.io/python4oceanographers/blog/2015/06/22/osm/

The example shows satellite imagery retrieved from NASA's Global Imagery
Browse Services for 5th Feb 2016. A true color MODIS image is shown on
the left, with the MODIS false color 'snow RGB' shown on the right.

"""
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
#import matplotlib.ticker as mticker
from  matplotlib.dates import YearLocator, MonthLocator, DateFormatter

import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from owslib.wmts import WebMapTileService

import geopandas as gpd
import pandas as pd
from pandas.plotting import table
import numpy as np
import shapely.wkt
import urllib
import json
import os

#plt.style.use('seaborn-white')
plt.rcParams['font.size'] = 14

def add_wmts_gibs_basemap(ax, date='2016-02-05'):
    """http://gibs.earthdata.nasa.gov/"""
    URL = 'http://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi'
    wmts = WebMapTileService(URL)

    # Layers for MODIS true color and snow RGB
    # NOTE: what other tiles available?: TONS!
    #https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-ReferenceLayers9Layers
    #layer = 'MODIS_Terra_SurfaceReflectance_Bands143'
    #layer = 'MODIS_Terra_CorrectedReflectance_Bands367'
    #layer = 'ASTER_GDEM_Greyscale_Shaded_Relief' #better zoomed in
    layer = 'SRTM_Color_Index' 
    #layer = 'BlueMarble_ShadedRelief' #static
    #layer = 'BlueMarble_NextGeneration'
    #layer = 'BlueMarble_ShadedRelief_Bathymetry'
    #layer = 'Reference_Labels'
    #layer = 'Reference_Features'

    ax.add_wmts(wmts, layer, wmts_kwargs={'time': date}) # alpha=0.5
    #NOTE: can access attributes:
    #wmts[layer].title
    return wmts


def add_xyz_tile(ax, url, zoom=6):
    """ Grab a map tile from the web """
    from cartopy.io.img_tiles import GoogleTiles
    # Not sure about how projections are handled...
    #url = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg'
    #url = 'http://tile.stamen.com/watercolor/{z}/{x}/{y}.png'
    #url = 'https://s3.amazonaws.com/elevation-tiles-prod/normal/{z}/{x}/{y}.png'
    tiler = GoogleTiles(url=url)
    ax.add_image(tiler, zoom)
    


def query_asf(snwe, sat='1A'):
    '''
    takes list of [south, north, west, east]
    '''
    miny, maxy, minx, maxx = snwe
    roi = shapely.geometry.box(minx, miny, maxx, maxy)
    #test.__geo_interface__ #geojson?
    polygonWKT = roi.to_wkt()

    polygon = urllib.parse.quote_plus(polygonWKT) #for URL
    base = 'https://api.daac.asf.alaska.edu/services/search/param?'
    #Note: 'intersects' gets entire global track! ALSO, better to use 'requests' librar?
    poly = 'intersectsWith={}'.format(polygon)
    plat = 'platform=Sentinel-{}'.format(sat) #1B
    proc = 'processingLevel=SLC'
    beam = 'beamMode=IW'
    #relativeOrbit=$ORBIT
    out =  'output=json > query{}.json'.format(sat)
    #out =  'output=kml > query.kml'
    querystr = '\\&'.join([base, poly, plat, beam, proc, out]) #escape backslash
    
    #poly= polygonWKT
    #url = '\&'.join([base, poly, plat, beam, proc, out]) #
    os.system('curl ' + querystr)
    print(querystr)
    
    

# query ASF catalog and get json back
def load_asf_json(jsonfile):
    ''' Convert JSON metadata from asf query to dataframe '''
    with open(jsonfile) as f:
        meta = json.load(f)[0] #list of scene dictionaries
    
    df = pd.DataFrame(meta)
    polygons = df.stringFootprint.apply(shapely.wkt.loads)
    gf = gpd.GeoDataFrame(df, 
                          crs={'init': 'epsg:4326'},
                          geometry=polygons)
    #gf['date'] = pd.to_datetime(df.sceneDate, format='%Y-%m-%d %H:%M:%S')
    #gf.to_file(outfile) #saves as shapefile
    return gf 



def main():
    # Plot setup
    plot_CRS = ccrs.Mercator()
    geodetic_CRS = ccrs.Geodetic()
    x0, y0 = plot_CRS.transform_point(30, -16, geodetic_CRS)
    x1, y1 = plot_CRS.transform_point(40, -6, geodetic_CRS)
    
    fig,ax = plt.subplots(figsize=(8,8), dpi=100, 
                          subplot_kw=dict(projection=plot_CRS))
    
    ax.set_xlim((x0, x1))
    ax.set_ylim((y0, y1))
    
    #wmts = add_wmts_gibs_basemap(ax)
    
    #url = 'http://tile.stamen.com/watercolor/{z}/{x}/{y}.png'
    #url = 'http://tile.stamen.com/terrain/{z}/{x}/{y}.png'
    #url = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg'
    #url = 'https://s3.amazonaws.com/elevation-tiles-prod/normal/{z}/{x}/{y}.png'
    #add_xyz_tile(ax, url, zoom=7) #higher number higher resultion #can't pass alpha argument...
    
    #ax.stock_img()
    
    # New 0.15 background method...
    os.environ["CARTOPY_USER_BACKGROUNDS"] = '/Users/scott/Data/cartopy_data'
    #ax.background_img("GrayEarth")
    ax.background_img("ne_shaded", "low")
    
    #ax.coastlines(resolution='50m', color='black')
    # http://scitools.org.uk/cartopy/docs/v0.15/matplotlib/feature_interface.html
    #https://gis.stackexchange.com/questions/88209/python-mapping-in-matplotlib-cartopy-color-one-country
    #borders = cfeature.BORDERS #low res
    borders = cfeature.NaturalEarthFeature(scale='10m',
                                 category='cultural',
                                 name='admin_0_boundary_lines_land')
    ax.add_feature(borders, facecolor='none', edgecolor='k')
    
    
    # Label Basemap
    '''
    txt = plt.text(-69, -20.5, 'MODIS' , fontsize=18,
                  color='wheat', transform=geodetic_CRS)
    txt.set_path_effects([PathEffects.withStroke(linewidth=5,
                                                 foreground='black')])
    '''
    
    # Add Rungwe
    ax.plot(-33.668, -9.135, 'k^', transform=geodetic_CRS)
    
    # Add deformation footprint
    # Add ROI gjson
    # Go directly from SNWE box: -11.75 -8.75 32.5 35.5
    snwe = [-11.75, -8.75, 32.5, 35.5]
    #query_asf(snwe, '1A') #downloads query.json
    #query_asf(snwe, '1B') #downloads query.json
    
    # Some basic coverage analysis
    gf = load_asf_json('query.json')
    orbits = gf.relativeOrbit.unique()
    #gf.to_file('query.geojson', driver='GeoJSON') #automatically rendered on github!
    #gf.iloc[0] #everything available
    #gf.plot()
    gf.groupby('relativeOrbit').fileName.count()
    gf.groupby('relativeOrbit').sceneDate.describe()
    gf.groupby('relativeOrbit').sceneDate.min()
    gf.groupby('relativeOrbit').sceneDate.max()
    
    # create a summary dataframe, insert as figure metadata
    #dfS = pd.DataFrame(index=gf.relativeOrbit.unique(), columns=['Count','Start','Stop'])
    #dfS['Count'] = gf.groupby('relativeOrbit').fileName.count()
    #dfS['Start'] = gf.groupby('relativeOrbit').sceneDate.min()
    #dfS['Stop'] = gf.groupby('relativeOrbit').sceneDate.max()
    # Above has weird dtypes...    
    # Add to plot! as a custom legend - see below for separate figure
    #table(ax, dfS, loc='upper right', zorder=10)
    
    # Add all scene footprints to the plot!
    '''
    ax.add_geometries(gf.geometry.values, ccrs.PlateCarree(), 
                  facecolor='none', 
                  edgecolor='k',
                  linewidth=1,
                  linestyle='-')
    '''
    
    # To unclutter, show cascaded union for each track in different colors
    #colors = ['c','m','b','y']
    colors = plt.cm.jet(np.linspace(0,1,orbits.size))
    #colors = plt.get_cmap('jet', orbits.size) #not iterable
    for orbit,color in zip(orbits, colors):
        df = gf.query('relativeOrbit == @orbit')
        poly = df.geometry.cascaded_union
        
        if df.flightDirection.iloc[0] == 'ASCENDING':
            linestyle = '--'
            #xpos, ypos = poly.bounds[0], poly.bounds[3] #upper left
            xpos,ypos = poly.centroid.x, poly.bounds[3] 
        else:
            linestyle = '-'
            #xpos, ypos = poly.bounds[2], poly.bounds[1] #lower right
            xpos,ypos = poly.centroid.x, poly.bounds[1] 
        
        
        ax.add_geometries([poly], 
                          ccrs.PlateCarree(), 
                          facecolor='none', 
                          edgecolor=color,
                          linewidth=2,
                          linestyle=linestyle)
        ax.text(xpos, ypos, orbit, color=color, fontsize=16, fontweight='bold', transform=geodetic_CRS)
    
    # Add some text labels
    #fs = 16
    #ax.text(-69.6, -20.6, 'A149', color='b', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    #ax.text(-67.0, -20.25, 'A76', color='y', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    #ax.text(-67.8, -24.9, 'D83', color='m', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    #ax.text(-69.6, -24.6, 'D156', color='c', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    
    
    # Warning, mixing PlateCarree and Mercator here... might not work
    # for large regions
    gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='-')
    gl.xlabels_top = False
    gl.ylabels_left = False
    #gl.xlines = False
    #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    #plt.show()
    plt.title('Sentinel-1 Coverage SEGMENT Project')
    plt.savefig('map_coverage.pdf', bbox_inches='tight')
    
    # Second figure for timeline
    '''
    df = pd.DataFrame(gf.relativeOrbit.astype('int'))
    #df['platform'] = gf.platform
    df['sceneDate'] = pd.to_datetime(gf.sceneDate)
    df['dateStr'] = df.sceneDate.apply(lambda x: x.strftime('%Y-%m-%d'))
    df['code'] = df.relativeOrbit.astype('category').cat.codes
    
    dfS = pd.DataFrame(index=df.relativeOrbit.unique())
    dfS['Count'] = df.groupby('relativeOrbit').dateStr.count()
    dfS['Start'] = df.groupby('relativeOrbit').dateStr.min()
    dfS['Stop'] =  df.groupby('relativeOrbit').dateStr.max()
    dfS.sort_index(inplace=True, ascending=False)
    dfS.index.name = 'Orbit'
    plot_timeline_table(df, dfS)
    '''
    
    #plot_timeline_table_new()


def plot_timeline(df):
    fig,ax = plt.subplots(figsize=(8,4))
    plt.scatter(df.sceneDate.values, df.code.values, c=df.code.values, cmap='viridis', s=40)
    plt.yticks(df.code.unique(), df.relativeOrbit.unique())
    plt.axvline('2016-04-22', color='gray', linestyle='dashed', label='Sentinel-1B launch')
    ax.xaxis.set_minor_locator(MonthLocator())
    ax.xaxis.set_major_locator(YearLocator())
    plt.legend()
    plt.ylabel('Orbit Number')
    fig.autofmt_xdate()
    plt.title('Sentinel-1 Coverage at Uturuncu')
    plt.savefig('timeline.pdf', bbox_inches='tight')


    
def plot_timeline_table(df, dfS):
    plt.rcParams['font.size'] = 14
    fig,ax = plt.subplots(figsize=(11,8.5))
    plt.scatter(df.sceneDate.values, df.code.values, c=df.code.values, cmap='viridis', s=40)
    plt.yticks(df.code.unique(), df.relativeOrbit.unique())
    plt.axvline('2016-04-22', color='gray', linestyle='dashed', label='Sentinel-1B launch')
    
    # Add to plot! as a custom legend
    table(ax, dfS, loc='top', zorder=10, fontsize=12,
          #colWidths=[0.2, 0.2, 0.2, 0.2],
          bbox=[0.1, 0.7, 0.8, 0.275] )#[left, bottom, width, height])
    
    ax.xaxis.set_minor_locator(MonthLocator())
    ax.xaxis.set_major_locator(YearLocator())
    plt.legend(loc='lower right')
    plt.ylim(-1,6)
    plt.ylabel('Orbit Number')
    fig.autofmt_xdate()
    plt.title('Sentinel-1 Coverage SEGMENT Project')
    plt.savefig('timeline_with_table.pdf', bbox_inches='tight')
    
    
def print_acquisitions(df, orbit=156):
    ''' Write acquistions to csv file '''
    dates = pd.to_datetime(df.query('relativeOrbit == @orbit').dateStr)
    tmp = pd.DataFrame(dates).drop_duplicates()
    tmp['dt'] = tmp.diff()
    filename = 'dates_{}.csv'.format(orbit)
    tmp.to_csv(filename, sep='\t', index=False)
    print('wrote {} acquisition dates to {}'.format(len(tmp), filename))


def get_inventory():
    gf = load_asf_json('query1A.json')
    dfA = pd.DataFrame(gf.relativeOrbit.astype('int'))
    dfA['platform'] = gf.platform
    dfA['sceneDate'] = pd.to_datetime(gf.sceneDate)
    dfA['dateStr'] = dfA.sceneDate.apply(lambda x: x.strftime('%Y-%m-%d'))
    dfA['code'] = dfA.relativeOrbit.astype('category').cat.codes

    gf = load_asf_json('query1B.json')
    dfB = pd.DataFrame(gf.relativeOrbit.astype('int'))
    dfB['platform'] = gf.platform
    dfB['sceneDate'] = pd.to_datetime(gf.sceneDate)
    dfB['dateStr'] = dfB.sceneDate.apply(lambda x: x.strftime('%Y-%m-%d'))
    dfB['code'] = dfB.relativeOrbit.astype('category').cat.codes

    return dfA, dfB

    
def plot_timeline_table_new():
    '''
    Improvements: use different symbols for 1A and 1B
    '''
    plt.rcParams['font.size'] = 16
    
    dfA,dfB = get_inventory()
    df = pd.concat([dfA,dfB])
    df.reset_index(inplace=True)
    dfS = pd.DataFrame(index=df.relativeOrbit.unique())
    dfS['Start'] = df.groupby('relativeOrbit').dateStr.min()
    dfS['Stop'] =  df.groupby('relativeOrbit').dateStr.max()
    dfS['#Dates'] = df.groupby('relativeOrbit').dateStr.nunique()
    #dfS['#Frames'] = df.groupby('relativeOrbit').dateStr.count()
    dfS.sort_index(inplace=True, ascending=False)
    dfS.index.name = 'Orbit'
    
    # Same colors as map
    orbits = df.relativeOrbit.unique()
    colors = plt.cm.jet(np.linspace(0,1, orbits.size))
    
    fig,ax = plt.subplots(figsize=(17,8))
    plt.scatter(dfA.sceneDate.values, dfA.code.values, c=dfA.code.values, cmap='jet', s=60, label='S1A')
    plt.scatter(dfB.sceneDate.values, dfB.code.values, c=dfB.code.values, cmap='jet', s=60, marker='d',label='S1B')
    plt.yticks(df.code.unique(), df.relativeOrbit.unique())
    plt.axvline('2016-04-22', color='gray', linestyle='dashed', label='Sentinel-1B launch')
    
    # Add to plot! as a custom legend
    table(ax, dfS, loc='top', zorder=10, fontsize=12,
          cellLoc = 'center', rowLoc = 'center',
          bbox=[0.1, 0.7, 0.6, 0.3] )#[left, bottom, width, height])
    
    ax.xaxis.set_minor_locator(MonthLocator())
    ax.xaxis.set_major_locator(YearLocator())
    plt.legend(loc='upper right')
    plt.ylim(-1,orbits.size+3)
    plt.ylabel('Orbit Number')
    fig.autofmt_xdate()
    plt.title('Sentinel-1 Coverage SEGMENT Project')
    plt.savefig('timeline_with_table.pdf', bbox_inches='tight')

if __name__ == '__main__':
    main()