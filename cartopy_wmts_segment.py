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
    layer = 'ASTER_GDEM_Greyscale_Shaded_Relief' #better zoomed in
    #layer = 'SRTM_Color_Index'

    ax.add_wmts(wmts, layer, wmts_kwargs={'time': date})
    #NOTE: can access attributes:
    #wmts[layer].title
    return wmts


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
    
    wmts = add_wmts_gibs_basemap(ax)
    
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
    #minx, miny, maxx, maxy
    roi = shapely.geometry.box(32.5, -11.75, 35.5, -8.75)
    #test.__geo_interface__ #geojson?
    polygonWKT = roi.to_wkt()
    '''
    df = gpd.read_file('/Users/scott/Dropbox/planet/map.geojson')
    ax.add_geometries(df.geometry.values, ccrs.PlateCarree(), 
                      facecolor='none', 
                      edgecolor='k',
                      linewidth=3,
                      linestyle='--')
    polygonWKT = df.geometry[0].to_wkt()
    '''
    polygon = urllib.parse.quote_plus(polygonWKT) #for URL
    base = 'https://api.daac.asf.alaska.edu/services/search/param?'
    #Note: 'intersects' gets entire global track! ALSO, better to use 'requests' librar?
    poly = 'intersectsWith={}'.format(polygon)
    plat = 'platform=Sentinel-{}'.format('1A') #1B
    proc = 'processingLevel=SLC'
    beam = 'beamMode=IW'
    #relativeOrbit=$ORBIT
    out =  'output=json > query.json'
    querystr = '\\&'.join([base, poly, plat, beam, proc, out])
    print(querystr)
    #os.system('curl ' + querystr) #does work from spyder
    
    # Some basic coverage analysis
    gf = load_asf_json('query.json')
    gf.relativeOrbit.unique()
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
    colors = ['c','m','b','y']
    for orbit,color in zip(gf.relativeOrbit.unique(),colors):
        df = gf.query('relativeOrbit == @orbit')
        ax.add_geometries([df.geometry.cascaded_union], 
                          ccrs.PlateCarree(), 
                          facecolor='none', 
                          edgecolor=color,
                          linewidth=2,
                          linestyle='-')
    
    # Add some text labels
    fs = 16
    ax.text(-69.6, -20.6, 'A149', color='b', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    ax.text(-67.0, -20.25, 'A76', color='y', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    ax.text(-67.8, -24.9, 'D83', color='m', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    ax.text(-69.6, -24.6, 'D156', color='c', fontsize=fs, fontweight='bold', transform=geodetic_CRS)
    
    
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
    plt.title('Sentinel-1 Coverage at Uturuncu')
    plt.savefig('map_coverage.pdf', bbox_inches='tight')
    
    # Second figure for timeline
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
    plt.title('Sentinel-1A Coverage at Uturuncu')
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
    gf = load_asf_json('/Users/scott/Dropbox/planet/query1A.json')
    dfA = pd.DataFrame(gf.relativeOrbit.astype('int'))
    dfA['platform'] = gf.platform
    dfA['sceneDate'] = pd.to_datetime(gf.sceneDate)
    dfA['dateStr'] = dfA.sceneDate.apply(lambda x: x.strftime('%Y-%m-%d'))
    dfA['code'] = dfA.relativeOrbit.astype('category').cat.codes

    gf = load_asf_json('/Users/scott/Dropbox/planet/query1B.json')
    dfB = pd.DataFrame(gf.relativeOrbit.astype('int'))
    dfB['platform'] = gf.platform
    dfB['sceneDate'] = pd.to_datetime(gf.sceneDate)
    dfB['dateStr'] = dfB.sceneDate.apply(lambda x: x.strftime('%Y-%m-%d'))
    dfB['code'] = dfB.relativeOrbit.astype('category').cat.codes
    
    df = pd.concat([dfA,dfB])
    df.reset_index(inplace=True)

    return df

    
def plot_timeline_table_new():
    '''
    Improvements: use different symbols for 1A and 1B
    '''
    plt.rcParams['font.size'] = 16
    
    df = get_inventory()
    
    dfS = pd.DataFrame(index=df.relativeOrbit.unique())
    dfS['Start'] = df.groupby('relativeOrbit').dateStr.min()
    dfS['Stop'] =  df.groupby('relativeOrbit').dateStr.max()
    dfS['#Dates'] = df.groupby('relativeOrbit').dateStr.nunique()
    #dfS['#Frames'] = df.groupby('relativeOrbit').dateStr.count()
    dfS.sort_index(inplace=True, ascending=False)
    dfS.index.name = 'Orbit'
    
    
    fig,ax = plt.subplots(figsize=(17,8))
    plt.scatter(dfA.sceneDate.values, dfA.code.values, c=dfA.code.values, cmap='viridis', s=60, label='S1A')
    plt.scatter(dfB.sceneDate.values, dfB.code.values, c=dfB.code.values, cmap='viridis', s=60, marker='d',label='S1B')
    plt.yticks(df.code.unique(), df.relativeOrbit.unique())
    plt.axvline('2016-04-22', color='gray', linestyle='dashed', label='Sentinel-1B launch')
    
    # Add to plot! as a custom legend
    table(ax, dfS, loc='top', zorder=10, fontsize=12,
          #colWidths=[0.2, 0.2, 0.2, 0.1],
          bbox=[0.1, 0.7, 0.8, 0.275] )#[left, bottom, width, height])
    
    ax.xaxis.set_minor_locator(MonthLocator())
    ax.xaxis.set_major_locator(YearLocator())
    plt.legend(loc='lower right')
    plt.ylim(-2,6)
    plt.ylabel('Orbit Number')
    fig.autofmt_xdate()
    plt.title('Sentinel-1 Coverage at Uturuncu')
    plt.savefig('timeline_with_table.pdf', bbox_inches='tight')

if __name__ == '__main__':
    main()