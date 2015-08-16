import numpy as np
from mpl_toolkits.basemap import Basemap

def pcolormesh(lon, lat, data, boundinglat=-55, lw=0.25,
                              ax=None, labels=[1,1,0,1], **kwargs):
    """Plot something in the southern ocean."""
    
    m = Basemap(projection='spstere',boundinglat=boundinglat,lon_0=180, ax=ax)
    #m.drawcoastlines()
    m.fillcontinents(color='0.5',lake_color='0.55')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.), linewidth=lw)
    m.drawmeridians(np.arange(-180.,181.,30.), labels=labels, linewidth=lw)
    m.drawmapboundary(fill_color='w')
    #x, y = m(nc.variables['TLON'][:], nc.variables['TLAT'][:])
    return m.pcolormesh(lon, lat, data, latlon=True, **kwargs), m

def contourf(lon, lat, data, clevs, boundinglat=-55, lw=0.25, 
                            ax=None, labels=[1,1,0,1], **kwargs):
    """Plot something in the southern ocean."""
    
    m = Basemap(projection='spstere',boundinglat=boundinglat,lon_0=180, ax=ax)
    #m.drawcoastlines()
    m.fillcontinents(color='0.5',lake_color='0.55')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.), linewidth=lw)
    m.drawmeridians(np.arange(-180.,181.,30.), labels=labels, linewidth=lw)
    m.drawmapboundary(fill_color='w')
    #x, y = m(nc.variables['TLON'][:], nc.variables['TLAT'][:])
    return m.contourf(lon, lat, data, clevs, latlon=True, **kwargs), m
