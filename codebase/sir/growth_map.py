import sys
sys.path.append('../../')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from codebase.sir.fetch import get_all_state_growth
import matplotlib as mpl
import matplotlib.cm as cm
from metpy.plots import USCOUNTIES


state_data,dates=get_all_state_growth(0.2,-10)

plot_name="./map_%s_%s.png"%(dates[0],dates[-1])

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.LambertConformal())

ax.set_extent([-125, -66.5, 20, 50], ccrs.Geodetic())

shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',
                                     category='cultural', name=shapename)

ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)

ax.set_title('Mean Doubling Time (%s to %s)'%(dates[0],dates[-1]))

colors=state_data.values()
max_val=10.0
min_val=0.0
norm=mpl.colors.Normalize(vmin=min_val,vmax=max_val)
cmap=cm.jet_r
#cmap=cm.terrain
m=cm.ScalarMappable(norm=norm,cmap=cmap)


#for state in shpreader.Reader(states_shp).geometries():
for astate in shpreader.Reader(states_shp).records():

    ### You want to replace the following code with code that sets the
    ### facecolor as a gradient based on the population density above
    #facecolor = [0.9375, 0.9375, 0.859375]
    name=astate.attributes['name']

    # `astate.geometry` is the polygon to plot
    ax.add_geometries([astate.geometry], 
                      ccrs.PlateCarree(),
                      facecolor=m.to_rgba(state_data[name]),
                      edgecolor='black')#m.to_rgba(state_data[name]))

#ax.add_feature(USCOUNTIES.with_scale("20m"))

cax=fig.add_axes([0.27,0.1,0.5,0.05])
cb=mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
cb.set_label("Days")

print("Saved %s"%plot_name)
plt.savefig(plot_name)
