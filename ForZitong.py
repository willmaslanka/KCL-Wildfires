from satpy import Scene
import numpy as np
from pyresample import create_area_def
from satpy.composites import GenericCompositor
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#%% UDF
def geolocating_VIIRS_data(filenames, bbox, bands=['M05','M04','M03']):
    """
    Master Function used for reprojecting and plotting the multi-band data 
    from VIIRS. Reprojects the RGB data to 1000m in Web Mercator Projection.

    Parameters
    ----------
    filenames : List of strings
        List of filenames for the reflectance and georeference files.
        Should be in the following order:
            Moderate-resolution Reflectance filename
            Moderate-resolution Geolocation filename
            Image-resolution Reflectance filename
            Image-resolution Geolocation filename
    bbox : List of strings
        List of strings that creates the bounding box.
        Should be +/- 2 degrees from the central plume identifier.
        Will look something like this
            ['lon_w', 'lat_s', 'lon_e', 'lat_n']
        Where lon_w is the western-most longitude
        and lat_s is the southern-most latitude
        and lon_e is the eastern-most longitude
        and lat_n is the norhtern-most latitude.
    bands: List of strings
        List of strings of the three bands used to make the colour composite.
        If undefined, set to RGB from VIIRS (M05, M04, M03)
        Can be either M or I bands.
        Should be in the order [Red-Band Green-Band Blue-Band]
        
        List of bands available: 

    Returns
    -------
    RGB : Array of projected RGB data from VIIRS. Geolocated to the projection
         stated in crs_custom.
    crs_custom : Custom produced CRS projection.
    """
    
    scn = Scene(reader='viirs_l1b', filenames=filenames)
    scn.load([bands[0], bands[1], bands[2], 'm_lon','m_lat'])
    
    dynamic_scene = resample_scene(bbox, scn, 1000) # 1000m
    
    composite_data = build_composite(dynamic_scene, 
                                     bands[0], bands[1], bands[2])
    
    RGB, crs_custom = get_composite_data(composite_data, scale='Percentage')
    
    return RGB, crs_custom

def area_def(bbox, resolution=1000):
    """
    Creates a new area definition variable, for a defined region (bounding box)
    and given resolution. Used to resample the VIIRS RGB data to a regular 
    grid, using a Web Mercator projection.

    Parameters
    ----------
    bbox : List of strings
        List of strings that creates the bounding box.
        Should be +/- 2 degrees from the central plume identifier.
        Will look something like this
            ['lon_w', 'lat_s', 'lon_e', 'lat_n']
        Where lon_w is the western-most longitude
        and lat_s is the southern-most latitude
        and lon_e is the eastern-most longitude
        and lat_n is the norhtern-most latitude.
    resolution : float
        Resolution wanted (in metres).
        If undefined, set to 1000m.

    Returns
    -------
    my_dynamic_area : Custom Area Definition Variable
    """
    
    my_dynamic_area = create_area_def('my_area', {'proj': 'webmerc', 
                                                  'lon_0': float(bbox[0]), 
                                                  'lat_0': float(bbox[1]), 
                                                  'lon_1': float(bbox[2]), 
                                                  'lat_1': float(bbox[3])},
                                                  resolution=resolution)
    return my_dynamic_area

def resample_scene(bbox, satpy_scene, resolution=1000):
    """
    Sub-function used to resample the scene.
        Creates the area definition, then resamples the data.

    Parameters
    ----------
    bbox : List of strings
        List of strings that creates the bounding box.
        Should be +/- 2 degrees from the central plume identifier.
        Will look something like this
            ['lon_w', 'lat_s', 'lon_e', 'lat_n']
        Where lon_w is the western-most longitude
        and lat_s is the southern-most latitude
        and lon_e is the eastern-most longitude
        and lat_n is the norhtern-most latitude.
    satpy_scene : Sat-Py Scene
        Sat-Py Scene variable, with the colour composite data from VIIRS
    resolution : float
        Resolution wanted (in metres).
        If undefined, set to 1000m.

    Returns
    -------
    dynamic_scn : Sat-Py Scene
        Resampled SatPy Scene, projected and resampled to the Web Meracator 
        projection, to the given resolution (default is 1000m). 

    """
    my_dynamic_area = area_def(bbox, resolution)
    dynamic_scn = satpy_scene.resample(my_dynamic_area)
    return dynamic_scn
	
def build_composite(dynamic_scene, Band1, Band2, Band3):
    """
    SatPy's Composite Builder Function. Builds the 3-band composite from the
    user defined bands.

    Parameters
    ----------
    satpy_scene : Sat-Py Scene
        Sat-Py Scene variable, with the colour composite data from VIIRS
    Band1 : String
        Band wanted in the "Red" position. Should be automated.
    Band2 : String
        Band wanted in the "Green" position. Should be automated.
    Band3 : String
        Band wanted in the "Blue" position. Should be automated.

    Returns
    -------
    composite_data : SatPy Composite Variable, containing the RGB-like swath 

    """
    compositor = GenericCompositor("overview")
    composite_data = compositor([dynamic_scene[Band1],
                                 dynamic_scene[Band2],
                                 dynamic_scene[Band3]])
    return composite_data

def get_composite_data(composite_data, scale='MinMax'):
    """
    Creates the SatPy RGB-like composite, from the SatPy scene variable. Scales
    it to a user defined scale, defaulting to a min-max scale.

    Parameters
    ----------
    composite_data : SatPy Composite Variable, containing the RGB-like swath 
    data
        DESCRIPTION.
    scale : String, optional
        User Defined Function for scaling the RGB data. 
        Choices are "MinMax" for scaling between the minimum and maximum values
        and "Percent"; which scales between the (pcnt)% and (100-pcnt)% values.
        The default is 'MinMax'.

    Returns
    -------
    image : TYPE
        DESCRIPTION.
    crs : TYPE
        DESCRIPTION.

    """
    image = np.asarray(composite_data).transpose(1,2,0)
    image = np.nan_to_num(image)
    if scale.lower() == 'percentage':
        image = np.interp(image, (np.percentile(image, 2),
                                  np.percentile(image,98)), (0, 1))
    elif scale.lower() == 'minmax':
        image = np.interp(image, (np.percentile(image, 0),
                                  np.percentile(image,100)), (0, 1))
    crs = composite_data.attrs["area"].to_cartopy_crs()
    
    return image, crs

def scaleMinMax(x):
    """
    Produces a Scaling Index between the smallest and largest values


    Parameters
    ----------
    x : Array
        Array to be scaled.

    Returns
    -------
    scle : Array
        Scaled Array.

    """
    scle = (x - np.nanmin(x))/(np.nanmax(x) - np.nanmin(x))
    return scle

def scalePercent(x, pcnt):
    """
    Produces a Scaling Index between the pcnt% and (100-pcnt)% values

    Parameters
    ----------
    x : Array to be scaled
        DESCRIPTION.
    pcnt : Integar
        How much to trim from the array.

    Returns
    -------
    scle : Array
        Scaled .

    """
    scle = (x - np.nanpercentile(x, pcnt)) / (np.nanpercentile(x, 100 - pcnt) - np.nanpercentile(x, pcnt))
    
    scle[scle >= 1] = 1
    scle[scle <= 0] = 0
    
    return scle	

#%%
data_dir = '/home/will/data/ROI_Siberia2020/'


mref_file = f'{data_dir}/CalRef/2020/158/VNP02MOD.A2020158.0112.001.2020158064930.nc'
mgeo_file = f'{data_dir}/GeoLoc/2020/158/VNP03MOD.A2020158.0112.001.2020158063732.nc'
iref_file = f'{data_dir}/CalRef/2020/158/VNP02IMG.A2020158.0112.001.2020158064930.nc'
igeo_file = f'{data_dir}/GeoLoc/2020/158/VNP03IMG.A2020158.0112.001.2020158063732.nc'

filenames=[mref_file, mgeo_file, iref_file, igeo_file]
bbox = ['125.69', '61.55', '129.69', '65.55'] #[lon_w, lat_s, lon_e, lat_n]

RGB, crs_custom = geolocating_VIIRS_data(filenames, bbox)

#%%
fig = plt.figure(figsize=(8,6))
ax = plt.subplot(1,1,1, projection = ccrs.Mercator())
ax.set_extent([float(bbox[0]), float(bbox[2]), float(bbox[1]), float(bbox[3])],
              crs=ccrs.PlateCarree())
ax.imshow(RGB, transform=crs_custom, extent=crs_custom.bounds, origin="upper")

# Plotting the Plume Outline
# p_outlon = plume outline longitudes
# p_outlat = plume outline latitudes

plt.plot(p_outlon, p_outlat, color='cyan', 
         zorder=2, transform=ccrs.PlateCarree())
