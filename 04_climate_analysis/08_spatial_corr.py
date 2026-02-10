import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import cftime
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Parámetros a modificar
nc_file = 'PLACEHOLDER/path/to/rocio.nc'
chronology_file_path = 'PLACEHOLDER/path/to/chronology.txt'
start_year = 1951
end_year = 2022
accumulated_days = 324
end_month = 7  # Mes de finalización (Julio)
end_day = 1    # Día de finalización

# Cargar el archivo netCDF
dataset = nc.Dataset(nc_file)

# Extraer variables
precipitation = dataset.variables['precipitation'][:, 0, :, :]
lon = dataset.variables['lon'][0, :, :]
lat = dataset.variables['lat'][0, :, :]
time = dataset.variables['time'][:]

# Convertir tiempo a fechas
dates = nc.num2date(time, units=dataset.variables['time'].units, only_use_cftime_datetimes=False)

# Calcular la precipitación acumulada para cada año
years = np.arange(start_year, end_year + 1)
accumulated_precipitation = np.zeros((len(years), precipitation.shape[1], precipitation.shape[2]))

for i, year in enumerate(years):
    start_date = cftime.DatetimeGregorian(year, end_month, end_day) - np.timedelta64(accumulated_days, 'D')
    end_date = cftime.DatetimeGregorian(year, end_month, end_day)
    mask = (dates >= start_date) & (dates < end_date)
    accumulated_precipitation[i, :, :] = np.sum(precipitation[mask, :, :], axis=0)

# Cargar archivo de cronología
chronology_data = np.genfromtxt(chronology_file_path, skip_header=1, missing_values='NA', filling_values=np.nan)

# Filtrar los datos de cronología para los años 1951-2022 y usar la columna 'res'
chronology_years = chronology_data[:, 0]
chronology_residuals = chronology_data[:, 2]
mask = (chronology_years >= start_year) & (chronology_years <= end_year)
filtered_chronology_residuals = chronology_residuals[mask]

# Calcular correlaciones espaciales y R
correlations = np.zeros((precipitation.shape[1], precipitation.shape[2]))
p_values = np.zeros((precipitation.shape[1], precipitation.shape[2]))

for i in range(precipitation.shape[1]):
    for j in range(precipitation.shape[2]):
        correlations[i, j], p_values[i, j] = pearsonr(filtered_chronology_residuals, accumulated_precipitation[:, i, j])

# Filtrar R significativos a nivel de significancia 0.05 (95%)
significant_r = np.where(p_values < 0.05, correlations, np.nan)

# Verificar si la latitud está en orden descendente y corregir si es necesario
if lat[0, 0] > lat[-1, 0]:  
    lat = lat[::-1, :]
    significant_r = significant_r[::-1, :]

# Verificar si la longitud está en orden ascendente y corregir si es necesario
if lon[0, 0] > lon[0, -1]:  
    lon = lon[:, ::-1]
    significant_r = significant_r[:, ::-1]

# Crear la nueva paleta de colores que va de #F1BC9C a #C40A0C pasando por tonos suaves de rojo
colors = ['#F1BC9C', '#E18466', '#D24734', '#C40A0C']
cmap_name = 'lightred_to_red'
lightred_to_red_cmap = LinearSegmentedColormap.from_list(cmap_name, colors)

# Plotear el mapa con la nueva paleta de colores y agregar el contorno de los países
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-', alpha=0.5)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linestyle='-', alpha=0.5)
contourf_plot = ax.contourf(lon, lat, significant_r, cmap=lightred_to_red_cmap)
cbar = plt.colorbar(contourf_plot, ax=ax)
cbar.set_label('R')
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.title('Correlación espacial con precipitación acumulada (R significativa a nivel 0.05)')

plt.show()