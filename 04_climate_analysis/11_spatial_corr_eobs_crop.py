import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import cftime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

# Parámetros a modificar
nc_file = 'PLACEHOLDER/path/to/eobs01.nc'
chronology_file_path = 'PLACEHOLDER/path/to/chronology.txt'
start_year = 1950
end_year = 2022
accumulated_days = 324
end_month = 7  # Mes de finalización (Julio)
end_day = 1    # Día de finalización

# Nombre para el archivo PDF de salida
output_pdf_filename = 'mapa_correlacion_precipitacion_recortado2.pdf' # Nombre descriptivo
# Opcional: definir un directorio de salida
output_directory = 'PLACEHOLDER/path/to/output_plots/'
# Crear el directorio si no existe
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
output_pdf_path = os.path.join(output_directory, output_pdf_filename)

# Definir los límites del mapa
map_lon_min = -10.0
map_lon_max = 19.0
map_lat_min = 30.0  # Sur de Libia (aproximado)
map_lat_max = 48.5  # Norte de Suiza (aproximado)
map_extent = [map_lon_min, map_lon_max, map_lat_min, map_lat_max]

# Cargar el archivo netCDF
dataset = nc.Dataset(nc_file)

# Extraer variables
rr_full = dataset.variables['rr'][:, :, :]
lon_full = dataset.variables['longitude'][:]
lat_full = dataset.variables['latitude'][:]
time = dataset.variables['time'][:]

# --- Recorte inicial de datos NetCDF para optimizar (opcional pero recomendado) ---
# Encontrar los índices para latitud y longitud que corresponden a la extensión del mapa
lat_indices = np.where((lat_full >= map_lat_min) & (lat_full <= map_lat_max))[0]
lon_indices = np.where((lon_full >= map_lon_min) & (lon_full <= map_lon_max))[0]

if len(lat_indices) == 0 or len(lon_indices) == 0:
    raise ValueError("La extensión del mapa definida no se solapa con los datos del NetCDF.")

# Asegurar que los índices son contiguos si es necesario o tomar min/max
lat_start_idx, lat_end_idx = lat_indices.min(), lat_indices.max()
lon_start_idx, lon_end_idx = lon_indices.min(), lon_indices.max()

# Extraer solo la porción relevante de los datos
lat = lat_full[lat_start_idx:lat_end_idx+1]
lon = lon_full[lon_start_idx:lon_end_idx+1]
rr = rr_full[:, lat_start_idx:lat_end_idx+1, lon_start_idx:lon_end_idx+1]
# ---------------------------------------------------------------------------------

# Convertir tiempo a fechas
dates = nc.num2date(time, units=dataset.variables['time'].units, only_use_cftime_datetimes=False)

# Calcular la precipitación acumulada para cada año
years = np.arange(start_year, end_year + 1)
accumulated_rr = np.zeros((len(years), len(lat), len(lon)))

for i, year in enumerate(years):
    end_date = cftime.DatetimeGregorian(year, end_month, end_day)
    start_date = end_date - np.timedelta64(accumulated_days, 'D')
    mask_time = (dates >= start_date) & (dates < end_date) # Renombrar para evitar confusión con mask_chronology
    
    # Sumar sobre la dimensión temporal de los datos 'rr' ya recortados espacialmente
    accumulated_rr[i, :, :] = np.sum(rr[mask_time, :, :], axis=0)


# Cargar archivo de cronología
chronology_data = np.genfromtxt(chronology_file_path, skip_header=1, missing_values='NA', filling_values=np.nan)

# Filtrar los datos de cronología para los años y usar la columna 'res'
chronology_years = chronology_data[:, 0]
chronology_residuals = chronology_data[:, 2]
mask_chronology = (chronology_years >= start_year) & (chronology_years <= end_year)
filtered_chronology_residuals = chronology_residuals[mask_chronology]

# Asegurarse de que accumulated_rr y filtered_chronology_residuals tengan la misma longitud
if len(filtered_chronology_residuals) != len(years):
    print(f"Advertencia: La longitud de la cronología filtrada ({len(filtered_chronology_residuals)}) no coincide con el número de años ({len(years)}). Alineando...")
    
    common_years_indices_in_rr = []
    common_years_indices_in_chronology = []
    chronology_years_filtered = chronology_years[mask_chronology]

    for idx_chron, yr_chron in enumerate(chronology_years_filtered):
        if yr_chron in years:
            try:
                common_years_indices_in_rr.append(list(years).index(yr_chron))
                common_years_indices_in_chronology.append(idx_chron)
            except ValueError:
                pass
            
    accumulated_rr_for_corr = accumulated_rr[common_years_indices_in_rr, :, :]
    filtered_chronology_residuals_for_corr = filtered_chronology_residuals[common_years_indices_in_chronology]
    
    if len(filtered_chronology_residuals_for_corr) != accumulated_rr_for_corr.shape[0]:
         raise ValueError("Error crítico: Desajuste en las longitudes de las series temporales para la correlación después del intento de alineación.")
else:
    accumulated_rr_for_corr = accumulated_rr
    filtered_chronology_residuals_for_corr = filtered_chronology_residuals

# Calcular correlaciones espaciales y R
correlations = np.zeros((len(lat), len(lon)))
p_values = np.zeros((len(lat), len(lon)))

for i in range(len(lat)):
    for j in range(len(lon)):
        valid_mask = ~np.isnan(filtered_chronology_residuals_for_corr) & ~np.isnan(accumulated_rr_for_corr[:, i, j])
        if np.sum(valid_mask) > 2:
            series1 = filtered_chronology_residuals_for_corr[valid_mask]
            series2 = accumulated_rr_for_corr[valid_mask, i, j]
            if np.var(series1) > 0 and np.var(series2) > 0:
                 correlations[i, j], p_values[i, j] = pearsonr(series1, series2)
            else:
                correlations[i, j] = np.nan
                p_values[i, j] = np.nan
        else:
            correlations[i, j] = np.nan
            p_values[i, j] = np.nan

# Filtrar R significativos a nivel de significancia 0.05 (95%) y solo valores positivos
significant_r = np.where((p_values < 0.05) & (correlations > 0), correlations, np.nan)

# --- Plotear el mapa ---
fig = plt.figure(figsize=(10, 8)) # Ajustar figsize si es necesario para la nueva extensión
ax = plt.axes(projection=ccrs.PlateCarree())

# Establecer la extensión del mapa a las coordenadas definidas
ax.set_extent(map_extent, crs=ccrs.PlateCarree())

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linestyle='-', alpha=0.7, edgecolor='black', linewidth=0.5)
# Opcional: añadir gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}


# Usar la paleta YlOrRd
# Asegúrate que `lon` y `lat` usados aquí son los recortados
contourf_plot = ax.contourf(lon, lat, significant_r, cmap='YlOrRd', transform=ccrs.PlateCarree(),
                                levels=np.linspace(0, np.nanmax(significant_r) if not np.all(np.isnan(significant_r)) else 1, 11), # 11 niveles, 10 colores
                                extend='max' # Muestra una flecha si los valores exceden el máximo nivel
                               )
cbar = plt.colorbar(contourf_plot, ax=ax, orientation='vertical', fraction=0.040, pad=0.06, aspect=30)
cbar.set_label('R (Coeficiente de correlación de Pearson)', fontsize=10)
cbar.ax.tick_params(labelsize=8)

plt.xlabel('Longitud', fontsize=10)
plt.ylabel('Latitud', fontsize=10)
plt.title(f'Correlación: Precip. ({accumulated_days} días) vs Cronología ({start_year}-{end_year})\n(R > 0, p < 0.05)', fontsize=12)

# --- Guardar la figura como PDF ---
plt.savefig(output_pdf_path, format='pdf', dpi=600, bbox_inches='tight')
print(f"Mapa guardado como: {output_pdf_path}")

# Mostrar el mapa
plt.show()

# Encontrar el valor más alto de R y sus coordenadas
if np.all(np.isnan(significant_r)):
    print("No se encontraron valores de R significativos y positivos en la región seleccionada.")
    max_r_value = np.nan
    max_r_latitude = np.nan
    max_r_longitude = np.nan
else:
    max_r_value = np.nanmax(significant_r)
    max_r_index = np.unravel_index(np.nanargmax(significant_r), significant_r.shape)
    # Las coordenadas de lat y lon deben corresponder a los arrays 'lat' y 'lon' recortados
    max_r_latitude = lat[max_r_index[0]]
    max_r_longitude = lon[max_r_index[1]]
    print(f"El valor más alto de R en la región seleccionada es {max_r_value:.4f} en las coordenadas (latitud, longitud): ({max_r_latitude:.2f}, {max_r_longitude:.2f})")

# Cerrar el dataset netCDF
dataset.close()