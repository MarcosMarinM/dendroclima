import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import cftime
# from matplotlib.colors import LinearSegmentedColormap # Ya no es necesario
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os # Para manejar rutas de archivo de forma más robusta

# Parámetros a modificar
nc_file = 'PLACEHOLDER/path/to/eobs01.nc'
chronology_file_path = 'PLACEHOLDER/path/to/chronology.txt'
start_year = 1950
end_year = 2022
accumulated_days = 324
end_month = 7  # Mes de finalización (Julio)
end_day = 1    # Día de finalización

# Nombre para el archivo PDF de salida
output_pdf_filename = 'mapa_correlacion_precipitacion.pdf'
# Opcional: definir un directorio de salida
output_directory = 'PLACEHOLDER/path/to/output_plots/' # Cambia esto a tu directorio deseado
# Crear el directorio si no existe
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
output_pdf_path = os.path.join(output_directory, output_pdf_filename)


# Cargar el archivo netCDF
dataset = nc.Dataset(nc_file)

# Extraer variables
rr = dataset.variables['rr'][:, :, :]
lon = dataset.variables['longitude'][:]
lat = dataset.variables['latitude'][:]
time = dataset.variables['time'][:]

# Convertir tiempo a fechas
dates = nc.num2date(time, units=dataset.variables['time'].units, only_use_cftime_datetimes=False)

# Calcular la precipitación acumulada para cada año
years = np.arange(start_year, end_year + 1)
accumulated_rr = np.zeros((len(years), len(lat), len(lon)))

for i, year in enumerate(years):
    end_date = cftime.DatetimeGregorian(year, end_month, end_day)
    start_date = end_date - np.timedelta64(accumulated_days, 'D')
    mask = (dates >= start_date) & (dates < end_date)
    accumulated_rr[i, :, :] = np.sum(rr[mask, :, :], axis=0)

# Cargar archivo de cronología
chronology_data = np.genfromtxt(chronology_file_path, skip_header=1, missing_values='NA', filling_values=np.nan)

# Filtrar los datos de cronología para los años y usar la columna 'res'
chronology_years = chronology_data[:, 0]
chronology_residuals = chronology_data[:, 2]
mask_chronology = (chronology_years >= start_year) & (chronology_years <= end_year)
filtered_chronology_residuals = chronology_residuals[mask_chronology]

# Asegurarse de que accumulated_rr y filtered_chronology_residuals tengan la misma longitud
if len(filtered_chronology_residuals) != len(years):
    print(f"Advertencia: La longitud de la cronología filtrada ({len(filtered_chronology_residuals)}) no coincide con el número de años ({len(years)}).")
    
    common_years_indices_in_rr = []
    common_years_indices_in_chronology = []
    chronology_years_filtered = chronology_years[mask_chronology]

    for idx_chron, yr_chron in enumerate(chronology_years_filtered):
        if yr_chron in years:
            try:
                common_years_indices_in_rr.append(list(years).index(yr_chron))
                common_years_indices_in_chronology.append(idx_chron)
            except ValueError: # El año de la cronología no está en 'years' (no debería pasar con el filtro actual)
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
            # Verificar que ambas series tengan varianza después de quitar NaNs
            series1 = filtered_chronology_residuals_for_corr[valid_mask]
            series2 = accumulated_rr_for_corr[valid_mask, i, j]
            if np.var(series1) > 0 and np.var(series2) > 0:
                 correlations[i, j], p_values[i, j] = pearsonr(series1, series2)
            else: # Si una serie es constante, la correlación no está bien definida o es 0/NaN
                correlations[i, j] = np.nan 
                p_values[i, j] = np.nan
        else:
            correlations[i, j] = np.nan
            p_values[i, j] = np.nan

# Filtrar R significativos a nivel de significancia 0.05 (95%) y solo valores positivos
significant_r = np.where((p_values < 0.05) & (correlations > 0), correlations, np.nan)

# --- Plotear el mapa ---
fig = plt.figure(figsize=(12, 8)) # Es buena práctica obtener el objeto Figure
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-', alpha=0.5, edgecolor='gray')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linestyle='-', alpha=0.5, edgecolor='gray')

# Usar la paleta YlOrRd
# Define los niveles para la barra de color si quieres más control
# Por ejemplo, si R varía entre 0 y 1:
# levels = np.linspace(np.nanmin(significant_r) if not np.all(np.isnan(significant_r)) else 0, 
#                      np.nanmax(significant_r) if not np.all(np.isnan(significant_r)) else 1, 
#                      11) # 11 niveles, 10 colores
# contourf_plot = ax.contourf(lon, lat, significant_r, levels=levels, cmap='YlOrRd', transform=ccrs.PlateCarree(), extend='max')

contourf_plot = ax.contourf(lon, lat, significant_r, cmap='YlOrRd', transform=ccrs.PlateCarree())
cbar = plt.colorbar(contourf_plot, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
cbar.set_label('R (Coeficiente de correlación de Pearson)')
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.title(f'Correlación: Precipitación ({accumulated_days} días) vs Cronología ({start_year}-{end_year})\n(R > 0, p < 0.05, Paleta YlOrRd)')

# --- Guardar la figura como PDF ---
# Se recomienda llamar a savefig ANTES de plt.show()
# bbox_inches='tight' recorta el espacio en blanco extra alrededor de la figura.
# dpi (dots per inch) es más relevante para formatos raster, pero es bueno incluirlo.
# format='pdf' asegura que se guarde como PDF.
plt.savefig(output_pdf_path, format='pdf', dpi=300, bbox_inches='tight')
print(f"Mapa guardado como: {output_pdf_path}")

# Mostrar el mapa
plt.show()

# Encontrar el valor más alto de R y sus coordenadas
if np.all(np.isnan(significant_r)):
    print("No se encontraron valores de R significativos y positivos.")
    max_r_value = np.nan
    max_r_latitude = np.nan
    max_r_longitude = np.nan
else:
    max_r_value = np.nanmax(significant_r)
    max_r_index = np.unravel_index(np.nanargmax(significant_r), significant_r.shape)
    max_r_latitude = lat[max_r_index[0]]
    max_r_longitude = lon[max_r_index[1]]
    print(f"El valor más alto de R es {max_r_value:.4f} en las coordenadas (latitud, longitud): ({max_r_latitude:.2f}, {max_r_longitude:.2f})")

# Cerrar el dataset netCDF
dataset.close()