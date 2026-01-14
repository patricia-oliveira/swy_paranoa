
// Anos de interesse
var startYear = 1985;
var endYear = 2024;

// Lista para armazenar imagens mensais por ano
var annualImages = {};

for (var year = startYear; year <= endYear; year++) {
  var monitoringYear = year.toString();
  var monthlyET = []; // Lista para imagens mensais deste ano
  
  for (var month = 1; month <= 12; month++) {
    // Configurar datas
    var startDate = monitoringYear + "-" + (month < 10 ? '0' : '') + month + '-01';
    var nextMonth = month === 12 ? 1 : month + 1;
    var endYearTemp = month === 12 ? (year + 1).toString() : monitoringYear;
    var endDate = ee.Date(endYearTemp + '-' + (nextMonth < 10 ? '0' : '') + nextMonth + '-01')
                  .advance(-1, 'day');
    
    // Obter dados mensais
    var dataset = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                  .filterDate(startDate, endDate.format('YYYY-MM-dd'))
                  .mean()
                  .clip(roi);
    
    var evapotranspiration = dataset.select('pet');
    monthlyET.push(evapotranspiration);
    
    // Exportar (mantido do seu código original)
    Export.image.toDrive({
      image: evapotranspiration,
      description: 'pet_' + year + '_' + month,
      folder: 'GEE_exports',
      fileNamePrefix: 'pet_' + year + '_' + month,
      region: roi,
      scale: 4638.3,
      crs: 'EPSG:31983',
      maxPixels: 1e13
    }); 
  } 
  
  // Calcular média anual e armazenar
  var annualMean = ee.ImageCollection(monthlyET).mean();
  annualImages[year] = annualMean;
  
  // Adicionar à visualização (opcional)
  Map.addLayer(annualMean, {min: 700, max: 1800, palette: ['blue', 'yellow', 'red']}, 
              'Média anual ' + year);
}

// Função para calcular e imprimir estatísticas anuais
function printAnnualMeans(year) {
  var image = annualImages[year];
  var stats = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 4638.3,
    maxPixels: 1e13
  });
  
  // Converter valor (TerraClimate PET está em 0.1 mm/mês)
  var petMean = ee.Number(stats.get('pet')).multiply(0.1);
  
  print('Ano ' + year + ': Média PET =', petMean, 'mm/mês');
  
  return petMean;
}

// Imprimir todas as médias anuais
print('--- MÉDIAS ANUAIS DE EVAPOTRANSPIRAÇÃO ---');
for (var year = startYear; year <= endYear; year++) {
  printAnnualMeans(year);
}
