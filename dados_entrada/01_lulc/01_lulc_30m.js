var Palettes = require('users/mapbiomas/modules:Palettes.js');
var palette = Palettes.get('classification9');
var vis = {'min': 0, 'max': 69, 'palette': palette, 'format': 'png'};

// Visualizar a região de interesse (ROI) no mapa
Map.addLayer(roi, {}, "ROI");
Map.centerObject(roi, 10);
var anosDesejados = [];
for (var i = 1985; i <= 2024; i++) {
  anosDesejados.push(i);
}

// Adiciona 2013 explicitamente
// anosDesejados.push(2013);
// Coleção MapBiomas
var colecao = ee.Image('projects/mapbiomas-public/assets/brazil/lulc/collection10/mapbiomas_brazil_collection10_integration_v2').clip(roi);
anosDesejados.forEach(function(ano) {
      // Visualizar a camada no mapa
      Map.addLayer(colecao.select('classification_' + ano), vis, 'MapBiomas LULC - ' + ano, true);
  
      // Exportação para o Google Drive
      Export.image.toDrive({
        image: colecao.select('classification_' + ano), 
        description: 'MapBiomas_30m_' + ano,           
        folder: 'mapbiomas_colecao_10',                          
        fileNamePrefix: 'land_use_30m_' + ano ,      
        region: roi,                                              
        scale: 30,                                                           
        crs: 'EPSG:4326',                                                    
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
          cloudOptimized: true                                                
        }
      });
    });
