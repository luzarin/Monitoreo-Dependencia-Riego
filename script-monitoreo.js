// ====================== FASE DE ENTRADA (BASE) ======================

// 1) Área de estudio (AOI)
var aoi = ee.FeatureCollection("projects/luzarin/assets/area_estudio");

// 2) Periodo temporal
var fechaInicio = '2023-09-01';
var fechaFin    = '2024-08-31';

// 3) Sensores
var S2_COLLECTION_ID = 'COPERNICUS/S2_SR_HARMONIZED';
var S1_COLLECTION_ID = 'COPERNICUS/S1_GRD';

// 4) Productos auxiliares
var DEM_ID = 'USGS/SRTMGL1_003';

// 5) Parámetros globalesp
var ESCALA = 10; // 10 m

// Verificación visual básica
Map.centerObject(aoi, 9);
Map.addLayer(aoi, {color: 'yellow'}, 'AOI_Asset', false);


// ====================== B1. SENTINEL-2 ==============================

// Máscara S2 por SCL (nubes/sombras/nieve) + escalar a reflectancia [0,1]
function mascaraNubesS2(img) {
  var scl = img.select('SCL');
  var good = scl.neq(1).and(scl.neq(2)).and(scl.neq(3)) // saturado, oscuro, sombra
                 .and(scl.neq(8)).and(scl.neq(9))       // nubes
                 .and(scl.neq(10)).and(scl.neq(11));    // cirrus, nieve
  return img.updateMask(good)
            .divide(10000)
            .copyProperties(img, ['system:time_start']);
}

// Índices S2 (NDVI, NDWI, NBR2)
function addS2Indices(img){
  var ndvi = img.select('B8').subtract(img.select('B4'))
                .divide(img.select('B8').add(img.select('B4'))).rename('NDVI');
  var ndwi = img.select('B3').subtract(img.select('B8'))
                .divide(img.select('B3').add(img.select('B8'))).rename('NDWI'); // McFeeters
  var nbr2 = img.select('B12').subtract(img.select('B11'))
                .divide(img.select('B12').add(img.select('B11'))).rename('NBR2');
  return img.addBands([ndvi, ndwi, nbr2]);
}

// Colección S2 filtrada, enmascarada, recortada, con índices
var S2_collection = ee.ImageCollection(S2_COLLECTION_ID)
  .filterBounds(aoi)
  .filterDate(fechaInicio, fechaFin)
  .map(mascaraNubesS2)
  .map(function(i){ return i.clip(aoi); })
  .map(addS2Indices);

// Estadísticos temporales S2 (median, max, min, stdDev)
var reducersS2 = ee.Reducer.median()
  .combine(ee.Reducer.max(),   '', true)
  .combine(ee.Reducer.min(),   '', true)
  .combine(ee.Reducer.stdDev(),'', true);

var S2_stats = S2_collection
  .select(['B2','B3','B4','B8','B11','B12','NDVI','NDWI','NBR2'])
  .reduce(reducersS2);

// Rango fenológico NDVI (max - min)
var S2_ndvi_range = S2_collection.select('NDVI')
  .reduce(ee.Reducer.max().combine(ee.Reducer.min(), '', true))
  .expression('b("NDVI_max") - b("NDVI_min")')
  .rename('NDVI_range');

// Salidas B1 (no se apilan con S1)
print('S2 imágenes válidas:', S2_collection.size());
Map.addLayer(S2_collection.median().select(['B4','B3','B2']), {min:0, max:0.3}, 'S2 RGB median', false);
Map.addLayer(S2_stats.select('NDVI_median'), {min:0, max:1}, 'S2 NDVI_median', false);


// ====================== B2. SENTINEL-1 ==============================

// Filtro básico S1 GRD IW, VV y VH (ambas órbitas)
var S1_collection = ee.ImageCollection(S1_COLLECTION_ID)
  .filterBounds(aoi)
  .filterDate(fechaInicio, fechaFin)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .map(function(i){ return i.clip(aoi); });

// dB a potencia lineal
function dbToLin(img, band){
  var x = img.select(band);
  return ee.Image(10).pow(x.divide(10)).rename(band + '_lin');
}

// Aproximación gamma0: gamma0 = sigma0 / cos(ángulo incidencia)
function toGamma0(img){
  var vv_lin = dbToLin(img, 'VV');
  var vh_lin = dbToLin(img, 'VH');
  var inc = img.select('angle').multiply(Math.PI/180);
  var cosi = inc.cos();
  var vv_g0 = vv_lin.divide(cosi).rename('VV_g0');
  var vh_g0 = vh_lin.divide(cosi).rename('VH_g0');
  return img.addBands([vv_g0, vh_g0], null, true);
}

// Suavizado Lee simple (3x3) en dominio lineal con reduceNeighborhood
function leeSimple(img, band){
  var nat = img.select(band);
  var kernel3 = ee.Kernel.square({radius: 1, units: 'pixels', normalize: false});
  var mean = nat.reduceNeighborhood({reducer: ee.Reducer.mean(), kernel: kernel3});
  var variance = nat.reduceNeighborhood({reducer: ee.Reducer.variance(), kernel: kernel3});
  var ci = variance.divide(mean.multiply(mean).add(1e-12));
  var W = ee.Image(1).divide(ee.Image(1).add(ci));
  var out = nat.multiply(W).add(mean.multiply(ee.Image(1).subtract(W)))
               .rename(band + '_lee');
  return out;
}

// Procesamiento S1: gamma0 + Lee + ratio VH/VV
var S1_proc = S1_collection
  .map(toGamma0)
  .map(function(i){
    var vv = leeSimple(i, 'VV_g0');
    var vh = leeSimple(i, 'VH_g0');
    var ratio = vh.divide(vv).rename('VVVH_ratio_lee'); // VH/VV
    return i.addBands([vv, vh, ratio], null, true);
  });

// Estadísticos temporales S1 (median, stdDev, p10/p50/p90) sobre VV/VH/ratio suavizados
var S1_stats = S1_proc
  .select(['VV_g0_lee','VH_g0_lee','VVVH_ratio_lee'])
  .reduce(ee.Reducer.median()
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.percentile([10,50,90]), '', true));


// Salidas B2
print('S1 imágenes válidas:', S1_collection.size());
Map.addLayer(S1_stats.select('VV_g0_lee_median'), {min:5e-4, max:2e-2}, 'S1 VV_median (gamma0, lee)', false);
Map.addLayer(S1_stats.select('VVVH_ratio_lee_median'), {min:0.2, max:0.8}, 'S1 ratio_median (VH/VV)', false);


// ====================== FASE C. CONSTRUCCIÓN DEL STACK ======================

// 1) DEM y pendiente
var dem = ee.Image(DEM_ID).clip(aoi);
var slope = ee.Terrain.slope(dem).rename('slope');

// 2) Definir CRS de salida (Zona 19 S para Chile central)
var CRS_SALIDA = 'EPSG:32719';

// 3) Combinar todas las capas en un solo raster multibanda
var stack = ee.Image.cat([
  S2_stats,
  S2_ndvi_range,
  S1_stats,
  dem.rename('elev'),
  slope
]);

// 4) Verificar estructura
print('Bandas en stack:', stack.bandNames());

// 5) Visualizaciones rápidas de control (opcional)
Map.addLayer(stack.select('NDVI_median').resample('bilinear'),
             {min:0, max:1}, 'NDVI_median', true);
Map.addLayer(stack.select('VV_g0_lee_median').resample('bilinear'),
             {min:0.02, max:0.3}, 'VV_median', false);
Map.addLayer(stack.select('elev'), {min:0, max:2000}, 'Elevación', false);


// ====================== FASE D: MUESTRAS DE ENTRENAMIENTO ======================

// Cargar puntos por clase (ya traen atributo 'clase' entero)
var c1 = ee.FeatureCollection('projects/luzarin/assets/C1');
var c2 = ee.FeatureCollection('projects/luzarin/assets/C2');
var c3 = ee.FeatureCollection('projects/luzarin/assets/C3');
var c4 = ee.FeatureCollection('projects/luzarin/assets/C4');
var c5 = ee.FeatureCollection('projects/luzarin/assets/C5');
var c6 = ee.FeatureCollection('projects/luzarin/assets/C6');
var c7 = ee.FeatureCollection('projects/luzarin/assets/C7');

// Fusionar
var training_fc = ee.FeatureCollection([c1, c2, c3, c4, c5, c6, c7]).flatten();

// Verificación mínima
print('Total puntos entrenamiento:', training_fc.size());
print('Distribución por clase:', training_fc.aggregate_histogram('clase'));
Map.addLayer(training_fc, {}, 'Muestras (todas)', false);

// Extraer valores del stack en cada punto
var muestras = stack.sampleRegions({
  collection: training_fc,
  properties: ['clase'],
  scale: ESCALA,
  tileScale: 4
});

// ====================== PREVISUALIZACIÓN LIGERA DE MUESTRAS ======================
var preview = muestras.limit(5);               // fuerza evaluación mínima
print('Preview filas:', preview);
print('Nº columnas (bandas):', stack.bandNames().size());

// ====================== BALANCEO DE CLASES PARA ENTRENAMIENTO ======================
var TARGET = 300;  // objetivo por clase

// Añade columna aleatoria estable
var withRand = muestras.randomColumn('rand', 1);

// Helper para recortar por clase
function capClass(code, n){
  var fc = withRand.filter(ee.Filter.eq('clase', code)).sort('rand');
  var size = fc.size();
  return ee.FeatureCollection(ee.Algorithms.If(size.lte(n), fc, fc.limit(n)));
}

// Concatenar clases
var trainBalanced = ee.FeatureCollection([
  capClass(1, TARGET),
  capClass(2, TARGET),
  capClass(3, TARGET),
  capClass(4, TARGET),
  capClass(5, TARGET),
  capClass(6, TARGET),
  capClass(7, TARGET)
]).flatten();

print('Tamaño trainBalanced:', trainBalanced.size());
print('Distribución balanceada:', trainBalanced.aggregate_histogram('clase'));

// ====================== SPLIT TRAIN/TEST (70/30) ======================
var wb = trainBalanced.randomColumn('split', 2);
var trainSet = wb.filter(ee.Filter.lt('split', 0.7));
var testSet  = wb.filter(ee.Filter.gte('split', 0.7));

print('Train size:', trainSet.size());
print('Test size:',  testSet.size());

// ====================== ENTRENnAR RF (sin variablesPerSplit) ======================
var bands = stack.bandNames();

var rf = ee.Classifier.smileRandomForest(300)  // 300, split raiz p por defecto
  .setOutputMode('CLASSIFICATION')
  .train({
    features: trainSet,
    classProperty: 'clase',
    inputProperties: bands
  });

// ====================== CLASIFICAR Y VALIDAR ======================
var mapa = stack.classify(rf).rename('clase');

var testPred = testSet.classify(rf);
var cm = testPred.errorMatrix('clase', 'classification');
print('Confusion matrix', cm);
print('Overall accuracy', cm.accuracy());
print('Kappa', cm.kappa());
print('Producers accuracy', cm.producersAccuracy());
print('Users accuracy', cm.consumersAccuracy());

// ====================== VISUALIZACIÓN ======================
var palette = [
  '#f2c74b', // 1 C1 agrícola perenne (amarillo)
  '#e67e22', // 2 C2 agrícola caducifolio (naranjo)
  '#1b5e20', // 3 C3 bosque esclerófilo (verde oscuro)
  '#b8860b', // 4 C4 matorral/espinal (mostaza/café)
  '#66bb6a', // 5 C5 bosque caducifolio (verde claro)
  '#9e9e9e', // 6 C6 urbano/suelo desnudo (gris)
  '#1f78b4'  // 7 C7 agua (azul)
];

Map.addLayer(mapa, {min:1, max:7, palette: palette}, 'Mapa clasificado (7 clases)', true);

// ====================== LEYENDA ======================

// Crear panel base
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Título
legend.add(ui.Label({
  value: 'Leyenda — Clasificación LCCS)',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}
}));

// Lista de clases
var classes = [
  {name: 'C1 Agrícola perenne de riego',       color: '#f2c74b'},
  {name: 'C2 Agrícola caducifolio de riego',   color: '#e67e22'},
  {name: 'C3 Bosque esclerófilo',              color: '#1b5e20'},
  {name: 'C4 Matorral o espinal',              color: '#b8860b'},
  {name: 'C5 Bosque caducifolio',              color: '#66bb6a'},
  {name: 'C6 Urbano / Suelo desnudo',          color: '#9e9e9e'},
  {name: 'C7 Agua',                            color: '#1f78b4'}
];


// Función para fila de leyenda
classes.forEach(function(item) {
  var row = ui.Panel({
    widgets: [
      ui.Label({style: {backgroundColor: item.color, padding: '8px', margin: '0 8px 0 0'}}),
      ui.Label({value: item.name, style: {fontSize: '12px'}})
    ],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
  legend.add(row);
});

// Agregar al mapa
Map.add(legend);

// ====================== EXPORT S2 RGB MEDIAN ======================
var s2_rgb_median = S2_collection
  .median()
  .select(['B4','B3','B2'])
  .rename(['R','G','B']);

Export.image.toDrive({
  image: s2_rgb_median,
  description: 'S2_RGB_median_2023-09_2024-08',
  folder: 'GEE_Exports',      // carpeta en Google Drive
  fileNamePrefix: 'S2_RGB_median_2023-09_2024-08',
  region: aoi.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});


// ====================== EXPORT LULC - LCCS ======================
Export.image.toDrive({
  image: mapa,
  description: 'LULC_Colchagua_LCCS_10m_2023-09_2024-08',
  folder: 'GEE_Exports',
  fileNamePrefix: 'LULC_Colchagua_LCCS_10m_2023-09_2024-08',
  region: aoi.geometry(),
  scale: 10,
  crs: 'EPSG:32719',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});