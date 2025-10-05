// Tiber River Water Area Monitoring Project
// ================================================================
// UMUT YILDIRIM 922797

// ======================================
// STEP 1: DEFINE STUDY AREA
// ======================================

// Tiber River region north of Rome
var studyArea = ee.Geometry.Polygon([[
  [12.30, 42.10], // Southwest
  [12.75, 42.10], // Southeast  
  [12.75, 42.55], // Northeast
  [12.30, 42.55], // Northwest
  [12.30, 42.10]  // Closing point
]]);

// Add study area to map
Map.centerObject(studyArea, 10);
Map.addLayer(studyArea, {color: 'red', fillColor: 'red'}, 'Study Area', false);

// Calculate area information dynamically
var studyAreaKm2 = studyArea.area().divide(1000000);
print('=== STUDY AREA INFORMATION ===');
print('Study Area Coordinates:', studyArea);
print('Study Area Center:', studyArea.centroid());
print('Study Area Size (km²):', studyAreaKm2);

// Set basemap to hybrid to see the river
Map.setOptions('HYBRID');

// ===============================================
// STEP 2: LOAD AND FILTER LANDSAT DATA
// ===============================================

// Define date ranges
var springSummer2022 = ee.DateRange('2022-03-01', '2022-08-31');
var springSummer2023 = ee.DateRange('2023-03-01', '2023-08-31');

// Load and filter Landsat 8 and 9 data
function loadLandsatData(dateRange, studyArea) {
  var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(dateRange.start(), dateRange.end())
    .filter(ee.Filter.lt('CLOUD_COVER', 20));
  
  var landsat9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(studyArea)
    .filterDate(dateRange.start(), dateRange.end())
    .filter(ee.Filter.lt('CLOUD_COVER', 20));
  
  return landsat8.merge(landsat9);
}

// Load 2022 and 2023 data seperately
var landsat2022 = loadLandsatData(springSummer2022, studyArea);
var landsat2023 = loadLandsatData(springSummer2023, studyArea);

// Check collection information
var imageCount2022 = landsat2022.size();
var imageCount2023 = landsat2023.size();
print('=== LANDSAT DATA STATISTICS ===');
print('2022 Spring-Summer Landsat Image Count:', imageCount2022);
print('2023 Spring-Summer Landsat Image Count:', imageCount2023);

// Cloud masking and scaling functions
function maskClouds(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBit = 1 << 3;
  var cloudShadowBit = 1 << 4;
  var mask = qa.bitwiseAnd(cloudBit).eq(0)
             .and(qa.bitwiseAnd(cloudShadowBit).eq(0));
  return image.updateMask(mask);
}

function scaleImage(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// Apply processing
var landsat2022Processed = landsat2022.map(maskClouds).map(scaleImage);
var landsat2023Processed = landsat2023.map(maskClouds).map(scaleImage);

// ==============================================
// STEP 3: CALCULATE WATER INDICES
// ==============================================

// Water indices calculation function (AWEI formulas corrected)
function calculateWaterIndices(image) {
  // NDWI (Normalized Difference Water Index)
  var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
  
  // MNDWI (Modified NDWI) - better for urban areas
  var mndwi = image.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');
  
  // AWEIsh (Automated Water Extraction Index - shadow)
  // Feyisa et al. (2014) formula
  var aweish = image.expression(
    'BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2', {
      'BLUE': image.select('SR_B2'),
      'GREEN': image.select('SR_B3'),
      'NIR': image.select('SR_B5'), 
      'SWIR1': image.select('SR_B6'),
      'SWIR2': image.select('SR_B7')
    }).rename('AWEIsh');
  
  // AWEInsh (AWEI - non-shadow)
  var aweinsh = image.expression(
    '4 * (GREEN - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)', {
      'GREEN': image.select('SR_B3'),
      'NIR': image.select('SR_B5'),
      'SWIR1': image.select('SR_B6'), 
      'SWIR2': image.select('SR_B7')
    }).rename('AWEInsh');
  
  return image.addBands([ndwi, mndwi, aweish, aweinsh]);
}

// Apply indices
var landsat2022WithIndices = landsat2022Processed.map(calculateWaterIndices);
var landsat2023WithIndices = landsat2023Processed.map(calculateWaterIndices);

// Create median composites (clipped with studyArea)
var composite2022 = landsat2022WithIndices.median().clip(studyArea);
var composite2023 = landsat2023WithIndices.median().clip(studyArea);

// ==============================================
// STEP 4: WATER AREA CLASSIFICATION
// ==============================================

// Water extraction with optimized threshold values
function extractWaterOptimized(composite, year) {
  // Optimized threshold for MNDWI (-0.1 for turbid water)
  var waterMNDWI = composite.select('MNDWI').gt(-0.1);
  
  // NDWI standard threshold
  var waterNDWI = composite.select('NDWI').gt(0);
  
  // AWEIsh threshold
  var waterAWEI = composite.select('AWEIsh').gt(0);
  
  // Combination approaches
  var waterCombined = waterMNDWI.and(waterNDWI); // More conservative
  var waterAggressive = waterMNDWI.or(waterNDWI.or(waterAWEI)); // More aggressive
  
  return {
    'MNDWI': waterMNDWI.rename('water_MNDWI_' + year),
    'NDWI': waterNDWI.rename('water_NDWI_' + year),
    'AWEI': waterAWEI.rename('water_AWEI_' + year),
    'Combined': waterCombined.rename('water_Combined_' + year),
    'Aggressive': waterAggressive.rename('water_Aggressive_' + year)
  };
}

// Extract water masks
var water2022 = extractWaterOptimized(composite2022, '2022');
var water2023 = extractWaterOptimized(composite2023, '2023');

// Visualization
var waterMaskVis = {
  min: 0,
  max: 1,
  palette: ['ffffff00', '0066ff'] // Transparent-blue
};

Map.addLayer(water2023.MNDWI, waterMaskVis, '2023 Water Mask (MNDWI)', true);
Map.addLayer(water2022.MNDWI, waterMaskVis, '2022 Water Mask (MNDWI)', false);
Map.addLayer(water2023.Combined, waterMaskVis, '2023 Water Mask (Combined)', false);

// ================================================
// STEP 5: DETAILED STATISTICS
// ================================================

// Enhanced statistics calculation function
function calculateWaterStats(waterMask, geometry, name) {
  var pixelArea = ee.Image.pixelArea();
  var waterArea = waterMask.multiply(pixelArea);
  
  var stats = waterArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geometry,
    scale: 30,
    maxPixels: 1e9,
    bestEffort: true
  });
  
  var areaM2 = ee.Number(stats.values().get(0));
  var areaKm2 = areaM2.divide(1000000);
  
  return {
    name: name,
    area_m2: areaM2,
    area_km2: areaKm2
  };
}

// Calculate area for each method
var stats2022_MNDWI = calculateWaterStats(water2022.MNDWI, studyArea, 'MNDWI_2022');
var stats2023_MNDWI = calculateWaterStats(water2023.MNDWI, studyArea, 'MNDWI_2023');
var stats2022_Combined = calculateWaterStats(water2022.Combined, studyArea, 'Combined_2022');
var stats2023_Combined = calculateWaterStats(water2023.Combined, studyArea, 'Combined_2023');

// Time series analysis
function calculateTimeSeriesStats(collection, geometry, method) {
  var timeSeriesData = collection.map(function(image) {
    // Create water mask based on method
    var waterMask;
    if (method === 'MNDWI') {
      waterMask = image.select('MNDWI').gt(-0.1);
    } else if (method === 'NDWI') {
      waterMask = image.select('NDWI').gt(0);
    } else {
      waterMask = image.select('MNDWI').gt(-0.1);
    }
    
    var pixelArea = ee.Image.pixelArea();
    var waterArea = waterMask.multiply(pixelArea);
    
    var totalArea = waterArea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 30,
      maxPixels: 1e9,
      bestEffort: true
    });
    
    return image.set({
      'water_area_m2': totalArea.values().get(0),
      'water_area_km2': ee.Number(totalArea.values().get(0)).divide(1000000),
      'date': image.date().format('YYYY-MM-dd'),
      'month': image.date().get('month')
    });
  });
  
  return timeSeriesData;
}

// Calculate time series data
var timeSeries2022 = calculateTimeSeriesStats(landsat2022WithIndices, studyArea, 'MNDWI');
var timeSeries2023 = calculateTimeSeriesStats(landsat2023WithIndices, studyArea, 'MNDWI');

// Extract statistics
var areas2022 = timeSeries2022.aggregate_array('water_area_km2');
var areas2023 = timeSeries2023.aggregate_array('water_area_km2');

// Comprehensive statistics
var detailedStats2022 = {
  mean: areas2022.reduce(ee.Reducer.mean()),
  median: areas2022.reduce(ee.Reducer.median()),
  min: areas2022.reduce(ee.Reducer.min()),
  max: areas2022.reduce(ee.Reducer.max()),
  stdDev: areas2022.reduce(ee.Reducer.stdDev()),
  count: areas2022.size()
};

var detailedStats2023 = {
  mean: areas2023.reduce(ee.Reducer.mean()),
  median: areas2023.reduce(ee.Reducer.median()),
  min: areas2023.reduce(ee.Reducer.min()),
  max: areas2023.reduce(ee.Reducer.max()),
  stdDev: areas2023.reduce(ee.Reducer.stdDev()),
  count: areas2023.size()
};

// ===============================
// STEP 6: RESULTS AND STATISTICS
// ===============================

// Change analysis
var meanDifference = ee.Number(detailedStats2023.mean).subtract(detailedStats2022.mean);
var percentChange = meanDifference.divide(detailedStats2022.mean).multiply(100);

// Dynamic summary statistics
var summaryStats = ee.Dictionary({
  'study_area_km2': studyAreaKm2,
  'analysis_period': 'March-August (Spring-Summer)',
  'landsat_images_2022': imageCount2022,
  'landsat_images_2023': imageCount2023,
  'mean_water_area_2022_km2': detailedStats2022.mean,
  'mean_water_area_2023_km2': detailedStats2023.mean,
  'median_water_area_2022_km2': detailedStats2022.median,
  'median_water_area_2023_km2': detailedStats2023.median,
  'min_water_area_2022_km2': detailedStats2022.min,
  'max_water_area_2022_km2': detailedStats2022.max,
  'min_water_area_2023_km2': detailedStats2023.min,
  'max_water_area_2023_km2': detailedStats2023.max,
  'std_dev_2022_km2': detailedStats2022.stdDev,
  'std_dev_2023_km2': detailedStats2023.stdDev,
  'percent_change': percentChange,
  'absolute_change_km2': meanDifference,
  'sample_count_2022': detailedStats2022.count,
  'sample_count_2023': detailedStats2023.count
});

// Print results
print('=== FINAL RESULTS ===');
print('Method: MNDWI-based water classification (threshold: -0.1)');
print('');
print('MAIN RESULTS:');
print('• Average water area change:', percentChange, '%');
print('• 2022 Average water area:', detailedStats2022.mean, 'km²');
print('• 2023 Average water area:', detailedStats2023.mean, 'km²');
print('• Absolute change:', meanDifference, 'km²');
print('');
print('All Statistics:', summaryStats);

// ===============================
// STEP 7: EXPORT AND VISUALIZATION
// ===============================

// Water change visualization
var waterChange = water2023.MNDWI.subtract(water2022.MNDWI);
var changeVis = {
  min: -1,
  max: 1,
  palette: ['red', 'white', 'blue'] // Red: water loss, Blue: water gain
};
Map.addLayer(waterChange, changeVis, 'Water Area Change (2023-2022)', false);

// Final visualization layers
Map.addLayer(composite2023, {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0, max: 0.3, gamma: 1.2
}, '2023 True Color', false);

// ===============================
// EXPORT FUNCTIONS
// ===============================

// 1. Export water masks
Export.image.toDrive({
  image: water2023.MNDWI,
  description: 'Tiber_Water_Mask_2023_MNDWI',
  folder: 'Tiber_Water_Analysis',
  scale: 30,
  region: studyArea,
  maxPixels: 1e9,
  crs: 'EPSG:4326'
});

Export.image.toDrive({
  image: water2022.MNDWI,
  description: 'Tiber_Water_Mask_2022_MNDWI',
  folder: 'Tiber_Water_Analysis',
  scale: 30,
  region: studyArea,
  maxPixels: 1e9,
  crs: 'EPSG:4326'
});

// 2. Export water change map
Export.image.toDrive({
  image: waterChange,
  description: 'Tiber_Water_Change_2023vs2022',
  folder: 'Tiber_Water_Analysis',
  scale: 30,
  region: studyArea,
  maxPixels: 1e9,
  crs: 'EPSG:4326'
});

// 3. Export RGB composites
Export.image.toDrive({
  image: composite2023.select(['SR_B4', 'SR_B3', 'SR_B2']),
  description: 'Tiber_RGB_2023',
  folder: 'Tiber_Water_Analysis',
  scale: 30,
  region: studyArea,
  maxPixels: 1e9,
  crs: 'EPSG:4326'
});

// 4. Export time series data as CSV
var timeSeriesFeatures2022 = ee.FeatureCollection(timeSeries2022.map(function(img) {
  return ee.Feature(null, {
    'date': img.get('date'),
    'month': img.get('month'),
    'water_area_km2': img.get('water_area_km2'),
    'year': 2022
  });
}));

var timeSeriesFeatures2023 = ee.FeatureCollection(timeSeries2023.map(function(img) {
  return ee.Feature(null, {
    'date': img.get('date'),
    'month': img.get('month'),
    'water_area_km2': img.get('water_area_km2'),
    'year': 2023
  });
}));

var combinedTimeSeries = timeSeriesFeatures2022.merge(timeSeriesFeatures2023);

Export.table.toDrive({
  collection: combinedTimeSeries,
  description: 'Tiber_Water_TimeSeries_Data',
  folder: 'Tiber_Water_Analysis',
  fileFormat: 'CSV'
});

// 5. Export summary statistics
var summaryFeature = ee.Feature(null, summaryStats);
var summaryCollection = ee.FeatureCollection([summaryFeature]);

Export.table.toDrive({
  collection: summaryCollection,
  description: 'Tiber_Water_Summary_Statistics',
  folder: 'Tiber_Water_Analysis',
  fileFormat: 'CSV'
});

// ===============================
// STEP 8: HISTOGRAM ANALYSIS
// ===============================

print('');
print('=== HISTOGRAM ANALYSIS ===');

// Create histograms for water areas
var histogram2022 = ui.Chart.array.values(areas2022, 0, null)
  .setChartType('ColumnChart')
  .setOptions({
    title: '2022 Water Area Distribution (km²)',
    hAxis: {title: 'Observation Number'},
    vAxis: {title: 'Water Area (km²)'},
    colors: ['#1f77b4'],
    legend: {position: 'none'}
  });

var histogram2023 = ui.Chart.array.values(areas2023, 0, null)
  .setChartType('ColumnChart')
  .setOptions({
    title: '2023 Water Area Distribution (km²)',
    hAxis: {title: 'Observation Number'},
    vAxis: {title: 'Water Area (km²)'},
    colors: ['#ff7f0e'],
    legend: {position: 'none'}
  });

// Display histograms
print('2022 Water Area Histogram:', histogram2022);
print('2023 Water Area Histogram:', histogram2023);

var combinedData = areas2022.map(function(area) {
  return [area, '2022'];
}).cat(areas2023.map(function(area) {
  return [area, '2023'];
}));

// Create MNDWI value histograms for water pixels (for comparison)
var waterIndexHistogram2022 = composite2022.select('MNDWI').updateMask(water2022.MNDWI);
var waterIndexHistogram2023 = composite2023.select('MNDWI').updateMask(water2023.MNDWI);

var histogram2022Index = ui.Chart.image.histogram({
  image: waterIndexHistogram2022,
  region: studyArea,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: '2022 Water Areas MNDWI Value Distribution',
  hAxis: {title: 'MNDWI Value'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#1f77b4']
});

var histogram2023Index = ui.Chart.image.histogram({
  image: waterIndexHistogram2023,
  region: studyArea,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: '2023 Water Areas MNDWI Value Distribution',
  hAxis: {title: 'MNDWI Value'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#ff7f0e']
});

print('2022 MNDWI Value Histogram:', histogram2022Index);
print('2023 MNDWI Value Histogram:', histogram2023Index);

// Monthly trend analysis
var monthlyStats2022 = timeSeries2022.aggregate_array('month');
var monthlyStats2023 = timeSeries2023.aggregate_array('month');

print('');
print('=== MONTHLY TREND ANALYSIS ===');
print('2022 Monthly Distribution:', monthlyStats2022);
print('2023 Monthly Distribution:', monthlyStats2023);

Map.centerObject(studyArea, 11);