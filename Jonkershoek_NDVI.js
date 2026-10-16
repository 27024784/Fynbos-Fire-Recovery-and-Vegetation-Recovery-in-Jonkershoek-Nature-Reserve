// ============================================================
// JONKERSHOEK FIRE RECOVERY - SENTINEL-2 NDVI TIME SERIES
// Pre-fire baseline + ~1-monthly post-fire observations
// Output: CSV for exponential recovery fitting (Wilson et al. 2015) in R
// ============================================================

// Polygon drawing:
var geometry = /* color: #d63000 */ee.Geometry.Polygon(
        [[[18.92659676134095, -33.972950202076085],
          [18.925223470325324, -33.98547687416176],
          [18.936896443958137, -33.9988555984173],
          [18.953719258899543, -34.00853251125362],
          [18.975691915149543, -34.01137845225214],
          [19.003157735462043, -34.01479345553949],
          [19.011054158801887, -34.007963311607845],
          [19.001097798938606, -33.99316278174168],
          [18.990798116321418, -33.98519219756833],
          [18.985304952258918, -33.97835967328018],
          [18.971572042102668, -33.973804351989386],
          [18.960242391223762, -33.97152659982028],
          [18.954749227161262, -33.9644082310517],
          [18.946852803821418, -33.96241498103158],
          [18.941359639758918, -33.9589978723321],
          [18.934149861926887, -33.959567399980486],
          [18.9293433433722, -33.96469297724166]]]);
          
// ------------------------------------------------------------
// 0. SETTINGS - everything you may want to change is here
// ------------------------------------------------------------
var FIRE_DATE = '2021-02-21';   // fire start/ignition date (YYYY-MM-DD)
var END_DATE  = '2026-01-01';   // "beginning of this year"
var STEP_MONTHS = 1;            // target sampling interval after the fire

var MAX_SCENE_CLOUD = 60;       // coarse scene-level filter (%). Kept loose on purpose;
                                // the real cloud screening is the pixel mask below.
var MIN_COVERAGE = 0.80;        // min. fraction of the burn scar that must be cloud-free
                                // for an image to be used (avoids biased partial means)

// Cloud masking method:
//  false = QA60 band (as in the methods text)
//  true  = Google Cloud Score+ (recommended; see note in section 3)
var USE_CLOUD_SCORE_PLUS = true;
var CS_THRESHOLD = 0.60;        // Cloud Score+ cs_cdf: higher = stricter (0.5-0.65 typical)

var SCALE = 10;                 // metres
var CRS = 'EPSG:32734';         // WGS84 / UTM 34S (covers Jonkershoek) - keeps pixel
                                // counts consistent between images

// ------------------------------------------------------------
// 1. STUDY AREA
// ------------------------------------------------------------
// Draw the burn-scar polygon with the geometry tools first. Earth Engine creates
// a variable called 'geometry'. If you named/imported it differently, change it here.
var burnScar = geometry;

Map.addLayer(burnScar, {color: 'red'}, 'Burn scar', true, 0.4);
Map.centerObject(burnScar, 13);
print('Burn scar area (ha):', burnScar.area(1).divide(10000));

var fireDate  = ee.Date(FIRE_DATE);
var startDate = fireDate.advance(-1, 'month');   // one month before the fire
var endDate   = ee.Date(END_DATE);

// ------------------------------------------------------------
// 2. LOAD SENTINEL-2 SURFACE REFLECTANCE
// ------------------------------------------------------------
// S2_SR_HARMONIZED is already surface reflectance and corrects the +1000 offset
// introduced in processing baseline 04.00 (Jan 2022), so pre- and post-2022 images
// are comparable. NDVI is a band ratio, so the 0.0001 scaling does not change it,
// but we apply it anyway so the reflectance values are correct if you reuse them.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(startDate, endDate)
  .filterBounds(burnScar)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD));

print('Sentinel-2 scenes after coarse filter:', s2.size());

// ------------------------------------------------------------
// 3. CLOUD MASK
// ------------------------------------------------------------
// NOTE: QA60 is not reliably populated for Sentinel-2 scenes acquired from
// roughly 2022 to early 2024, so a QA60 mask lets clouds through for much of your
// time series and would bias NDVI low at random dates. Cloud Score+ works for the whole
// archive. If you must stay with QA60, set USE_CLOUD_SCORE_PLUS = false and
// check the 'coverage' column to see how many bad dates slip through.

function maskQA60(image) {
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;   // opaque clouds
  var cirrusBitMask = 1 << 11;   // cirrus
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask);
}

function maskCloudScorePlus(image) {
  return image.updateMask(image.select('cs_cdf').gte(CS_THRESHOLD));
}

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var prepared = USE_CLOUD_SCORE_PLUS
  ? s2.linkCollection(csPlus, ['cs_cdf']).map(maskCloudScorePlus)
  : s2.map(maskQA60);

// ------------------------------------------------------------
// 4. NDVI (scaled to reflectance, then NDVI = (NIR - Red) / (NIR + Red))
// ------------------------------------------------------------
function calculateNDVI(image) {
  var refl = image.select(['B8', 'B4']).multiply(0.0001);
  return refl.normalizedDifference(['B8', 'B4'])
    .rename('NDVI')
    .copyProperties(image, ['system:time_start']);
}

var ndviScenes = prepared.map(calculateNDVI);

// ------------------------------------------------------------
// 5. MOSAIC SCENES ACQUIRED ON THE SAME DAY
// ------------------------------------------------------------
// Neighbouring tiles / overlapping orbits give several images for one date.
// Mosaicking them gives one observation per day, so the burn scar is not split
// across two partial images.
var dayList = ee.List(
  ndviScenes.map(function (img) {
    return ee.Feature(null, {day: img.date().format('YYYY-MM-dd')});
  }).aggregate_array('day')
).distinct();

var ndviDaily = ee.ImageCollection(dayList.map(function (d) {
  d = ee.String(d);
  var dayStart = ee.Date(d);
  return ndviScenes
    .filterDate(dayStart, dayStart.advance(1, 'day'))
    .mosaic()
    .set({'system:time_start': dayStart.millis(), 'date': d});
}));

print('Daily NDVI observations:', ndviDaily.size());

// ------------------------------------------------------------
// 6. MEAN NDVI OVER THE BURN SCAR FOR EACH DATE
// ------------------------------------------------------------
// Also records 'coverage' = fraction of burn-scar pixels that are cloud-free,
// so you can see how trustworthy each mean is.
var totalPixels = ee.Number(
  ee.Image.constant(1).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: burnScar,
    crs: CRS,
    scale: SCALE,
    maxPixels: 1e9
  }).get('constant')
);

var meanAndCount = ee.Reducer.mean().combine({
  reducer2: ee.Reducer.count(),
  sharedInputs: true
});

function calculateMeanNDVI(image) {
  var stats = image.select('NDVI').reduceRegion({
    reducer: meanAndCount,
    geometry: burnScar,
    crs: CRS,
    scale: SCALE,
    maxPixels: 1e9
  });
  var date = ee.Date(image.get('system:time_start'));
  var days = date.difference(fireDate, 'day');
  return ee.Feature(null, {
    'date': image.get('date'),
    'time_start': date.millis(),
    'days_since_fire': days,
    'years_since_fire': days.divide(365.25),
    'NDVI': stats.get('NDVI_mean'),
    'valid_pixels': stats.get('NDVI_count'),
    'coverage': ee.Number(stats.get('NDVI_count')).divide(totalPixels),
    'period': ee.Algorithms.If(days.lt(0), 'pre-fire', 'post-fire')
  });
}

var allObs = ndviDaily.map(calculateMeanNDVI)
  .filter(ee.Filter.notNull(['NDVI']))
  .filter(ee.Filter.gte('coverage', MIN_COVERAGE))
  .sort('time_start');

print('Usable observations (coverage >= ' + MIN_COVERAGE + '):', allObs.size());

// ------------------------------------------------------------
// 7. SELECT ONE OBSERVATION PER ~1-MONTH WINDOW
// ------------------------------------------------------------
// - One pre-fire baseline window (one month before fire up to the fire date).
// - Then consecutive STEP_MONTHS windows starting at the fire date.
// In each window we keep the image with the largest cloud-free coverage.
// Windows with no usable image simply produce no row.

function bestInWindow(winStart, winEnd) {
  return allObs
    .filter(ee.Filter.and(
      ee.Filter.gte('time_start', ee.Date(winStart).millis()),
      ee.Filter.lt('time_start', ee.Date(winEnd).millis())))
    .sort('coverage', false)
    .limit(1);
}

var preFire = bestInWindow(startDate, fireDate);

var nWindows = endDate.difference(fireDate, 'month').divide(STEP_MONTHS).ceil();
var postFireList = ee.List.sequence(0, nWindows.subtract(1)).map(function (i) {
  var winStart = fireDate.advance(ee.Number(i).multiply(STEP_MONTHS), 'month');
  var winEnd   = fireDate.advance(ee.Number(i).add(1).multiply(STEP_MONTHS), 'month');
  return bestInWindow(winStart, winEnd);
});

var series = preFire
  .merge(ee.FeatureCollection(postFireList).flatten())
  .sort('time_start');

print('Selected ~1-monthly observations:', series.size());
print('Time series:', series);

// ------------------------------------------------------------
// 8. CHART: mean NDVI vs years since fire
// ------------------------------------------------------------
// The Wilson et al. (2015) exponential model is fitted in R (nls), so no
// trendline is drawn here.
var chart = ui.Chart.feature.byFeature({
  features: series,
  xProperty: 'years_since_fire',
  yProperties: ['NDVI']
})
  .setChartType('ScatterChart')
  .setOptions({
    title: 'Mean NDVI of the Jonkershoek burn scar after the fire',
    hAxis: {title: 'Years since fire'},
    vAxis: {title: 'Mean NDVI'},
    pointSize: 5,
    legend: {position: 'none'}
  });
print(chart);

// ------------------------------------------------------------
// 9. MAP LAYERS FOR A VISUAL SANITY CHECK
// ------------------------------------------------------------
var ndviVis = {min: -0.1, max: 0.9, palette: ['brown', 'yellow', 'green']};

Map.addLayer(
  ndviDaily.filterDate(startDate, fireDate).median().clip(burnScar),
  ndviVis, 'NDVI pre-fire (median)', false);

Map.addLayer(
  ndviDaily.filterDate(endDate.advance(-12, 'month'), endDate).median().clip(burnScar),
  ndviVis, 'NDVI most recent 12 months (median)', true);

