/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var R = ee.ImageCollection("users/FJX/China_yearly_R");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
/*****************************************************************
 * Loess Plateau statistic BoxPlot
*****************************************************************/

var fensha = ee.Image("users/joeyqmf83/fensha_500m"),
            nianli = ee.Image("users/joeyqmf83/nianli_500m"),
            xisha = ee.Image("users/joeyqmf83/xisha_500m"),
            youjizhi = ee.Image("users/joeyqmf83/youjizhi_500m");


var list = ee.List([])

for(var i = 0; i < 29; i++){
  list = list.add(ee.Image('users/FJX/LS_China/LS_' + i))
}

var col = ee.ImageCollection.fromImages(list)
var merge = col.mosaic()


/***** End of imports. If edited, may not auto-convert in the playground. *****/

/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var ProvinceList = ['Anhui','Beijing','Chongqing',
'Fujian','Gansu','Guangdong','Guangxi','Guizhou','Hainan',
'Hebei','Heilongjiang','Henan','Hubei','Hunan','Jiangsu',
'Jiangxi','Jilin','Liaoning','Neimenggu','Ningxia',
'Qinghai','Shan_xi','Shandong','Shanghai','Shanxi','Sichuan',
'Tianjin', 'Xinjiang','Xizang','Yunnan',
'Zhejiang']

var blank_provinces = ['Aomen', 'Taiwan', 'Xianggang']
var slope_province = ['Neimenggu', 'Xinjiang', 'Xizang']



/*****************************************************************
 * process Landsat images
*****************************************************************/
// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              // .addBands(thermalBands, null, true);
}

// reomove cloud for Landsat-8
function cloudRemoval(image) { 
  var cloudShadowBitMask = (1 << 4); 
  var cloudsBitMask = (1 << 3); 
  var qa = image.select('QA_PIXEL'); 
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0) 
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0)); 
  var mask2 = image.select("blue").gt(0.2);                 
  return image.updateMask(mask).updateMask(mask2.not()).toDouble()
              .copyProperties(image)
              .copyProperties(image, ["system:time_start"]);
} 


// Assign a common name to the sensor-specific bands.
var LC9_BANDS = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL']; //Landsat 9
var LC8_BANDS = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL']; //Landsat 8
var LC7_BANDS = ['SR_B1','SR_B2','SR_B3', 'SR_B4','SR_B5','SR_B7','QA_PIXEL']; //Landsat 7
var LC5_BANDS = ['SR_B1','SR_B2','SR_B3', 'SR_B4','SR_B5','SR_B7','QA_PIXEL']; //Landsat 5
var S2_BANDS  = ['B2','B3','B4','B8','B11','B12']; // Sentinel-2
var STD_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2','QA_PIXEL'];


/*****************************************************************
 * calculate the K factor
*****************************************************************/
function cal_K(SAN, SIL, CLA, C, SN1) {
  var A = ((((SAN.multiply((SIL.divide(-100)).add(1))).multiply(-0.0256)).exp()).multiply(0.3)).add(0.2)
  var B = SIL.divide(CLA.add(SIL))
  var C_ = (C.multiply(0.25)).divide(C.add(((C.multiply(-2.95)).add(3.72)).exp()))
  var D = (SN1.multiply(0.7)).divide(SN1.add(((SN1.multiply(22.9)).add(-5.51)).exp()))
  var k_factor = SAN.expression("A*(B**0.3)*(1-C)*(1-D)", {
      "A": A,
      "B": B,
      "C": C_,
      "D": D,
  }).rename('K_Factor').multiply(0.1317);
  return k_factor;
}

/*****************************************************************
 * calculate the FVC and NDVI
*****************************************************************/
//*****NDVI植被管理因子计算*****
function FVC(img){
  var ndvi = img.select("NDVI");
 var fvc = img.expression(
   "108.49 * NDVI + 0.717",
   {
     "NDVI": ndvi,
   }
 ).rename('FVC').divide(100);
 return fvc;
}

function NDVI(img) {
 var nir = img.select("nir");
 var red = img.select("red");
 var ndvi = img.expression(
   "(nir - red)/(nir + red)",
   {
     "nir": nir,
     "red": red
   }
 ).rename("NDVI");
 return ndvi;
}


function cal_C_Factor(landCover, FVC, NDVI){
  var C = landCover.expression(
    'land_cover == 1 ? Cropland_cal : ' +
    'land_cover == 2 ? Forest_cal : ' +
    'land_cover == 3 ? Shrub_cal : ' +
    'land_cover == 4 ? Grass_cal : ' +
    'land_cover == 7 ? Barren_cal : ' +
    'default_formula',
    {
      'land_cover': landCover,
      'Cropland_cal': FVC.expression('1 - fvc', {'fvc': FVC.select('FVC')}),
      'Forest_cal': FVC.expression('((1-ndvi)/2)**(1+ndvi)', {'ndvi': NDVI.select('NDVI')}),
      'Shrub_cal': FVC.expression('min+(max - min) * (1 - fvc)', {'min':0.01, 'max': 0.15, 'fvc': FVC.select('FVC')}),
      'Grass_cal': FVC.expression('(1 - ndvi) / 2', {'ndvi': NDVI.select('NDVI')}),
      'Barren_cal': FVC.expression('min+(max - min) * (1 - fvc)', {'min':0.1, 'max': 0.5, 'fvc': FVC.select('FVC')}),
      'default_formula': ee.Image(0)
    }
  );
  
  return C.rename('C_Factor')
}


/***********calculate the area ***********************************/
var calculateMeanArea = function(image, feature, loc, scale) {
    var areas = image.reduceRegion({
      reducer: ee.Reducer.mean().group({
      groupField: 1,
      groupName: 'class',
    }),
      geometry: feature.geometry(),
      scale: scale,
      tileScale:4,
      maxPixels: 1e13
    })
 
    var classAreas = ee.List(areas.get('groups'))
    var classAreaLists = classAreas.map(function(item) {
      var areaDict = ee.Dictionary(item)
      var classNumber = ee.Number(
        areaDict.get('class')).format()
      var area = ee.Number(
        areaDict.get('mean'))
      return ee.List([classNumber, area])
    })
 
    var result = ee.Dictionary(classAreaLists.flatten())
    return result
}

var calculateMedianArea = function(image, feature, loc, scale) {
    var areas = image.reduceRegion({
      reducer: ee.Reducer.median().group({
      groupField: 1,
      groupName: 'class',
    }),
      geometry: feature.geometry(),
      scale: scale,
      tileScale:4,
      maxPixels: 1e13
    })
 
    var classAreas = ee.List(areas.get('groups'))
    var classAreaLists = classAreas.map(function(item) {
      var areaDict = ee.Dictionary(item)
      var classNumber = ee.Number(
        areaDict.get('class')).format()
      var area = ee.Number(
        areaDict.get('median'))
      return ee.List([classNumber, area])
    })
 
    var result = ee.Dictionary(classAreaLists.flatten())
    return result
}

var calculateP25Area = function(image, feature, loc, scale) {
    var areas = image.reduceRegion({
      reducer: ee.Reducer.percentile([25]).group({
      groupField: 1,
      groupName: 'class',
    }),
      geometry: feature.geometry(),
      scale: scale,
      tileScale:4,
      maxPixels: 1e13
    })
 
    var classAreas = ee.List(areas.get('groups'))
    var classAreaLists = classAreas.map(function(item) {
      var areaDict = ee.Dictionary(item)
      var classNumber = ee.Number(
        areaDict.get('class')).format()
      var area = ee.Number(
        areaDict.get('p25'))
      return ee.List([classNumber, area])
    })
 
    var result = ee.Dictionary(classAreaLists.flatten())
    return result
}

var calculateP75Area = function(image, feature, loc, scale) {
    var areas = image.reduceRegion({
      reducer: ee.Reducer.percentile([75]).group({
      groupField: 1,
      groupName: 'class',
    }),
      geometry: feature.geometry(),
      scale: scale,
      tileScale:4,
      maxPixels: 1e13
    })
 
    var classAreas = ee.List(areas.get('groups'))
    var classAreaLists = classAreas.map(function(item) {
      var areaDict = ee.Dictionary(item)
      var classNumber = ee.Number(
        areaDict.get('class')).format()
      var area = ee.Number(
        areaDict.get('p75'))
      return ee.List([classNumber, area])
    })
 
    var result = ee.Dictionary(classAreaLists.flatten())
    return result
}

var calculateMinArea = function(image, feature, loc, scale) {
    var areas = image.reduceRegion({
      reducer: ee.Reducer.min().group({
      groupField: 1,
      groupName: 'class',
    }),
      geometry: feature.geometry(),
      scale: scale,
      tileScale:4,
      maxPixels: 1e13
    })
 
    var classAreas = ee.List(areas.get('groups'))
    var classAreaLists = classAreas.map(function(item) {
      var areaDict = ee.Dictionary(item)
      var classNumber = ee.Number(
        areaDict.get('class')).format()
      var area = ee.Number(
        areaDict.get('min'))
      return ee.List([classNumber, area])
    })
 
    var result = ee.Dictionary(classAreaLists.flatten())
    return result
}


//切割实验区
function gridSplitMN(roi, n, m) {
  var coordinateList = roi.bounds().coordinates().flatten();
  var xmin = coordinateList.get(0).getInfo();
  var ymin = coordinateList.get(1).getInfo();
  var xmax = coordinateList.get(4).getInfo();
  var ymax = coordinateList.get(5).getInfo();
  // Get dx & dy.
  var dx = (xmax-xmin)/n;
  var dy = (ymax-ymin)/m;
  // Generate all the rects.
  var gridList=ee.List([]);
  var fid = 0;
  for(var i=0;i<n;i++)
  {
    for(var j=0;j<m;j++)
    { 
      if(i==n-1)
      {
        if(j == m-1)
        { 
          var coords=[xmin+i*dx,ymin+j*dy,xmax,ymax];
        }
        else {
          var coords=[xmin+i*dx,ymin+j*dy,xmax,ymin+(j+1)*dy]
        }
      }
      else if (j==m-1)
      {
        var coords=[xmin+i*dx,ymin+j*dy,xmin+(i+1)*dx,ymax]
      }
      else
      {
        var coords=[xmin+i*dx,ymin+j*dy,xmin+(i+1)*dx,ymin+(j+1)*dy]
      }
      var rect = ee.Feature(ee.Algorithms.GeometryConstructors.Rectangle(ee.List(coords)));
      
      // Filter.
      if(rect.intersects(roi).getInfo()) {
        var intersect = rect.intersection(roi);
        intersect = intersect.set('grid_id', fid);
        fid++;
        gridList=gridList.add(intersect);
      }
    }
  }
  // Return.
  return ee.FeatureCollection(gridList);
}

/*****************************************************************
        * provincal loop calculate
*****************************************************************/
var RList = []
var meanList = ee.List([]),
      minList = ee.List([]),
      maxList = ee.List([]),
      per_25_list = ee.List([]),
      per_75_list = ee.List([]),
      medianList = ee.List([]);
for(var year=2015; year<=2015 ;year++){
  
  
  for(var j=0;j<ProvinceList.length;j++){
    
    var loc = ProvinceList[j]
    
    
    var table = ee.FeatureCollection("users/joeyqmf83/"+loc);
    
    /*** split the roi ***/
    if(loc == 'Neimenggu'){
      var grid = gridSplitMN(table.geometry(), 4, 4)
      var gridlist = grid.getInfo()['features']
      var flag = true
    } else {
      var gridlist = [table]
      var flag = false
    }
    for(var k = 0; k < gridlist.length; k++){
      if(flag){
        var table = ee.FeatureCollection(ee.Feature(gridlist[k]))
      } else {
        var table = gridlist[k]
      }
      print('Current process is '+ loc + year+ k)
      var landCover = ee.Image('users/FJX/CLCD/CLCD_v01_'+ year +'_albert').clip(table)
        
      
      /*****************************************************************
          * R factor
      *****************************************************************/
      var L9Col = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                      .filterBounds(table)
                      .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                      .map(applyScaleFactors)
                      .select(LC9_BANDS,STD_NAMES)
                      .map(cloudRemoval);
      var L8Col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                      .filterBounds(table)
                      .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                      .map(applyScaleFactors)
                      .select(LC8_BANDS,STD_NAMES)
                      .map(cloudRemoval);
      var L7Col = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                      .filterBounds(table)
                      .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                      .map(applyScaleFactors)
                      .select(LC7_BANDS,STD_NAMES)
                      .map(cloudRemoval);
      var L5Col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                      .filterBounds(table)
                      .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                      .map(applyScaleFactors)
                      .select(LC5_BANDS,STD_NAMES)
                      .map(cloudRemoval);
                      
      var LandsatCol = ee.ImageCollection(L9Col.merge(L8Col).merge(L7Col).merge(L5Col))
      
      /** caculate the NDVI **/
      var final_image = LandsatCol.mean().clip(table);
      var image_ndvi = NDVI(final_image)
      var ndvi = image_ndvi.select('NDVI').clip(table)
      var ndvi = ndvi.where(ndvi.lt(-1), -1);
      var ndvi = ndvi.where(ndvi.gt(1), 1);
      var fvc = FVC(ndvi).clip(table)
      
      var C_Factor = cal_C_Factor(landCover, fvc, ndvi).clip(table)
      
      /*****************************************************************
          * R factor
      *****************************************************************/
      var R_Factor = R.filter(ee.Filter.date(year+'-01-01', year+'-12-31')).first().clip(table)
      
      var scale = 30
      if(scale < 11132){
        var R_Factor = R_Factor.resample('bilinear').reproject({
          crs: R_Factor.projection(),
          scale: scale
        });
      }
      R_Factor = R_Factor.clip(table)
      
      /*****************************************************************
          * K factor
      *****************************************************************/
      var SAN = xisha.clip(table)
      var SIL = fensha.clip(table)
      var CLA = nianli.clip(table)
      var C = youjizhi.clip(table)
      var SN1 = (SAN.multiply(-1).add(1)).divide(100)
      
      var K_Factor = cal_K(SAN, SIL, CLA, C, SN1).clip(table)
      
      if(scale < 500){
        var K_Factor = K_Factor.resample('bilinear').reproject({
          crs: K_Factor.projection(),
          scale: scale
        });
      }
      K_Factor = K_Factor.clip(table)
      
      /*****************************************************************
          * LS factor
      *****************************************************************/
      var LS_Factor = merge.clip(table)
    
      /*****************************************************************
          * P factor
      *****************************************************************/
      if(slope_province.indexOf(loc) > -1){
        var slope_0 = ee.Image('users/FJX/China_Slope/'+loc+'_slope_deg0');
        var slope_1 = ee.Image('users/FJX/China_Slope/'+loc+'_slope_deg1');
        var slope = ee.ImageCollection([slope_0, slope_1]).mosaic().clip(table)
      }else{
        var slope = ee.Image('users/FJX/China_Slope/'+loc+'_slope_deg').clip(table);
      }
      
      landCover = landCover.rename('lulc')
      slope = slope.rename('slope')
      var lulc_slope = landCover.addBands(slope)
      
      var P_Factor = lulc_slope.expression(
  
        "(b('lulc') == 2) ? 1" +
        ": (b('lulc') == 7) ? 1" +
        ": (b('lulc') == 3) ? 0.29" +
        ": (b('lulc') == 4) ? 0.41" +
        ": (b('slope') < 5) and (b('lulc')==1) ? 0.49" +
        ": (b('slope') < 7) and (b('lulc')==1) ? 0.59" +
        ": (b('slope') < 9) and (b('lulc')==1) ? 0.65" +
        ": (b('slope') < 12) and (b('lulc')==1) ? 0.70" +
        ": (b('slope') < 20) and (b('lulc')==1) ? 0.81" +
        ": (b('slope') < 24) and (b('lulc')==1) ? 0.95" +
        ": (b('slope') > 24) and (b('lulc')==1) ? 1" +
        ": 0"
      ).rename('P').clip(table);
    
      /*****************************************************************
          * Soil Erosion A
      *****************************************************************/
      var A = C_Factor.multiply(R_Factor).multiply(LS_Factor).multiply(K_Factor).multiply(P_Factor).clip(table).rename('A')
      
      print(loc, scale)
      
      
      var A_1 = A.where(A.lt(0), 0)
      var A_1 = A_1.where(A_1.gte(1e9), 0)
      
      var image = A_1.addBands(ee.Image(year)).clip(table)
      
      // Map.addLayer(image1)
      
      var p90 = A_1.reduceRegion({
        reducer:ee.Reducer.percentile([90]),
        geometry: table,
        scale: scale,
        maxPixels:1e13,
        tileScale:4
      })
      
      var mask = A_1.lte(ee.Number(p90.get("A")))
      
      var image1 = image.updateMask(mask)
      
      var areaMin = calculateMinArea(image1, table, loc, scale)
      var areaMedian = calculateMedianArea(image1, table, loc, scale)
      var areaP25 = calculateP25Area(image1, table, loc, scale)
      var areaMean = calculateMeanArea(image1, table, loc, scale)
      var areaP75 = calculateP75Area(image1, table, loc, scale)
      
      var mean = areaMean.get(year+"")
      var min = areaMin.get(year+"")
      var max = ee.Number(p90.get("A"))
      var p25 = areaP25.get(year+"")
      var p75 = areaP75.get(year+"")
      var median = areaMedian.get(year+"")
      
      var res = ee.Dictionary({'min': min ,
                          'per_25': p25,
                          'median':median,
                          'mean':mean,
                          'per_75': p75,
                          'max': max,
      })
      
      
      Export.table.toDrive({
        collection: ee.FeatureCollection([
          ee.Feature(null, res)
        ]),
        fileNamePrefix: year+"_"+loc+'_part'+k+'_China',
        description: year+"_"+loc+'_part'+k+'_China',
        folder:'BoxPlot_2012',
        fileFormat: 'CSV'
      });  
    }
  }
  
}
  
print('calculate finish! start export the csv')


