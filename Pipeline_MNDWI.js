var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA"),
    ROI = 
    /* color: #ffc82d */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-114.82199490319702, 36.032990821501826],
          [-114.71214020803126, 35.9907817947358],
          [-114.6352418165249, 36.09293294705135],
          [-114.44298315666965, 36.053027492700956],
          [-114.41003965321518, 36.07961634439728],
          [-114.29469217487626, 35.99744778664728],
          [-114.20441353059697, 35.99942569924499],
          [-114.11411707917281, 36.013000856384444],
          [-114.05301658215707, 36.12732034796362],
          [-114.1244218122217, 36.16169626753615],
          [-114.22877934993814, 36.041873894439554],
          [-114.36335140766104, 36.17721863347228],
          [-114.32284546617387, 36.41018167031547],
          [-114.32284110959833, 36.53938181974587],
          [-114.3947369290206, 36.423847484509714],
          [-114.47095138196966, 36.19368723484828],
          [-114.66740510425369, 36.17165999925907],
          [-114.85132363645926, 36.12162176669618]]]),
    las_vegas_lake = 
    /* color: #98ff00 */
    /* shown: false */
    ee.Geometry.Point([-114.40624286651715, 36.24827139434781]);

function water_extension(image_collection,area,roi,start_time,stop_time){
  /* FUNCTION DESCRIPTION
    This function, using MNDWI, compute the square kilometers covered by water
    in an area of interest  in a given interval.
  
    Parameters:
    - image_collection, a Sentine2 image collection.
    - area, geometry point useful to filter collection by location (filterBounds)
    - roi (region of interest) must be a Polygon, area in which the MNDWI is to be calculated.
    - start_time must be a string in format "YYYY-MM-DD",   interval in which function start to collect satellites images.
    - end_time must be a string in format "YYYY-MM-DD",  interval in which function end to collect satellites images.
  */

  
  
  // Filter the image collection by a time interval and region of interst
  // Also sort the images of the filtered collection from the less cloudy
  // image to the most cloudy image
  var L8_filtered = image_collection.filterDate(start_time,stop_time)
                                    .filterBounds(area)
                                    .sort("CLOUD_COVER");
                              
                                    
  // Get the first image, since we sorted by cloud cover the first one will be the one with less clouds
  var image = L8_filtered.first();
  
  
  // Compute the MNDWI on the less cloudy image by using the GREEN(B3)
  // and SWIR1(Shortwave infrared 1,B6) bands of the less cloudy image 
  var mndwi = image.normalizedDifference(["B3", "B6"]);


  // Create a mask for areas (pixels) covered by water by 
  // setting a threshold on the MNDWI  
  var mask_MNDWI = mndwi.gt(0.0); 
  
  
  // Use the mask to identify and extract the areas (pixels) covered by water
  // from the original image
  var masked_image = image.updateMask(mask_MNDWI);


  // Find the edges of the water areas
  var edge = ee.Algorithms.CannyEdgeDetector(mndwi, 0.99);

  
  // Retrieve the nominal spatial resolution of the 
  // original less cloudy image 
  var nominal_spatial_resolution = image.select("B3").projection().nominalScale();


  // Count the number of pixels covered by  water
  var n_pixels_from_reducer = masked_image.reduceRegion({
                                reducer: ee.Reducer.count(),
                                maxPixels: 1e19,
                                scale: nominal_spatial_resolution,
                                bestEffort: false,
                                geometry: roi ,
                                tileScale:1 //default is 1
                              });

  
  
  
  // Estimate the extension of the area covered by water
  var area_m2 =ee.Number(n_pixels_from_reducer.get("B3")).multiply(nominal_spatial_resolution).multiply(nominal_spatial_resolution);
  var area_km2 = area_m2.divide(1000).divide(1000); 


  /*---------------*/
  // This part is optional, km^2 of water already measured,
  //althougth it helps to visualize the results obtained
  
  // String that helps to format the title
  var interval = start_time + " /  "+ stop_time;
  
  Map.addLayer(image, {"bands": ["B4", "B3", "B2"],min: 0, max: 0.3}, interval + " landsat 8 true color image (RGB");
  //Map.addLayer(mndwi, {min:-1, max:1}, interval + " MNDWI"); 
  Map.centerObject(roi ,12);
  //Map.addLayer(image, {"bands": ["B6", "B5", "B3"],min: 0, max: 0.3}, interval + " landsat 8 false color image (SWIR1-NIR-G");
  Map.addLayer(mask_MNDWI ,{}, interval + " mask_MNDWI");
 // Map.addLayer(masked_image, {"bands": ["B4", "B3", "B2"], min: 0,max: 0.3}, interval + " masked image (RGB)");
 // Map.addLayer(edge.updateMask(edge), {palette:["ffffff"]},interval +" reservoir edge");
  // Map.addLayer(image_bounding_box , {}, interval +" image bounding box", false) ;
  /*---------------*/                
                        
                      
  return area_km2
}


// Get the differece between 2013 and 2023
var first_detection=water_extension(L8,las_vegas_lake,ROI,"2013-05-01","2013-07-01");
// print(first_detection);       

var second_detection=water_extension(L8,las_vegas_lake,ROI,"2023-05-01","2023-07-01");
// print(second_detection);

var difference = first_detection.subtract(second_detection);
print("Differenza between 2013 and 2023");
print(difference);



// Get the time series
var values = []; 
for (var index = 13; index <= 23; index = index +1){
  var detection=water_extension(L8,las_vegas_lake,ROI,"20"+index+"-05-01","20"+index+"-07-01");
  print("20"+index+"-05-01","20"+index+"-07-01");
  print(detection);
  values[index] = detection;
}


print(values);

