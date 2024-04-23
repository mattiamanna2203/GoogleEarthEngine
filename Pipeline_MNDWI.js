var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA"),
    ROI = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[-114.83163475962702, 36.03578956834516],
          [-114.72177614345583, 35.99358057614947],
          [-114.64487505515038, 36.09573174654227],
          [-114.45260954581576, 36.05582631532806],
          [-114.41966487123592, 36.082415167826085],
          [-114.30431329794959, 36.00024660559123],
          [-114.15600412683114, 35.991358408312294],
          [-114.06262906983675, 36.13011912072826],
          [-114.11755675998677, 36.146753314995294],
          [-114.23839812399031, 36.04467270572757],
          [-114.37297495930814, 36.180017456452326],
          [-114.29607400038307, 36.57803335569396],
          [-114.35100327394889, 36.597880628431874],
          [-114.43938178648162, 36.48298157280254],
          [-114.4805787791427, 36.19648605606854],
          [-114.85684456695964, 36.123311196494164]]]);


function water_extension(image_collection,roi,start_time,stop_time){
  /* FUNCTION DESCRIPTION
    This function, using MNDWI, compute the square kilometers covered by water
    in an area of interest  in a given interval.
  
    Parameters:
    - image_collection, a Sentine2 image collection.
    - roi (region of interest) must be a Polygon, area in which the MNDWI is to be calculated.
    - start_time must be a string in format "YYYY-MM-DD",   interval in which function start to collect satellites images.
    - end_time must be a string in format "YYYY-MM-DD",  interval in which function end to collect satellites images.
  */

  
  
  // Filter the image collection by a time interval and region of interst
  // Also sort the images of the filtered collection from the less cloudy
  // image to the most cloudy image
  var L8_filtered = image_collection.filterDate(start_time,stop_time)
                                    .filterBounds(roi)
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



  // Retrieve the bounding box of the original less cloudy image using the geometry and bounds methods
  // var image_bounding_box = image.geometry().bounds(nominal_spatial_resolution);


                        
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
  Map.addLayer(mndwi, {min:-1, max:1}, interval + " MNDWI"); 
  Map.centerObject(roi ,12);
  Map.addLayer(image, {"bands": ["B6", "B5", "B3"],min: 0, max: 0.3}, interval + " landsat 8 false color image (SWIR1-NIR-G");
  Map.addLayer(mask_MNDWI ,{}, interval + " mask_MNDWI");
  Map.addLayer(masked_image, {"bands": ["B4", "B3", "B2"], min: 0,max: 0.3}, interval + " masked image (RGB)");
  Map.addLayer(edge.updateMask(edge), {palette:["ffffff"]},interval +" reservoir edge");
  Map.addLayer(image_bounding_box , {}, interval +" image bounding box", false) ;
  /*---------------*/                
                        
                      
  return area_km2
}



var first_detection=water_extension(L8,ROI,"2013-05-01","2013-06-01");
print(first_detection);       

var second_detection=water_extension(L8,ROI,"2023-05-01","2023-06-01");
print(second_detection);





var difference = first_detection.subtract(second_detection)
print(difference);

