/*Below is the JavaScript code representing the current imports. 
To transfer them to another script, paste this code into the editor 
and click "Convert" in the suggestion tooltip.*/
var S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    Points = /* color: #98ff00 */ee.Geometry.MultiPoint(
        [[11.38068168053687, 46.32613597439522],
         [11.539732737336621, 46.42933275388789],
         [11.427087313021094, 46.41464512605198],
         [11.931246837967011, 46.70352629487342]]),
    Areacolpita = /* color: #000000 */ee.Geometry.Polygon(
        [[[11.570912634706584, 46.60254905870851],
          [11.210900464512692, 46.11256008059441],
          [11.812401929356442, 46.24093368117597],
          [11.927758374668942, 46.72507960483141]]]);




function maskS2clouds_and_radiometric_scaling(image) {
  /* FUNCTION DESCRIPTION
  Function to mask the clouds and apply radiometric scaling.
  This function will be used in the "compute_NDVI" function.
  */
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).multiply(0.0001)
}



function compute_NDVI(image_collection,roi,start_time,end_time,ndvi_threshold){
  /* FUNCTION DESCRIPTION
  This function, using NDVI, compute the square kilometers of healthy vegetations
  in an area of interest  in a given interval.
  
  Parameters:
  - image_collection, a Sentine2 image collection.
  - roi (region of interest) must be a Polygon, area in which the NDVI is to be calculated.
  - start_time must be a string in format "YYYY-MM-DD",   interval in which function start to collect satellites images.
  - end_time must be a string in format "YYYY-MM-DD",  interval in which function end to collect satellites images.
  - ndvi_threshold, in order to change it faster, this allows to play with this threhold in order to find the better.
  */

  /*---------------*/
  // Filter collection for time interval, area and clouds.
  var filtered_collection = image_collection.filterDate(start_time,end_time) // Filter by time
                                     .filterBounds(roi)       // Filter by region of interest
                                     .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",50));


  // Cloud masking and radiometric scaling
  var filtered_collection = filtered_collection.map(maskS2clouds_and_radiometric_scaling);

  
  // Aggregate images and select only Red, Green, Blue and Near-Infrared bands
  var filtered_collection_aggregated = filtered_collection.median().select("B4","B3","B2","B8");
  
  
  // Compute NDVI band and rename it, remember that B8 is the near infrared band and B4 is the red band, at least in this image collection.
  // Quickly calculate NDVI using a built in function from google API
  var NDVI = filtered_collection_aggregated.normalizedDifference(["B8","B4"]).rename("NDVI"); // rename the band as NDVI
    
    
  // Now highlight the healty vegetation, that is, the one with the NDVI > ndvi threshold.
  var healty_vegetation_mask = NDVI.gt(ndvi_threshold);


  // Mask (remove) the pixels with unhealthy vegetation (NDVI < ndvi threshold), 
  // this will allow us to count all the pixels in the image
  // that have healthy vegation, and so have a measure of the healty vegetation.
  var masked_aggregated_image = filtered_collection_aggregated.updateMask(healty_vegetation_mask);
  
  // Retrieve the nominal spatial resolution of filtered collection
  var scale_of_analysis = filtered_collection.first().select("B8").projection().nominalScale(); 


  // Count the number of pixels with healthy vegetations using reduceRegion
  var unmasked_pixel_count = masked_aggregated_image.reduceRegion({
    
        // counting the number of unmasked pixels
        reducer: ee.Reducer.count(),
        
        // inside my ROI (region of interest)
        geometry: roi,
        
        // at the scale of analysis, that is the resolution of the satellite image.
        scale: scale_of_analysis,
        
        maxPixels: 1e9 // Force Google Earth Engine to process a great numberf of  pixels
  });


  // Estimate the extension of the area covered by healthy vegetation
  var extensions_of_healty_vegetation_area_m2 = ee.Number(unmasked_pixel_count.get("B8"))
                                                          .multiply(scale_of_analysis)
                                                          .multiply(scale_of_analysis);
  var extensions_of_healty_vegetation_area_km2 = extensions_of_healty_vegetation_area_m2.divide(1000)
                                                                                        .divide(1000);


  /*---------------*/
  // This part is optional, km^2 of healty vegetation already measured,
  //althougth it helps to visualize the results obtained
  
  
  // Centering the map in the region of interest
  Map.centerObject(roi);
  
  
  // Visualize the image in RGB
  var title = start_time + " / " +  end_time + " Aggregated image RGB (true color composite)";
  Map.addLayer(filtered_collection_aggregated,{bands:["B4","B3","B2"],min:0,max:0.3},title);
  
  
  // Visualize the NDVI
  var title_ndvi = start_time + " / " +  end_time + " NDVI";
  Map.addLayer(NDVI,{min:-1,max:1},title_ndvi); 


  // Visualize FALSE COLOR COMPOSITE,  pixels that contains vegetations will be in a distinct red color
  var title_false_color_composite = start_time + " / " +  end_time + " FALSE COLOR COMPOSITE per la vegetazione NIR Green Blue channels";
  Map.addLayer(filtered_collection_aggregated,{bands:["B8","B4","B3"],min:0,max:0.3},title_false_color_composite);
  
  
  // Layer to highlight the healty vegetation, that is, the one with the NDVI > ndvi threshold.
  var title_highlight_ndvi = start_time + " / " +  end_time + " healty_vegetation_mask";
  Map.addLayer(healty_vegetation_mask,{},title_highlight_ndvi);
      
  
  // Show the masking for the unhealthy vegation
  var title_masked_aggregated_image = start_time + " / " +  end_time + " masked_aggregated_image";
  Map.addLayer(masked_aggregated_image,{bands:["B4","B3","B2"],min:0,max:0.3},title_masked_aggregated_image);
  /*---------------------------*/

  // Count the number of total pixels that have been object of the study 
  var total_pixel_count = filtered_collection_aggregated.reduceRegion({reducer: ee.Reducer.count(),geometry: roi,scale: scale_of_analysis, maxPixels: 1e9 });
  var pixel_count_result = total_pixel_count.get("B8");
  var extensions_of_roi_m2 = ee.Number(pixel_count_result).multiply(scale_of_analysis).multiply(scale_of_analysis);
  var extensions_of_roi_km2 = extensions_of_roi_m2.divide(1000).divide(1000);
  print("Area of interest");
  print(extensions_of_roi_km2);
  print("Green area");
  print(extensions_of_healty_vegetation_area_km2);
  return extensions_of_healty_vegetation_area_km2

}

var ndvi_threshold = 0.5;

// Apply the compute_NDVI function and get the km^2 oh healty vegetation before the storm
var green_area_before_storm = compute_NDVI(S2,Areacolpita,"2018-06-01","2018-07-30",ndvi_threshold);


print("");

// Apply the compute_NDVI function and get the km^2 oh healty vegetation after the storm
var green_area_after_storm = compute_NDVI(S2,Areacolpita,"2019-06-01","2019-07-30",ndvi_threshold);

print("Total of green area damaged")

// Subtraction to get the km^2 of vegation area damaged
print(green_area_before_storm.subtract(green_area_after_storm));
