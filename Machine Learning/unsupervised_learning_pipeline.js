/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var L7_SR_coll = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
    countries = ee.FeatureCollection("FAO/GAUL/2015/level0"),
    training_area = 
    /* color: #cfd5d6 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[56.10205724798542, 24.793825297323743],
          [56.354641529846205, 24.99732365669551],
          [56.266945969128415, 25.483217339156557],
          [56.03347041128055, 25.783359879706946],
          [55.486916672253415, 25.473299447035775],
          [55.212258469128415, 25.14306520023253],
          [55.11499852794632, 25.044166315515152],
          [56.069192062878415, 25.003750719397363]]]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
function kmeans_cluster_analysis(image_collection,training_area,roi,n_cluster,start_time,stop_time,cluster_palette){
  /* FUNCTION DESCRIPTION
    This function perform  unsupervised learning.
    In particular it perform the  kmeans cluster analysis in a given region of interest.
  
    Parameters:
    - image_collection, specify a LANDSAT image collection.
    - area, geometry point useful to filter collection by location (filterBounds)
    - roi (region of interest), must be a rectangle or polygon  area in which the cluster must be calculated.
    - n_cluster, number of cluster wanted
    - start_time must be a string in format "YYYY-MM-DD",   interval in which function start to collect satellites images.
    - end_time must be a string in format "YYYY-MM-DD",  interval in which function end to collect satellites images.
    
    Optional parameters:
    
    - cluster_palette must be a list of colors, will set up the colors of the clusters.
  */
  
  
  
  
  // Center the map in the zone of interest
  //Map.centerObject(training_area);
  
  
  
  /*-------------------------------------------------------------------------------------*/
 // Get the collection filter it and extract the median image.
  
  
  
  /*------------------------------*/
  //3) Take the LANDSAT Collection, filter it filter as necessary;
  //Filter by Date, Zone, and clouds percentage
  var L7_SR_coll_filtered = L7_SR_coll.filterDate(start_time,stop_time)
                                      .filterBounds(roi)
                                      .filter(ee.Filter.lt('CLOUD_COVER',50));
  /*------------------------------*/
  //4) Apply radiometric scaling and cloud masking to the collection;

  // Define the function that applies radiometric scaling factors.
  function applyRadiometricScaleFactors(image){
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
    return image.addBands(opticalBands, null, true)
                .addBands(thermalBand, null, true);
  }
  
  
  // Define the function that applies cloud masking
  function cloud_maskL7sr(image) {
    // Bits 3 and 4 are cloud shadow and cloud, respectively
    var cloudShadowBitMask = 1 << 3;
    var cloudsBitMask = 1 << 4;
    // Get the pixel QA band.
    var qa = image.select('QA_PIXEL');
    // Both flags should be set to zero, indicating clear conditions
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    // Return the masked image without the QA bands
    return image.updateMask(mask);
  }
  
  
  // Apply radiometric scaling and cloud masking
  L7_SR_coll_filtered= L7_SR_coll_filtered.map(applyRadiometricScaleFactors)
                                          .map(cloud_maskL7sr);
    
  /*-------------------------------------------------------------------------------------*/                                   
  // Get the spatial resolution of the image.
  // The spatial resolution ill be crucial in the training phase, since it will be required in order to
  // select a training set trought the SAMPLE function.
  var L7_spatial_resolution = L7_SR_coll_filtered.first().select('SR_B1').projection().nominalScale();

  
  // As always select the median image also use the .clip() function to select only
  // the region of interest from the image.
  // This became very useful since the landsat image are very big and maybe sometimes we want to restrict
  // our analysis. it also allow to avoid useless and heavy computation.
  
  
  /*------------------------------*/
  //5)  Select the median image, cut it out in the roi;
  var L7_SR_median_image = L7_SR_coll_filtered.median().clip(roi).select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6']);

  // Add the median image 
  Map.addLayer(L7_SR_median_image, {bands:['SR_B3', 'SR_B2', 'SR_B1'], min:0, max:0.3}, 'L7_SR_median_image RGB', true);
  /*-------------------------------------------------------------------------------------*/
  
  //print("L7_SR_median_image",L7_SR_median_image);
  
  /*-------------------------------------------------------------------------------------*/
  // Perform the cluster analysis
  
  /*------------------------------*/
  // 6) Sample from the latter image in the specified area and define a traning set
  // Get the training set trough the sample function, 
  // The sample function extract 5000 pixels from the trai
  var training_dataset = L7_SR_median_image.sample(
    {
      numPixels: 5000,
      region: training_area,
      // the spatial scale at which we select our 5000 pixels
      // to create our training dataset
      scale: L7_spatial_resolution, // should be 30 m
      geometries: true
    }
    
    );
    

  // Add the training set to the map
  //print("training_dataset",training_dataset);
  Map.addLayer(training_dataset, {}, 'training_dataset', false);
  
  
  

  /*------------------------------*/
  // 7) Define the cluster algorithm (kmeans) and train it on the traning set;
  
  // Define the cluster algorithm and specify the number of clusters.
  // Since the kmeans is an non hierarchical cluster algorithm the number of clusters must be specified.
  
  var clusterer = ee.Clusterer.wekaKMeans(n_cluster);
  
  
  
  // Train the cluster algorithm on the training data
  var trained_clusterer = clusterer.train(training_dataset);
  
  
  
  /*------------------------------*/
  //8) Apply the trained algorithm on the roi;
  // Use the trained algorithm on the whole image
  var clustered_L7_SR_median_image = L7_SR_median_image.cluster(trained_clusterer);
  

  // Add the clustered image layer
  if (cluster_palette === undefined || cluster_palette === null){
    Map.addLayer(clustered_L7_SR_median_image, {min:0, max: (n_cluster -1)}, 'clustered_L7_SR_median_image');
  
  }
  else{
    Map.addLayer(clustered_L7_SR_median_image, {palette: cluster_palette, min:0, max: (n_cluster -1) }, 'clustered_L7_SR_median_image');
  }
  



  
  var total_pixel_count = L7_SR_median_image.reduceRegion({reducer: ee.Reducer.count(),geometry: roi,scale: L7_spatial_resolution, maxPixels: 1e9 });
  var pixel_count_result = total_pixel_count.get("SR_B3");
  var extensions_of_roi_m2 = ee.Number(pixel_count_result).multiply(L7_spatial_resolution).multiply(L7_spatial_resolution);
  var extensions_of_roi_km2 = extensions_of_roi_m2.divide(1000).divide(1000);
  print("Extension area of interest",extensions_of_roi_km2);


  var total_pixel_count_training = L7_SR_median_image.reduceRegion({reducer: ee.Reducer.count(),geometry: training_area ,scale: L7_spatial_resolution, maxPixels: 1e9 });
  var pixel_count_result_training = total_pixel_count_training.get("SR_B3");
  var extensions_of_roi_m2_training = ee.Number(pixel_count_result_training).multiply(L7_spatial_resolution).multiply(L7_spatial_resolution);
  var extensions_of_roi_km2_training = extensions_of_roi_m2_training.divide(1000).divide(1000);
  print("Extension area of training",extensions_of_roi_km2_training);



  var total_pixel_count_cluster = clustered_L7_SR_median_image.reduceRegion({reducer: ee.Reducer.count(),geometry: roi,scale: L7_spatial_resolution, maxPixels: 1e9 });
  var pixel_count_result_cluster = total_pixel_count_cluster.get("cluster");
  var extensions_of_roi_m2_cluster = ee.Number(pixel_count_result_cluster).multiply(L7_spatial_resolution).multiply(L7_spatial_resolution);
  var extensions_of_roi_km2_cluster = extensions_of_roi_m2_cluster.divide(1000).divide(1000);
  print("Extension area of cluster",extensions_of_roi_km2_cluster);




  // 9) Visualize the result
  var somma = 0;
  for (var index = 0; index <= (n_cluster-1); index = index +1){
    var mask = clustered_L7_SR_median_image.eq(index);
    var count = mask.reduceRegion({reducer: ee.Reducer.sum(),
                                   geometry: roi ,
                                   scale: L7_spatial_resolution,
                                   maxPixels: 1e9}).get('cluster');
                                   
    var title = "Cluster " + index;                                  
    var area_m2  = ee.Number(count).multiply(L7_spatial_resolution).multiply(L7_spatial_resolution);
    var area_km2 = area_m2.divide(1000).divide(1000); 
    
    var somma = ee.Number(somma).add(area_km2);
    var extension = title + " extension: ";
    
    Map.addLayer(mask,{},title,false);
    print(extension,area_km2);
  }
 
  print("Sum of the area of clusters  (roi): ",somma);
}




// Abu Dhabi and Dubai clusters
// Training is the zone near to abu dhabi
// Then a zone including Abu Dhabi, Dubai, Al Ain, Al Jirab, Al Mirfa....
var palette = ['yellow',    // Sand
               'brown',     // Mountains
               'black'      // Urban area 
              ];


// Define a geometry with the United Arab Emirates borders
var country = countries.filter(ee.Filter.eq('ADM0_CODE',255));

Map.addLayer(country,['green'],'Emirati Arabi Uniti',false);




kmeans_cluster_analysis(L7_SR_coll,       // Image collection
                        training_area,    // Training area
                        country,          // Region of interest
                        3,                // cluster
                        '2023-01-01',     // Start time 
                        '2023-12-31',     // Stop time
                        palette           // Color palette
                        );






/*
Estendere l'area a tutti gli emirati arabi tranne la zona montuosa e fare una analisi di 
quanto deserto e quanta civiltà c'è.
*/
