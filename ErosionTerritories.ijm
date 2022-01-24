/*
* Part of FISH Macros based on initial iVision scripts written by Paul Perry and translated to IJ1 by Matthew Pearson in 2015/16                                                                                                                                                                                                
*                                                                                                                                                                                                                                                                                                                                             
* Erosion Territories                                                                                                                                                                                                                                                                                                                                      
* Nuclei segmented using gaussian blur of radius 2 and huang threshold and the 'watershed irregular features' from BioVoxxel
* For each one, area is measured and divided into a number of bins selected by user                                                                                                                                                                                                                                                                                                                                       
* Distance map is created from nuclear map                                                                                                                                                                                                                                                                                                                                        
* Different threshold values are iteratively gone through so that the number of bins are roughly equal area
* Intensity mean is measured in each of these bins for the remaining 2 or 3 channels                                                        
* Output is a spreadsheet with area and intensity means per channel and each nuclei with bins drawn on as an RGB tiff
* 
* Bins go from outside in so bin 1 is the outermost one.
* 
* Note: needs Biovoxxel update site
*                                                                                                                     Laura Murphy                                                                                                                                                                                                                                                                                                                              IGMM Advanced Imaging Resource
*                                                                                                                     November 2019                                                                                                                                                                                                
*/


//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 0: Preparation steps: get directories from user and set up arrays to store results that are being added to
//--------------------------------//-----------------------------------------------------------------------------------

// Get user input and ouput folder
inputFolder = getDirectory("Parent directory of your image folders");
outputFolder = getDirectory("Choose the folder where you want to save your results");

Image_Name = newArray();
Nucleus_No = newArray();
Nucleus_Area = newArray();
Bin_No = newArray()
Bin_Area = newArray()
DAPI_MeanIntensity = newArray();
FITC_MeanIntensity = newArray();
TxRd_MeanIntensity = newArray();
Cy5_MeanIntensity = newArray();

colArray = newArray("Blue", "Green", "Red", "Cyan");
roiManager("Reset");

// get list of images in nested folders with extension tif - function at bottom of macro
dirList = newArray();
dirList = getFileTree(inputFolder, dirList);

for (img = 0; img < dirList.length; img++) {
	path = dirList[img];
	run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Grayscale concatenate_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	getDimensions(width, height, channels, slices, frames);
	
	ori = getTitle();
	imgName = File.nameWithoutExtension;
	
	if (img == 0) {
		Dialog.create("User input");
		
		Dialog.addMessage("Enter the number of concentric territories you want the nucleus split into");
		Dialog.addNumber("Bins", 4);
		
		Dialog.addMessage("Enter the order of your channels below");
		Dialog.addNumber("DAPI", 4);
		Dialog.addNumber("FITC", 1);
		Dialog.addNumber("TxRd", 2);
		Dialog.addNumber("Cy5", 3);
		
		Dialog.show;
	
		//Store user input
		bins = Dialog.getNumber;
		
		channelOrder = newArray(channels);
	
		for (ch = 0; ch < channels; ch++){
			channelOrder[ch] = Dialog.getNumber;
		}                              
		
	}              
	
	getVoxelSize(width, height, depth, unit);
	if (unit == "pixels"){
		waitForUser("Image has no scaling and will be skipped");
		run("Close");
		continue;
	}
	
	for (lut = 0; lut < channelOrder.length; lut++){
		selectWindow(ori);
		Stack.setChannel(channelOrder[lut]);
		run(colArray[lut]);
	}              
	
	selectWindow(ori);
	run("Duplicate...", "title=Display duplicate");
	
	selectWindow(ori);
	Stack.setChannel(channelOrder[0]);
	run("Duplicate...", "title=Nuclear_Mask");
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold("Huang dark");
	run("Convert to Mask"); 
	
	// These lines may need changed between users/ images.
	run("Watershed Irregular Features", "erosion=1 convexity_threshold=0.99 separator_size=0-Infinity");
	run("Analyze Particles...", "size=50-700 circularity=0.7-1.00 exclude include add");
	
	nNuclei = roiManager("Count");
	
	for (nuclei = 0; nuclei < nNuclei; nuclei++){
		selectWindow("Nuclear_Mask");
		roiManager("Select", nuclei);
		title = "Nuclei_" + nuclei + 1;
		run("Duplicate...", "title=" + title);
	
		selectWindow(ori);
		roiManager("Select", nuclei);
		run("Duplicate...", "duplicate");
		rename("Measure" + title);
		
		selectWindow("Display");
		roiManager("Select", nuclei);
		run("Duplicate...", "duplicate");
		rename("Display" + title);
	}
	
	for (nucImg = 0; nucImg < nNuclei; nucImg++){
		title = "Nuclei_" + nucImg + 1;
		roiManager("Reset");
		selectWindow(title);
		
		run("Distance Map");
		run("Clear Outside");
		
		getStatistics(wholeArea, mean, min, wholeMax, std, histogram);
		nucArea = wholeArea;
		run("Select None");
		
		fraction = wholeArea/bins;
		areas = newArray(bins);
		
		for (area = 0; area< areas.length; area++){
			addtoArray = wholeArea;
			areas[area] = addtoArray;
			wholeArea = wholeArea-fraction;
		}
		
		for (b = 0; b < areas.length; b++) {
			currentArea = areas[b];
		
			for (t = 1; t < wholeMax; t++) {
				getStatistics(area, mean, min, max, std, histogram);
				run("Select None");
				setThreshold(t, wholeMax);
				run("Create Selection");
				getStatistics(area, mean, min, max, std, histogram);
		
				if (area <= currentArea){
					roiManager("Add");
						break;
						continue;
				}
		
			}
			k = max;
		}
	
		nROIs = roiManager("Count");
		toDelete = newArray(nROIs);
		
		for (bins = 0; bins < nROIs; bins++){
			toDelete[bins] = bins;
		}
		
		for (r = 0; r < nROIs-1; r++){
			roiManager("select", newArray(r, r+1));
			roiManager("XOR");
			Roi.setName("BIN_" + r+1);
			roiManager("Add");
			roiManager("select", newArray(r, r+1));
		}
		
		roiManager("select", nROIs-1);
		Roi.setName("BIN_" + nROIs);
		roiManager("Add");
		
		roiManager("select", toDelete);
		roiManager("Delete");
	
		resetThreshold();
		
		for (territories = 0; territories < roiManager("Count"); territories++){                        	
			Image_Name = Array.concat(Image_Name, imgName);
			Nucleus_No = Array.concat(Nucleus_No, nucImg+1);
			Nucleus_Area = Array.concat(Nucleus_Area, nucArea);
			Bin_No = Array.concat(Bin_No, territories+1);    
			
			selectWindow("Measure" + title);
			roiManager("Select", territories);
			getStatistics(binArea, mean, min, max, std, histogram);
			
			Bin_Area = Array.concat(Bin_Area, binArea);

			Stack.setChannel(channelOrder[0]);
			roiManager("Select", territories);
			getStatistics(binArea, DAPImean, min, max, std, histogram);
			DAPI_MeanIntensity = Array.concat(DAPI_MeanIntensity, DAPImean);
			
			Stack.setChannel(channelOrder[1]);
			roiManager("Select", territories);
			getStatistics(binArea, FITCmean, min, max, std, histogram);
			FITC_MeanIntensity = Array.concat(FITC_MeanIntensity, FITCmean);
			
			Stack.setChannel(channelOrder[2]);
			roiManager("Select", territories);
			getStatistics(binArea, TxRdmean, min, max, std, histogram);
			TxRd_MeanIntensity = Array.concat(TxRd_MeanIntensity, TxRdmean);

			if (channelOrder.length == 4){
				Stack.setChannel(channelOrder[3]);
				roiManager("Select", territories);
				getStatistics(binArea, Cy5mean, min, max, std, histogram);
				Cy5_MeanIntensity = Array.concat(Cy5_MeanIntensity, Cy5mean);
			}
		}
		
		selectWindow("Display" + title);
		run("Make Composite");
		run("RGB Color");
		roiManager("Show all without labels");                  
		run("Flatten");
		saveAs("Tiff", outputFolder + imgName + "_" + title);
		
	}
	
	run("Close All");
	roiManager("Reset");
	run("Clear Results");
	
	Array.show(Image_Name, Nucleus_No, Nucleus_Area, Bin_No, Bin_Area, DAPI_MeanIntensity, FITC_MeanIntensity, TxRd_MeanIntensity, Cy5_MeanIntensity);
	
}

saveAs("Results", outputFolder + File.separator + "Results.csv");
run("Close");


//--------------------------------//-----------------------------------------------------------------------
//-- Part 5: Let user know the macro is complete
//--------------------------------//-----------------------------------------------------------------------

Dialog.create("Progress"); 
Dialog.addMessage("Macro Complete!");
Dialog.show;

//--------------------------------//-----------------------------------------------------------------------
//-- Functions
//--------------------------------//-----------------------------------------------------------------------

function getFileTree(dir , fileTree){
list = getFileList(dir);

for(f = 0; f < list.length; f++){
	if (matches(list[f], "(?i).*\\.(tif|tiff|nd2|lif|ndpi|mvd2|ims|oib)$"))
		fileTree = Array.concat(fileTree, dir + list[f]);
	if(File.isDirectory(dir + File.separator + list[f]))
		fileTree = getFileTree(dir + list[f],fileTree);
}
return fileTree;
