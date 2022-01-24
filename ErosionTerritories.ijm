/*                                                                                                                                                                                                                                                                             
* Erosion Territories                                                                                                                                                                                                                                                                                                                                      
* 
* Expected Input:
* The expected input for this script are 2D non-timelapse images in tiff/tif/ND2 format that are either 3 or 4 channel 
* One channel is expect to be a DAPI or other nuclear stain
* At start of macro, user selects folder with images in (nested folders dealth with fine), sets channel order and chooses desired number of bins
* 
* Description of script workflow:
* The DAPI channel is used to segment the nuclei (using gaussian blur and huang threshold and the 'watershed irregular features' from BioVoxxel - see installation notes)
* Each nucleus area is duplicated from original image to isolate nuclei for further processing
* For each nuclei image, the area is measured and the expected area of the bins is calculated based on the number requsted
* Using the expected areas of these bins and iteratively thresholding on a distance map of the nucleus ROI bins are produced
* Mean intensity of these ROIs is measured in all the channels                                                        
* 
* Output:
* Each nucleus found in each image will be saved as an RGB image with the outline of created bins overlaid on the image
* A .csv file will be saved with the areas of the nuclei, the bins and the mean intensity of each bin for each channel in the image
* Bin/ territory naming in the results go from outside in so bin 1 is the outermost one.
* 
* Installation Notes: 
* Download FIJI here - https://imagej.net/software/fiji/downloads
* Details of the Biovoxxel Update site - https://imagej.net/plugins/biovoxxel-toolbox
* How to add an update site - https://imagej.net/update-sites/                                                        
* How to use an ImageJ macro from Github - https://github.com/IGC-Advanced-Imaging-Resource/ImageJ_Macros
* This script is hosted here with an example image - https://github.com/IGC-Advanced-Imaging-Resource/Purshouse2022_paper       
* 																														Written by Laura Murphy (laura.murphy@ed.ac.uk)                                                                                                                                                                                                                                                                                                                             
* 																													  	IGC Advanced Imaging Resource (https://www.ed.ac.uk/institute-genetics-cancer/facilities/advanced-imaging-resource)
*                                                                                                                     	First written: November 2019  Last updated January 2022                                                                                                        
*/


//--------------------------------//-------------------------------------------------------------------------------------------------------
//-- Part 0: Preparation steps: get directories from user, set up arrays to store results that are being added to and some other setup
//--------------------------------//-------------------------------------------------------------------------------------------------------

// Get user to select folder with images for input and the output folder where they want results saved
input_folder = getDirectory("Parent directory of your image folders");
output_folder = getDirectory("Choose the folder where you want to save your results");

// Set up empty arrays that will be filled during script and displayed and saved at the end
Image_Name = newArray();
Nucleus_No = newArray();
Nucleus_Area = newArray();
Bin_No = newArray()
Bin_Area = newArray()
DAPI_MeanIntensity = newArray();
FITC_MeanIntensity = newArray();
TxRd_MeanIntensity = newArray();
Cy5_MeanIntensity = newArray();

// Magic numbers for nuclei segmentation
filter_sigma = 2; // Sigma for Gaussian blur
erosion_cycle = 1; // Erosion cycle for watershed
watershed_convexity = 0.99; // Convexity threshold for watershed
min_particle_size = 50; // Minimum particle size for Analyze Particles
max_particle_size = 700; // Maxmimum particle size for Analyze Particles
min_particle_circularity = 0.7; // Minimum circularity for Analyze Particles
max_particle_circularity = 1.00; // Maxmimum circularity for Analyze Particles

// Array of channel colours
colour_array = newArray("Blue", "Green", "Red", "Cyan");
roiManager("Reset");

// Get list of images in nested folders with certain extensions  - function at bottom of macro
directory_list = newArray();
directory_list = get_file_tree(input_folder, directory_list);

//--------------------------------//--------------------------------------------------------------------------------------------
//-- Part 1: Image loop starts, each image goes through nuclear segmentation and each nuclei is duplicated from original image
//--------------------------------//--------------------------------------------------------------------------------------------

// Start looping through images in input folders
for (img = 0; img < directory_list.length; img++) {
	path = directory_list[img];
	run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	getDimensions(width, height, channels, slices, frames);

	original_image = getTitle();
	image_name = File.nameWithoutExtension;

	// When first image is open, get user to input desired number of bins and assign channels in order of the stack
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
		
		channel_order = newArray(channels);
	
		for (ch = 0; ch < channels; ch++){
			channel_order[ch] = Dialog.getNumber;
		}                              
		
	}              

	// Get pixel sizing, if image isn't scaled then skip image
	getVoxelSize(width, height, depth, unit);
	
	if (unit == "pixels"){
		waitForUser("Image has no scaling and will be skipped");
		run("Close");
		continue;
	}

	// Set channel colours 
	for (lut = 0; lut < channel_order.length; lut++){
		selectWindow(original_image);
		Stack.setChannel(channel_order[lut]);
		run(colour_array[lut]);
	}              

	// Duplicate image to preserve for display use 
	selectWindow(original_image);
	run("Duplicate...", "title=Display duplicate");

	// Segment nuclei in image
	selectWindow(original_image);
	Stack.setChannel(channel_order[0]);
	run("Duplicate...", "title=Nuclear_Mask");
	run("Gaussian Blur...", "sigma=" + filter_sigma);
	setAutoThreshold("Huang dark");
	run("Convert to Mask"); 
	run("Watershed Irregular Features", "erosion=" + erosion_cycle + " convexity_threshold=" + watershed_convexity + " separator_size=0-Infinity");
	run("Analyze Particles...", "size=" + min_particle_size + "-" + max_particle_size + " circularity=" + min_particle_circularity + "-" + max_particle_circularity + " exclude include add");

	// Store number of nuclei
	number_of_nuclei = roiManager("Count");

	// Loop through nuclei, duplicate image of each one in image, naming them so they can be processed further	
	for (nuclei = 0; nuclei < number_of_nuclei; nuclei++){
		selectWindow("Nuclear_Mask");
		roiManager("Select", nuclei);
		title = "Nuclei_" + nuclei + 1;
		run("Duplicate...", "title=" + title);
	
		selectWindow(original_image);
		roiManager("Select", nuclei);
		run("Duplicate...", "duplicate");
		rename("Measure" + title);
		
		selectWindow("Display");
		roiManager("Select", nuclei);
		run("Duplicate...", "duplicate");
		rename("Display" + title);
	}

//--------------------------------//--------------------------------------------------------------------------------------------
//-- Part 2: Looping through isolated nuclei image starts, creation of bins
//--------------------------------//--------------------------------------------------------------------------------------------

	// Start of loop for going through the duplicated images of the solo nuclei
	for (image_of_nucleus = 0; image_of_nucleus < number_of_nuclei; image_of_nucleus++){
		title = "Nuclei_" + image_of_nucleus + 1;
		roiManager("Reset");
		selectWindow(title);

		// Create distance map to use for erosion
		run("Distance Map");
		run("Clear Outside");

		// Get area of whole nuclei and maximum value of distance map (further distance in nuclei mask from edge)
		getStatistics(whole_area, mean, min, whole_max, std, histogram);
		nuclear_area = whole_area;
		run("Select None");

		// Work out the desired area each bin should be by dividing total nuclear area by number of bins
		area_fraction = whole_area/bins;
		prebins = newArray(bins);

		// Work out the expected areas of the shells (i.e. the smaller circles that share centroid with nucleus)
		for (areas_of_prebins = 0; areas_of_prebins < prebins.length; areas_of_prebins++){
			addtoArray = whole_area;
			prebins[areas_of_prebins] = addtoArray;
			whole_area = whole_area-area_fraction;
		}

		// For each of the chosen bin number, iteratively reduce threshold of distance map until area of thresholded area is lower/equal to expected
		// Once that is reached, add to ROI Manager.
		for (bin = 0; bin < prebins.length; bin++) {
			current_area_aim = prebins[bin];
		
			for (threshold = 1; threshold < whole_max; threshold++) {
				getStatistics(area, mean, min, max, std, histogram);
				run("Select None");
				setThreshold(threshold, whole_max);
				run("Create Selection");
				getStatistics(area, mean, min, max, std, histogram);
		
				if (area <= current_area_aim){
					roiManager("Add");
						break;
						continue;
				}
		
			}
		}
				
		// Use prebins to create bins by using XOR function between the pre-bins in the ROI Manager
		number_of_rois = roiManager("Count");
		to_be_deleted = newArray(number_of_rois);

		// Add current ROIs to array to delete later
		for (bins = 0; bins < number_of_rois; bins++){
			to_be_deleted[bins] = bins;
		}


		// Loop through rois, selecting pairs to get XOR area between them, add to ROI and name them.
		for (roi = 0; roi < number_of_rois-1; roi++){
			roiManager("select", newArray(roi, roi+1));
			roiManager("XOR");
			Roi.setName("BIN_" + roi+1);
			roiManager("Add");
			roiManager("select", newArray(roi, roi+1));
		}
		
		roiManager("select", number_of_rois-1);
		Roi.setName("BIN_" + number_of_rois);
		roiManager("Add");

		// Delete the shells from earlier
		roiManager("select", to_be_deleted);
		roiManager("Delete");
	
		resetThreshold();

//--------------------------------//--------------------------------------------------------------------------------------------
//-- Part 3: Time to save results, for each bin, the image name, nucleus and nucleus area is saved to keep arrays same length
		//   Then for each of the bins (also called territories), the area is measured and then the intensity of each of the channels
		//   All things are saved to arrays which are displayed once complete. On the last image, the arrays are saved to the output folder
//--------------------------------//--------------------------------------------------------------------------------------------

		// Loop through final shapes, adding nuclei results to arrays to keep arrays even
		for (territories = 0; territories < roiManager("Count"); territories++){                        	
			Image_Name = Array.concat(Image_Name, image_name);
			Nucleus_No = Array.concat(Nucleus_No, image_of_nucleus+1);
			Nucleus_Area = Array.concat(Nucleus_Area, nuclear_area);
			Bin_No = Array.concat(Bin_No, territories+1);    
			
			selectWindow("Measure" + title);
			roiManager("Select", territories);
			getStatistics(bin_area, mean, min, max, std, histogram);
			Bin_Area = Array.concat(Bin_Area, bin_area);

			Stack.setChannel(channel_order[0]);
			roiManager("Select", territories);
			getStatistics(bin_area, mean_DAPI, min, max, std, histogram);
			DAPI_MeanIntensity = Array.concat(DAPI_MeanIntensity, mean_DAPI);
			
			Stack.setChannel(channel_order[1]);
			roiManager("Select", territories);
			getStatistics(bin_area, mean_FITC, min, max, std, histogram);
			FITC_MeanIntensity = Array.concat(FITC_MeanIntensity, mean_FITC);
			
			Stack.setChannel(channel_order[2]);
			roiManager("Select", territories);
			getStatistics(binArea, mean_TXRD, min, max, std, histogram);
			TxRd_MeanIntensity = Array.concat(TxRd_MeanIntensity, mean_TXRD);

			// If the images are 4 channel, then save the Cy5 channel
			if (channel_order.length > 3){
				Stack.setChannel(channel_order[3]);
				roiManager("Select", territories);
				getStatistics(binArea, mean_Cy5, min, max, std, histogram);
				Cy5_MeanIntensity = Array.concat(Cy5_MeanIntensity, mean_Cy5);
			}
		}

		// Take original image, convert to RGB and save with bin boundaries overlayed on image for checking.
		selectWindow("Display" + title);
		run("Make Composite");
		run("RGB Color");
		roiManager("Show all without labels");                  
		run("Flatten");
		saveAs("Tiff", output_folder + image_name + "_" + title);
		
	}

	// Set up for next image coming up by closing open images and reset ROI manager and results table
	run("Close All");
	roiManager("Reset");
	run("Clear Results");
	
	// Show updated array of results each time an image is completed. If image was 3 channel, the Cy5 array will be empty.
	Array.show(Image_Name, Nucleus_No, Nucleus_Area, Bin_No, Bin_Area, DAPI_MeanIntensity, FITC_MeanIntensity, TxRd_MeanIntensity, Cy5_MeanIntensity);
	
}

// Save table of arrays to user directed output folder now all images are processed
saveAs("Results", output_folder + File.separator + "Results.csv");
run("Close");


//--------------------------------//-----------------------------------------------------------------------
//-- Part 4: Let user know the macro is complete
//--------------------------------//-----------------------------------------------------------------------

Dialog.create("Progress"); 
Dialog.addMessage("Macro Complete!");
Dialog.show;

//--------------------------------//-----------------------------------------------------------------------
//-- Functions
//--------------------------------//-----------------------------------------------------------------------

function get_file_tree(directory, file_tree){
file_list = getFileList(directory);

for(f = 0; f < file_list.length; f++){
	if (matches(file_list[f], "(?i).*\\.(tif|tiff|nd2)$"))
		file_tree = Array.concat(file_tree, directory + file_list[f]);
	if(File.isDirectory(directory + File.separator + file_list[f]))
		file_tree = get_file_tree(directory + file_list[f], file_tree);
}
return file_tree;