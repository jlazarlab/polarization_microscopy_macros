var composite_img_id = 0;
var w=0; // image width
var h=0; // image height

var bleedthrough = 0.13; //0.05 = 5% - NH; 0.18 - Prague
var pol_modulation = "columns";
var first_pixel_pol = "horizontal";
var pixel_shape = "square";
var two_pol_color_scheme = "red-green";
var multipol_color_scheme = "white-red";
var orig_image_close = false;
var range_lock = false;
var locked_rg_range = 1.3;
var make_16_bit = false;
var process_whole_stack = false;
var ratiometric_fitting = true;
var look_for_pol_angles = 2;		        //number of phases that should be determined by fitting LD data
var phase = 0;                              //polarization angle error
var phase_choice = "Identical, found from fit";
var phase_fixed = false;
var phase1_deg = 0;
var phase2_deg = 0;
var polarization_direction = "horizontal";  //polarization used in single-polarization images
var polarization_purity = 100;              //for single polarization only, in %
var polarization_purity_hor = 100;          //for single polarization only, in %
var polarization_purity_ver = 100;          //for single polarization only, in %
var mixed_image = true;

var do_segmentation = false;
var membranes_or_filaments = "Analyze a membrane";
var look_for_vesicles = false;
var	determine_rmax_alpha0_sigma = "none";
var cutoff_distance_from_snake = 0.75;
var save_files = true;

// special character definitions
var thetachar     = fromCharCode(0x03b8);
var alpha0char    = fromCharCode(0x03b1)+fromCharCode(0x2080);
var sigmachar     = fromCharCode(0x03c3);
var muchar        = fromCharCode(0x03bc);
var degreechar    = fromCharCode(0x00b0);
var plusminuschar = fromCharCode(0x00b1);
var subtwochar    = fromCharCode(0x2082);
var rmaxchar      = "r"+fromCharCode(0x2098)+fromCharCode(0x2090)+fromCharCode(0x2093);
var suptwochar    = fromCharCode(0x00b2);
var backspacechar = fromCharCode(0x0008);
var rightarrowchar     = fromCharCode(0x2192);


//'t' macro variables
var rmax_1PPM = 0;
var B_1P_fit = 0;
var B_2P_fit = 0;
var C_2P_fit = 0;
var alpha0_interval = 2;                 // sampling interval, in degrees
var sigma_interval = 2;                  // sampling interval, in degrees
var composition_percentage_interval = 5; // sampling interval, in degrees
var	min_alpha0_major = 0;
var	max_alpha0_major = 90;
var	min_sigma_major = 0;
var	max_sigma_major = 90;
var	min_alpha0_minor = 0;
var	max_alpha0_minor = 90;
var	min_sigma_minor = 0;
var	max_sigma_minor = 90;

var property_kind_1P = "log2(rmax)";
var property_kind_2P = "log2(rmax)";
var property_mean_1P = 0;
var property_mean_2P = 0;
var interval_value_1P = 10;
var interval_value_2P = 10;
var interval_kind_1P = "%";
var interval_kind_2P = "%";

var rmax_2P = 0;
var log2rmax_1P = 0;
var log2rmax_2P = 0;
var interval_1P = 0;
var interval_2P = 0;
var interval_1P_perc = 10;
var interval_2P_perc = 10;

// multipolarization image processing variables
var starting_pol_angle = 0;
var bkgd = 0;
var pol_angle_increment = 45;
var F_cutoff_for_LD_display = 0.2;
var F_cutoff_for_azimuth_display = 0.2;
var azimuth_pixel_bin_size = 10;
var show_LD = "Crude/fast";
var show_azimuths = "Yes";
var azimuth_color = "direction";
var azimuth_scaling_factor = 1.0;
var azimuth_size_lock = false;
var use_fitting_for_LD = false;
var show_additional_plots = false;




macro "Process a polarization image (mixed image) [g]"{

	//operating_system = getInfo("os.name");
	number_of_polarizations = 2;
	mixed_image = true;
	
	raw_image_name = getTitle();
	raw_image_id = getImageID();
	number_of_slices = nSlices;
	file_name_base = substring(raw_image_name,0,lengthOf(raw_image_name)-4);

	w = getWidth;
	h = getHeight;
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	composite_img_id = 0;

	read_params_from_config_file();   // reads parameters from a config file; stores them in global vars
	color_scheme = two_pol_color_scheme;
	nonguessable_params = newArray(pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, polarization_purity_hor, polarization_purity_ver);



//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------- process a single slice ----------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------


	setBatchMode(true);

	selectImage(raw_image_id);
	run("Duplicate...", " ");  // duplicates an image or a single slice from a stack
	rename(raw_image_name);  
	raw_slice_image_id = getImageID();


	raw_composite_img_id = turn_RPM_image_into_a_composite(raw_slice_image_id, bleedthrough, pol_modulation, first_pixel_pol, pixel_shape);

	guessed_params = guess_params(raw_composite_img_id); // returns hor_bkgd, ver_bkgd, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj

	processed_composite_img_id = make_processed_composite(raw_composite_img_id, guessed_params); //creates a composite image (slice 1 = Fh, slice 2 = Fv, by using the guessed processing parameters

	processing_params = Array.concat(nonguessable_params, guessed_params);  //pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, polarization_purity_hor, polarization_purity_ver, composite_img_id, hor_min, ver_min, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj

	color_scheme = processing_params[4];
	bright_adj = processing_params[9];
	rg_range = processing_params[10];
	rg_image_id = make_rg_image(processed_composite_img_id, bright_adj, rg_range, color_scheme);  //creates a two-color (RGB) image showing fluorescence intensity as brightness, LD as hue

	setBatchMode("show");  // shows initial RGB image

	selectImage(rg_image_id);

	do{// keep showing the parameter setting dialog until the user is happy with the resulting RGB image
		param_modified = 0;
		setBatchMode(true);

		new_processing_params = get_manual_params(processing_params, rg_image_id);  //displays image processing dialog, collects values of image processing parameters

		//Array.print(new_processing_params);
		//print(color_scheme);


		//if pol_modulation || first_pixel_pol || bleedthrough || polarization_purity_hor || polarization_purity_ver are changed, get rid of the raw composite image and make a new one
		if((new_processing_params[0] != processing_params[0]) || (new_processing_params[1] != processing_params[1]) || (new_processing_params[2] != processing_params[2]) || (new_processing_params[5] != processing_params[5]) || (new_processing_params[6] != processing_params[6])){  
			selectImage(raw_composite_img_id);
			close();
			raw_composite_img_id = turn_RPM_image_into_a_composite(raw_slice_image_id, new_processing_params[2], new_processing_params[0], new_processing_params[1], new_processing_params[3]);
		}	

		// check if parameters have been changed by user, put newly entered parameters into parameter array
		for(param_number = 0; param_number < processing_params.length; param_number++){
			if(new_processing_params[param_number] != processing_params[param_number]){
				param_modified = 1;
				processing_params[param_number] = new_processing_params[param_number];
			}
		}


		nonguessable_params[0] = processing_params[0];
		nonguessable_params[1] = processing_params[1];
		nonguessable_params[2] = processing_params[2];
		nonguessable_params[3] = processing_params[3];
		nonguessable_params[4] = processing_params[4];
		nonguessable_params[5] = processing_params[5];
		nonguessable_params[6] = processing_params[6];


		if(param_modified == 1){
			selectImage(rg_image_id);
			close();
			selectImage(processed_composite_img_id);
			close();

			guessed_params = Array.slice(processing_params, 7);
			processed_composite_img_id = make_processed_composite(raw_composite_img_id, guessed_params);

			processing_params = Array.concat(nonguessable_params, guessed_params);  //pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, composite_img_id, hor_min, ver_min, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj
			color_scheme = processing_params[4];
			bright_adj = processing_params[9];
			rg_range = processing_params[10];
			rg_image_id = make_rg_image(processed_composite_img_id, bright_adj, rg_range, color_scheme);  
			
		}
		
		setBatchMode("show");

	
	}while (param_modified == 1)


	two_pol_color_scheme = color_scheme;
	save_params_to_config_file();

	selectImage(raw_composite_img_id);
	close();

	
	if(matches(determine_rmax_alpha0_sigma, "none") == false){
		selectImage(rg_image_id);
		run("Duplicate...", "title="+file_name_base+"_SEG.tif");
		segmentation_image_id = getImageID();
		segmentation_image_title = getTitle();
		setBatchMode("show");
		
		getPixelSize(pixel_size_unit, pixel_width, pixel_height);
		if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um" || pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm" || pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
			if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
			if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
			if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}
		}
		else{cutoff_distance_from_snake_in_pixels = cutoff_distance_from_snake;}
		run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 
			
		
		if(membranes_or_filaments == "Analyze a membrane"){
				
			// make a snake, manually reject some pixels
			make_membrane_snake(segmentation_image_id);
	
			check_with_user = true;
			while(check_with_user == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				check_with_user = false;
				run("Copy");
				setTool("multipoint");
				waitForUser("Modify selection or make a new selection, such as by using the polygon or segmented line tool, \n\nfollowed by a spline fit (pressing 'f'). For precise shape fitting of round vesicles, select \n3 points belonging to the vesicle outline. For fitting the shape of a cell, select 1 point\n approximately in the center of the cell. Erase unwanted pixels by the paintbrush tool.\n \nTo restore the selection made by the macro, press Ctrl+Shift+E.\n \nWhen done, click OK.");
				getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
				//print("Number of snake points/spline nodes:", snake_points_x.length);
		
				if(selectionType() == 6){ // (poly)line
					makeSelection("polygon", snake_points_x, snake_points_y);
				}
		
				if(snake_points_x.length == 1){  // a single point
					check_with_user = true;		
					center_x = snake_points_x[0];
					center_y = snake_points_y[0];
					radius = (h+w)/4-1;
					look_for_vesicles = false;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, false); //'false' = do not look for a circle by my algorithm
				}
				if(selectionType() == 10 && snake_points_x.length == 3){ // three points
					check_with_user = true;		
					run("Fit Circle");		
					Roi.getBounds(x, y, width, height)
					center_x = x + width/2;
					center_y = y + height/2;
					radius   = width/2;
					look_for_vesicles = true;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, true); //'true' = look for a circle by my algorithm
				}
				if(selectionType() == 1){  // circle/oval
					run("Interpolate");
					run("Fit Spline");
				}
				//print("selection type="+selectionType()+" snake_points_x.length="+snake_points_x.length);
			}
			Overlay.addSelection("orange");	
		}	
		
		if(membranes_or_filaments == "Analyze a filament"){
				
			// make a snake, manually reject some pixels
			//make_filament_snake(segmentation_image_id);
			run("Select None");
			run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 
			setTool("multipoint");

			check_with_user1 = true;		
			check_with_user2 = true;		
						
			while(check_with_user1 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Select a filament", "Mark a few points on (or close to) the filament by the 'multipoint' selection\n tool (to use the filament automatic tracing algorithm), or \nselect a filament using the 'segmented line' selection tool..\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user2 = true; 		// a single point or no selection
				}	
				else{
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){	// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user1 = true;
					} 
					else{
						if(selectionType() == 10){ trace_a_filament(segmentation_image_id); } 	// multipoint
						check_with_user1 = false;
					}
				}
			}
				
			while(check_with_user2 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Adjust/confirm selection", "If desired, modify selection or make a new selection. \nYou may also erase unwanted pixels by the paintbrush tool.\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user2 = true; // a single point or no selection
				}	
				else{	
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){		// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user2 = true;
					} 
					else{
						if(selectionType() == 10){     // multipoint
							trace_a_filament(segmentation_image_id); 
							check_with_user2 = true;
						} 
						else{
							check_with_user2 = false;
						}
					}	
				}				
			}
			makeSelection("polyline", snake_points_x, snake_points_y);
			
			Overlay.addSelection("orange");	
		}	
		


		setBatchMode(true);
		selectImage(segmentation_image_id);
		setBatchMode("hide");
	
		// segmentation_image is the background-subtracted raw image, with some pixels manually removed by user
		// Now that I have a snake (defined by snake_points_x, snake_points_y) from user selection, I should combine
		// the segmentation_image and the snake to make a _BIN image (pixels close to the snake, not manually rejected = 1, others = NaN)
	
	
		BIN_image_id = make_BIN_image(segmentation_image_id, rg_image_id, cutoff_distance_from_snake_in_pixels);	//the segmentation image has a selection representing the snake
		
	
		selectImage(BIN_image_id);
		rename(file_name_base+"_BIN.tif");
	
		selectImage(segmentation_image_id);
		close();
	

		
		if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){LD_quantitation_results = perform_LD_quantitation_1PPM(BIN_image_id, processed_composite_img_id);}  // the BIN image carries a selection (-> snake coordinates)}
		else{LD_quantitation_results = perform_LD_quantitation_2PPM(BIN_image_id, processed_composite_img_id);}

		plot_image_id =       parseInt(LD_quantitation_results[0]);
		hyperstack_image_id = parseInt(LD_quantitation_results[1]);
		theta_image_id =      parseInt(LD_quantitation_results[2]);
		fitting_data_text =  	       LD_quantitation_results[3];
		fitting_results_text =         LD_quantitation_results[4];
		selectImage(plot_image_id);
		setBatchMode("show");
		selectImage(hyperstack_image_id);
		setBatchMode("show");

		save_params_to_config_file();
	

		if(save_files == false){
			selectImage(theta_image_id);
			close();
			selectImage(BIN_image_id);
			close();
			//selectImage(processed_composite_img_id);
			//close();
		}
		else{
			selectImage(rg_image_id);
			saveAs("png");
		
			//print(File.directory);
			file_name_base = replace(File.name, "_range_.+","");
			
			selectImage(theta_image_id);
			save(File.directory+File.separator+file_name_base+"_THT.tif");
			close();
	
			selectImage(BIN_image_id);
			save(File.directory+File.separator+file_name_base+"_BIN.tif");
			close();

			selectImage(processed_composite_img_id);
			setSlice(1);
			run("Duplicate...", " ");
			save(File.directory+File.separator+file_name_base+"_HOR.tif");
			close();

			selectImage(processed_composite_img_id);
			setSlice(2);
			run("Duplicate...", " ");
			save(File.directory+File.separator+file_name_base+"_VER.tif");
			close();

			
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_1P_ALL.txt");}
			else{File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_2P_ALL.txt");}
	
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_1P_FIT.txt");}
			else{File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_2P_FIT.txt");}
	
			selectImage(plot_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_FIT");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_1P_FIT.tif");
				saveAs(File.directory+File.separator+file_name_base+"_1P_FIT.png");
			}
			else{
				rename(file_name_base+"_2P_FIT");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_2P_FIT.tif");
				saveAs(File.directory+File.separator+file_name_base+"_2P_FIT.png");
			}
	
			selectImage(hyperstack_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_1P_GFIT.tif");
			}
			else{
				rename(file_name_base+"_2P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_2P_GFIT.tif");
			}
		}		
	}

	selectImage(raw_slice_image_id);
	close();
	setTool("polygon");


	if(determine_rmax_alpha0_sigma == "none" && save_files == true){
		selectImage(rg_image_id);
		saveAs("png");			
	}	

	selectImage(processed_composite_img_id);
	close();



//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------ multislice image ------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------


	if(number_of_slices > 1 && process_whole_stack == true){


		selectImage(rg_image_id);
		close();

		selectImage(raw_image_id);
		w = getWidth();
		h = getHeight();

		final_stack_width = w;
		final_stack_height = h;
		if(pixel_shape == "rectangle" && pol_modulation = "columns"){final_stack_width  = w/2;}
		if(pixel_shape == "rectangle" && pol_modulation = "rows" . ){final_stack_height = h/2;}

		newImage("RGB_stack", "RGB Black", final_stack_width, final_stack_height, number_of_slices);
		rgb_stack_image_id = getImageID();		

		if(make_16_bit){
			newImage("hor_stack", "16-bit", final_stack_width , final_stack_height, number_of_slices);
			hor_stack_id = getImageID();
			newImage("ver_stack", "16-bit", final_stack_width , final_stack_height, number_of_slices);
			ver_stack_id = getImageID();
		}

		for (slice_number = 1; slice_number <= number_of_slices; slice_number++){
			selectImage(raw_image_id);
			setSlice(slice_number);
			run("Duplicate...", "title =[temp_slice.tif]");
			raw_slice_image_id = getImageID();
			raw_composite_img_id = turn_RPM_image_into_a_composite(raw_slice_image_id, bleedthrough, pol_modulation, first_pixel_pol, pixel_shape);
			selectImage(raw_composite_img_id);
			rename(slice_number);

			guessed_params = Array.slice(processing_params, 7);
			processed_composite_img_id = make_processed_composite(raw_composite_img_id, guessed_params);

			bright_adj = processing_params[9];
			rg_range = processing_params[10];
			rg_image_id = make_rg_image(processed_composite_img_id, bright_adj, rg_range, color_scheme);  


			if(make_16_bit){
				selectImage(processed_composite_img_id);
				setSlice(1);
				run("Select All");
				run("Copy");
				selectImage(hor_stack_id);
				setSlice(slice_number);
				run("Paste");

				selectImage(processed_composite_img_id);
				setSlice(2);
				run("Select All");
				run("Copy");
				selectImage(ver_stack_id);
				setSlice(slice_number);
				run("Paste");
			}

			selectImage(rg_image_id);
			run("Select All");
			run("Copy");

			selectImage(rgb_stack_image_id);
			setSlice(slice_number);
			run("Paste");

			selectImage(raw_composite_img_id);
			close();
			selectImage(processed_composite_img_id);
			close();
			selectImage(rg_image_id);
			close();
			selectImage(raw_slice_image_id);
			close();


		}


		if(make_16_bit){
			selectImage(hor_stack_id);
			setSlice(1);
			run("Enhance Contrast", "saturated=0.35");
			hor_stack_name = replace(raw_image_name, ".oib", ".tif");
			hor_stack_name = replace(hor_stack_name, ".tif", "_HOR.tif");
			rename(hor_stack_name);

			selectImage(ver_stack_id);
			setSlice(1);
			run("Enhance Contrast", "saturated=0.35");
			ver_stack_name = replace(raw_image_name, ".oib", ".tif");
			ver_stack_name = replace(ver_stack_name, ".tif", "_VER.tif");
			rename(ver_stack_name);
		}


		selectImage(rgb_stack_image_id);
		rg_image_name = replace(raw_image_name, ".tif", "_range_") + processing_params[11] + "_RG.tiff";
		rename(rg_image_name);


		setBatchMode("show");
		selectImage(rgb_stack_image_id);
	
	}


	if(orig_image_close){
		selectImage(raw_image_id);
		close();
	}

	setTool("polygon");
	setBatchMode("exit and display");

}








macro "Process a polarization image (composite image) [z]"{

	number_of_polarizations = 2;
	mixed_image = false;

	raw_image_name = getTitle();
	raw_image_id = getImageID();
	number_of_slices = nSlices;

	file_name_base = substring(raw_image_name,0,lengthOf(raw_image_name)-4);

	w = getWidth;
	h = getHeight;
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	composite_img_id = 0;

	read_params_from_config_file();   // reads parameters from a config file; stores them in global vars
	color_scheme = two_pol_color_scheme;
	nonguessable_params = newArray("N/A", first_pixel_pol, 0, "N/A", color_scheme, polarization_purity_hor, polarization_purity_ver);



	setBatchMode(true);

	selectImage(raw_image_id);
	run("Duplicate...", "duplicate");  // duplicates the whole stack

	if(endsWith(raw_image_name, ".tif") || endsWith(raw_image_name, ".oib")){
		raw_composite_name = replace(raw_image_name, ".tif", "_rcomp.tif");
		raw_composite_name = replace(raw_image_name, ".oib", "_rcomp.tif");
	}
	else{
		raw_composite_name = raw_image_name+"_rcomp.tif";		
	}
	rename(raw_composite_name);  
	raw_composite_img_id = getImageID();

	correct_negative_values(raw_composite_img_id);

	
	raw_slice_image_id = getImageID();


	guessed_params = guess_params(raw_composite_img_id); // returns hor_bkgd, ver_bkgd, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj


	processed_composite_img_id = make_processed_composite(raw_composite_img_id, guessed_params);
	//creates a composite image (slice 1 = Fh, slice 2 = Fv, by using the guessed processing parameters

	processing_params = Array.concat(nonguessable_params, guessed_params);  //pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, polarization_purity_hor, polarization_purity_ver, composite_img_id, hor_min, ver_min, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj


	bright_adj = processing_params[9];
	rg_range = processing_params[10];
	rg_image_id = make_rg_image(processed_composite_img_id, bright_adj, rg_range, color_scheme);  //creates a two-color (RGB) image showing fluorescence intensity as brightness, LD as hue

	setBatchMode("show");
   // shows initial RGB image

	selectImage(rg_image_id);

	do{
	// keep showing the parameter setting dialog until the user is happy with the resulting RGB image
		param_modified = 0;
		setBatchMode(true);

		new_processing_params = get_manual_params(processing_params, rg_image_id);  //displays image processing dialog, collects values of image processing parameters

		
		// check if parameters have been changed by user, put newly entered parameters into parameter array
		for(param_number = 0; param_number < processing_params.length; param_number++){
			if(new_processing_params[param_number] != processing_params[param_number]){
				param_modified = 1;
				processing_params[param_number] = new_processing_params[param_number];
			}
		}


		nonguessable_params[0] = processing_params[0];
		nonguessable_params[1] = processing_params[1];
		nonguessable_params[2] = processing_params[2];
		nonguessable_params[3] = processing_params[3];
		nonguessable_params[4] = processing_params[4];
		nonguessable_params[5] = processing_params[5];
		nonguessable_params[6] = processing_params[6];


		if(param_modified == 1){
			selectImage(rg_image_id);
			close();
			selectImage(processed_composite_img_id);
			close();

			guessed_params = Array.slice(processing_params, 7);
			processed_composite_img_id = make_processed_composite(raw_composite_img_id, guessed_params);

			processing_params = Array.concat(nonguessable_params, guessed_params);  //pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, composite_img_id, hor_min, ver_min, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj
			bright_adj = processing_params[9];
			rg_range = processing_params[10];
			rg_image_id = make_rg_image(processed_composite_img_id, bright_adj, rg_range, color_scheme);  
			
		}
		
		setBatchMode("show");

	
	}while (param_modified == 1)

	save_params_to_config_file();

	selectImage(raw_composite_img_id);
	close();

	
	if(determine_rmax_alpha0_sigma != "none"){
		selectImage(rg_image_id);
	
	run("Duplicate...", "title="+file_name_base+"_SEG.tif");
		segmentation_image_id = getImageID();
		segmentation_image_title = getTitle();
		setBatchMode("show");
		
		getPixelSize(pixel_size_unit, pixel_width, pixel_height);
		if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um" || pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm" || pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
			if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
			if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
			if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}
		}
		else{cutoff_distance_from_snake_in_pixels = cutoff_distance_from_snake;}
		run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 
			

		if(membranes_or_filaments == "Analyze a membrane"){
		
			// make a snake, manually reject some pixels
			make_membrane_snake(segmentation_image_id);
	
			check_with_user = true;
			while(check_with_user == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				check_with_user = false;
				run("Copy");
				setTool("multipoint");
				waitForUser("Modify selection or make a new selection, such as by using the polygon or segmented line tool, \n\nfollowed by a spline fit (pressing 'f'). For precise shape fitting of round vesicles, select \n3 points belonging to the vesicle outline. For fitting the shape of a cell, select 1 point\n approximately in the center of the cell. Erase unwanted pixels by the paintbrush tool.\n \nWhen done, click OK.");
				getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
		
				if(selectionType() == 6){ // (poly)line
					makeSelection("polygon", snake_points_x, snake_points_y);
				}
		
				if(snake_points_x.length == 1){  // a single point
					check_with_user = true;		
					center_x = snake_points_x[0];
					center_y = snake_points_y[0];
					radius = (h+w)/4-1;
					look_for_vesicles = false;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, false); //'false' = do not look for a circle by my algorithm
				}
				if(selectionType() == 10 && snake_points_x.length == 3){ // three points
					check_with_user = true;		
					run("Fit Circle");		
					Roi.getBounds(x, y, width, height)
					center_x = x + width/2;
					center_y = y + height/2;
					radius   = width/2;
					look_for_vesicles = true;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, true); //'true' = look for a circle by my algorithm
				}
				if(selectionType() == 1){  // circle/oval
					run("Interpolate");
					run("Fit Spline");
				}
			}
			Overlay.addSelection("orange");	
		}


		if(membranes_or_filaments == "Analyze a filament"){
				
			// make a snake, manually reject some pixels
			//make_filament_snake(segmentation_image_id);
			run("Select None");
			setTool("multipoint");
			run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

			check_with_user1 = true;		
			check_with_user2 = true;		
			
			while(check_with_user1 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Select a filament", "Mark a few points on (or close to) the filament by the 'multipoint' selection\n tool (to use the filament automatic tracing algorithm), or \nselect a filament using the 'segmented line' selection tool..\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user2 = true;} 		// a single point or no selection
				}	
				else{
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){	// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user1 = true;
					} 
					else{
						if(selectionType() == 10){ trace_a_filament(segmentation_image_id); } 	// multipoint
						check_with_user1 = false;
					}
				}
			}
				
			while(check_with_user2 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Adjust/confirm selection", "If desired, modify selection or make a new selection. \nYou may also erase unwanted pixels by the paintbrush tool.\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user2 = true; // a single point or no selection
				}	
				else{	
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){		// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user2 = true;
					} 
					else{
						if(selectionType() == 10){     // multipoint
							trace_a_filament(segmentation_image_id); 
							check_with_user2 = true;
						} 
						else{
							check_with_user2 = false;
						}
					}	
				}				
			}
			makeSelection("polyline", snake_points_x, snake_points_y);
		
			Overlay.addSelection("orange");	
		}	




		setBatchMode(true);
		selectImage(segmentation_image_id);
		setBatchMode("hide");
		selectImage(processed_composite_img_id);
		setBatchMode("hide");
	
		// segmentation_image is the background-subtracted raw image, with some pixels manually removed by user
		// Now that I have a snake (defined by snake_points_x, snake_points_y) from user selection, I should combine
		// the segmentation_image and the snake to make a _BIN image (pixels close to the snake, not manually rejected = 1, others = NaN)
	
	
		BIN_image_id = make_BIN_image(segmentation_image_id, rg_image_id, cutoff_distance_from_snake_in_pixels);	//the segmentation image has a selection representing the snake
		
	
		selectImage(BIN_image_id);
		rename(file_name_base+"_BIN.tif");
	
		selectImage(segmentation_image_id);
		close();
	

		
		if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){LD_quantitation_results = perform_LD_quantitation_1PPM(BIN_image_id, processed_composite_img_id);}  // the BIN image carries a selection (-> snake coordinates)}
		else{LD_quantitation_results = perform_LD_quantitation_2PPM(BIN_image_id, processed_composite_img_id);}

		plot_image_id =       parseInt(LD_quantitation_results[0]);
		hyperstack_image_id = parseInt(LD_quantitation_results[1]);
		theta_image_id =      parseInt(LD_quantitation_results[2]);
		fitting_data_text =  	       LD_quantitation_results[3];
		fitting_results_text =         LD_quantitation_results[4];
		selectImage(plot_image_id);
		setBatchMode("show");
		selectImage(hyperstack_image_id);
		setBatchMode("show");

		save_params_to_config_file();
	

		if(save_files == false){
			selectImage(theta_image_id);
			close();
			selectImage(BIN_image_id);
			close();
			selectImage(processed_composite_img_id);
			close();
		}
		else{
			selectImage(rg_image_id);
			saveAs("png");
		
			file_name_base = replace(File.name, "_range_.+","");
			
			selectImage(theta_image_id);
			save(File.directory+File.separator+file_name_base+"_THT.tif");
			close();
	
			selectImage(BIN_image_id);
			save(File.directory+File.separator+file_name_base+"_BIN.tif");
			close();

			selectImage(processed_composite_img_id);
			setSlice(1);
			run("Duplicate...", " ");
			save(File.directory+File.separator+file_name_base+"_HOR.tif");
			close();

			selectImage(processed_composite_img_id);
			setSlice(2);
			run("Duplicate...", " ");
			save(File.directory+File.separator+file_name_base+"_VER.tif");
			close();

			selectImage(processed_composite_img_id);
			close();
			
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_1P_ALL.txt");}
			else{File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_2P_ALL.txt");}
	
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_1P_FIT.txt");}
			else{File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_2P_FIT.txt");}
	
			selectImage(plot_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_FIT");
				saveAs(File.directory+File.separator+file_name_base+"_1P_FIT.png");
			}
			else{
				rename(file_name_base+"_1P_FIT");
				saveAs(File.directory+File.separator+file_name_base+"_2P_FIT.png");
			}
	
			selectImage(hyperstack_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_1P_GFIT.tif");
			}
			else{
				rename(file_name_base+"_2P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_2P_GFIT.tif");
			}
		}		
	}

	setTool("polygon");


	if(determine_rmax_alpha0_sigma == "none" && save_files == true){
		selectImage(rg_image_id);
		saveAs("png");			
	}	


	if(orig_image_close){
		selectImage(raw_image_id);
		close();
	}

	setTool("polygon");
	setBatchMode("exit and display");

}












macro "Process a polarization image (multiple polarizations) [m]"{

	raw_multipol_stack_image_id = getImageID();
	raw_image_name = getTitle();
	file_name_base = substring(raw_image_name,0,lengthOf(raw_image_name)-4);
	number_of_polarizations = nSlices;

	if(number_of_polarizations < 4){
		showMessage("A stack with at least 4 slices is required");
		exit;
	}
	
	w = getWidth;
	h = getHeight;
	hor_exists = true;
	ver_exists = true;
	use_fitting_for_LD = false;
	log2r_azimuth_r_squared_image_id = 0;
	
	//guessing initial values
	pol_angle_increment = 180/number_of_polarizations;

	read_params_from_config_file();   // reads parameters from a config file; stores them in global vars
	color_scheme = multipol_color_scheme;

	if(show_LD == "Accurate/slow"){show_LD = "Crude/fast";}

	setBatchMode(true);
	//print("\\Clear");


	guessed_multipol_params = guess_multipol_params(raw_multipol_stack_image_id); // returns bkgd, bright_adj, rg_range
	// Array.print(guessed_multipol_params);

	processing_params = newArray(guessed_multipol_params[0], guessed_multipol_params[1], show_LD, color_scheme, guessed_multipol_params[2], F_cutoff_for_LD_display, show_azimuths, azimuth_color, starting_pol_angle, pol_angle_increment, azimuth_pixel_bin_size, azimuth_scaling_factor, F_cutoff_for_azimuth_display, use_fitting_for_LD);
	// first_pol, pol_angle_increment, color_scheme, composite_img_id, bkgd, bright_adj, rg_range, show_azimuths, azimuth_color, azimuth_bin_size, azimuth_cutoff


	make_new_log2r_azimuth_image = true;
	processed_multipol_images = make_processed_multipol_img(raw_multipol_stack_image_id, processing_params, log2r_azimuth_r_squared_image_id, make_new_log2r_azimuth_image);
	processed_multipol_img_id = processed_multipol_images[0];
	log2r_azimuth_r_squared_image_id = processed_multipol_images[1];
	
	selectImage(processed_multipol_img_id);
	setBatchMode("show");


	// keep showing the parameter setting dialog until the user is happy with the resulting RGB image	
	do{
		param_modified = 0;
		setBatchMode(true);

		//print("processing_params:");
		//Array.print(processing_params);
		new_processing_params = get_manual_multipol_params(processing_params);  //displays image processing dialog, collects values of image processing parameters

		bkgd 							 = new_processing_params[0];
		bright_adj      				 = new_processing_params[1];
		show_LD  						 = new_processing_params[2];
		color_scheme 					 = new_processing_params[3];
		rg_range 						 = new_processing_params[4];
		F_cutoff_for_LD_display			 = new_processing_params[5];
		show_azimuths        			 = new_processing_params[6];
		azimuth_color					 = new_processing_params[7];
		starting_pol_angle				 = new_processing_params[8];
		pol_angle_increment				 = new_processing_params[9];
		azimuth_pixel_bin_size			 = new_processing_params[10];
		azimuth_scaling_factor           = new_processing_params[11];
		F_cutoff_for_azimuth_display	 = new_processing_params[12];		
		use_fitting_for_LD				 = new_processing_params[13];


		//show_LD = "None";
	
		//print("new_processing_params:");
		//Array.print(new_processing_params);

		if(processing_params[5] < processing_params[12]){
			lower_of_old_LD_azimuth_cutoffs = processing_params[5];
		}
		else{
			lower_of_old_LD_azimuth_cutoffs = processing_params[5];
		}
		
		if(new_processing_params[5] < new_processing_params[12]){
			lower_of_new_LD_azimuth_cutoffs = new_processing_params[5];
		}
		else{
			lower_of_new_LD_azimuth_cutoffs = new_processing_params[5];
		}

		make_new_log2r_azimuth_image = false;

		
		// check if parameters have been changed by user, put newly entered parameters into parameter array
		for(param_number = 0; param_number < processing_params.length; param_number++){
			if(toString(new_processing_params[param_number]) != toString(processing_params[param_number])){
				
				param_modified = 1;


				if(param_number == 0 || param_number == 1 || param_number == 2 || param_number == 8 || param_number == 9 || param_number == 10 || lower_of_new_LD_azimuth_cutoffs < lower_of_old_LD_azimuth_cutoffs){
					make_new_log2r_azimuth_image = true;
				}

				
				processing_params[param_number] = new_processing_params[param_number];
			}
		}

		//print("processing_params:");


		if(param_modified == 1){
			
			selectImage(processed_multipol_img_id);
			close();

			processed_multipol_images = make_processed_multipol_img(raw_multipol_stack_image_id, processing_params, log2r_azimuth_r_squared_image_id, make_new_log2r_azimuth_image); // returns an RGB image with an azimuth overlay
			processed_multipol_img_id = processed_multipol_images[0];
			log2r_azimuth_r_squared_image_id = processed_multipol_images[1];
			selectImage(processed_multipol_img_id);
			setBatchMode("show");
			
		}

	
	}while (param_modified == 1)

	if(isOpen(log2r_azimuth_r_squared_image_id)){
		selectImage(log2r_azimuth_r_squared_image_id);
		close();
	}



////////////////////////////////// LD quantitation //////////////////////////////
////////////////////////////////// LD quantitation //////////////////////////////
////////////////////////////////// LD quantitation //////////////////////////////
////////////////////////////////// LD quantitation //////////////////////////////



	if(determine_rmax_alpha0_sigma != "none"){    // = if LD quantitation: yes
		selectImage(raw_multipol_stack_image_id);
		run("Z Project...", "projection=[Average Intensity]");
		multipol_stack_z_average_image_id = getImageID();
		
		selectImage(processed_multipol_img_id);		 
	
		run("Duplicate...", "title="+file_name_base+"_SEG.tif");
		//run("Duplicate...", "title=_SEG.tif");
		segmentation_image_id = getImageID();
		segmentation_image_title = getTitle();
		//Overlay.remove;
		setBatchMode("show");
		
		getPixelSize(pixel_size_unit, pixel_width, pixel_height);
		//print(pixel_size_unit+", "+pixel_width+", "+pixel_height);
		if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um" || pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm" || pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
			if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
			if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
			if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}
		}
		else{cutoff_distance_from_snake_in_pixels = cutoff_distance_from_snake;}
		run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

		
		//print("snake width in pixels: "+cutoff_distance_from_snake_in_pixels);
		//print(pixel_size_unit+", "+pixel_width+", "+pixel_height);

		snake_points_x = newArray();
		snake_points_y = newArray();


		if(membranes_or_filaments == "Analyze a membrane"){
	
				
			// make a snake, manually reject some pixels
			// make_membrane_snake(segmentation_image_id);
			make_membrane_snake(segmentation_image_id);
			setBatchMode("show");
	
			check_with_user = true;
			while(check_with_user == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				check_with_user = false;
				run("Copy");
				setTool("multipoint");

					waitForUser("Modify selection or make a new selection, such as by using the polygon or segmented line tool, \n\nfollowed by a spline fit (pressing 'f'). For precise shape fitting of round vesicles, select \n3 points belonging to the vesicle outline. For fitting the shape of a cell, select 1 point\n approximately in the center of the cell. Erase unwanted pixels by the paintbrush tool.\n \nTo restore the selection made by the macro, press Ctrl+Shift+E.\n \nWhen done, click OK.");
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline

/*
				while(snake_points_x.length < 2){			
					waitForUser("Modify selection or make a new selection, such as by using the polygon or segmented line tool, \n\nfollowed by a spline fit (pressing 'f'). For precise shape fitting of round vesicles, select \n3 points belonging to the vesicle outline. For fitting the shape of a cell, select 1 point\n approximately in the center of the cell. Erase unwanted pixels by the paintbrush tool.\n \nTo restore the selection made by the macro, press Ctrl+Shift+E.\n \nWhen done, click OK.");
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
	
					if(snake_points_x.length < 2){
						showMessage("Please select at least two points");
					}
				}

*/				
				//print("Number of snake points/spline nodes:", snake_points_x.length);
	
				//print("Selection type: "+selectionType()+", num. of snake points: "+snake_points_x.length);
		
				if(selectionType() == 6){ // (poly)line
					// makeSelection("polygon", snake_points_x, snake_points_y);
				}
		
				if(snake_points_x.length == 1){  // a single point
					//print("a single point selection");
					check_with_user = true;		
					center_x = snake_points_x[0];
					center_y = snake_points_y[0];
					radius = (h+w)/4-1;
					look_for_vesicles = false;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, false); //'false' = do not look for a circle by my algorithm; makes a polyline fitted by a spline
				}
				if(selectionType() == 10 && snake_points_x.length == 3){ // three points
					check_with_user = true;		
					run("Fit Circle");		
					Roi.getBounds(x, y, width, height)
					center_x = x + width/2;
					center_y = y + height/2;
					radius   = width/2;
					look_for_vesicles = true;
					fit_outline_by_spline(segmentation_image_id, center_x, center_y, radius, number_of_polarizations, true); //'true' = look for a circle by my algorithm; makes a polyline fitted by a spline
				}
				if(selectionType() == 1){  // circle/oval
					run("Interpolate");
					run("Fit Spline");
				}
				//print("selection type="+selectionType()+" snake_points_x.length="+snake_points_x.length);
			}
	
			if(selectionType() == 1 || selectionType() == 6){  // circle/oval
				makeSelection("polyline", snake_points_x, snake_points_y);
			}
	
			// HSV color definition
			make_HSV_snake_overlay(segmentation_image_id, snake_points_x, snake_points_y);
			make_HSV_snake_overlay(processed_multipol_img_id, snake_points_x, snake_points_y);
	
	
			setBatchMode(true);
		}

		if(membranes_or_filaments == "Analyze a filament"){
				
			// make a snake, manually reject some pixels
			// make_filament_snake(segmentation_image_id);
			
			run("Select None");
			setTool("multipoint");
			run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

			check_with_user1 = true;		
			check_with_user2 = true;		
			
			while(check_with_user1 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Select a filament", "Mark a few points on (or close to) the filament by the 'multipoint' selection\n tool (to use the filament automatic tracing algorithm), or \nselect a filament using the 'segmented line' selection tool..\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user1 = true; 		// a single point or no selection
				}	
				else{
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){	// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user1 = true;
					} 
					else{
						if(selectionType() == 10){  // multipoint
							selectImage(multipol_stack_z_average_image_id);
							makeSelection("polyline", snake_points_x, snake_points_y);
							trace_a_filament(multipol_stack_z_average_image_id); 
							run("Copy");
							selectImage(segmentation_image_id);
							run("Restore Selection");
							check_with_user1 = false;
						}
					}
				}
			}
				
			while(check_with_user2 == true){		// keep making a new snake if the user selects one point (= center of a cell) or three points (= defining a circular vesicle outline)
				waitForUser("Adjust/confirm selection", "If desired, modify selection or make a new selection. \nYou may also erase unwanted pixels by the paintbrush tool.\n \nWhen done, click OK.");
				if(selectionType() == -1){
					showMessage("At least 2 points are required. Please adjust your selection.");
					check_with_user2 = true; // a single point or no selection
				}	
				else{	
					getSelectionCoordinates(snake_points_x, snake_points_y);	// x, y coordinates of the spline
					if(snake_points_x.length < 2){		// a single point or no selection
						showMessage("At least 2 points are required. Please adjust your selection.");
						check_with_user2 = true;
					} 
					else{
						if(selectionType() == 10){     // multipoint
							selectImage(multipol_stack_z_average_image_id);
							makeSelection("polyline", snake_points_x, snake_points_y);
							trace_a_filament(multipol_stack_z_average_image_id); 
							run("Copy");
							selectImage(segmentation_image_id);
							run("Restore Selection");
							check_with_user2 = true;
						} 
						else{
							check_with_user2 = false;
						}
					}	
				}				
			}
			makeSelection("polyline", snake_points_x, snake_points_y);
			
			//Overlay.addSelection("orange");	
			make_HSV_snake_overlay(segmentation_image_id, snake_points_x, snake_points_y);
			make_HSV_snake_overlay(processed_multipol_img_id, snake_points_x, snake_points_y);

		}	
	

	
		// segmentation_image is the background-subtracted RG image, with some pixels manually removed by user
		// Now that I have a snake (defined by snake_points_x, snake_points_y) from user selection, I should combine
		// the segmentation_image and the snake to make a _BIN image (pixels close to the snake, not manually rejected = 1, others = NaN)
	
	
		BIN_image_id = make_BIN_image(segmentation_image_id, processed_multipol_img_id, cutoff_distance_from_snake_in_pixels);	//the segmentation image has a selection representing the snake
		


	
		selectImage(BIN_image_id);
		rename(file_name_base+"_BIN.tif");
		//rename("_BIN.tif");

		selectImage(multipol_stack_z_average_image_id);
		close();
	
		selectImage(segmentation_image_id);
		close();
	

		
		if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){ // the BIN image carries a selection (-> snake coordinates)}
			selectImage(raw_multipol_stack_image_id);
			run("Duplicate...", "title=bkgd_subbed_multipol_stack duplicate");
			run("Subtract...", "value="+bkgd+" stack");
			bkgd_subbed_multipol_stack_image_id = getImageID();
			
			//LD_quantitation_results = perform_LD_quantitation_multipol_1PPM(BIN_image_id, raw_multipol_stack_image_id);
			LD_quantitation_results = perform_LD_quantitation_1PPM(BIN_image_id, bkgd_subbed_multipol_stack_image_id);  // the BIN image contains values of 1 and NaN, and carries a selection (= snake points)
			if(show_additional_plots == true){
				additional_plot_ids = show_LD_azimuths_along_the_snake(BIN_image_id, bkgd_subbed_multipol_stack_image_id);
				log2rmax_r_squared_plot_image_id = additional_plot_ids[0];
				azimuths_thetas_plot_image_id = additional_plot_ids[1];
			}
		}  
		else{  // 2PPM
			selectImage(raw_multipol_stack_image_id);
			run("Duplicate...", "title=bkgd_subbed_multipol_stack duplicate");
			run("Subtract...", "value="+bkgd+" stack");
			bkgd_subbed_multipol_stack_image_id = getImageID();

			LD_quantitation_results = perform_LD_quantitation_2PPM(BIN_image_id, bkgd_subbed_multipol_stack_image_id);
			if(show_additional_plots == true){
				additional_plot_ids = show_LD_azimuths_along_the_snake(BIN_image_id, bkgd_subbed_multipol_stack_image_id);
				log2rmax_r_squared_plot_image_id = additional_plot_ids[0];
				azimuths_thetas_plot_image_id = additional_plot_ids[1];
			}
		}

		plot_image_id =       parseInt(LD_quantitation_results[0]);
		hyperstack_image_id = parseInt(LD_quantitation_results[1]);
		theta_image_id =      parseInt(LD_quantitation_results[2]);
		fitting_data_text =  	       LD_quantitation_results[3];
		fitting_results_text =         LD_quantitation_results[4];

		selectImage(plot_image_id);
		rename("Plot image");
		setBatchMode("show");

		selectImage(hyperstack_image_id);
		rename("Hyperstack image");
		setBatchMode("show");

		selectImage(theta_image_id);
		rename("Theta image");
		//setBatchMode("show");

		print("Output of perform_LD_quantitation_1PPM");
		print("--------------------------------------");
		print("Fitting data text:");
		print(fitting_data_text);
		print("");
		print("Fitting results text:");
		print(fitting_results_text);
		
		//setBatchMode("exit and display");

		
		save_params_to_config_file();
	

		if(save_files == false){
			selectImage(theta_image_id);
			close();
			selectImage(BIN_image_id);
			close();
			selectImage(bkgd_subbed_multipol_stack_image_id);
			close();
			//selectImage(processed_composite_img_id);
			//close();
		}
		else{
			selectImage(processed_multipol_img_id);
			saveAs("tiff");
			png_file_name = replace(File.name, ".tif",".png");
			save(File.directory+File.separator+png_file_name);
		
			//print(File.directory);
			file_name_base = replace(File.name, "_range_.+","");
			file_name_base = replace(file_name_base, ".png","");

			
			selectImage(theta_image_id);
			save(File.directory+File.separator+file_name_base+"_THT.tif");
			close();
	
			selectImage(BIN_image_id);
			save(File.directory+File.separator+file_name_base+"_BIN.tif");
			close();

			selectImage(bkgd_subbed_multipol_stack_image_id);
			save(File.directory+File.separator+file_name_base+"_B-SUBBED.tif");
			close();

			
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_1P_ALL.txt");}
			else{File.saveString(fitting_data_text, File.directory+File.separator+file_name_base+"_2P_ALL.txt");}
	
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_1P_FIT.txt");}
			else{File.saveString(fitting_results_text, File.directory+File.separator+file_name_base+"_2P_FIT.txt");}
	
			selectImage(plot_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_FIT");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_1P_FIT.tif");
				saveAs(File.directory+File.separator+file_name_base+"_1P_FIT.png");
			}
			else{
				rename(file_name_base+"_2P_FIT");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_2P_FIT.tif");
				saveAs(File.directory+File.separator+file_name_base+"_2P_FIT.png");
			}
	
			selectImage(hyperstack_image_id);
			if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){
				rename(file_name_base+"_1P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_1P_GFIT.tif");
			}
			else{
				rename(file_name_base+"_2P_GFIT");
				saveAs(File.directory+File.separator+file_name_base+"_2P_GFIT.tif");
			}

			if(show_additional_plots == true){
				selectImage(log2rmax_r_squared_plot_image_id);
				rename(file_name_base+"_L2R_trace");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_L2R_trace.tif");
				saveAs(File.directory+File.separator+file_name_base+"_L2R_trace.png");
				
				selectImage(azimuths_thetas_plot_image_id);
				rename(file_name_base+"_azimuth_trace");
				saveAs("Tiff", File.directory+File.separator+file_name_base+"_azimuth_trace.tif");
				saveAs(File.directory+File.separator+file_name_base+"_azimuth_trace.png");
			}		
		}		
	}



	if(determine_rmax_alpha0_sigma == "none" && save_files == true){
		selectImage(processed_multipol_img_id);
		saveAs("tiff");
		png_file_name = replace(File.name, ".tif",".png");
		save(File.directory+File.separator+png_file_name);
			
	}	


	
	
	
	
////////////////////////////////// end of LD quantitation //////////////////////////////
////////////////////////////////// end of LD quantitation //////////////////////////////
////////////////////////////////// end of LD quantitation //////////////////////////////
////////////////////////////////// end of LD quantitation //////////////////////////////
	
	
	setBatchMode("exit and display");
	
	save_params_to_config_file();


	if(orig_image_close){
		selectImage(raw_multipol_stack_image_id);
		close();
	}


	//print("Done!");
	

}








macro "Combine 1P, 2P data [c]"{

	setBatchMode(true);

	operating_system = getInfo("os.name");
	print("\\Clear");

	image_of_1P_data_image_id = 0;
	image_of_1P_counts_image_id = 0;
	image_of_2P_data_image_id = 0;
	image_of_2P_counts_image_id = 0;
	number_of_theta_bins = 181;
	number_of_1P_FIT_files = 0;
	number_of_2P_FIT_files = 0;

	binned_1P_log2r_values = newArray(181);
	binned_1P_log2r_counts = newArray(181);
	binned_2P_log2r_values = newArray(181);
	binned_2P_log2r_counts = newArray(181);
	combined_fitting_1P_2P_parsing_summary = "";
	combined_fitting_1P_2P_SUMMARY_text = "";

	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();
	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();

	run("Image Expression Parser (Macro)", "expression=[(4*B - 7*A)/(-10 + A + 2*B)] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("1PPM_B.tif");
	B_1PPM_image_id = getImageID();


	num_of_alpha0s = getWidth;
	num_of_sigmas  = getHeight;


//######################compiling/plotting 1P data##########################################################################

	list_of_1P_ALL_files = newArray();
	list_of_1P_thetas = newArray();
	list_of_1P_log_r_values = newArray();
	list_of_1P_B1P_values = newArray();
	list_of_1P_log2rmax_values = newArray();
	list_of_1P_thetas_from_one_ALL_file = newArray();
	list_of_1P_log_r_values_from_one_ALL_file = newArray();
	list_of_1P_fitting_results = "";
	path_1P = File.openDialog("Select a file in a directory containing 1P data");
	directory_1P = File.getParent(path_1P);
	combined_fitting_1P_2P_SUMMARY_text = combined_fitting_1P_2P_SUMMARY_text + "1P data directory: "+directory_1P+"\r\n";

	//directory_1P = getDirectory("Select a directory containing 1PPM data");
	file_list_1P = getFileList(directory_1P);

	color_list = newArray("gray", "lightGray", "blue", "cyan",  "green", "yellow", "orange", "red", "pink", "magenta", "black");
	color_count = 0;
	plot_legend = "";

	Plot.create("Graph of 1P/2P data" , "Membrane orientation ("+thetachar+")", "");
	Plot.setLimits(-90,270,NaN,NaN);
	Plot.setFrameSize(450, 200);
	Plot.setLineWidth(2);
	Plot.setColor(color_list[color_count]);

	// check whether 1P_FIT.txt files were made for 1 polarization or for 2 polarizations
	number_of_1_polarization_FIT_files = 0;
	number_of_2_polarization_FIT_files = 0;
	number_of_polarizations_1P = 0;
	
	for (i=0; i<file_list_1P.length; i++) {
		if (endsWith(file_list_1P[i], "_1P_FIT.txt")){
			contents_of_1P_FIT_file = File.openAsString(directory_1P +File.separator+ file_list_1P[i]);
			contents_of_1P_FIT_file_rows = split(contents_of_1P_FIT_file,"\n");
	
			for(j=0; j<contents_of_1P_FIT_file_rows.length; j++){
				key_value_pair = split(contents_of_1P_FIT_file_rows[j]);
				if(key_value_pair.length > 1){
					if(key_value_pair[0] == "number_of_polarizations"){
						number_of_polarizations_1P = key_value_pair[1];
					}				
				}
			}
		

			if(number_of_polarizations_1P == 1){
				number_of_1_polarization_FIT_files = number_of_1_polarization_FIT_files + 1;
			}
			else{
				number_of_2_polarization_FIT_files = number_of_2_polarization_FIT_files + 1;
			}
		}
	}

	if(number_of_1_polarization_FIT_files > 0 && number_of_2_polarization_FIT_files > 0){
		if(number_of_1_polarization_FIT_files > number_of_1_polarization_FIT_files){
			number_of_polarizations_1P = 1;
		}
		else{
			number_of_polarizations_1P = 2;
		}
		showMessage("1PPM images acquired with one ("+number_of_1_polarization_FIT_files+") and two ("+number_of_2_polarization_FIT_files+") polarizations \nfound. Only files acquired with "+number_of_polarizations_1P+" polarization will be used for further analysis.");			
	}
	else{
		if(number_of_1_polarization_FIT_files > 0){
			number_of_polarizations_1P = 1;			
		}
		if(number_of_2_polarization_FIT_files > 0){
			number_of_polarizations_1P = 2;			
		}	
	}




	for (i=0; i<file_list_1P.length; i++) {

		list_of_1P_thetas_from_one_ALL_file = newArray();
		list_of_1P_thetas_from_one_ALL_file_deg = newArray();
		list_of_1P_Fhs_from_one_ALL_file = newArray();
		list_of_1P_Fvs_from_one_ALL_file = newArray();
		list_of_1P_log_r_values_from_one_ALL_file = newArray();
		offset = 0;
		normalization_factor = 0;
		phase = 0;
		phase2 = 0;
		B1P = 0;
		polarization_direction = "";
		file_parsed = false;

		if (endsWith(file_list_1P[i], "_1P_FIT.txt")){
			contents_of_1P_FIT_file = File.openAsString(directory_1P +File.separator+ file_list_1P[i]);

			print(file_list_1P[i]+"num of pols "+number_of_polarizations_1P+" index of n_of_p = 1 "+indexOf(contents_of_1P_FIT_file, "number_of_polarizations = 1"));
			combined_fitting_1P_2P_parsing_summary = combined_fitting_1P_2P_parsing_summary + "################################\r\n"+file_list_1P[i]+"\r\nnumber of polarizations "+number_of_polarizations_1P+"\r\nindex of number of polarizations = "+indexOf(contents_of_1P_FIT_file, "\r\nnumber_of_polarizations = ")+"\r\n\r\n";

			// one polarization
			if(  number_of_polarizations_1P == 1 && indexOf(contents_of_1P_FIT_file, "number_of_polarizations = 1") > 0  ) {
				lines = split(contents_of_1P_FIT_file,"\n");
				for (j=0; j<lines.length; j++){
					//print(lines[j]);
					if(startsWith(lines[j],"normalization_factor =")){
						normalization_factor_line = split(lines[j]," ");
						normalization_factor = parseFloat(normalization_factor_line[2]);
					}
					if(startsWith(lines[j],"offset =")){
						offset_line = split(lines[j]," ");
						offset = parseFloat(offset_line[2]);
					}
					if(startsWith(lines[j],"phase =")){
						phase_line = split(lines[j]," ");
						phase = parseFloat(phase_line[2]);
					}
					if(startsWith(lines[j],"a =")){
						a_line = split(lines[j]," ");
						a = parseFloat(a_line[2]);
					}
					if(startsWith(lines[j],"B1P =")){
						B1P_line = split(lines[j]," ");
						B1P = parseFloat(B1P_line[2]);
						list_of_1P_B1P_values = Array.concat(list_of_1P_B1P_values, B1P);
						log2rmax = log( (1 + B1P)/(1 - B1P) )/log(2);
						list_of_1P_log2rmax_values = Array.concat(list_of_1P_log2rmax_values, log2rmax);
					}
					if(startsWith(lines[j],"polarization_direction =")){
						polarization_direction_line = split(lines[j]," ");
						polarization_direction = polarization_direction_line[2];
					}
					if(startsWith(lines[j],"r(max) ")){
						list_of_1P_fitting_results = list_of_1P_fitting_results + "\n" + lines[j+1];
					}
				}
				FIT_file_parsed = true;	

				name_of_1P_ALL_file = replace(file_list_1P[i],"_1P_FIT.txt","_1P_ALL.txt");
				contents_of_1P_ALL_file = File.openAsString(directory_1P +File.separator+ name_of_1P_ALL_file);

				lines = split(contents_of_1P_ALL_file,"\n");
				for(j=1; j<lines.length; j++){    //skip the 1st line (headers)
					values = split(lines[j],"\t");
					theta = parseFloat(values[0]) - phase;
					if(polarization_direction == "vertical"){
						theta = theta - PI/2;
					}
					while(theta < -PI/2){
						theta = theta + PI;
					}
					while(theta > PI/2){
						theta = theta - PI;
					}
					r = parseFloat(values[5]);
					if(r > 0 && indexOf(r, "Infinity") < 0){
						log2r = r / normalization_factor;
						list_of_1P_thetas_from_one_ALL_file = Array.concat(list_of_1P_thetas_from_one_ALL_file, theta);
						list_of_1P_thetas_from_one_ALL_file_deg = Array.concat(list_of_1P_thetas_from_one_ALL_file_deg, theta/3.14159*180);
						list_of_1P_log_r_values_from_one_ALL_file = Array.concat(list_of_1P_log_r_values_from_one_ALL_file, log2r);
					}
				}
	
				list_of_1P_thetas = Array.concat(list_of_1P_thetas, list_of_1P_thetas_from_one_ALL_file);
				list_of_1P_log_r_values = Array.concat(list_of_1P_log_r_values, list_of_1P_log_r_values_from_one_ALL_file);						 
				 
			}
			// two polarizations -------------------------------------
			else{ 
			    if( number_of_polarizations_1P > 1 && indexOf(contents_of_1P_FIT_file, "number_of_polarizations = ") > 0)  {
					lines = split(contents_of_1P_FIT_file,"\n");
					for (j=0; j<lines.length; j++){
						if(startsWith(lines[j],"normalization_factor =")){
							normalization_factor_line = split(lines[j]," ");
							normalization_factor = parseFloat(normalization_factor_line[2]);
						}
						if(startsWith(lines[j],"offset =")){
							offset_line = split(lines[j]," ");
							offset = parseFloat(offset_line[2]);
						}
						if(startsWith(lines[j],"a =")){
							a_line = split(lines[j]," ");
							a = parseFloat(a_line[2]);
						}
						if(startsWith(lines[j],"phase =")){
							phase_line = split(lines[j]," ");
							phase = parseFloat(phase_line[2]); 
						}
						if(startsWith(lines[j],"phase2 =")){
							phase2_line = split(lines[j]," ");
							phase2 = parseFloat(phase2_line[2]); 
						}
						if(startsWith(lines[j],"B1P =")){
							B1P_line = split(lines[j]," ");
							B1P = parseFloat(B1P_line[2]);
//B1P = - B1P;							
							list_of_1P_B1P_values = Array.concat(list_of_1P_B1P_values, B1P);
							log2rmax = log( (1 + B1P)/(1 - B1P) )/log(2);
							list_of_1P_log2rmax_values = Array.concat(list_of_1P_log2rmax_values, log2rmax);
						}
						if(startsWith(lines[j],"polarization_direction =")){
							polarization_direction_line = split(lines[j]," ");
							polarization_direction = parseFloat(polarization_direction_line[2]);
						}
						if(startsWith(lines[j],"r(max) ")){
							list_of_1P_fitting_results = list_of_1P_fitting_results + "\n" + lines[j+1];
						}
					}
					FIT_file_parsed = true;				 
				}
				if( number_of_polarizations_1P > 1 && indexOf(contents_of_1P_FIT_file, "number_of_polarizations = ") == -1  ){
					lines = split(contents_of_1P_FIT_file,"\n");
					B1P_line = split(lines[2]," ");
					offset_line = split(lines[4]," ");
					phase_line  = split(lines[3]," ");
					B1P = parseFloat(B1P_line[2]);

//B1P = - B1P;
					
					offset = parseFloat(offset_line[2]);
					phase  = parseFloat(phase_line[2]);

					list_of_1P_B1P_values = Array.concat(list_of_1P_B1P_values, B1P);
					log2rmax = log( (1 + B1P)/(1 - B1P) )/log(2);
					list_of_1P_log2rmax_values = Array.concat(list_of_1P_log2rmax_values, log2rmax);

					FIT_file_parsed = true;
				}
				if(FIT_file_parsed == true){
					name_of_1P_ALL_file = replace(file_list_1P[i],"_1P_FIT.txt","_1P_ALL.txt");
					contents_of_1P_ALL_file = File.openAsString(directory_1P +File.separator+ name_of_1P_ALL_file);

					lines = split(contents_of_1P_ALL_file,"\n");
					for(j=1; j<lines.length; j++){    //skip the 1st line (headers)
						values = split(lines[j],"\t");
						theta = parseFloat(values[0]);

//theta = theta + PI/2;

						while(theta < -PI/2){theta = theta + PI;}
						while(theta >  PI/2){theta = theta - PI;}

						Fh = parseFloat(values[1]);
						Fv = parseFloat(values[2]);
						r  = parseFloat(values[5]);
						num_of_pixels = parseFloat(values[6]);


						if(r > 0 && indexOf(r, "Infinity") < 0 && num_of_pixels > 0 && indexOf(r, "NaN") < 0){
							log2r = log(r)/log(2)-offset;
							//log2r = log(r)/log(2)-0.28;
							list_of_1P_thetas_from_one_ALL_file = Array.concat(list_of_1P_thetas_from_one_ALL_file, theta);
							list_of_1P_thetas_from_one_ALL_file_deg = Array.concat(list_of_1P_thetas_from_one_ALL_file_deg, theta/3.14159*180);
							list_of_1P_Fhs_from_one_ALL_file = Array.concat(list_of_1P_Fhs_from_one_ALL_file,Fh/num_of_pixels);
							list_of_1P_Fvs_from_one_ALL_file = Array.concat(list_of_1P_Fvs_from_one_ALL_file,Fv/num_of_pixels);
							list_of_1P_log_r_values_from_one_ALL_file = Array.concat(list_of_1P_log_r_values_from_one_ALL_file, log2r);
						}
					}			
				}
				
				if(phase2 == 0){
					phase2 = phase;
				}	
				for(k=0; k<list_of_1P_log_r_values_from_one_ALL_file.length; k++){
					Fh = list_of_1P_Fhs_from_one_ALL_file[k];
					Fv = list_of_1P_Fvs_from_one_ALL_file[k];
					theta = list_of_1P_thetas_from_one_ALL_file[k];
					r = Fh/Fv;

					log2r_experiment = log(Fh/Fv)/log(2);
					log2r_skewed_fit = offset + log(  (1 + B1P*cos(2*theta-2*phase)) / (1 - B1P*cos(2*theta-2*phase2))  )/log(2);
					delta_log2r = log2r_experiment - log2r_skewed_fit;

					log2r_straight_fit = log(  (1 + B1P*cos(2*theta)) / (1 - B1P*cos(2*theta))  )/log(2);
					log2r_straightened_experiment = log2r_straight_fit + delta_log2r;

					list_of_1P_log_r_values_from_one_ALL_file[k] = log2r_straightened_experiment;

				}
			
				list_of_1P_log_r_values = Array.concat(list_of_1P_log_r_values, list_of_1P_log_r_values_from_one_ALL_file);
				list_of_1P_thetas = Array.concat(list_of_1P_thetas, list_of_1P_thetas_from_one_ALL_file);

				
			}

			print("parsed "+FIT_file_parsed+" offset "+offset+" norm f "+normalization_factor+" B1P "+B1P+" phase "+phase+" phase2 "+phase2+" pol direction "+polarization_direction);

			

			//make/expand an image of 1P data, 1P data counts
			number_of_1P_FIT_files = number_of_1P_FIT_files + 1;
			
			if(number_of_1P_FIT_files == 1){
				newImage("1P_data", "32-bit black", number_of_theta_bins, number_of_1P_FIT_files, 1);
				image_of_1P_data_image_id = getImageID();
				newImage("1P_counts", "16-bit black", number_of_theta_bins, number_of_1P_FIT_files, 1);
				image_of_1P_counts_image_id = getImageID();				
			}
			else{
				selectImage(image_of_1P_data_image_id);
				run("Canvas Size...", "width="+number_of_theta_bins+" height="+number_of_1P_FIT_files+" position=Top-Center zero");	
				selectImage(image_of_1P_counts_image_id);
				run("Canvas Size...", "width="+number_of_theta_bins+" height="+number_of_1P_FIT_files+" position=Top-Center zero");	
			}
			for(j = 0; j < list_of_1P_thetas_from_one_ALL_file.length; j++){
				theta = list_of_1P_thetas_from_one_ALL_file[j];
				theta = theta + PI/2;	// theta is now between 0 and PI
				theta = round(theta/PI*180);   // theta is now in degrees
				selectImage(image_of_1P_data_image_id);
				setPixel(theta, number_of_1P_FIT_files-1, getPixel(theta, number_of_1P_FIT_files-1)+list_of_1P_log_r_values_from_one_ALL_file[j]);
				selectImage(image_of_1P_counts_image_id);
				setPixel(theta, number_of_1P_FIT_files-1, getPixel(theta, number_of_1P_FIT_files-1)+1);
			}




			plot_legend = plot_legend + name_of_1P_ALL_file + "\n";
			Plot.setColor(color_list[color_count]);

Plot.setColor("gray");
			
			
			Plot.setLineWidth(1);
			Plot.add("dot", list_of_1P_thetas_from_one_ALL_file_deg, list_of_1P_log_r_values_from_one_ALL_file);

			print("thetas line 1834");
			Array.print(list_of_1P_thetas_from_one_ALL_file_deg);
			print("log2rs");
			Array.print(list_of_1P_log_r_values_from_one_ALL_file);
			
			color_count = color_count + 1;
			if(color_count > 10){
				color_count = color_count - 11;
			}		
		}
	}


	if(image_of_1P_data_image_id < 0 && image_of_1P_counts_image_id < 0){
		imageCalculator("Divide", image_of_1P_data_image_id, image_of_1P_counts_image_id);
		selectImage(image_of_1P_counts_image_id);
		close();	
	}

	Array.getStatistics(list_of_1P_B1P_values, list_of_1P_B1P_values_min, list_of_1P_B1P_values_max, list_of_1P_B1P_values_mean, list_of_1P_B1P_values_stdDev);
	list_of_1P_B1P_values_2x_stdErr = 2* (list_of_1P_B1P_values_stdDev/sqrt(list_of_1P_B1P_values.length));
	Array.getStatistics(list_of_1P_log2rmax_values, list_of_1P_log2rmax_values_min, list_of_1P_log2rmax_values_max, list_of_1P_log2rmax_values_mean, list_of_1P_log2rmax_values_stdDev);
	list_of_1P_log2rmax_values_2x_stdErr = 2* (list_of_1P_log2rmax_values_stdDev/sqrt(list_of_1P_log2rmax_values.length));


	Array.getStatistics(list_of_1P_log_r_values, list_of_1P_log_r_values_min, list_of_1P_log_r_values_max, list_of_1P_log_r_values_mean, list_of_1P_log_r_values_stdev); 

/*
print("list of theta values");
Array.print(list_of_1P_thetas);
print("list_of_1P_log_r_values");
Array.print(list_of_1P_log_r_values);
print(list_of_1P_log_r_values_min, list_of_1P_log_r_values_max, list_of_1P_log_r_values_mean, list_of_1P_log_r_values_stdev);
brk
*/

//######################compiling/plotting 2P data##########################################################################

	list_of_2P_ALL_files = newArray();
	list_of_2P_thetas = newArray();
	list_of_2P_log_r_values = newArray();
	list_of_2P_alpha0s = newArray();
	list_of_2P_sigmas = newArray();
	list_of_2P_thetas_from_one_ALL_file = newArray();
	list_of_2P_log_r_values_from_one_ALL_file = newArray();
	list_of_2P_B2P_values = newArray();
	list_of_2P_C2P_values = newArray();
	list_of_2P_log2rmax_values = newArray();
	list_of_2P_fitting_results = "";
	
	path_2P = File.openDialog("Select a file in a directory containing 2P data");
	directory_2P = File.getParent(path_2P);
	file_list_2P = getFileList(directory_2P);

	color_count = 0;

	Plot.setLineWidth(1);
	Plot.setColor(color_list[color_count]);
	
	// check whether 2P_FIT.txt files were made for 1 polarization or for 2 polarizations
	number_of_1_polarization_2P_FIT_files = 0;
	number_of_2_polarization_2P_FIT_files = 0;
	number_of_polarizations_2P = 0;
	
	for (i=0; i<file_list_2P.length; i++) {
		if (endsWith(file_list_2P[i], "_2P_FIT.txt")){

			contents_of_2P_FIT_file = File.openAsString(directory_2P +File.separator+ file_list_2P[i]);
			contents_of_2P_FIT_file_rows = split(contents_of_2P_FIT_file,"\n");
	
			for(j=0; j<contents_of_2P_FIT_file_rows.length; j++){
				key_value_pair = split(contents_of_2P_FIT_file_rows[j]);
				if(key_value_pair.length > 1){
					if(key_value_pair[0] == "number_of_polarizations"){
						number_of_polarizations_2P = key_value_pair[1];
					}				
				}
			}
			

			if(number_of_polarizations_2P == 1){
				number_of_1_polarization_2P_FIT_files = number_of_1_polarization_2P_FIT_files + 1;
			}
			else{
				number_of_2_polarization_2P_FIT_files = number_of_2_polarization_2P_FIT_files + 1;
			}
		}
	}

	if(number_of_1_polarization_2P_FIT_files > 0 && number_of_2_polarization_2P_FIT_files > 0){
		if(number_of_1_polarization_2P_FIT_files > number_of_1_polarization_2P_FIT_files){
			number_of_polarizations_2P = 1;
		}
		else{
			number_of_polarizations_2P = 2;
		}
		showMessage("2PPM images acquired with one ("+number_of_1_polarization_2P_FIT_files+") and two ("+number_of_2_polarization_2P_FIT_files+") polarizations \nfound. Only files acquired with "+number_of_polarizations_2P+" polarization will be used for further analysis.");			
	}
	else{
		if(number_of_1_polarization_2P_FIT_files > 0){
			number_of_polarizations_2P = 1;			
		}
		if(number_of_2_polarization_2P_FIT_files > 0){
			number_of_polarizations_2P = 2;			
		}	
	}


	for (i=0; i<file_list_2P.length; i++) {

		list_of_2P_thetas_from_one_ALL_file = newArray();
		list_of_2P_thetas_from_one_ALL_file_deg = newArray();
		list_of_2P_Fhs_from_one_ALL_file = newArray();
		list_of_2P_Fvs_from_one_ALL_file = newArray();
		list_of_2P_log_r_values_from_one_ALL_file = newArray();




		offset = 0;
		normalization_factor = 0;
		phase = 0;
		phase2 = 0;
		B2P = 0;
		C2P = 0;
		polarization_direction = "";
		FIT_file_parsed = false;
		alpha0 = -1;
		sigma = -1;


		if (endsWith(file_list_2P[i], "_2P_FIT.txt")){
			contents_of_2P_FIT_file = File.openAsString(directory_2P +File.separator+ file_list_2P[i]);

			print(file_list_2P[i]+"num of pols "+number_of_polarizations_2P+" index of n_of_p = 1 "+indexOf(contents_of_2P_FIT_file, "number_of_polarizations = 1"));
			combined_fitting_1P_2P_parsing_summary = combined_fitting_1P_2P_parsing_summary + "################################\r\n"+file_list_2P[i]+"\r\nnumber of polarizations "+number_of_polarizations_2P+"\r\nindex of number of polarizations = "+indexOf(contents_of_2P_FIT_file, "\r\nnumber_of_polarizations = ")+"\r\n\r\n";


			// one polarization
			if(  number_of_polarizations_2P == 1 && indexOf(contents_of_2P_FIT_file, "number_of_polarizations = 1") > 0  ) {
				lines = split(contents_of_2P_FIT_file,"\n");
				for (j=0; j<lines.length; j++){
					if(startsWith(lines[j],"normalization_factor =")){
						normalization_factor_line = split(lines[j]," ");
						normalization_factor = parseFloat(normalization_factor_line[2]);
					}
					if(startsWith(lines[j],"offset =")){
						offset_line = split(lines[j]," ");
						offset = parseFloat(offset_line[2]);
					}
					if(startsWith(lines[j],"phase =")){
						phase_line = split(lines[j]," ");
						phase = parseFloat(phase_line[2]);
					}
					if(startsWith(lines[j],"B2P =")){
						B2P_line = split(lines[j]," ");
						B2P = parseFloat(B2P_line[2]);
					}
					if(startsWith(lines[j],"C2P =")){
						C2P_line = split(lines[j]," ");
						C2P = parseFloat(C2P_line[2]);
					}
					if(startsWith(lines[j],"polarization_direction =")){
						polarization_direction_line = split(lines[j]," ");
						polarization_direction = polarization_direction_line[2];
					}
					if(startsWith(lines[j],"alpha0 =")){
						alpha0_line = split(lines[j]," ");
						alpha0 = alpha0_line[2];
						alpha0 = replace(alpha0, degreechar, "");
						alpha0 = parseFloat(alpha0);
					}
					if(startsWith(lines[j],"sigma =")){
						sigma_line = split(lines[j]," ");
						sigma = sigma_line[2];
						sigma = replace(sigma, degreechar, "");
						sigma = parseFloat(sigma);
					}
					if(startsWith(lines[j],"r(max) ")){
						list_of_2P_fitting_results = list_of_2P_fitting_results + "\n" + lines[j+1];
					}
				}
				list_of_2P_fitting_results = list_of_2P_fitting_results + "\t" + alpha0 + "\t" + sigma;
				FIT_file_parsed = true;	

				list_of_2P_B2P_values = Array.concat(list_of_2P_B2P_values, B2P);
				list_of_2P_C2P_values = Array.concat(list_of_2P_C2P_values, C2P);
				log2rmax = log( (1 + B2P + C2P)/(1 - B2P + C2P) )/log(2);
				list_of_2P_log2rmax_values = Array.concat(list_of_2P_log2rmax_values, log2rmax);


				name_of_2P_ALL_file = replace(file_list_2P[i],"_2P_FIT.txt","_2P_ALL.txt");
				contents_of_2P_ALL_file = File.openAsString(directory_2P +File.separator+ name_of_2P_ALL_file);

				lines = split(contents_of_2P_ALL_file,"\n");
				for(j=1; j<lines.length; j++){    //skip the 1st line (headers)
					values = split(lines[j],"\t");
					theta = parseFloat(values[0]) - phase;
					if(polarization_direction == "vertical"){
						theta = theta - 3.14159/2;
					}
					while(theta < -3.14159/2){
						theta = theta + 3.14159;
					}
					while(theta > 3.14159/2){
						theta = theta - 3.14159;
					}
					r = parseFloat(values[5]);
					if(r > 0 && indexOf(r, "Infinity") < 0 && indexOf(r, "NaN") < 0){
						log2r = r / normalization_factor;
						list_of_2P_thetas_from_one_ALL_file = Array.concat(list_of_2P_thetas_from_one_ALL_file, theta + 3.14159);
						list_of_2P_thetas_from_one_ALL_file_deg = Array.concat(list_of_2P_thetas_from_one_ALL_file_deg, theta/3.14159*180 + 180);
						list_of_2P_log_r_values_from_one_ALL_file = Array.concat(list_of_2P_log_r_values_from_one_ALL_file, log2r);
					}
				}
	
				list_of_2P_thetas = Array.concat(list_of_2P_thetas, list_of_2P_thetas_from_one_ALL_file);
				list_of_2P_log_r_values = Array.concat(list_of_2P_log_r_values, list_of_2P_log_r_values_from_one_ALL_file);						 

				plot_legend = plot_legend + name_of_2P_ALL_file + "\n";
				Plot.setColor(color_list[color_count]);

				Plot.setLineWidth(1);
				Plot.add("dot", list_of_2P_thetas_from_one_ALL_file_deg, list_of_2P_log_r_values_from_one_ALL_file);
				color_count = color_count + 1;
				if(color_count > 10){
					color_count = color_count - 11;
				}	

				fitting_equation_2P = "y = a * (1 + b * cos(2*x) + c * cos(4*x) )";



			}
			else{ // two polarizations
			    if( number_of_polarizations_2P > 1 && indexOf(contents_of_2P_FIT_file, "number_of_polarizations = ") > 0)  {
					lines = split(contents_of_2P_FIT_file,"\n");
					for (j=0; j<lines.length; j++){
						if(startsWith(lines[j],"normalization_factor =")){
							normalization_factor_line = split(lines[j]," ");
							normalization_factor = parseFloat(normalization_factor_line[2]);
						}
						if(startsWith(lines[j],"offset =")){
							offset_line = split(lines[j]," ");
							offset = parseFloat(offset_line[2]);
						}
						if(startsWith(lines[j],"phase =")){
							phase_line = split(lines[j]," ");
							phase = parseFloat(phase_line[2]);
						}
						if(startsWith(lines[j],"phase2 =")){
							phase2_line = split(lines[j]," ");
							phase2 = parseFloat(phase2_line[2]);
						}
						if(startsWith(lines[j],"B2P =")){
							B2P_line = split(lines[j]," ");
							B2P = parseFloat(B2P_line[2]);
						}
						if(startsWith(lines[j],"C2P =")){
							C2P_line = split(lines[j]," ");
							C2P = parseFloat(C2P_line[2]);
						}
						if(startsWith(lines[j],"polarization_direction =")){
							polarization_direction_line = split(lines[j]," ");
							polarization_direction = parseFloat(polarization_direction_line[2]);
						}
						if(startsWith(lines[j],"alpha0 =")){
							alpha0_line = split(lines[j]," ");
							alpha0 = alpha0_line[2];
							alpha0 = replace(alpha0, degreechar, " ");
							alpha0 = parseFloat(alpha0);
						}
						if(startsWith(lines[j],"sigma =")){
							sigma_line = split(lines[j]," ");
							sigma = sigma_line[2];
							sigma = replace(sigma, degreechar, " ");
							sigma = parseFloat(sigma);
						}
						if(startsWith(lines[j],"r(max) ")){
							list_of_2P_fitting_results = list_of_2P_fitting_results + "\n" + lines[j+1];
						}
					}
					list_of_2P_fitting_results = list_of_2P_fitting_results + "\t" + alpha0 + "\t" + sigma;
					FIT_file_parsed = true;				 
				}

				
				if( number_of_polarizations_2P > 1 && indexOf(contents_of_2P_FIT_file, "number_of_polarizations = ") == -1  ){
					lines = split(contents_of_2P_FIT_file,"\n");
					offset_line = split(lines[6]," ");
					phase_line  = split(lines[5]," ");
					offset = parseFloat(offset_line[2]);
					phase  = parseFloat(phase_line[2]);
					B_C_line = split(lines[9]," ");
					B2P = parseFloat(replace(B_C_line[4],",","0"));
					C2P = parseFloat(B_C_line[7]);

					FIT_file_parsed = true;
				}
				list_of_2P_B2P_values = Array.concat(list_of_2P_B2P_values, B2P);
				list_of_2P_C2P_values = Array.concat(list_of_2P_C2P_values, C2P);
				log2rmax = log( (1 + B2P + C2P)/(1 - B2P + C2P) )/log(2);
				list_of_2P_log2rmax_values = Array.concat(list_of_2P_log2rmax_values, log2rmax);
				if(alpha0 > -1){
					list_of_2P_alpha0s = Array.concat(list_of_2P_alpha0s, alpha0);
					list_of_2P_sigmas  = Array.concat(list_of_2P_sigmas, sigma);			
				}


				if(FIT_file_parsed == true){   // parse the 2P_ALL.txt file
					name_of_2P_ALL_file = replace(file_list_2P[i],"_2P_FIT.txt","_2P_ALL.txt");
					contents_of_2P_ALL_file = File.openAsString(directory_2P +File.separator+ name_of_2P_ALL_file);

					lines = split(contents_of_2P_ALL_file,"\n");
					for(j=1; j<lines.length; j++){    //skip the 1st line (headers)
						values = split(lines[j],"\t");
						theta = parseFloat(values[0]);

//theta = theta + PI/2;						

						while(theta < -PI/2){theta = theta + PI;}
						while(theta >  PI/2){theta = theta - PI;}

						Fh = parseFloat(values[1]);
						Fv = parseFloat(values[2]);
						r = parseFloat(values[5]);
						if(r > 0 && indexOf(r, "Infinity") < 0 && indexOf(r, "NaN") < 0){
							log2r = log(r)/log(2)-offset;

							list_of_2P_thetas_from_one_ALL_file = Array.concat(list_of_2P_thetas_from_one_ALL_file, theta + 3.14159); 
							list_of_2P_thetas_from_one_ALL_file_deg = Array.concat(list_of_2P_thetas_from_one_ALL_file_deg, theta/3.14159*180 + 180); 
							list_of_2P_Fhs_from_one_ALL_file = Array.concat(list_of_2P_Fhs_from_one_ALL_file,Fh);
							list_of_2P_Fvs_from_one_ALL_file = Array.concat(list_of_2P_Fvs_from_one_ALL_file,Fv);
							list_of_2P_log_r_values_from_one_ALL_file = Array.concat(list_of_2P_log_r_values_from_one_ALL_file, log2r);

						}
					}


					if(phase2 == 0){
						phase2 = phase;
					}
					for(k=0; k<list_of_2P_log_r_values_from_one_ALL_file.length; k++){
						Fh = list_of_2P_Fhs_from_one_ALL_file[k];
						Fv = list_of_2P_Fvs_from_one_ALL_file[k];
						theta = list_of_2P_thetas_from_one_ALL_file[k];

						log2r_experiment = log(Fh/Fv)/log(2);
						log2r_skewed_fit = offset + log(  (1 + B2P*cos(2*theta-2*phase) + C2P*cos(4*theta-4*phase)) / (1 - B2P*cos(2*theta-2*phase2) + C2P*cos(4*theta-4*phase))  )/log(2);

						delta_log2r = log2r_experiment - log2r_skewed_fit;

						log2r_straight_fit = log(  (1 + B2P*cos(2*theta) + C2P*cos(4*theta)) / (1 - B2P*cos(2*theta) + C2P*cos(4*theta))  )/log(2);
						log2r_straightened_experiment = log2r_straight_fit + delta_log2r;

						list_of_2P_log_r_values_from_one_ALL_file[k] = log2r_straightened_experiment;

						
						if(isNaN(log2r_straightened_experiment)){}
						else{
							list_of_2P_thetas = Array.concat(list_of_2P_thetas, list_of_2P_thetas_from_one_ALL_file[k]);
							list_of_2P_log_r_values = Array.concat(list_of_2P_log_r_values, list_of_2P_log_r_values_from_one_ALL_file[k]);
						}
						
					}


				
					plot_legend = plot_legend + name_of_2P_ALL_file + "\n";
					Plot.setColor(color_list[color_count]);

Plot.setColor("gray");
					
					Plot.setLineWidth(1);
					Plot.add("dot", list_of_2P_thetas_from_one_ALL_file_deg, list_of_2P_log_r_values_from_one_ALL_file);
					color_count = color_count + 1;
					if(color_count > 10){
						color_count = color_count - 11;
					}					
				}

				fitting_equation_2P = "y = a + log( (1 + b * cos(2*x) + c * cos(4*x))/(1 - b * cos(2*x) + c * cos(4*x) ) )/log(2)";

			}

			print("parsed "+FIT_file_parsed+" offset "+offset+" norm f "+normalization_factor+" B2P "+B2P+" C2P "+C2P+" phase "+phase+" phase2 "+phase2+" pol direction "+polarization_direction);

			if(alpha0 == -1){
				showStatus("Calculating individual "+alpha0char+", "+sigmachar+" values");
				newImage("2P_SSD_image", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
				SSD_2P_image_id = getImageID();

				for(data_point_2P_index = 0; data_point_2P_index < list_of_2P_thetas_from_one_ALL_file.length; data_point_2P_index++){
					theta = list_of_2P_thetas_from_one_ALL_file[data_point_2P_index];
					log2r = list_of_2P_log_r_values_from_one_ALL_file[data_point_2P_index];
				
					if (number_of_polarizations_2P > 1){
						expectation_expression_2P = "log( (1 + A * cos(" + toString(2*theta) + ") + B * cos(" + toString(4*theta) + ") ) / (1 - A * cos(" + toString(2*theta) + ") + B * cos(" + toString(4*theta) + ") ) ) / log(2)";	
					}
					else{
						expectation_expression_2P = "1 + A * cos(" + toString(2*theta) + ") + B * cos(" + toString(4*theta) + ")";
					}
		
					run("Image Expression Parser (Macro)", "expression=["+ expectation_expression_2P +"] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");	
					showStatus("Calculating individual "+alpha0char+", "+sigmachar+" values");
					LD_prediction_image_id = getImageID();
					run("Subtract...", "value="+log2r);
					run("Square");
					imageCalculator("Add", SSD_2P_image_id, LD_prediction_image_id);
					showStatus("Calculating individual "+alpha0char+", "+sigmachar+" values");

					selectImage(LD_prediction_image_id);
					close();
				}	

				selectImage(SSD_2P_image_id);
				getRawStatistics(nPixels, mean, min, max);
				run("Macro...", "code=[if(isNaN(v)) v="+max+"]");  // replace 'NaN' by max sqdev found
				getRawStatistics(nPixels, mean, min, max);
	
				// find alpha0, sigma with minimum SSD		
				run("Find Maxima...", "noise="+max+" output=[Point Selection] light");
				getSelectionBounds(alpha0, sigma_row_number, w, h);
		
				list_of_2P_alpha0s = Array.concat(list_of_2P_alpha0s, alpha0);		
				list_of_2P_sigmas = Array.concat(list_of_2P_sigmas, num_of_sigmas - sigma_row_number - 1);		

				selectImage(SSD_2P_image_id);
				close();
									
			}


			//make/expand an image of 2P data, 2P data counts
			number_of_2P_FIT_files = number_of_2P_FIT_files + 1;

			if(number_of_2P_FIT_files == 1){
				newImage("2P_data", "32-bit black", number_of_theta_bins, number_of_2P_FIT_files, 1);
				image_of_2P_data_image_id = getImageID();
				newImage("2P_counts", "16-bit black", number_of_theta_bins, number_of_2P_FIT_files, 1);
				image_of_2P_counts_image_id = getImageID();				
			}
			else{
				selectImage(image_of_2P_data_image_id);
				run("Canvas Size...", "width="+number_of_theta_bins+" height="+number_of_2P_FIT_files+" position=Top-Center zero");	
	
				selectImage(image_of_2P_counts_image_id);
				run("Canvas Size...", "width="+number_of_theta_bins+" height="+number_of_2P_FIT_files+" position=Top-Center zero");	
			}
			for(j = 0; j < list_of_2P_thetas_from_one_ALL_file.length; j++){
				theta = list_of_2P_thetas_from_one_ALL_file[j];
				theta = theta - PI/2;	// theta is now between 0 and PI
				theta = round(theta/PI*180);   // theta is now in degrees
				selectImage(image_of_2P_data_image_id);
				pixel_log2r_init_val = getPixel(theta, number_of_2P_FIT_files-1);
				setPixel(theta, number_of_2P_FIT_files-1, getPixel(theta, number_of_2P_FIT_files-1)+list_of_2P_log_r_values_from_one_ALL_file[j]);
				pixel_log2r_final_val = getPixel(theta, number_of_2P_FIT_files-1);

				selectImage(image_of_2P_counts_image_id);
				pixel_count_init_val = getPixel(theta, number_of_2P_FIT_files-1);
				setPixel(theta, number_of_2P_FIT_files-1, getPixel(theta, number_of_2P_FIT_files-1)+1);
				pixel_count_final_val = getPixel(theta, number_of_2P_FIT_files-1);
			}
		}
	}			

	if(image_of_2P_data_image_id < 0 && image_of_2P_counts_image_id < 0){
		imageCalculator("Divide", image_of_2P_data_image_id, image_of_2P_counts_image_id);
		selectImage(image_of_2P_counts_image_id);
		close();	
	}
	

	Array.getStatistics(list_of_2P_log_r_values, list_of_2P_log_r_values_min, list_of_2P_log_r_values_max, list_of_2P_log_r_values_mean, list_of_2P_log_r_values_stdev); 

	print("\n\nList of "+alpha0char+" values (from individual 2PPM files):");
	Array.print(list_of_2P_alpha0s);

	print("List of "+sigmachar+" values (from individual 2PPM files):");
	Array.print(list_of_2P_sigmas);
	print("\n");




	// 2P statistics
	// If alpha0, sigma values are not present in the individual 2P_FIT.txt files, finds alpha0, sigma for each 2PPM image.
	// If alpha0, sigma values are listed in the 2P_FIT.txt files (which they should be, if the 2P_FIT.txt files were
	// created by macro v13.88 or later), then it just pools them together.
	// Once the list of alpha0, sigma values is available, it calculates mean +- 2*StdErr




	if(list_of_2P_alpha0s.length < list_of_2P_log2rmax_values.length && image_of_2P_data_image_id < 0){   // find alpha0, sigma values for individual 2PPM images, if the values cannot be found in the 2P_FIT.txt files
		list_of_2P_alpha0s = Array.trim(list_of_2P_alpha0s, 0);
		list_of_2P_sigmas = Array.trim(list_of_2P_sigmas, 0);
		
		
		newImage("2P_LD_predictions", "32-bit black", num_of_alpha0s+1, num_of_sigmas+1, number_of_theta_bins);
		LD_prediction_stack_image_id = getImageID();


		// this creates a stack of 2P LD predictions for all alpha0, sigma combinations and all (181) values of theta bins
		for(slice_number = 1; slice_number <= number_of_theta_bins; slice_number++){
			showStatus("Preparing 2P analysis");
			showProgress(slice_number, number_of_theta_bins);
			theta = (slice_number-1)/(number_of_theta_bins-1)*PI;
	
			if (number_of_polarizations_2P > 1){
				expectation_expression_2P = "log( (1 + A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ") ) / (1 - A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ") ) ) / log(2)";	
			}
			else{
				expectation_expression_2P = "1 + A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ")";
			}

			run("Image Expression Parser (Macro)", "expression=["+ expectation_expression_2P +"] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");	
			LD_prediction_image_id = getImageID();
			
			showStatus("Preparing 2P analysis");
	
			run("Select All");
			run("Copy");
	
			selectImage(LD_prediction_stack_image_id);
			setSlice(slice_number);
			run("Paste");
	
			selectImage(LD_prediction_image_id);
			close();
	
		}

		// this compares the LD prediction stack to individual 2P LD datasets
		for(data_set_number = 0; data_set_number < number_of_2P_FIT_files; data_set_number++){    
			showProgress(data_set_number, number_of_2P_FIT_files);

			selectImage(LD_prediction_stack_image_id);
			run("Duplicate...", "title=rmsd_stack duplicate");
			rmsd_stack_image_id = getImageID();
			for(theta = 0; theta < number_of_theta_bins; theta++){	// from each slice of the stack of 2P LD prediction subtract the observed value
				selectImage(image_of_2P_data_image_id);
				ld_value = getPixel(theta, data_set_number);
				selectImage(rmsd_stack_image_id);
				setSlice(theta+1);
				if( isNaN(ld_value) ){
					run("Set...", "value=0 slice");
				}
				else{
					run("Subtract...", "value="+ld_value+" slice");				
				}
				showStatus(" Obtaining "+alpha0char+", "+sigmachar+" values for each 2P image: "+toString(data_set_number)+"/"+toString(number_of_2P_FIT_files));
			}

			
			run("Square", "stack");		// calculate a square of the deviations (for each slice of the 2P LD prediction stack)
			run("Z Project...", "projection=[Sum Slices]");  // sum the square deviations
			z_projection_image_id = getImageID();


			// replace any 'NaN' pixels with the max sqdev found
			getRawStatistics(nPixels, mean, min, max);

			run("Macro...", "code=[if(isNaN(v)) v="+max+"]");  // replace 'NaN' by max sqdev found

			// find alpha0, sigma with minimum SSD		
			run("Find Maxima...", "noise="+max+" output=[Point Selection] light");
			getSelectionBounds(alpha0, sigma_row_number, w, h);
	
			list_of_2P_alpha0s = Array.concat(list_of_2P_alpha0s, alpha0);		
			list_of_2P_sigmas = Array.concat(list_of_2P_sigmas, num_of_sigmas - sigma_row_number - 1);		

			selectImage(rmsd_stack_image_id);
			close();
			selectImage(z_projection_image_id);
			close();
						
		}

		
		selectImage(LD_prediction_stack_image_id);
		close();
		
	}





//######################################################################################################################################
//##################data fitting########################################################################################################
//######################################################################################################################################






	setBatchMode(true);


	newImage("combined_1P_r_squared_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	r_squared_1P_image_id = getImageID();		

	newImage("combined_2P_r_squared_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	r_squared_2P_image_id = getImageID();		

	newImage("combined_1P_RMSD_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	RMSD_1P_image_id = getImageID();		

	newImage("combined_2P_RMSD_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	RMSD_2P_image_id = getImageID();		

	newImage("combined_1P_chi_sq.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	chi_sq_1P_image_id = getImageID();

	newImage("combined_2P_chi_sq.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	chi_sq_2P_image_id = getImageID();

	showStatus("Finding "+alpha0char+", "+sigmachar);


	min_RMSD_from_alpha0_sigma_search_1P = 100000;
	min_RMSD_from_alpha0_sigma_search_2P = 100000;
	min_RMSD_from_alpha0_sigma_search_1P_2P = 100000;
	min_alpha0_1P = 0;
	min_sigma_1P = 0;
	min_alpha0_2P = 1;
	min_sigma_2P = 0;
	min_alpha0_1P_2P = 0;
	min_sigma_1P_2P = 0;
	min_sigma_1P_row_number = -1;
	min_sigma_2P_row_number = -1;
	min_B_1P = 0;
	rmax_1P = 0;
	min_B_2P = 0;
	min_C_2P = 0;
	min_B_1P_2P = 0;
	min_C_1P_2P = 0;
	

	for(i = 0; i < list_of_1P_thetas.length; i++){
		showProgress(i, list_of_1P_thetas.length);
		showStatus("Calculating "+alpha0char+", "+sigmachar+" from combined 1P data");

		if(number_of_polarizations_1P == 1){
			expectation_expression = " 1 + A * cos(2 * " + toString(list_of_1P_thetas[i]) + " )";					
		}
		else{
			expectation_expression = "log( (1 + A * cos(2 * " + toString(list_of_1P_thetas[i]) + ") ) / (1 - A * cos(2 * " + toString(list_of_1P_thetas[i]) + ") ) ) / log(2)";		
		}
		run("Image Expression Parser (Macro)", "expression=["+ expectation_expression +"] a=1PPM_B.tif b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		expectation_image_id = getImageID();
		rename("expectation image");

		//RMSD
		run("Duplicate...", "title=sq_dev");
		sq_dev_image_id = getImageID();
		run("Subtract...", "value="+toString(list_of_1P_log_r_values[i]));
		run("Square");
		imageCalculator("Add 32-bit", RMSD_1P_image_id, sq_dev_image_id); 
		showStatus("Calculating "+alpha0char+", "+sigmachar+" from combined 1P data");


		//Chi-squared
		selectImage(expectation_image_id);
		run("Abs");
		run("Add...", "value=0.01");
		imageCalculator("Divide 32-bit", sq_dev_image_id, expectation_image_id);
		imageCalculator("Add 32-bit", chi_sq_1P_image_id, sq_dev_image_id);
		showStatus("Calculating "+alpha0char+", "+"sigmachar");
		
		selectImage(sq_dev_image_id);
		close();
		
		selectImage(expectation_image_id);
		close();
	}


	if(list_of_1P_thetas.length > 0){
		selectImage(RMSD_1P_image_id);
		run("Divide...", "value="+toString(list_of_1P_thetas.length));
		run("Square Root");
		getRawStatistics(nPixels, mean, min_RMSD_from_alpha0_sigma_search_1P, max_RMSD_from_alpha0_sigma_search_1P);
		run("Select All");
		run("Copy");


		selectImage(r_squared_1P_image_id);
		run("Paste");
		
		run("Divide...", "value="+list_of_1P_log_r_values_stdev);
		run("Square");
		run("Multiply...", "value=-1");
		run("Add...", "value=1");


		getRawStatistics(nPixels, mean, min_r_squared_from_alpha0_sigma_search_1P, max_r_squared_from_alpha0_sigma_search_1P);

		run("Find Maxima...", "noise="+max_r_squared_from_alpha0_sigma_search_1P+" output=[Point Selection]");
		getSelectionBounds(min_alpha0_1P, min_sigma_1P_row_number, selection_width, selection_height);

		if(min_alpha0_1P==0 && min_sigma_1P_row_number==0 && selection_width==num_of_alpha0s && selection_height==num_of_sigmas){
			for (col = 0; col <= num_of_alpha0s; col++) {
				for (row = 0; row <= num_of_sigmas; row++){
					if(getPixel(col,row)==max_r_squared_from_alpha0_sigma_search_1P){
						min_alpha0_1P = col;
						min_sigma_1P_row_number = row;
					}
				}
			}
		}


		min_sigma_1P = num_of_sigmas - 1 - min_sigma_1P_row_number;
		selectImage(B_1PPM_image_id);
		min_B_1P = getPixel(min_alpha0_1P, min_sigma_1P_row_number);
		rmax_1P = (1 + min_B_1P)/(1 - min_B_1P);
		selectImage(chi_sq_1P_image_id);
		run("Divide...", "value="+toString(list_of_1P_thetas.length*1.2));
		getRawStatistics(nPixels, mean, min_chi_sq, max_chi_sq);
		for (col = 0; col <= num_of_alpha0s; col++) {
			showProgress(col, num_of_alpha0s);
			for (row = 0; row <= num_of_sigmas; row++){
				if(getPixel(col,row) != getPixel(col,row)){
					setPixel(col,row,max_chi_sq);
				}
			}
		}
	
		run("phase");
		setMinAndMax(0, 1);
	
	}

print("list_of_2P_log_r_values:");
Array.print(list_of_2P_log_r_values);
	for(i = 0; i < list_of_2P_thetas.length; i++){
		showProgress(i, list_of_2P_thetas.length);
		showStatus("Calculating "+alpha0char+", "+sigmachar+" from combined 2P data");

		if(number_of_polarizations_2P == 1){
			expectation_expression = " (1 + B * cos(" + toString(2*list_of_2P_thetas[i]) + ") + C * cos(" + toString(4*list_of_2P_thetas[i]) + "))";					
		}
		else{
			expectation_expression = toString(offset) + toString("+ log( (1 + B * cos(" + toString(2*list_of_2P_thetas[i]) + ") + C * cos(" + toString(4*list_of_2P_thetas[i]) + ")) / (1 - B * cos(" + toString(2*list_of_2P_thetas[i]) + ") + C * cos(" + toString(4*list_of_2P_thetas[i]) + ") ) ) / log(2)");
		}



		run("Image Expression Parser (Macro)", "expression=["+ expectation_expression +"] b=2PPM_B.tif c=2PPM_C.tif d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		expectation_image_id = getImageID();
		rename("expectation image");


		//RMSD
		run("Duplicate...", "title=sq_dev");
		sq_dev_image_id = getImageID();

		run("Subtract...", "value="+toString(list_of_2P_log_r_values[i]));
		
		run("Square");
		imageCalculator("Add 32-bit", RMSD_2P_image_id, sq_dev_image_id); 
		showStatus("Calculating "+alpha0char+", "+sigmachar+" from combined 2P data");

		//Chi-squared
		selectImage(expectation_image_id);
		run("Abs");
		run("Add...", "value=0.01");
		imageCalculator("Divide 32-bit", sq_dev_image_id, expectation_image_id);
		imageCalculator("Add 32-bit", chi_sq_2P_image_id, sq_dev_image_id);
		showStatus("Calculating "+alpha0char+", "+sigmachar);



		
		selectImage(sq_dev_image_id);
		close();
		
		selectImage(expectation_image_id);
		close();
	}

	if(list_of_2P_thetas.length > 0){

		selectImage(RMSD_2P_image_id);
		run("Divide...", "value="+toString(list_of_2P_thetas.length));
		run("Square Root");
		getRawStatistics(nPixels, mean, min_RMSD_from_alpha0_sigma_search_2P, max_RMSD_from_alpha0_sigma_search_2P);
		run("Select All");
		run("Copy");

		selectImage(r_squared_2P_image_id);
		run("Paste");
		
		run("Divide...", "value="+list_of_2P_log_r_values_stdev);
		run("Square");
		run("Multiply...", "value=-1");
		run("Add...", "value=1");



		getRawStatistics(nPixels, mean, min_r_squared_from_alpha0_sigma_search_2P, max_r_squared_from_alpha0_sigma_search_2P);
		// eliminate invalid pixels from r-squared image
		for (col = 0; col <= num_of_alpha0s; col++) {
			showProgress(col, num_of_alpha0s);
			for (row = 0; row <= num_of_sigmas; row++){
				if(getPixel(col,row) != getPixel(col,row)){
					setPixel(col,row,min_r_squared_from_alpha0_sigma_search_2P);
				}
			}
		}



		run("Find Maxima...", "noise="+max_r_squared_from_alpha0_sigma_search_2P+" output=[Point Selection]");
		getSelectionBounds(min_alpha0_2P, min_sigma_2P_row_number, selection_width, selection_height);

		if(min_alpha0_2P==0 && min_sigma_2P_row_number==0 && selection_width==num_of_alpha0s && selection_height==num_of_sigmas){
			for (col = 0; col <= num_of_alpha0s; col++) {
				for (row = 0; row <= num_of_sigmas; row++){
					if(getPixel(col,row)==max_r_squared_from_alpha0_sigma_search_2P){
						min_alpha0_2P = col;
						min_sigma_2P_row_number = row;
					}
				}
			}
		}

		
		min_sigma_2P = num_of_sigmas - 1 - min_sigma_2P_row_number;
		selectImage(B_2PPM_image_id);
		min_B_2P = getPixel(min_alpha0_2P, min_sigma_2P_row_number);
		selectImage(C_2PPM_image_id);
		min_C_2P = getPixel(min_alpha0_2P, min_sigma_2P_row_number);
		rmax_2P = (1 + min_B_2P + min_C_2P)/(1 - min_B_2P + min_C_2P);


		selectImage(chi_sq_2P_image_id);
		run("Divide...", "value="+toString(list_of_2P_thetas.length*1.2));
		getRawStatistics(nPixels, mean, min_chi_sq, max_chi_sq);
		for (col = 0; col <= num_of_alpha0s; col++) {
			showProgress(col, num_of_alpha0s);
			for (row = 0; row <= num_of_sigmas; row++){
				if(getPixel(col,row) != getPixel(col,row)){
					setPixel(col,row,max_chi_sq);
				}
			}
		}
	
		run("phase");
		setMinAndMax(0, 1);

	
	}


	imageCalculator("Add create 32-bit", RMSD_1P_image_id, RMSD_2P_image_id);
	rename("combined_1P_2P_RMSD_a0_s.tif");
	RMSD_1P_2P_image_id = getImageID();
	getRawStatistics(nPixels, mean, min_RMSD_from_alpha0_sigma_search_1P_2P, max_RMSD_from_alpha0_sigma_search_1P_2P);

	if(list_of_1P_thetas.length > 0 && list_of_2P_thetas.length > 0){
		imageCalculator("Average create 32-bit", r_squared_1P_image_id, r_squared_2P_image_id);	
		rename("combined_1P_2P_r_squared_a0_s.tif");
	}
	else{
		if(list_of_1P_thetas.length > 0){
			selectImage(r_squared_1P_image_id);
			run("Duplicate...", "title=combined_1P_2P_r_squared_a0_s.tif");
		}
		if(list_of_2P_thetas.length > 0){
			selectImage(r_squared_2P_image_id);
			run("Duplicate...", "title=combined_1P_2P_r_squared_a0_s.tif");
		}	
	}
	r_squared_1P_2P_image_id = getImageID();



	getRawStatistics(nPixels, mean, min_r_squared_from_alpha0_sigma_search_1P_2P, max_r_squared_from_alpha0_sigma_search_1P_2P);
	run("Find Maxima...", "noise="+max_r_squared_from_alpha0_sigma_search_1P_2P+" output=[Point Selection]");
	getSelectionBounds(min_alpha0_1P_2P, min_sigma_1P_2P_row_number, selection_width, selection_height);

	if(min_alpha0_1P_2P==0 && min_sigma_1P_2P_row_number==0 && selection_width==num_of_alpha0s && selection_height==num_of_sigmas){
		for (col = 0; col <= num_of_alpha0s; col++) {
			for (row = 0; row <= num_of_sigmas; row++){
				if(getPixel(col,row)==max_r_squared_from_alpha0_sigma_search_1P_2P){
					min_alpha0_1P_2P = col;
					min_sigma_1P_2P_row_number = row;
				}
			}
		}
	}

	
	
	min_sigma_1P_2P = num_of_sigmas - 1 - min_sigma_1P_2P_row_number;
	selectImage(B_2PPM_image_id);
	min_B_1P_2P = getPixel(min_alpha0_1P_2P, min_sigma_1P_2P_row_number);
	selectImage(C_2PPM_image_id);
	min_C_1P_2P = getPixel(min_alpha0_1P_2P, min_sigma_1P_2P_row_number);


	imageCalculator("Add create 32-bit", chi_sq_1P_image_id, chi_sq_2P_image_id);
	rename("combined_1P_2P_chi_sq_a0_s.tif");
	chi_sq_1P_2P_image_id = getImageID();


////////////////////////////////////////////////////////////////////////////////////
//////////////////////// make hyperstack of results ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



	directory_1P_parts = split(directory_1P, "/\\");
	subdirectory_1P_name = directory_1P_parts[lengthOf(directory_1P_parts)-1];
	directory_2P_parts = split(directory_2P, "/\\");
	subdirectory_2P_name = directory_2P_parts[lengthOf(directory_2P_parts)-1];

	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	year = toString(year);
	month = month + 1;
	if(month < 10){month = "0"+toString(month);}
	else{month = toString(month);}
	if(dayOfMonth < 10){dayOfMonth = "0"+toString(dayOfMonth);}
	else{dayOfMonth = toString(dayOfMonth);}
	if(hour < 10){hour = "0"+toString(hour);}
	else{hour = toString(hour);}
	if(minute < 10){minute = "0"+toString(minute);}
	else{minute = toString(minute);}
	if(second < 10){second = "0"+toString(second);}
	else{second = toString(second);}


	hyperstack_name = "goodness_of_fit_"+subdirectory_1P_name+"_"+subdirectory_2P_name+"_"+year+"_"+month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second;
	hyperstack_info_text = year+"/"+month+"/"+dayOfMonth+", "+hour+":"+minute+":"+second+"\n";
	hyperstack_info_text = hyperstack_info_text + "1P data directory: "+directory_1P+"\n2P data directory: "+directory_2P+"\n";
	
	newImage(hyperstack_name, "32-bit color-mode", num_of_alpha0s, num_of_sigmas, 4, 1, 3);
	hyperstack_image_id = getImageID();
	setColor(100,100,255);  //  blue


	// 1P r^2
	selectImage(r_squared_1P_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(1);
	Overlay.drawString("r"+suptwochar+" 1P", 5, 20);
	Overlay.setPosition(1,1,1);
	Overlay.show;	
	run("Set Label...", "label=[r"+suptwochar+" 1P]");
	run("Paste");
	run("Cyan Hot");
	setMinAndMax(0, 1);
	getRawStatistics(nPixels, mean_r_square, min_r_square, max_r_square);

	
	// 1P RMSD
	selectImage(RMSD_1P_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(2);
	Overlay.drawString("RMSD 1P", 5, 20);
	Overlay.setPosition(2,1,1);
	Overlay.show;	
	run("Set Label...", "label=[RMSD 1P]");
	run("Paste");
	run("phase");
	setMinAndMax(min_RMSD_from_alpha0_sigma_search_1P, 5*min_RMSD_from_alpha0_sigma_search_1P);


	// 1P Chi^2	
	selectImage(chi_sq_1P_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(3);
	Overlay.drawString("Chi"+suptwochar+" 1P", 5, 20);
	Overlay.setPosition(3,1,1);
	Overlay.show;
	run("Set Label...", "label=[Chi"+suptwochar+" 1P]");
	run("Paste");
	run("phase");
	setMinAndMax(0, 1);

	// 1P alpha0, sigma
	selectImage(B_1PPM_image_id);
	run("Duplicate...", "title=alpha0s_sigmas_within_95CI");
	alpha0s_sigmas_within_95CI_image_id = getImageID();
	run("Macro...", "code=[if(v > "+toString(list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr)+" && v <"+toString(list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr)+") v=1; else v = 0]");
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(4);
	run("Set Label...", "label=["+alpha0char+", "+sigmachar+" 1P]");
	Overlay.drawString(alpha0char+", "+sigmachar+" 1P", 5, 20);
	Overlay.setPosition(4,1,1);
	Overlay.show;
	run("Paste");
	run("Grays");
	selectImage(alpha0s_sigmas_within_95CI_image_id);
	close();

	// 2P r^2

	setColor(255,100,100);//light red
	selectImage(r_squared_2P_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(5);
	Overlay.drawString("r"+suptwochar+" 2P", 5, 20);
	Overlay.setPosition(1,1,2);
	Overlay.show;	
	run("Set Label...", "label=[r"+suptwochar+" 2P]");
	run("Paste");
	run("Cyan Hot");
	setMinAndMax(0, 1);


	// 2P RMSD
	selectImage(RMSD_2P_image_id);
	run("Copy");
	run("Specify...", "width=1 height=1 x="+min_alpha0_2P+" y="+min_sigma_2P_row_number);
	selectImage(hyperstack_image_id);
	setSlice(6);
	Overlay.drawString("RMSD 2P", 5, 20);
	Overlay.setPosition(2,1,2);
	Overlay.show;
	run("Set Label...", "label=[RMSD 2P]");
	run("Paste");
	run("phase");
	setMinAndMax(min_RMSD_from_alpha0_sigma_search_2P, 5*min_RMSD_from_alpha0_sigma_search_2P);


	// 2P Chi^2
	selectImage(chi_sq_2P_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(7);
	Overlay.drawString("Chi"+suptwochar+" 2P", 5, 20);
	Overlay.setPosition(3,1,2);
	Overlay.show;
	run("Set Label...", "label=[Chi"+suptwochar+" 2P]");
	run("Paste");
	run("phase");
	setMinAndMax(0, 1);


	// 2P alpha0, sigma
	selectImage(hyperstack_image_id);
	setSlice(8);
	Overlay.drawString(alpha0char+", "+sigmachar+" 2P", 5, 20);
	Overlay.setPosition(4,1,2);
	Overlay.show;
	run("Set Label...", "label=["+alpha0char+", "+sigmachar+" 2P]");
	for (i = 0; i < list_of_2P_alpha0s.length; i++) {
		alpha0 = list_of_2P_alpha0s[i];
		sigma = list_of_2P_sigmas[i];
		sigma_row_number = num_of_sigmas - 1 - sigma;
		setPixel(alpha0, sigma_row_number, getPixel(alpha0, sigma_row_number)+1);
	}
	setPixel(min_alpha0_2P, min_sigma_2P_row_number, 255);
	run("Grays");



	// 1P, 2P r^2
	setColor(255,100,255);//light magenta
	selectImage(r_squared_1P_2P_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(9);
	Overlay.drawString("r"+suptwochar+" 1P, 2P", 5, 20);
	Overlay.setPosition(1,1,3);
	Overlay.show;	
	run("Set Label...", "label=[r"+suptwochar+" 1P, 2P]");
	run("Paste");
	run("Cyan Hot");
	setMinAndMax(0, 1);



	// 1P, 2P RMSD
	selectImage(RMSD_1P_2P_image_id);
	run("Copy");
	run("Specify...", "width=1 height=1 x="+min_alpha0_1P_2P+" y="+min_sigma_1P_2P_row_number);
	selectImage(hyperstack_image_id);
	setSlice(10);
	Overlay.drawString("RMSD 1P, 2P", 5, 20);
	Overlay.setPosition(2,1,3);
	Overlay.show;
	run("Set Label...", "label=[RMSD 1P, 2P]");
	run("Paste");

	run("phase");
	setMinAndMax(min_RMSD_from_alpha0_sigma_search_1P_2P, 5*min_RMSD_from_alpha0_sigma_search_1P_2P);


	// 1P, 2P Chi^2
	selectImage(chi_sq_1P_2P_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(11);
	Overlay.drawString("Chi"+suptwochar+" 1P, 2P", 5, 20);
	Overlay.setPosition(3,1,3);
	Overlay.show;
	run("Set Label...", "label=[Chi"+suptwochar+" 1P, 2P]");
	run("Paste");

	run("phase");
	setMinAndMax(0, 1);

	// 1P, 2P alpha0, sigma
	selectImage(hyperstack_image_id);
	setSlice(12);
	Overlay.drawString(alpha0char+", "+sigmachar+" 1P, 2P", 5, 20);
	Overlay.setPosition(4,1,3);
	Overlay.show;
	run("Set Label...", "label=["+alpha0char+", "+sigmachar+" 1P, 2P]");
	setPixel(min_alpha0_1P_2P, min_sigma_1P_2P_row_number, 255);
	run("Grays");




	// do 1PPM statistics


	selectImage(B_1PPM_image_id);
	run("Select None");
	setColor(50,50,255);  //  blue
	for(row = 1; row < num_of_sigmas; row++){
		for(col = 1; col < num_of_alpha0s; col++){
			mark_pixels = false;
			if( ( getPixel(col-1, row) < (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) > (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col-1, row) < (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) > (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col-1, row) > (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) < (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col-1, row) > (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) < (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if(mark_pixels == true){
				selectImage(hyperstack_image_id);
				setSlice(1);
				Overlay.drawLine(col-1, row, col, row);
				selectImage(B_1PPM_image_id);
				mark_pixels = false;
			}
			if( ( getPixel(col, row-1) < (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) > (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col, row-1) < (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) > (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col, row-1) > (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) < (list_of_1P_B1P_values_mean - list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if( ( getPixel(col, row-1) > (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) && ( getPixel(col, row) < (list_of_1P_B1P_values_mean + list_of_1P_B1P_values_2x_stdErr) ) ){mark_pixels = true;}
			if(mark_pixels == true){
				selectImage(hyperstack_image_id);
				setSlice(1);
				Overlay.drawLine(col, row-1, col, row);
				selectImage(B_1PPM_image_id);
				mark_pixels = false;
			}
		}
	}

	selectImage(hyperstack_image_id);
	Overlay.show;

	if(list_of_1P_thetas.length > 0){  // print 1P data statistics
		statistics_1P_text = "1P data statistics:\n";
		statistics_1P_text = statistics_1P_text + "log"+subtwochar+"(rmax) = "+toString(round(list_of_1P_log2rmax_values_mean*100)/100)+plusminuschar+toString(round(list_of_1P_log2rmax_values_2x_stdErr*100)/100)+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+list_of_1P_log2rmax_values_min+", max = "+list_of_1P_log2rmax_values_max+", stdev = "+list_of_1P_log2rmax_values_stdDev+")\n";	
		statistics_1P_text_few_special_chars = "1P data statistics:\n";
		statistics_1P_text_few_special_chars = statistics_1P_text_few_special_chars + "log2(rmax) = "+toString(round(list_of_1P_log2rmax_values_mean*100)/100)+plusminuschar+toString(round(list_of_1P_log2rmax_values_2x_stdErr*100)/100)+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+list_of_1P_log2rmax_values_min+", max = "+list_of_1P_log2rmax_values_max+", stdev = "+list_of_1P_log2rmax_values_stdDev+")\n";	

		print("\n"+statistics_1P_text);
	}
	



	setSlice(1);
	run("Select None");



	if(list_of_2P_alpha0s.length > 0){  // indicate the alpha0, sigma variance among individual 2P data sets
		Array.getStatistics(list_of_2P_alpha0s, alpha0_2P_min, alpha0_2P_max, alpha0_2P_mean, alpha0_2P_stdev);
		Array.getStatistics(list_of_2P_sigmas, sigma_2P_min, sigma_2P_max, sigma_2P_mean, sigma_2P_stdev);
		Array.getStatistics(list_of_2P_log2rmax_values, list_of_2P_log2rmax_min, list_of_2P_log2rmax_max, list_of_2P_log2rmax_mean, list_of_2P_log2rmax_stdev);
		statistics_2P_text = "2P data statistics:\n";
		statistics_2P_text = statistics_2P_text + alpha0char+" = "+toString(round(alpha0_2P_mean*10)/10)+degreechar+plusminuschar+toString(round((2*alpha0_2P_stdev/sqrt(list_of_2P_alpha0s.length))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+alpha0_2P_min+degreechar+", max = "+alpha0_2P_max+degreechar+", stdev = "+alpha0_2P_stdev+degreechar+")\n";
		statistics_2P_text = statistics_2P_text + sigmachar+" = "+toString(round(sigma_2P_mean*10)/10)+degreechar+plusminuschar+toString(round((2*sigma_2P_stdev/sqrt(list_of_2P_sigmas.length))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+sigma_2P_min+degreechar+", max = "+sigma_2P_max+degreechar+", stdev = "+sigma_2P_stdev+degreechar+")\n";
		statistics_2P_text = statistics_2P_text + "log"+subtwochar+"(rmax) = "+toString(round(list_of_2P_log2rmax_mean*100)/100)+plusminuschar+toString(round((2*list_of_2P_log2rmax_stdev/sqrt(list_of_2P_log2rmax_values.length))*100)/100)+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+list_of_2P_log2rmax_min+", max = "+list_of_2P_log2rmax_max+", stdev = "+list_of_2P_log2rmax_stdev+")\n";
		statistics_2P_text_few_special_chars = "2P data statistics:\n";
		statistics_2P_text_few_special_chars = statistics_2P_text_few_special_chars + "alpha0 = "+toString(round(alpha0_2P_mean*10)/10)+degreechar+plusminuschar+toString(round((2*alpha0_2P_stdev/sqrt(list_of_2P_alpha0s.length))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+alpha0_2P_min+degreechar+", max = "+alpha0_2P_max+degreechar+", stdev = "+alpha0_2P_stdev+degreechar+")\n";
		statistics_2P_text_few_special_chars = statistics_2P_text_few_special_chars + "sigma = "+toString(round(sigma_2P_mean*10)/10)+degreechar+plusminuschar+toString(round((2*sigma_2P_stdev/sqrt(list_of_2P_sigmas.length))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+sigma_2P_min+degreechar+", max = "+sigma_2P_max+degreechar+", stdev = "+sigma_2P_stdev+degreechar+")\n";
		statistics_2P_text_few_special_chars = statistics_2P_text_few_special_chars + "log2(rmax) = "+toString(round(list_of_2P_log2rmax_mean*100)/100)+plusminuschar+toString(round((2*list_of_2P_log2rmax_stdev/sqrt(list_of_2P_log2rmax_values.length))*100)/100)+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+list_of_2P_log2rmax_min+", max = "+list_of_2P_log2rmax_max+", stdev = "+list_of_2P_log2rmax_stdev+")\n";
		print(statistics_2P_text);
		

		selectImage(hyperstack_image_id);	// created on lines #

		for(i = 0; i < list_of_2P_alpha0s.length; i++){
			setPixel(list_of_2P_alpha0s[i], num_of_sigmas - list_of_2P_sigmas[i] - 1, getPixel(list_of_2P_alpha0s[i], num_of_sigmas - list_of_2P_sigmas[i] - 1)+1);
		}

		setColor(255,0,0);  // variance of alpha0, sigma values from individual 2P traces
		Overlay.drawRect(round(alpha0_2P_mean-2*alpha0_2P_stdev), round(num_of_sigmas-sigma_2P_mean-1 - 2*sigma_2P_stdev), 4*alpha0_2P_stdev, 4*sigma_2P_stdev);
		Overlay.add;
		Overlay.drawLine(alpha0_2P_mean, num_of_sigmas-sigma_2P_mean-1, alpha0_2P_mean, num_of_sigmas-sigma_2P_mean-1);
		Overlay.add;
		Overlay.show;		
	}

	if(list_of_2P_thetas.length > 0){  // mark alpha0, sigma from the fit of all 2P data
		setColor(255,150,150);  // lighter red than used for the alpha0, sigma interval rectangle
		Overlay.drawLine(min_alpha0_2P, min_sigma_2P_row_number, min_alpha0_2P, min_sigma_2P_row_number);
		Overlay.show;
	}


	if(list_of_2P_thetas.length > 0){  // mark alpha0, sigma from the fit of all 1P, 2P data
		setColor(255,0,255);  // lighter magenta than for the 1P, 2P data rectangle
		Overlay.drawLine(min_alpha0_1P_2P, min_sigma_1P_2P_row_number, min_alpha0_1P_2P, min_sigma_1P_2P_row_number);
		Overlay.show;
	}

	
	run("Properties...", "channels=4 slices=1 frames=3 unit="+degreechar+" pixel_width=1.0000 pixel_height=-1 voxel_depth=1.0000 frame=[0.5 sec] origin=0,90");

	setBatchMode("show");



	selectImage(r_squared_1P_image_id);
	close();
	selectImage(RMSD_1P_image_id);
	close();
	selectImage(chi_sq_1P_image_id);
	close();
	selectImage(r_squared_2P_image_id);
	close();
	selectImage(RMSD_2P_image_id);
	close();
	selectImage(chi_sq_2P_image_id);
	close();
	selectImage(r_squared_1P_2P_image_id);
	close();
	selectImage(RMSD_1P_2P_image_id);
	close();
	selectImage(chi_sq_1P_2P_image_id);
	close();



	showStatus("Plotting...");

	

	if(list_of_1P_thetas.length > 0){
		x_1P = newArray(101);
		x_1P_deg = newArray(101);
		fit_1P = newArray(101);
		for(i = 0; i < 101; i++){
			showProgress(i, 100);
			x_1P_deg[i] = -90+180/100*i;
			x_1P[i] = x_1P_deg[i]/180*PI;
			if(number_of_polarizations_1P > 1){
				fit_1P[i] = log((1 + min_B_1P*cos(2*x_1P[i]))/(1 - min_B_1P*cos(2*x_1P[i])))/log(2);
			}
			else{
				fit_1P[i] = 1 + min_B_1P*cos(2*x_1P[i]);				
			}
		} 
		Plot.setColor("blue");
		Plot.setLineWidth(2);
		Plot.add("line", x_1P_deg, fit_1P);
		Plot.addText("--- 1PPM fit: r(max) = "+rmax_1P+" ("+toString(1/rmax_1P)+"), B(1P) = "+min_B_1P+", r"+suptwochar+" = "+max_r_squared_from_alpha0_sigma_search_1P, 0, 0.1); 
		Plot.addText("1PPM compiled data: log2r(max) = "+list_of_1P_log2rmax_values_mean+plusminuschar+list_of_1P_log2rmax_values_2x_stdErr+", B(1P) = "+list_of_1P_B1P_values_mean+plusminuschar+list_of_1P_B1P_values_2x_stdErr, 0, 0.175); 

	}		

	if(list_of_2P_thetas.length > 0){
		initialGuesses = newArray(1, min_B_2P, min_C_2P);
		Fit.doFit(fitting_equation_2P, list_of_2P_thetas, list_of_2P_log_r_values, initialGuesses);
		min_B_2P_non_gaussian = Fit.p(1);
		min_C_2P_non_gaussian = Fit.p(2);
		r_squared_2P_non_gaussian_fit = Fit.rSquared;
		rmax_2P_non_gaussian = (1 + min_B_2P_non_gaussian + min_C_2P_non_gaussian)/(1 - min_B_2P_non_gaussian + min_C_2P_non_gaussian);
		log2rmax_2P_non_gaussian = log(rmax_2P_non_gaussian)/log(2);
		
		x_2P_a = newArray(101);
		x_2P_a_deg = newArray(101);
		fit_2P_a = newArray(101);
		fit_2P_a_non_gaussian = newArray(101);
		x_2P_b = newArray(101);
		x_2P_b_deg = newArray(101);
		fit_2P_b = newArray(101);
		fit_2P_b_non_gaussian = newArray(101);
		min_B_1P = (-7*min_B_2P + 4*min_C_2P)/(-10 + min_B_2P + 2*min_C_2P);
		min_B_1P_non_gaussian = (-7*min_B_2P_non_gaussian + 4*min_C_2P_non_gaussian)/(-10 + min_B_2P_non_gaussian + 2*min_C_2P_non_gaussian);
		for(i = 0; i < 101; i++){
			showProgress(i, 100);
			x_2P_a_deg[i] = -90+180/100*i;
			x_2P_a[i] = x_2P_a_deg[i]/180*PI;
			if(number_of_polarizations_1P == 1 || (number_of_polarizations_1P == 0 && number_of_polarizations_2P == 1)){
				fit_2P_a[i] = 1 + min_B_1P*cos(2*x_2P_a[i]);				
				fit_2P_a_non_gaussian[i] = 1 + min_B_1P_non_gaussian*cos(2*x_2P_a[i]);				
			}
			else{
				fit_2P_a[i] = log((1 + min_B_1P*cos(2*x_2P_a[i]))/(1 - min_B_1P*cos(2*x_2P_a[i])))/log(2);			
				fit_2P_a_non_gaussian[i] = log((1 + min_B_1P_non_gaussian*cos(2*x_2P_a[i]))/(1 - min_B_1P_non_gaussian*cos(2*x_2P_a[i])))/log(2);			
			}
		} 
		for(i = 0; i < 101; i++){
			showProgress(i, 100);
			x_2P_b_deg[i] = 90+180/100*i;
			x_2P_b[i] = x_2P_b_deg[i]/180*PI;
			if(number_of_polarizations_2P == 1){
				fit_2P_b[i] = 1 + min_B_2P*cos(2*x_2P_b[i]) + min_C_2P*cos(4*x_2P_b[i]);
				fit_2P_b_non_gaussian[i] = 1 + min_B_2P_non_gaussian*cos(2*x_2P_b[i]) + min_C_2P_non_gaussian*cos(4*x_2P_b[i]);
			}
			else{
				fit_2P_b[i] = log((1 + min_B_2P*cos(2*x_2P_b[i]) + min_C_2P*cos(4*x_2P_b[i]))/(1 - min_B_2P*cos(2*x_2P_b[i]) + min_C_2P*cos(4*x_2P_b[i])))/log(2);
				fit_2P_b_non_gaussian[i] = log((1 + min_B_2P_non_gaussian*cos(2*x_2P_b[i]) + min_C_2P_non_gaussian*cos(4*x_2P_b[i]))/(1 - min_B_2P_non_gaussian*cos(2*x_2P_b[i]) + min_C_2P_non_gaussian*cos(4*x_2P_b[i])))/log(2);
			}
		} 

		Plot.setColor("#8b0015"); //dark red/marroon
		Plot.setLineWidth(2);
		Plot.add("line", x_2P_a_deg, fit_2P_a_non_gaussian);
		Plot.add("line", x_2P_b_deg, fit_2P_b_non_gaussian);
		Plot.addText("--- 2PPM fit (non-gaussian): r(max) = "+rmax_2P_non_gaussian+", B(2P) = "+min_B_2P_non_gaussian+", C(2P) = "+min_C_2P_non_gaussian+", r"+suptwochar+" = "+r_squared_2P_non_gaussian_fit, 0, 0.25); 
		
		Plot.setColor("red");
		Plot.setLineWidth(2);
		Plot.add("line", x_2P_a_deg, fit_2P_a);
		Plot.add("line", x_2P_b_deg, fit_2P_b);
		Plot.addText("--- 2PPM fit: "+alpha0char+" = "+min_alpha0_2P+degreechar+", "+sigmachar+" = "+min_sigma_2P+degreechar+", r(max) = "+rmax_2P+", B(2P) = "+min_B_2P+", C(2P) = "+min_C_2P+", r"+suptwochar+" = "+max_r_squared_from_alpha0_sigma_search_2P, 0, 0.325); 
		
	
	}


	if(list_of_2P_thetas.length > 0 && list_of_1P_thetas.length > 0){
		x_1P_2P_a = newArray(101);
		x_1P_2P_a_deg = newArray(101);
		fit_1P_2P_a = newArray(101);
		x_1P_2P_b = newArray(101);
		x_1P_2P_b_deg = newArray(101);
		fit_1P_2P_b = newArray(101);
		min_B_1P = (-7*min_B_1P_2P + 4*min_C_1P_2P)/(-10 + min_B_1P_2P + 2*min_C_1P_2P);
		for(i = 0; i < 101; i++){
			showProgress(i, 100);
			x_1P_2P_a_deg[i] = -90+180/100*i;
			x_1P_2P_a[i] = x_1P_2P_a_deg[i]/180*PI;
			if(number_of_polarizations_1P == 1){
				fit_1P_2P_a[i] = 1 + min_B_1P*cos(2*x_1P_2P_a[i]);				
			}
			else{
				fit_1P_2P_a[i] = log((1 + min_B_1P*cos(2*x_1P_2P_a[i]))/(1 - min_B_1P*cos(2*x_1P_2P_a[i])))/log(2);			
			}
		} 
		for(i = 0; i < 101; i++){
			showProgress(i, 100);
			x_1P_2P_b_deg[i] = 90+180/100*i;
			x_1P_2P_b[i] = x_1P_2P_a_deg[i]/180*PI;
			if(number_of_polarizations_2P == 1){
				fit_1P_2P_b[i] = 1 + min_B_1P_2P*cos(2*x_1P_2P_b[i]) + min_C_1P_2P*cos(4*x_1P_2P_b[i]);
			}
			else{
				fit_1P_2P_b[i] = log((1 + min_B_1P_2P*cos(2*x_1P_2P_b[i]) + min_C_1P_2P*cos(4*x_1P_2P_b[i]))/(1 - min_B_1P_2P*cos(2*x_1P_2P_b[i]) + min_C_1P_2P*cos(4*x_1P_2P_b[i])))/log(2);				
			}
		} 
		Plot.setColor("#b900ff");	//magenta with a bit extra blue in it	
		Plot.setLineWidth(2);
		Plot.add("line", x_1P_2P_a_deg, fit_1P_2P_a);
		Plot.add("line", x_1P_2P_b_deg, fit_1P_2P_b);
		Plot.setColor("#b900ff");	//magenta with a bit extra blue in it	
		Plot.addText("--- Combined 1PPM/2PPM fit: "+alpha0char+" = "+min_alpha0_1P_2P+degreechar+", "+sigmachar+" = "+min_sigma_1P_2P+degreechar+", B(2P) = "+min_B_1P_2P+", C(2P) = "+min_C_1P_2P+", r"+suptwochar+" = "+max_r_squared_from_alpha0_sigma_search_1P_2P, 0, 0.4); 
	}



	Plot.setColor("black");
	if(number_of_polarizations_1P == 1 || (number_of_polarizations_1P == 0 && number_of_polarizations_2P == 1)){
		Plot.addText("   1P data (normalized fluor. intensity)", 0, 0);		 
	}
	else{
		Plot.addText("                         1P data (log"+subtwochar+"(r))", 0, 0);
	}
	Plot.setJustification("right");
	if(number_of_polarizations_2P == 1){
		Plot.addText("2P data  (normalized fluor. intensity)", 1, 0);
	}
	else{
		Plot.addText("2P data (log"+subtwochar+"(r))                       ", 1, 0);
	}
	Plot.addLegend(plot_legend, "Auto");
	Plot.addLegend(plot_legend, "No Legend Transparent");
	if(number_of_polarizations_1P > 1 || number_of_polarizations_2P > 1){
		Plot.setLimits(-90,270,NaN,NaN);
	}
	else{
		Plot.setLimits(-90,270,0,NaN);	
	}
	Plot.setFrameSize(450, 200);
	Plot.show();

	saveAs("Tiff", directory_1P+File.separator+"plot_of_1PPM_2PPM_data.tif");
	saveAs("png", directory_1P+File.separator+"plot_of_1PPM_2PPM_data.png");
	saveAs("Tiff", directory_2P+File.separator+"plot_of_1PPM_2PPM_data.tif");
	saveAs("png", directory_2P+File.separator+"plot_of_1PPM_2PPM_data.png");




	setBatchMode("show");

	fitting_summary = File.open(directory_1P+File.separator+"combined_1P_2P_data_fitting_log_"+year+"_"+month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second+".txt");
	print(fitting_summary,combined_fitting_1P_2P_parsing_summary);
	print(fitting_summary,combined_fitting_1P_2P_SUMMARY_text);
	print(fitting_summary,"rmax\t1/rmax\tlog2(rmax)\tphase/phase2\tr^2");
	print(fitting_summary,list_of_1P_fitting_results);
	print(fitting_summary,"\r\n\r\n2P data directory: "+directory_2P+"\r\n");
	print(fitting_summary,"rmax\t1/rmax\tlog2(rmax)\tphase/phase2\tr^2\talpha0(deg)\tsigma(deg)");
	print(fitting_summary,list_of_2P_fitting_results);
	File.close(fitting_summary);

	if(directory_1P != directory_2P){
		fitting_summary = File.open(directory_1P+File.separator+"combined_1P_2P_data_fitting_log_"+year+"_"+month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second+".txt");
		print(fitting_summary,combined_fitting_1P_2P_parsing_summary);
		print(fitting_summary,combined_fitting_1P_2P_SUMMARY_text);
		print(fitting_summary,"rmax\t1/rmax\tlog2(rmax)\tphase/phase2\tr^2");
		print(fitting_summary,list_of_1P_fitting_results);
		print(fitting_summary,"\r\n\r\n2P data directory: "+directory_2P+"\r\n");
		print(fitting_summary,"rmax\t1/rmax\tlog2(rmax)\tphase/phase2\tr^2\talpha0(deg)\tsigma(deg)");
		print(fitting_summary,list_of_2P_fitting_results);
		File.close(fitting_summary);
	}



	selectImage(hyperstack_image_id);


	setBatchMode("show");



//


////////////////////////////////////////////////////////////////////////////////////
//////////////////  do bootstrapping statistics on 1P, 2P data  ////////////////////
////////////////////////////////////////////////////////////////////////////////////


	
	Dialog.create("1P, 2P Boostrapping Statistics");
	Dialog.addCheckbox("Calculate statistics", false);
	Dialog.addNumber("Number of synthetic data sets", 100);
	Dialog.addMessage("Statistical analysis on the combined data is performed by \nbootstrapping: generating and analyzing data sets similar \nto the experimental data set. The higher the number of data sets, \nthe more accurate the analysis. \n\nAt least 100 synthetic data sets should be used.");

	Dialog.show();
	run_stats = Dialog.getCheckbox();
	number_of_subsets = Dialog.getNumber();
	


	if(run_stats == true && list_of_1P_thetas.length > 0 && list_of_2P_thetas.length > 0){

		setBatchMode(true);
		prediction_stack_image_id = 0;
		

		

		// combined 1P, 2P statistics	
		// parse 1P and 2P data, generate a synthetic data set (100 samples of real data sets)
		list_of_1P_2P_alpha0s = newArray(number_of_subsets);
		list_of_1P_2P_sigmas  = newArray(number_of_subsets);


		
		if(image_of_1P_data_image_id < -1 && image_of_2P_data_image_id < -1){
			showStatus("Preparing synthetic data sets");
			showProgress(0.1);
	
			// 1P subsets:
			// 1) Convert the image data into a stack, with each data set in one slice; replace 'NaN' with '0'
			// 2) Make another stack, but holding 'counts', by replacing 'NaN' by '0', everything else by '1'
			// 3) Generate a stack of randomly chosen data sets by selecting randomly slices from the data stack
			// 4) Combine (sum up) slabs of this stack of thickness the same as the number of original data sets, thereby making a stack of synthetic data sets
			// 5) Convert this stack to an image (z -> y)
			// 6) Divide this image by a similarly made image of counts
			selectImage(image_of_1P_data_image_id);
			run("Duplicate...", "1P_counts");
			showStatus("Preparing synthetic data sets: 1PPM");
			rename("1P_counts");
			showStatus("Preparing synthetic data sets: 1PPM");
			image_of_1P_data_counts_image_id = getImageID();
			run("Macro...", "code=v=(v==v)");            // replace 'NaN' pixels by '0', everything else by '1' 
			showStatus("Preparing synthetic data sets: 1PPM");
			selectImage(image_of_1P_data_image_id);
			run("Macro...", "code=[if(isNaN(v)) v=0]");  // replace 'NaN' by '0'
			showStatus("Preparing synthetic data sets: 1PPM");
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			showStatus("Preparing synthetic data sets: 1PPM");
			image_of_1P_data_stack_id = getImageID();
			selectImage(image_of_1P_data_counts_image_id);
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			showStatus("Preparing synthetic data sets: 1PPM");
			image_of_1P_data_counts_stack_id = getImageID();
			rename("1P counts, stack");
			showStatus("Preparing synthetic data sets: 1PPM");
	
			list_of_1P_slice_numbers = toString(1+floor(random * number_of_1P_FIT_files));
			for(j = 1; j < number_of_subsets * number_of_1P_FIT_files; j++){
				showProgress(j, number_of_subsets * number_of_1P_FIT_files);
				list_of_1P_slice_numbers = list_of_1P_slice_numbers + ","+toString(1+floor(random * number_of_1P_FIT_files));
			}
	
	
			selectImage(image_of_1P_data_stack_id);
			run("Make Substack...", "  slices="+list_of_1P_slice_numbers);
			showStatus("Preparing synthetic data sets: 1PPM");
			substack_1P_data_before_grouped_z_project_image_id = getImageID();
			rename("1P data sets generated by Make Substack");
			run("Grouped Z Project...", "projection=[Sum Slices] group="+number_of_1P_FIT_files);
			showStatus("Preparing synthetic data sets: 1PPM");
			rename("1P data grouped Z projection");
			grouped_z_projection_image_id = getImageID();
			run("Make Montage...", "columns=1 rows="+number_of_subsets+" scale=1");
			showStatus("Preparing synthetic data sets: 1PPM");
			rename("1P data Make Montage");
			showStatus("Preparing synthetic data sets: 1PPM");
			combined_subset_1P_data_image_id = getImageID();
			selectImage(substack_1P_data_before_grouped_z_project_image_id);
			close();
			selectImage(grouped_z_projection_image_id);
			close();
			selectImage(image_of_1P_data_stack_id);
			close();
	
	
			selectImage(image_of_1P_data_counts_stack_id);
			run("Make Substack...", "  slices="+list_of_1P_slice_numbers);
			showStatus("Preparing synthetic data sets: 1PPM");
			substack_1P_counts_before_grouped_z_project_image_id = getImageID();
			run("Grouped Z Project...", "projection=[Sum Slices] group="+number_of_1P_FIT_files);
			showStatus("Preparing synthetic data sets: 1PPM");
			grouped_z_projection_image_id = getImageID();
			run("Make Montage...", "columns=1 rows="+number_of_subsets+" scale=1");
			showStatus("Preparing synthetic data sets: 1PPM");
			combined_subset_1P_count_image_id = getImageID();
			selectImage(substack_1P_counts_before_grouped_z_project_image_id);
			close();
			selectImage(grouped_z_projection_image_id);
			close();
			selectImage(image_of_1P_data_counts_stack_id);
			close();
	
			imageCalculator("Divide create 32-bit", combined_subset_1P_data_image_id, combined_subset_1P_count_image_id);
			showStatus("Preparing synthetic data sets: 1PPM");
			normalized_combined_subset_1P_data_image_id = getImageID();
			rename("normalized_1P_subsets");
	
			selectImage(combined_subset_1P_count_image_id);
			close();


			// 2P subsets
			showProgress(0.2);
			selectImage(image_of_2P_data_image_id);
			run("Duplicate...", "2P_counts");
			showStatus("Preparing synthetic data sets: 2PPM");
			rename("2P_counts");
			image_of_2P_data_counts_image_id = getImageID();
			run("Macro...", "code=v=(v==v)");            // replace 'NaN' pixels by '0', everything else by '1' 
			showStatus("Preparing synthetic data sets: 2PPM");
			selectImage(image_of_2P_data_image_id);
			run("Macro...", "code=[if(isNaN(v)) v=0]");  // replace 'NaN' by '0'
			showStatus("Preparing synthetic data sets: 2PPM");
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			showStatus("Preparing synthetic data sets: 2PPM");
			image_of_2P_data_stack_id = getImageID();
			selectImage(image_of_2P_data_counts_image_id);
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			showStatus("Preparing synthetic data sets: 2PPM");
			image_of_2P_data_counts_stack_id = getImageID();
			rename("2P counts, stack");
	
			list_of_2P_slice_numbers = toString(1+floor(random * number_of_2P_FIT_files));
			for(j = 1; j < number_of_subsets * number_of_2P_FIT_files; j++){
				showProgress(j, number_of_subsets * number_of_1P_FIT_files);
				list_of_2P_slice_numbers = list_of_2P_slice_numbers + ","+toString(1+floor(random * number_of_2P_FIT_files));
			}
	
	
			selectImage(image_of_2P_data_stack_id);
			run("Make Substack...", "  slices="+list_of_2P_slice_numbers);
			showStatus("Preparing synthetic data sets: 2PPM");
			substack_2P_data_before_grouped_z_project_image_id = getImageID();
			rename("2P data sets generated by Make Substack");
			run("Grouped Z Project...", "projection=[Sum Slices] group="+number_of_2P_FIT_files);
			showStatus("Preparing synthetic data sets: 2PPM");
			rename("2P data grouped Z projection");
			grouped_z_projection_image_id = getImageID();
			run("Make Montage...", "columns=1 rows="+number_of_subsets+" scale=1");
			showStatus("Preparing synthetic data sets: 2PPM");
			rename("2P data Make Montage");
			combined_subset_2P_data_image_id = getImageID();
			selectImage(substack_2P_data_before_grouped_z_project_image_id);
			close();
			selectImage(grouped_z_projection_image_id);
			close();
			selectImage(image_of_2P_data_stack_id);
			close();
	
	
			selectImage(image_of_2P_data_counts_stack_id);
			run("Make Substack...", "  slices="+list_of_2P_slice_numbers);
			showStatus("Preparing synthetic data sets: 2PPM");
			substack_2P_counts_before_grouped_z_project_image_id = getImageID();
			run("Grouped Z Project...", "projection=[Sum Slices] group="+number_of_2P_FIT_files);
			showStatus("Preparing synthetic data sets: 2PPM");
			grouped_z_projection_image_id = getImageID();
			run("Make Montage...", "columns=1 rows="+number_of_subsets+" scale=1");
			showStatus("Preparing synthetic data sets: 2PPM");
			combined_subset_2P_count_image_id = getImageID();
			selectImage(substack_2P_counts_before_grouped_z_project_image_id);
			close();
			selectImage(grouped_z_projection_image_id);
			close();
			selectImage(image_of_2P_data_counts_stack_id);
			close();
	
			imageCalculator("Divide create 32-bit", combined_subset_2P_data_image_id, combined_subset_2P_count_image_id);
			showStatus("Preparing synthetic data sets: 2PPM");
			normalized_combined_subset_2P_data_image_id = getImageID();
			rename("normalized_2P_subsets");
			showStatus("Preparing synthetic data sets: 2PPM");
	
			selectImage(combined_subset_2P_count_image_id);
			close();
	
	
			// combine 1P, 2P subsets
			run("Combine...", "stack1=normalized_1P_subsets stack2=normalized_2P_subsets");
			all_subset_image_id = getImageID();
			rename("Combined 1P, 2P subsets");
	
	
			selectImage(combined_subset_1P_data_image_id);
			close();
			selectImage(combined_subset_2P_data_image_id);
			close();
	
	
			selectImage(all_subset_image_id);
			run("Reslice [/]...", "output=1.000 start=Left rotate avoid");
			subsets_turned_into_stack_image_id = getImageID();		
			run("Size...", "width="+num_of_alpha0s+" height="+toString(number_of_subsets * num_of_sigmas)+" depth="+toString(2 * number_of_theta_bins)+" constrain interpolation=None");
	
						
			showStatus("Synthetic data sets made!");
	
			selectImage(all_subset_image_id);
			close();
	
			showProgress(0.3);
		
		
	
			// make a stack of 1P, 2P LD predictions (x, y, z) = (alpha0, sigma, theta), so value stored in
			// a pixel of coordinates x, y, z is the log2(r) for a particular combination of alpha0, sigma, theta
			
			// make the first 181 slices: of 1P LD for individual thetas 	
	
			showStatus("Preparing analysis: 1PPM");
			selectImage(B_1PPM_image_id);	
			substack_string = "  slices=1";
			for(j = 1; j < number_of_theta_bins; j++){substack_string = substack_string + ",1";}
			run("Make Substack...", substack_string);  // makes a stack of 181 B_1PPM.tif images
			rename("Prediction stack 1P, 2P");
			showStatus("Preparing analysis: 1PPM");
			prediction_stack_1P_2P_image_id = getImageID();
	
			// turns the stack of B_1PPM.tif images into a stack of 1P log2(r) predictions
			if (number_of_polarizations_1P > 1){run("Macro...", "code=[v = log(   (1 + v*cos(2*(z-1)/180*3.14159 - 3.14159))  /  (1 - v*cos(2*(z-1)/180*3.14159 - 3.14159))      )   /log(2)] stack");}
			else{run("Macro...", "code=[v = (1 + v*cos(2*(z-1)/180*3.14159 - 3.14159))] stack");}
			showStatus("Preparing analysis: 1PPM");
	
			showProgress(0.4);
	
	
			// make the second 181 slices: of 2P LD for individual thetas 
			selectImage(prediction_stack_1P_2P_image_id);
			setSlice(number_of_theta_bins);
			for(slice_number = 0; slice_number < number_of_theta_bins; slice_number++){
				run("Add Slice");
			}
			for(slice_number = 0; slice_number < number_of_theta_bins; slice_number++){
				showStatus("Preparing analysis: 2PPM");
				showProgress(slice_number, number_of_theta_bins);
	
				// make an image of 2P log2(r) prediction
				theta = slice_number/(number_of_theta_bins-1)*PI;
				if (number_of_polarizations_2P > 1){expectation_expression_2P = "log( (1 + A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ") ) / (1 - A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ") ) ) / log(2)";	}
				else{expectation_expression_2P = "1 + A * cos(" + toString(2*theta - 2*PI/2) + ") + B * cos(" + toString(4*theta - 4*PI/2) + ")";}
				run("Image Expression Parser (Macro)", "expression=["+ expectation_expression_2P +"] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	
				showStatus("Preparing analysis: 2PPM");
		
				expectation_image_id = getImageID();
				rename("expectation image");
		
				run("Select All");
				run("Copy");
				showStatus("Preparing analysis: 2PPM");
		
				selectImage(prediction_stack_1P_2P_image_id);
				setSlice(slice_number + number_of_theta_bins + 1);
				run("Paste");
				showStatus("Preparing analysis: 2PPM");
		
				selectImage(expectation_image_id);
				close();
		
			}
	
			showProgress(0.5);
	
			selectImage(prediction_stack_1P_2P_image_id);
			
			prediction_stack_height = getHeight();
			orig_name = getTitle();
			suffix = toString(abs(prediction_stack_1P_2P_image_id));
			suffix = substring(suffix, lengthOf(suffix)-4);
			temp_name = "temp"+suffix;
			rename(temp_name);
	
			showStatus("Analyzing combined 1P, 2P data");
	
			prediction_stack_height = getHeight();
			while(prediction_stack_height < num_of_sigmas * number_of_subsets){
				run("Combine...", "stack1="+temp_name+" stack2="+temp_name+" combine");
				showStatus("Analyzing combined 1P, 2P data");
				prediction_stack_1P_2P_image_id = getImageID();
				prediction_stack_height = getHeight();		
				rename(temp_name);
			}
			rename(orig_name);
			run("Canvas Size...", "width="+num_of_alpha0s+" height="+toString(num_of_sigmas * number_of_subsets)+" position=Top-Center");
			showStatus("Analyzing combined 1P, 2P data");
	
			showProgress(0.6);
		
			// make an image that will store the alpha0, sigma combinations for the individual data sets
			newImage("alpha0_sigma_cloud_for_subsets", "16-bit black", num_of_alpha0s, num_of_sigmas, 1);
			alpha0_sigma_cloud_image_id = getImageID();
		
			imageCalculator("Subtract 32-bit stack", prediction_stack_1P_2P_image_id, subsets_turned_into_stack_image_id);
			selectImage(subsets_turned_into_stack_image_id);
			close();
			selectImage(prediction_stack_1P_2P_image_id);
			run("Square", "stack");
			run("Macro...", "code=[if(isNaN(v)) v=0] stack");  // replace 'NaN' by '0'

			showProgress(0.7);

			showStatus("Analyzing combined 1P, 2P data");
			run("Z Project...", "projection=[Sum Slices]");
			showStatus("Analyzing combined 1P, 2P data");
			subset_SSD_image = getImageID();
			selectImage(prediction_stack_1P_2P_image_id);
			close();
	
			selectImage(subset_SSD_image);
			getRawStatistics(nPixels, mean, min, max);
			run("Macro...", "code=[if(isNaN(v)) v="+max+"] stack");  // replace 'NaN' by max SSD value
			getRawStatistics(nPixels, mean, min, max);
			showStatus("Analyzing combined 1P, 2P data");
	
		
			for(data_set_number = 0; data_set_number < number_of_subsets; data_set_number++){
				showProgress(data_set_number, number_of_subsets);
				selectImage(subset_SSD_image);
				makeRectangle(0, data_set_number * num_of_sigmas, num_of_alpha0s, num_of_sigmas);
				run("Duplicate...", " ");
				run("Find Maxima...", "noise="+max+" output=[Point Selection] light");
				getSelectionBounds(alpha0, sigma_row_number, w, h);
				close();
				selectImage(alpha0_sigma_cloud_image_id);
				setPixel(alpha0, sigma_row_number, getPixel(alpha0, sigma_row_number)+1);
	
				list_of_1P_2P_alpha0s[data_set_number] = alpha0;
				list_of_1P_2P_sigmas[data_set_number] = num_of_sigmas - sigma_row_number - 1;
			
			}		
	
			selectImage(subset_SSD_image);
			close();
				
			showStatus("Analyzing combined 1P, 2P data");
		
			Array.getStatistics(list_of_1P_2P_alpha0s, alpha0_min, alpha0_max, alpha0_mean, alpha0_stdev);
			Array.getStatistics(list_of_1P_2P_sigmas, sigma_min, sigma_max, sigma_mean, sigma_stdev);

			statistics_1P_2P_text = "Combined 1P, 2P statistics (bootstrapping):\n";
			statistics_1P_2P_text = statistics_1P_2P_text + alpha0char+" = "+toString(round(alpha0_mean*10)/10)+degreechar+plusminuschar+toString(round((2*alpha0_stdev/sqrt(number_of_subsets))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+alpha0_min+degreechar+", max = "+alpha0_max+degreechar+", stdev = "+alpha0_stdev+degreechar+")\n";
			statistics_1P_2P_text = statistics_1P_2P_text + sigmachar+" = "+toString(round(sigma_mean*10)/10)+degreechar+plusminuschar+toString(round((2*sigma_stdev/sqrt(number_of_subsets))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+sigma_min+degreechar+", max = "+sigma_max+degreechar+", stdev = "+sigma_stdev+degreechar+")\n";
			statistics_1P_2P_text_few_special_chars = "Combined 1P, 2P statistics (bootstrapping):\n";
			statistics_1P_2P_text_few_special_chars = statistics_1P_2P_text_few_special_chars + "alpha0 = "+toString(round(alpha0_mean*10)/10)+degreechar+plusminuschar+toString(round((2*alpha0_stdev/sqrt(number_of_subsets))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+alpha0_min+degreechar+", max = "+alpha0_max+degreechar+", stdev = "+alpha0_stdev+degreechar+")\n";
			statistics_1P_2P_text_few_special_chars = statistics_1P_2P_text_few_special_chars + "sigma = "+toString(round(sigma_mean*10)/10)+degreechar+plusminuschar+toString(round((2*sigma_stdev/sqrt(number_of_subsets))*10)/10)+degreechar+" (mean "+plusminuschar+" 2*stderr (95% CI); min = "+sigma_min+degreechar+", max = "+sigma_max+degreechar+", stdev = "+sigma_stdev+degreechar+")\n";
			print(statistics_1P_2P_text);

		


			selectImage(alpha0_sigma_cloud_image_id);
			run("Select All");
			run("Copy");
			showStatus("Analyzing combined 1P, 2P data");
		
			selectImage(hyperstack_image_id);	// created on lines #3891
			setSlice(12);						// 1P, 2P alpha0, sigma
			run("Paste");
			setColor(255,100,255);  // magenta
			Overlay.show;
		
			selectImage(hyperstack_image_id);	// created on lines #3891
			setSlice(1);						// 
			run("Select None");

			hyperstack_info_text = hyperstack_info_text + "\n" + statistics_1P_text_few_special_chars + "\n" + statistics_2P_text_few_special_chars + "\n" + statistics_1P_2P_text_few_special_chars+"\n\n";
			setMetadata("Info", hyperstack_info_text);
			
			save(directory_2P+File.separator+hyperstack_name+".tif");
			save(directory_1P+File.separator+hyperstack_name+".tif");

			log_file_name = "1P_2P_combined_fit_"+subdirectory_1P_name+"_"+subdirectory_2P_name+"_"+year+"_"+month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second+".log";

			output_log = File.open(directory_1P+File.separator+log_file_name);
			print(output_log,hyperstack_info_text);
			File.close(output_log);

			if(directory_1P != directory_2P){
				output_log = File.open(directory_2P+File.separator+log_file_name);
				print(output_log,hyperstack_info_text);
				File.close(output_log);			
			}


			selectImage(alpha0_sigma_cloud_image_id);
			close();
			selectImage(image_of_1P_data_counts_image_id);
			close();
			selectImage(image_of_2P_data_counts_image_id);
			close();

		
		}

	}

	if(number_of_1P_FIT_files > 0){
		selectImage(image_of_1P_data_image_id);
		close();
	}
	if(number_of_2P_FIT_files > 0){
		selectImage(image_of_2P_data_image_id);
		close();
	}

	selectImage(B_1PPM_image_id);
	close();
	selectImage(B_2PPM_image_id);
	close();
	selectImage(C_2PPM_image_id);
	close();


	selectImage(hyperstack_image_id);	

	setBatchMode("exit and display");
	showStatus("Done!");


}







macro "Fit by two Gaussian distributions [t]"{
	setBatchMode(true);

	predicted_B_2PPM = 0;
	predicted_C_2PPM = 0;
	predicted_b_1PPM = 0;
	predicted_rmax_2PPM = 0;
	predicted_rmax_1PPM = 0;


	rmax_1PPM = 0;
	log2rmax_1PPM = 0;
	rmax_2PPM = 0;
	log2rmax_2PPM = 0;




	Dialog.create("Fitting by two gaussian distributions");
	Dialog.addMessage("Enter the parameters provided\nby the 'g' or 'c' macros. Entering\na 1P fit parameter B is optional.\n ");
	Dialog.addNumber("1P fit parameter B:", B_1P_fit);
	Dialog.addNumber("2P fit parameter B:", B_2P_fit);
	Dialog.addNumber("2P fit parameter C:", C_2P_fit);
	Dialog.addMessage("Search intervals (major distribution):");

	Dialog.addNumber(alpha0char+", from ", min_alpha0_major, 0, 2, degreechar);
	Dialog.addToSameRow();
	Dialog.addNumber("to:", max_alpha0_major, 0, 2,degreechar);
	Dialog.addNumber(sigmachar+", from ", min_sigma_major, 0, 2, degreechar);
	Dialog.addToSameRow();
	Dialog.addNumber("to:", max_sigma_major, 0, 2,degreechar);
	Dialog.addMessage("Search intervals (minor distribution):");
	Dialog.addNumber(alpha0char+", from ", min_alpha0_minor, 0, 2, degreechar);
	Dialog.addToSameRow();
	Dialog.addNumber("to:", max_alpha0_minor, 0, 2, degreechar);
	Dialog.addNumber(sigmachar+", from ", min_sigma_minor, 0, 2, degreechar);
	Dialog.addToSameRow();
	Dialog.addNumber("to:", max_sigma_minor, 0, 2, degreechar);

	Dialog.addMessage("Sampling intervals:");
	Dialog.addNumber(alpha0char+" (in "+degreechar+"):",alpha0_interval);
	Dialog.addNumber(sigmachar+" (in "+degreechar+"):",sigma_interval);
	Dialog.addNumber("Composition (in %):",composition_percentage_interval);
	
	Dialog.show();
	B_1P_fit = Dialog.getNumber();
	B_2P_fit = Dialog.getNumber();
	C_2P_fit = Dialog.getNumber();
	
	min_alpha0_major = Dialog.getNumber();
	max_alpha0_major = Dialog.getNumber();
	min_sigma_major = Dialog.getNumber();
	max_sigma_major = Dialog.getNumber();
	
	min_alpha0_minor = Dialog.getNumber();
	max_alpha0_minor = Dialog.getNumber();
	min_sigma_minor = Dialog.getNumber();
	max_sigma_minor = Dialog.getNumber();

	
	alpha0_interval = Dialog.getNumber();
	sigma_interval = Dialog.getNumber();
	composition_percentage_interval = Dialog.getNumber();
		
	if(B_1P_fit != 0){
		rmax_1PPM = (1 + B_1P_fit)/(1 - B_1P_fit);
		log2rmax_1PPM = log(rmax_1PPM)/log(2);
	}

	rmax_2PPM = (1 + B_2P_fit + C_2P_fit)/(1 - B_2P_fit + C_2P_fit);
	log2rmax_2PPM = log(rmax_2PPM)/log(2);
	
	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();
	run("Duplicate...", "title=B_C_b_stack");
	B_C_b_stack_image_id = getImageID();
	
	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();
	run("Select All");
	run("Copy");

	selectImage(B_C_b_stack_image_id);
	run("Add Slice");
	setSlice(2);
	run("Paste");
	
	
	B_1PPM_image_id = 0;
	if(rmax_1PPM > 0){
		run("Image Expression Parser (Macro)", "expression=[(4*B - 7*A)/(-10 + A + 2*B)] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		rename("1PPM_B.tif");
		B_1PPM_image_id = getImageID();
		run("Select All");
		run("Copy");

		selectImage(B_C_b_stack_image_id);
		run("Add Slice");
		setSlice(3);
		run("Paste");
	}


	
	selectImage(B_C_b_stack_image_id);	
	run("Specify...", "width="+toString(max_alpha0_major - min_alpha0_major + 1)+" height="+toString(max_sigma_major - min_sigma_major +1)+" x="+min_alpha0_major+" y="+toString(90-max_sigma_major));	
	run("Duplicate...", "title=B_C_b_major duplicate");
	B_C_b_major_stack_id = getImageID();
	run("Size...", "width="+toString(floor((getWidth-1)/alpha0_interval)+1)+" height="+toString(floor((getHeight-1)/sigma_interval)+1)+" interpolation=Bicubic");
	num_of_alpha0s_major = getWidth;
	num_of_sigmas_major  = getHeight;

	selectImage(B_C_b_stack_image_id);
	run("Specify...", "width="+toString(max_alpha0_minor - min_alpha0_minor + 1)+" height="+toString(max_sigma_minor - min_sigma_minor +1)+" x="+min_alpha0_minor+" y="+toString(90-max_sigma_minor));	
	run("Duplicate...", "title=B_C_b_minor duplicate");
	B_C_b_minor_stack_id = getImageID();
	run("Size...", "width="+toString(floor((getWidth-1)/alpha0_interval)+1)+" height="+toString(floor((getHeight-1)/sigma_interval)+1)+" interpolation=Bicubic");
	num_of_alpha0s_minor = getWidth;
	num_of_sigmas_minor  = getHeight;



	showStatus("Preparations for fitting...");


	B_C_b_tiled_stack_id = t_make_a_tiled_stack(B_C_b_major_stack_id, num_of_alpha0s_minor, num_of_sigmas_minor);   //major fraction
	
	B_C_b_stretched_stack_id = t_make_a_stretched_stack(B_C_b_minor_stack_id, num_of_alpha0s_major, num_of_sigmas_major);  //minor fraction

	selectImage(B_2PPM_image_id);
	close();
	selectImage(C_2PPM_image_id);
	close();
	if(rmax_1PPM > 0){
		selectImage(B_1PPM_image_id);
		close();
	}
	selectImage(B_C_b_stack_image_id);
	close();


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	print(" ");
	print("Values entered by user:");
	print("B(1P): "+B_1P_fit+"   B(2P): "+B_2P_fit+"   C(2P): "+C_2P_fit+"   [rmax(1P): "+toString((1+B_1P_fit)/(1-B_1P_fit))+", rmax(2P): "+toString((1+B_2P_fit+C_2P_fit)/(1-B_2P_fit+C_2P_fit))+"]   "+alpha0char+" interval: "+alpha0_interval+degreechar+"   "+sigmachar+" interval: "+sigma_interval+degreechar);
	print(" ");
	print("  #1%    "+alpha0char+", "+sigmachar+"           #2%   "+alpha0char+", "+sigmachar+"           RMSD           entropy   rmax(1P) rmax(2P)	B(1P)         B(2P)         C(2P) ");


	for(minor_fraction_percentage = 0; minor_fraction_percentage < 0.501; minor_fraction_percentage = minor_fraction_percentage + composition_percentage_interval/100){
		
		major_fraction_percentage = 1 - minor_fraction_percentage;
		showStatus("Analyzing "+toString(minor_fraction_percentage*100)+":"+toString(major_fraction_percentage*100));
		showProgress(minor_fraction_percentage, 0.5);

		sum_of_B_C_b_major_minor_values_stack_id = t_make_a_sum_of_B_C_b_major_minor_values_stack(B_C_b_tiled_stack_id, B_C_b_stretched_stack_id, minor_fraction_percentage);
		rename("B_C_b_values_"+toString(100 - 100*minor_fraction_percentage)+":"+toString(100*minor_fraction_percentage)+"_mix");


		if(rmax_1PPM > 0){    	//use rmax(1P), rmax(2P) to find the smallest RMSD
			rmax_1PPM_image_id = t_make_rmax_1PPM_image(sum_of_B_C_b_major_minor_values_stack_id);

			selectImage(rmax_1PPM_image_id);
			run("Log");																					// Ln(Rmax_1P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_1P)/Ln(2) = Log2(Rmax_1P)
			run("Subtract...", "value="+log2rmax_1PPM);
			run("Square");																				// square deviation from log2 of 1P_rmax calculated from entered B2P, C2P values

			rmax_2PPM_image_id = t_make_rmax_2PPM_image(sum_of_B_C_b_major_minor_values_stack_id);

			selectImage(rmax_2PPM_image_id);
			run("Log");																					// Ln(Rmax_2P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_2P)/Ln(2) = Log2(Rmax_2P)
			run("Subtract...", "value="+log2rmax_2PPM);
			run("Square");																				// square deviation from log2 of 2P_rmax calculated from entered B2P, C2P values

			imageCalculator("Add 32-bit", rmax_2PPM_image_id, rmax_1PPM_image_id);				
			selectImage(rmax_1PPM_image_id);
			close();

			selectImage(rmax_2PPM_image_id);
			rename("SSD_1P_2P_"+toString(100 - 100*minor_fraction_percentage)+":"+toString(100*minor_fraction_percentage));
			SSD_image_id = getImageID();
		

//--------------------theta = 15 degrees----------------------//

			log2r_1PPM_15_deg = log( (1+B_1P_fit*cos(2*15/180*PI))/(1-B_1P_fit*cos(2*15/180*PI)) )/log(2);
			log2r_2PPM_15_deg = log( (1+B_2P_fit*cos(2*15/180*PI)+C_2P_fit*cos(4*15/180*PI))/(1-B_2P_fit*cos(2*15/180*PI)+C_2P_fit*cos(4*15/180*PI)) )/log(2);

			rmax_1PPM_image_id_15 = t_make_rmax_1PPM_image_15(sum_of_B_C_b_major_minor_values_stack_id);
			selectImage(rmax_1PPM_image_id_15);
			run("Log");																					// Ln(Rmax_1P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_1P)/Ln(2) = Log2(Rmax_1P)
			run("Subtract...", "value="+log2r_1PPM_15_deg);
			run("Square");																				// square deviation from log2 of 1P_rmax calculated from entered B2P, C2P values

			rmax_2PPM_image_id_15 = t_make_rmax_2PPM_image_15(sum_of_B_C_b_major_minor_values_stack_id);
			selectImage(rmax_2PPM_image_id_15);
			run("Log");																					// Ln(Rmax_2P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_2P)/Ln(2) = Log2(Rmax_2P)
			run("Subtract...", "value="+log2r_2PPM_15_deg);
			run("Square");																				// square deviation from log2 of 2P_rmax calculated from entered B2P, C2P values

			imageCalculator("Add 32-bit", SSD_image_id, rmax_1PPM_image_id_15);				
			selectImage(rmax_1PPM_image_id_15);
			close();

			imageCalculator("Add 32-bit", SSD_image_id, rmax_2PPM_image_id_15);				
			selectImage(rmax_2PPM_image_id_15);
			close();


//--------------------theta = 15 degrees----------------------//


//--------------------theta = 30 degrees----------------------//

			log2r_1PPM_30_deg = log( (1+B_1P_fit*cos(2*30/180*PI))/(1-B_1P_fit*cos(2*30/180*PI)) )/log(2);
			log2r_2PPM_30_deg = log( (1+B_2P_fit*cos(2*30/180*PI)+C_2P_fit*cos(4*30/180*PI))/(1-B_2P_fit*cos(2*30/180*PI)+C_2P_fit*cos(4*30/180*PI)) )/log(2);

			rmax_1PPM_image_id_30 = t_make_rmax_1PPM_image_30(sum_of_B_C_b_major_minor_values_stack_id);
			selectImage(rmax_1PPM_image_id_30);
			run("Log");																					// Ln(Rmax_1P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_1P)/Ln(2) = Log2(Rmax_1P)
			run("Subtract...", "value="+log2r_1PPM_30_deg);
			run("Square");																				// square deviation from log2 of 1P_rmax calculated from entered B2P, C2P values

			rmax_2PPM_image_id_30 = t_make_rmax_2PPM_image_30(sum_of_B_C_b_major_minor_values_stack_id);
			selectImage(rmax_2PPM_image_id_30);
			run("Log");																					// Ln(Rmax_2P)
			run("Divide...", "value=0.693147");															// Ln(Rmax_2P)/Ln(2) = Log2(Rmax_2P)
			run("Subtract...", "value="+log2r_2PPM_30_deg);
			run("Square");																				// square deviation from log2 of 2P_rmax calculated from entered B2P, C2P values

			imageCalculator("Add 32-bit", SSD_image_id, rmax_1PPM_image_id_30);				
			selectImage(rmax_1PPM_image_id_30);
			close();

			imageCalculator("Add 32-bit", SSD_image_id, rmax_2PPM_image_id_30);				
			selectImage(rmax_2PPM_image_id_30);
			close();


//--------------------theta = 30 degrees----------------------//

			selectImage(SSD_image_id);
			run("Divide...", "value=8");   // I am using 1P, 2P values at 0o, 15o, 30o, but the values at 45o are 0, so I can say I add them too, so I should divide by 8




		
		}
		else{   			//use B2P, C2P to find the smallest RMSD
			setSlice(1);
			run("Subtract...", "value="+B_2P_fit+" slice");
			setSlice(2);
			run("Subtract...", "value="+C_2P_fit+" slice");
			run("Multiply...", "value=4 slice");				//account for the 4x smaller range of C2P (compared to B2P)
			run("Square","stack");
			run("Z Project...", "projection=[Sum Slices]");
			rename("SSD_2P_"+toString(100 - 100*minor_fraction_percentage)+":"+toString(100*minor_fraction_percentage));
			SSD_image_id = getImageID();

		}

		selectImage(sum_of_B_C_b_major_minor_values_stack_id);
		close();

		selectImage(SSD_image_id);
		run("Multiply...", "value=1000000");
		getRawStatistics(nPixels, mean, min, max);
		for(j=0; j<getWidth(); j++){
			for(k=0; k<getHeight(); k++){
				if(getPixel(j, k) == min){
					min_SSD_x_coord = j;
					min_SSD_y_coord = k;
				}
			}
		}

		showStatus("Analyzing "+toString(minor_fraction_percentage*100)+":"+toString(major_fraction_percentage*100));
		showProgress(minor_fraction_percentage, 0.5);
		

		alpha0_major = min_alpha0_major + (min_SSD_x_coord % num_of_alpha0s_major) * alpha0_interval;
		sigma_major  = min_sigma_major  + (num_of_sigmas_major - 1 - min_SSD_y_coord % num_of_sigmas_major) * sigma_interval;
		alpha0_minor = min_alpha0_minor + floor(min_SSD_x_coord / num_of_alpha0s_major * alpha0_interval);
		sigma_minor  = min_sigma_minor  + (num_of_sigmas_minor - 1 - floor(min_SSD_y_coord/num_of_sigmas_major)) * sigma_interval;

		
		min_SSD = getPixel(min_SSD_x_coord, min_SSD_y_coord)/1000000;
		min_RMSD = sqrt(min_SSD);
		entropy = (1 + log(2*3.14159))/2 + minor_fraction_percentage*log(sigma_minor) + major_fraction_percentage*log(sigma_major);


		selectImage(B_C_b_tiled_stack_id);
		setZCoordinate(0);
		major_B_2PPM = getPixel(min_SSD_x_coord, min_SSD_y_coord);
		setZCoordinate(1);
		major_C_2PPM = getPixel(min_SSD_x_coord, min_SSD_y_coord);
		
		selectImage(B_C_b_stretched_stack_id);
		setZCoordinate(0);
		minor_B_2PPM = getPixel(min_SSD_x_coord, min_SSD_y_coord);
		setZCoordinate(1);
		minor_C_2PPM = getPixel(min_SSD_x_coord, min_SSD_y_coord);
		
		predicted_B_2PPM = major_fraction_percentage * major_B_2PPM + minor_fraction_percentage * minor_B_2PPM;
		predicted_C_2PPM = major_fraction_percentage * major_C_2PPM + minor_fraction_percentage * minor_C_2PPM;
		predicted_b_1PPM = major_fraction_percentage * (4 * major_C_2PPM - 7 * major_B_2PPM)/(-10 + major_B_2PPM + 2 * major_C_2PPM) + minor_fraction_percentage * (4 * minor_C_2PPM - 7 * minor_B_2PPM)/(-10 + minor_B_2PPM + 2 * minor_C_2PPM);
		
		predicted_rmax_1PPM = (1 + predicted_b_1PPM)/(1 - predicted_b_1PPM);
		predicted_rmax_2PPM = (1 + predicted_B_2PPM +  predicted_C_2PPM)/(1 - predicted_B_2PPM +  predicted_C_2PPM);

		formatted_output = t_format_output(major_fraction_percentage, alpha0_major, sigma_major, minor_fraction_percentage, alpha0_minor, sigma_minor, min_RMSD, entropy, predicted_rmax_1PPM, predicted_rmax_2PPM, predicted_b_1PPM, predicted_B_2PPM, predicted_C_2PPM);
		print(formatted_output);

		selectImage(SSD_image_id);
		close();
				
		
	}


	selectImage(B_C_b_minor_stack_id);
	close();
	selectImage(B_C_b_major_stack_id);
	close();
	
	selectImage(B_C_b_tiled_stack_id);
	close();
	selectImage(B_C_b_stretched_stack_id);
	close();
		
	
	showStatus("Done! Check the Log window for results.");

	
	setBatchMode("exit and display");

}





macro "Determine alpha0, sigma from rmax, log2(rmax) [d]"{
	setBatchMode(true);
	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();
	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();

	run("Image Expression Parser (Macro)", "expression=[(4*B - 7*A)/(-10 + A + 2*B)] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("1PPM_B.tif");
	B_1PPM_image_id = getImageID();

	run("Image Expression Parser (Macro)", "expression=[log((1 + A)/(1 - A))/log(2)] a=1PPM_B.tif b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("1PPM_log2rmax.tif");
	log2rmax_1PPM_image_id = getImageID();
	setPixel(0,90,1000);
	setPixel(90,90,-1000);

	run("Image Expression Parser (Macro)", "expression=[log((1 + A + B)/(1 - A + B))/log(2)] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("2PPM_log2rmax.tif");
	log2rmax_2PPM_image_id = getImageID();
	setPixel(0,90,1000);
	setPixel(90,90,-1000);

	alpha0_sigma_1P_image_id = 0;
	alpha0_sigma_2P_image_id = 0;
	rgb_alpha0_sigma_1P_2P_image_id = 0;
	
	parameters_changed = true;
	while(parameters_changed == true){
		Dialog.create("Finding "+alpha0char+", "+sigmachar+" from rmax or log2(rmax) values");
		Dialog.addMessage("Enter values for 1PPM, 2PPM, or both:");


		Dialog.addChoice("1PPM: ", newArray("log2(rmax)", "rmax"), property_kind_1P);
		Dialog.addToSameRow();
		Dialog.addNumber("=", property_mean_1P);
		Dialog.addToSameRow();
		Dialog.addNumber(plusminuschar, interval_value_1P);
		Dialog.addToSameRow();
		Dialog.addChoice("", newArray("%", ""), interval_kind_1P);

		Dialog.addChoice("2PPM: ", newArray("log2(rmax)", "rmax"), property_kind_2P);
		Dialog.addToSameRow();
		Dialog.addNumber("=", property_mean_2P);
		Dialog.addToSameRow();
		Dialog.addNumber(plusminuschar, interval_value_2P);
		Dialog.addToSameRow();
		Dialog.addChoice("", newArray("%", ""), interval_kind_2P);


		Dialog.show();

		property_kind_1P_new = Dialog.getChoice();		
		property_mean_1P_new = Dialog.getNumber();		
		interval_value_1P_new = Dialog.getNumber();
		interval_kind_1P_new = Dialog.getChoice();

		property_kind_2P_new = Dialog.getChoice();		
		property_mean_2P_new = Dialog.getNumber();		
		interval_value_2P_new = Dialog.getNumber();
		interval_kind_2P_new = Dialog.getChoice();



		if(alpha0_sigma_1P_image_id < 0 && property_kind_1P_new == property_kind_1P && property_kind_2P_new == property_kind_2P && property_mean_1P_new == property_mean_1P && property_mean_2P_new == property_mean_2P && interval_value_1P_new == interval_value_1P && interval_value_2P_new == interval_value_2P && interval_kind_1P_new == interval_kind_1P && interval_kind_2P_new == interval_kind_2P){		
			parameters_changed = false;
		}
		else{
			property_kind_1P = property_kind_1P_new;
			property_kind_2P = property_kind_2P_new;
			property_mean_1P = property_mean_1P_new;
			property_mean_2P = property_mean_2P_new;
			interval_value_1P = interval_value_1P_new;
			interval_value_2P = interval_value_2P_new;
			interval_kind_1P = interval_kind_1P_new;
			interval_kind_2P = interval_kind_2P_new;
			
			if(property_kind_1P == "log2(rmax)"){log2rmax_1P = property_mean_1P;}
			else{log2rmax_1P = log(property_mean_1P)/log(2);}

			if(property_kind_2P == "log2(rmax)"){log2rmax_2P = property_mean_2P;}
			else{log2rmax_2P = log(property_mean_2P)/log(2);}

			if(interval_kind_1P == "%"){interval_1P = abs(log2rmax_1P*interval_value_1P/100);}
			else{interval_1P = interval_value_1P;}

			if(interval_kind_2P == "%"){interval_2P = abs(log2rmax_2P*interval_value_2P/100);}
			else{interval_2P = interval_value_2P;}

		

			if(rgb_alpha0_sigma_1P_2P_image_id < 0){
				selectImage(rgb_alpha0_sigma_1P_2P_image_id );
				close();
			}
	
			if(log2rmax_1P != 0){
				selectImage(log2rmax_1PPM_image_id);
				run("Scale...", "x=5 y=5 width=451 height=451 interpolation=Bicubic average create");
				rename("alpha0_sigma_1P");
				run("Subtract...", "value="+log2rmax_1P);
				run("Square");
				run("Square Root");			

				alpha0_sigma_1P_image_id = getImageID();			
				setMinAndMax(0, abs(interval_1P));
				run("16-bit");
				run("8-bit");
			}
			else{
				selectImage(log2rmax_1PPM_image_id);
				run("Scale...", "x=5 y=5 width=451 height=451 interpolation=Bicubic average create");
				rename("alpha0_sigma_1P");
				alpha0_sigma_1P_image_id = getImageID();			
				run("16-bit");
				run("8-bit");
				run("Set...", "value=255");
			}

			

			if(log2rmax_2P != 0){
				selectImage(log2rmax_2PPM_image_id);
				run("Scale...", "x=5 y=5 width=451 height=451 interpolation=Bicubic average create");
				rename("alpha0_sigma_2P");
				run("Subtract...", "value="+log2rmax_2P);
				run("Square");
				run("Square Root");
				alpha0_sigma_2P_image_id = getImageID();
				
				setMinAndMax(0, abs(interval_2P));
				run("16-bit");
				run("8-bit");
			}
			else{
				selectImage(log2rmax_2PPM_image_id);
				run("Scale...", "x=5 y=5 width=451 height=451 interpolation=Bicubic average create");
				rename("alpha0_sigma_2P");
				alpha0_sigma_2P_image_id = getImageID();			
				run("16-bit");
				run("8-bit");
				run("Set...", "value=255");
				
			}

			newImage("test", "RGB white", 451, 451, 1);
			rgb_alpha0_sigma_1P_2P_image_id = getImageID();
			setLineWidth(1);
			for(row=0; row < 455; row++){
				for(column = 0; column < 455; column++){
					selectImage(alpha0_sigma_1P_image_id);
					SSD_1P_pixel_value = getPixel(column, row);
					selectImage(alpha0_sigma_2P_image_id);
					SSD_2P_pixel_value = getPixel(column, row);
					selectImage(rgb_alpha0_sigma_1P_2P_image_id);
					if(SSD_1P_pixel_value == 255 && SSD_2P_pixel_value == 255){}
					if(SSD_1P_pixel_value < 255 && SSD_2P_pixel_value == 255){
						// I need to make blue
						setColor(SSD_1P_pixel_value, SSD_1P_pixel_value, 255);
						drawLine(column, row, column, row);
					}
					if(SSD_1P_pixel_value == 255 && SSD_2P_pixel_value < 255){
						// I need to make red
						setColor(255, SSD_2P_pixel_value, SSD_2P_pixel_value);
						drawLine(column, row, column, row);
					}
					if(SSD_1P_pixel_value < 255 && SSD_2P_pixel_value < 255){
						// I need to make magenta
						setColor( (255+SSD_1P_pixel_value)/2, 255-(255-SSD_1P_pixel_value) - (255-SSD_2P_pixel_value), (255+SSD_2P_pixel_value)/2);
						drawLine(column, row, column, row);
					}
				}
			}


			rename("alpha0, sigma pairs compatible with data");
			run("Properties...", "channels=1 slices=1 frames=1 unit=° pixel_width=0.2000000 pixel_height=-0.2000000 voxel_depth=1.0000000 origin=0,454");
			rename("1P log2rmax="+log2rmax_1P+plusminuschar+interval_value_1P+interval_kind_1P+", 2P log2rmax="+log2rmax_2P+plusminuschar+interval_value_2P+interval_kind_2P);

			setFont("SansSerif", 16);
			setColor(0,0,255);  //  blue
			Overlay.drawString("1PPM", 10, 30);
			setColor(255,0,0);  //  red
			Overlay.drawString("2PPM", 10, 50);
			Overlay.show;	
			run("Set Label...", "label=[x="+alpha0char+", y="+sigmachar+"]");

			
			setBatchMode("show");
			rgb_alpha0_sigma_1P_2P_image_id = getImageID();


			selectImage(alpha0_sigma_1P_image_id);
			close();
			selectImage(alpha0_sigma_2P_image_id);
			close();



		}
		
	}

	selectImage(B_1PPM_image_id);
	close();
	selectImage(B_2PPM_image_id);
	close();
	selectImage(C_2PPM_image_id);
	close();
	selectImage(log2rmax_1PPM_image_id);
	close();
	selectImage(log2rmax_2PPM_image_id);
	close();
	setBatchMode("exit and display");
	


}





macro "Predict LD [p]"{
	setBatchMode(true);


	initial_alpha0 = 0;
	initial_sigma = 0;
	initial_num_of_polarizations = "2";
	initial_modify_plot = false;
	
	macro_path = getDirectory("macros");
	if(File.exists(macro_path+"LD_prediction_config.txt")){
		config_text = File.openAsString(macro_path+"LD_prediction_config.txt");
		config_rows = split(config_text,"\n");
	
		for(i=0; i<config_rows.length; i++){
			key_value_pair = split(config_rows[i]);
			if(key_value_pair[0] == "last_alpha0"){
				initial_alpha0 = key_value_pair[1];
			}
			if(key_value_pair[0] == "last_sigma"){	
				initial_sigma = key_value_pair[1];
			}
			if(key_value_pair[0] == "last_num_of_polarizations"){	
				initial_num_of_polarizations = key_value_pair[1];
			}
			if(key_value_pair[0] == "last_modify_plot"){
				initial_modify_plot = key_value_pair[1];
			}
		}
	}


	Dialog.create("LD prediction");
	Dialog.addMessage("Enter the parameters ("+alpha0char+", "+sigmachar+")\nof the orientational distribution \nof fluorophores:\n ");
	Dialog.addNumber("            "+alpha0char+" (in "+degreechar+"):", initial_alpha0);
	Dialog.addNumber("            "+sigmachar+" (in "+degreechar+"):", initial_sigma);
	Dialog.addRadioButtonGroup("Number of polarizations:", newArray("1","2"), 2, 1, initial_num_of_polarizations)
	Dialog.addCheckbox("Modify an existing plot", initial_modify_plot);
	Dialog.show()

	alpha0 = Dialog.getNumber();
	sigma = Dialog.getNumber();
	number_of_polarizations = Dialog.getRadioButton();
	use_existing_plot = Dialog.getCheckbox();
	number_of_polarizations = parseInt(number_of_polarizations);

	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();
	B_2PPM = getPixel(alpha0, 90-sigma);

	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();
	C_2PPM = getPixel(alpha0, 90-sigma);

	run("Image Expression Parser (Macro)", "expression=[(4*B - 7*A)/(-10 + A + 2*B)] a=2PPM_B.tif b=2PPM_C.tif c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("1PPM_B.tif");
	B_1PPM_image_id = getImageID();
	B_1PPM = getPixel(alpha0, 90-sigma);

	selectImage(B_1PPM_image_id);
	close();
	selectImage(B_2PPM_image_id);
	close();
	selectImage(C_2PPM_image_id);
	close();
	

	number_of_samples = 181; //-90o to 90o
	fit_x_values_deg = newArray(number_of_samples);
	fit_1P_LD_values = newArray(number_of_samples);
	fit_2P_LD_values = newArray(number_of_samples);

	if(number_of_polarizations == 2){
		for(i = 0; i < number_of_samples; i++){
			theta_deg = 180/(number_of_samples-1)*i - 90;
			fit_x_values_deg[i] = theta_deg;
			fit_1P_LD_values[i] = log((1 + B_1PPM * cos(2*theta_deg/180*PI))/(1 - B_1PPM * cos(2*theta_deg/180*PI)))/log(2);
			fit_2P_LD_values[i] = log((1 + B_2PPM * cos(2*theta_deg/180*PI) + C_2PPM * cos(4*theta_deg/180*PI))/(1 - B_2PPM * cos(2*theta_deg/180*PI) + C_2PPM * cos(4*theta_deg/180*PI)))/log(2);	
		}	
	}
	else{  // number_of_polarizations = 1
		F_max_1P = 0;
		F_max_2P = 0;
		for(i = 0; i < number_of_samples; i++){
			theta_deg = 180/(number_of_samples-1)*i - 90;
			fit_x_values_deg[i] = theta_deg;
			fit_1P_LD_values[i] = 1 + B_1PPM * cos(2*theta_deg/180*PI);
			fit_2P_LD_values[i] = 1 + B_2PPM * cos(2*theta_deg/180*PI) + C_2PPM * cos(4*theta_deg/180*PI);		
			if(fit_1P_LD_values[i] > F_max_1P){
				F_max_1P = fit_1P_LD_values[i];
			}
			if(fit_2P_LD_values[i] > F_max_2P){
				F_max_2P = fit_2P_LD_values[i];
			}	
		}
		// normalize the fluorescence intensities
		for(i = 0; i < number_of_samples; i++){
			fit_1P_LD_values[i] = fit_1P_LD_values[i]/F_max_1P;
			fit_2P_LD_values[i] = fit_2P_LD_values[i]/F_max_2P;
		}				
	}
	
	rmax_1P = (1 + B_1PPM)/(1 - B_1PPM);
	rmax_2P = (1 + B_2PPM + C_2PPM)/(1 - B_2PPM + C_2PPM);
	log2_rmax_1P = log(rmax_1P)/log(2);
	log2_rmax_2P = log(rmax_2P)/log(2);

	legend_1P = alpha0char+" = "+alpha0+degreechar+", "+sigmachar+" = "+sigma+degreechar+", 1P: r(max) = "+rmax_1P+" ("+1/rmax_1P+"), log"+subtwochar+"(r(max)) = "+log2_rmax_1P+", B1P = "+B_1PPM;
	legend_2P = alpha0char+" = "+alpha0+degreechar+", "+sigmachar+" = "+sigma+degreechar+", 2P: r(max) = "+rmax_2P+" ("+1/rmax_2P+"), log"+subtwochar+"(r(max)) = "+log2_rmax_2P+", B2P = "+B_2PPM+", C2P = "+C_2PPM;
	legend = legend_1P + "\n" + legend_2P;


	if(number_of_polarizations == 2){
		Plot.create("LD prediction "+alpha0char+"="+alpha0+degreechar+", "+sigmachar+"="+sigma+degreechar, "Membrane orientation, "+thetachar+" ["+degreechar+"]", "Log"+subtwochar+"(r)");	
		Plot.setLimits(-90, 90, -abs(log2_rmax_2P)*1.2, abs(log2_rmax_2P)*1.2);		
	}
	else{
		Plot.create("LD prediction "+alpha0char+"="+alpha0+degreechar+", "+sigmachar+"="+sigma+degreechar, "Membrane orientation, "+thetachar+" ["+degreechar+"]", "Normalized fluorescence intensity");	
		Plot.setLimits(-90, 90, -0.1, 1.4);				
	}

	Plot.setLineWidth(1);
	Plot.setColor("Blue");
	Plot.add("line", fit_x_values_deg, fit_1P_LD_values);

	Plot.setColor("Red");
	Plot.add("line", fit_x_values_deg, fit_2P_LD_values);
	Plot.addLegend(legend, "Auto");


	setBatchMode("exit and display");


	if(use_existing_plot == false){
		Plot.show();	
	}
	else{
		//look for previously made plot window
		image_titles = getList("image.titles");	
		previous_plot_title = "LD prediction "+alpha0char+"="+initial_alpha0+degreechar+", "+sigmachar+"="+initial_sigma+degreechar;
		for (i = 0; i<image_titles.length; i++){
			if(startsWith(image_titles[i], previous_plot_title)){
				selectImage(image_titles[i]);
			}
		}
		//update it!
		Plot.update();
		rename("LD prediction "+alpha0char+"="+alpha0+degreechar+", "+sigmachar+"="+sigma+degreechar);	
			
	}


	config_text2 = File.open(macro_path+"LD_prediction_config.txt");
	print(config_text2, "last_alpha0\t"+alpha0+"\r\n");
	print(config_text2, "last_sigma\t"+sigma+"\r\n");
	print(config_text2, "last_num_of_polarizations\t"+number_of_polarizations+"\r\n");
	print(config_text2, "last_modify_plot\t"+use_existing_plot+"\r\n");
	File.close(config_text2);

	operating_system = getInfo("os.name");
	if(matches(operating_system, "Windows")){
		run("Main Window [enter]");	
	}
	else{
		run("Main Window [return]");	
	}

	
}




macro "Accessory: Fit spline [f]" {
	run("Fit Spline");
} 




macro "Accessory: Align stack [a]"{
	List.setCommands;
	plugin_list = List.getList;
	if(indexOf(plugin_list, "StackReg ") > 0){
		run("StackReg ", "transformation=[Translation]");
	}
	else{
		if(indexOf(plugin_list, "StackReg") > 0){
			run("StackReg", "transformation=[Translation]");
		}
		else{
			showMessage("Please install the EPFL StackReg plugin to be able to align slices of image stacks.");			
		}
	}

	showStatus("Stack alignment finished!");
	
}


macro "Accessory: Correct bleaching within a stack [b]"{

	// calculating background for each slice
	
	setBatchMode(true);
	stack_image_id = getImageID();
	//run("32-bit");
	getDimensions(width, height, channels, slices, frames);
	Stack.getPosition(initial_channel, initial_slice, initial_frame);


	// if nothing selected, find a 5% x 5% rectangle the lowest avg intensity
	getSelectionBounds(x, y, sel_width, sel_height);
	//print(x+", "+y+", "+sel_width+", "+sel_height);
	if(x==0 && y==0 && sel_width == width && sel_height == height){
		min_mean = 100000;
		run("Select None");
		run("Z Project...", "projection=[Average Intensity]");
	
		sum_image_id = getImageID();

		sel_width = width/10;
		sel_height = height/10;
		
		for(x = width/20; x < width - width/10; x = x + width/20){
			for(y = height/20; y < height - height/10; y = y + height/20){
				makeRectangle(x, y, width/10, height/10);
				getStatistics(area, mean);
				if(mean < min_mean){
					min_mean = mean;
					min_x = x;
					min_y = y;

					//print(min_x+", "+min_y);
				}
			}
		}
		close();
	}
	else {
		min_x = x;
		min_y = y;
	}


	background_values = newArray(frames);
	mean_intensities  = newArray(frames);
	bkgd_subbed_intensities = newArray(frames);
	slice_numbers     = newArray(frames);

	for (frame = 1; frame <= frames; frame++) {
		slice_numbers[frame-1] = frame;
		makeRectangle(min_x, min_y, sel_width, sel_height);

		Stack.setPosition(initial_channel, initial_slice, frame);
		getStatistics(area, mean);
		background_values[frame-1] = mean;
		bkgd_current_frame = mean;
		if(frame == 1){bkgd_first_frame = mean;}
		
		run("Select None");
		//run("Subtract...", "value="+toString(floor(mean))+" slice");	
		getStatistics(area, mean);
		mean_intensities[frame-1] = mean;	
		mean_intensity_current_frame = mean;
		bkgd_subbed_intensities[frame-1] = mean_intensities[frame-1] - background_values[frame-1];
		
		//if(frame == 1){mean_intensity_first_frame = mean;}

		//run("Macro...", "code=[v = ("+mean_intensity_first_frame+" - "+bkgd_first_frame+")/("+mean_intensity_current_frame+" - "+bkgd_current_frame+")*(v-"+bkgd_current_frame+") + "+bkgd_current_frame+"]");
	}

	fitting_equation = "y = a - b * x";
	Fit.doFit(fitting_equation, slice_numbers, bkgd_subbed_intensities);
	a = Fit.p(0);
	b = Fit.p(1);

	for (frame = 1; frame <= frames; frame++) {
		Stack.setPosition(initial_channel, initial_slice, frame);
		run("Macro...", "code=[v = ("+toString(a / (a - b * frame))+"*(v-"+background_values[frame-1]+") + "+background_values[frame-1]+")]");		
	}



	Stack.setPosition(initial_channel, initial_slice, initial_frame);
	run("Select None");
	setBatchMode("exit and display");

	showStatus("Bleaching compensation finished!");

}




macro "Accessory: Show current two-polarization color scheme [I]"{
	setBatchMode(true);

	w = 30;
	h = 250;

	read_params_from_config_file();
	bkgd = 0;

	newImage("Intensity gradienat image", "32-bit black", w, h, 1);
	run("Macro...", "code=(1-x/w)*1.1");
	intensity_gradient_image_id = getImageID();
	

	//run("Macro...", "code=(h/2-y)/h*4");  			// makes a vertical gradient from +2 to -2
	//run("Macro...", "code=pow(2,(h/2-y)/h*4)");  		// makes a vertical gradient from +2 to 1/2

	newImage("Fh, Fv composite", "32-bit black", w, h, 2);
	Fh_Fv_composite_image_id = getImageID();
	setSlice(1);
	//run("Macro...", "code=sqrt(pow(2,(h/2-y)/h*4))");  		// makes a vertical gradient from +2 to +1/2
	run("Macro...", "code=sqrt(pow(2,(h/2-y)/h*4))*100");  		// makes a vertical gradient from 200 to 50
	setSlice(2);
	//run("Macro...", "code=1/sqrt(pow(2,(h/2-y)/h*4))");  		// makes a vertical gradient from +1/2 to +2
	run("Macro...", "code=1/sqrt(pow(2,(h/2-y)/h*4))*100");  	// makes a vertical gradient from 50 to 200

	imageCalculator("Multiply stack", Fh_Fv_composite_image_id, intensity_gradient_image_id);
		
	rg_image_id = make_rg_image(Fh_Fv_composite_image_id, 125, 4, two_pol_color_scheme);

	selectImage(Fh_Fv_composite_image_id);
	close();
	
	selectImage(intensity_gradient_image_id);
	close();

	selectImage(rg_image_id);


	setBatchMode("exit and display");
}



macro "Accessory: Show current multi-polarization color scheme [i]"{
	setBatchMode(true);

	w = 30;
	h = 250;

	read_params_from_config_file();

	bkgd = 0;
	

	newImage("Intensity gradient image", "32-bit black", w, h, 1);
	run("Macro...", "code=(1-x/w)");
	intensity_gradient_image_id = getImageID();
	

	//run("Macro...", "code=(h/2-y)/h*4");  			// makes a vertical gradient from +2 to -2
	//run("Macro...", "code=pow(2,(h/2-y)/h*4)");  		// makes a vertical gradient from +2 to 1/2

	//run("Macro...", "code=(h-y)/h*2");  			// makes a vertical gradient from 0 to 2

	newImage("Fh, Fv composite", "32-bit black", w, h, 2);
	Fh_Fv_composite_image_id = getImageID();
	setSlice(2);
	run("Macro...", "code=sqrt(pow(2,(h-y)/h*2))");  		// makes a vertical gradient from 100 to 50
	setSlice(1);
	run("Macro...", "code=sqrt(pow(2,-(h-y)/h*2))");  		// makes a vertical gradient from 100 to 200

	run("Z Project...", "projection=[Average Intensity]");
	rename("Avg composite [i] body");
	avg_composite_image = getImageID();

	imageCalculator("Divide stack", Fh_Fv_composite_image_id, avg_composite_image);

	imageCalculator("Multiply stack", Fh_Fv_composite_image_id, intensity_gradient_image_id);
		
	rg_image_id = make_rg_image(Fh_Fv_composite_image_id, 1, 3.85, multipol_color_scheme);  // composite image id, bright_adj, rg_range, color scheme

	selectImage(Fh_Fv_composite_image_id);
	close();
	
	selectImage(intensity_gradient_image_id);
	close();
	selectImage(avg_composite_image);
	close();

	selectImage(rg_image_id);


	setBatchMode("exit and display");
}



macro "Accessory: Show azimuth color wheel [o]"{
	setBatchMode(true);
	image_size = 200;
	newImage("Azimuth color wheel overlay", "16-bit black", image_size, image_size, 1);
	color_wheel_overlay_image_id = getImageID();
	run("Macro...", "  code=[if(d>w/5 && d<4/10*w){v=a/3.14159*360 +  1;}]");


	for (col = 0; col < image_size; col++) {
		for (row = 0; row < image_size; row++) {
			selectImage(color_wheel_overlay_image_id);
			azimuth = getPixel(col, row);
			if(azimuth > 0){
				azimuth = azimuth - 1;  // compenasting for using a range 1-361 earlier
				azimuth = azimuth + 300; // rotating the color wheel to make the colors consistent with Oldenbourg & al.
				hue = azimuth;
				if(hue > 360){hue = hue - 360;}
				if(hue > 360){hue = hue - 360;}
				if(hue > 360){hue = hue - 360;}

				hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
				if(hue >= 0 && hue < 60){
					red = 255;
					green = round(hsv_x * 255);
					blue = 0;
				}
		
				if(hue >= 60 && hue < 120){
					red = round(hsv_x * 255);
					green = 255;
					blue = 0;
				}
		
				if(hue >= 120 && hue < 180){
					red = 0;
					green = 255;
					blue = round(hsv_x * 255);
				}
				if(hue >= 180 && hue < 240){
					red = 0;
					green = round(hsv_x * 255);
					blue = 255;
				}
		
				if(hue >= 240 && hue < 300){
					red = round(hsv_x * 255);
					green = 0;
					blue = 255;
				}
				if(hue >= 300 && hue < 360){
					red = 255;
					green = 0;
					blue = round(hsv_x * 255);
				}
				setColor(red, green, blue);
				setLineWidth(1);
				
	
				Overlay.drawLine(col, row, col, row);
				Overlay.show;			
				
			}
		}
	}
	run("Flatten");
	color_wheel_RGB_image_id = getImageID();
	
	selectImage(color_wheel_overlay_image_id);
	close();

	selectImage(color_wheel_RGB_image_id);
	rename("Azimuth color wheel");
	
	setBatchMode("exit and display");
}





macro "Accessory: Show color coding for macro 'c' [e]"{
	if(isOpen("Color coding for macro 'c'")==true){
		selectImage("Color coding for macro 'c'");
		close();
	}

	newImage("Color coding for macro 'c'", "32-bit black", 500, 150, 1);

	setColor(50, 50, 255);
	setLineWidth(1);
	Overlay.drawLine(20, 42, 60, 42);		
	Overlay.show
		Overlay.drawString(alpha0char+", "+sigmachar+" values consistent with 1P data", 80, 50);
	Overlay.show;	
	

	setColor(255, 0, 0);
	setLineWidth(1);
	Overlay.drawLine(40, 62, 40, 62);		
	Overlay.show
	Overlay.drawRect(35, 57, 10, 10);	
	Overlay.show	
	Overlay.drawString("Mean, 95% CI of "+alpha0char+", "+sigmachar+" values from individual 2P experiments", 80, 70);
	Overlay.show;	


	setColor(255, 150, 150);
	setLineWidth(1);
	Overlay.drawLine(40, 82, 40, 82);		
	Overlay.show
	Overlay.drawString("Best match of pooled 2P data", 80, 90);
	Overlay.show;	


	setColor(255, 0, 255);
	setLineWidth(1);
	Overlay.drawLine(40, 102, 40, 102);		
	Overlay.show
	Overlay.drawString("Best match of pooled 1P, 2P data", 80, 110);
	Overlay.show;	


}




macro "Accessory: Close window [w]" {
	run("Close"); 
} 




macro "Accessory: Reload macros [r]" {
	run("Install...", getInfo("macro.filepath"));
}







//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------functions------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------



function read_params_from_config_file(){
	macro_path = getDirectory("macros");
	if(File.exists(macro_path+"pol_deconvolution_config.txt")){
		config_text = File.openAsString(macro_path+"pol_deconvolution_config.txt");
		config_rows = split(config_text,"\n");
	
		for(i=0; i<config_rows.length; i++){
			key_value_pair = split(config_rows[i]);
			if(key_value_pair[0] == "bleedthrough"               ){bleedthrough = key_value_pair[1];}
			if(key_value_pair[0] == "pol_modulation"             ){pol_modulation = key_value_pair[1];}
			if(key_value_pair[0] == "first_pixel_pol"            ){first_pixel_pol = key_value_pair[1];}
			if(key_value_pair[0] == "polarization_direction"     ){polarization_direction = key_value_pair[1];}			
			if(key_value_pair[0] == "pixel_shape"                ){pixel_shape = key_value_pair[1];}
			//if(key_value_pair[0] == "color_scheme"               ){color_scheme = key_value_pair[1];}
			if(key_value_pair[0] == "two_pol_color_scheme"       ){two_pol_color_scheme = key_value_pair[1];}
			if(key_value_pair[0] == "multipol_color_scheme"      ){multipol_color_scheme = key_value_pair[1];}
			if(key_value_pair[0] == "polarization_purity_hor"    ){polarization_purity_hor = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "polarization_purity_ver"    ){polarization_purity_ver = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "polarization_purity"        ){polarization_purity = parseFloat(key_value_pair[1]);}	
			if(key_value_pair[0] == "look_for_pol_angles"        ){look_for_pol_angles = key_value_pair[1];}	
			if(key_value_pair[0] == "ratiometric_fitting"        ){ratiometric_fitting = key_value_pair[1];}	
			if(key_value_pair[0] == "phase_fixed"                ){phase_fixed = key_value_pair[1];}	
			if(key_value_pair[0] == "phase1_deg"                 ){phase1_deg = parseFloat(key_value_pair[1]);}	
			if(key_value_pair[0] == "phase2_deg"                 ){phase2_deg = parseFloat(key_value_pair[1]);}	
			if(key_value_pair[0] == "orig_image_close"           ){orig_image_close = key_value_pair[1];}
			if(key_value_pair[0] == "make_16_bit"                ){make_16_bit = key_value_pair[1];}
			if(key_value_pair[0] == "do_segmentation"            ){do_segmentation = key_value_pair[1];}
			if(key_value_pair[0] == "look_for_vesicles"          ){look_for_vesicles = key_value_pair[1];}
			if(key_value_pair[0] == "membranes_or_filaments"     ){membranes_or_filaments = key_value_pair[3];}
				if(membranes_or_filaments == "membrane"){membranes_or_filaments = "Analyze a membrane";}
				if(membranes_or_filaments == "filament"){membranes_or_filaments = "Analyze a filament";}
				
			if(key_value_pair[0] == "save_files"                 ){save_files = key_value_pair[1];}
			if(key_value_pair[0] == "cutoff_distance_from_snake" ){cutoff_distance_from_snake = key_value_pair[1];}
			if(key_value_pair[0] == "determine_rmax_alpha0_sigma"){determine_rmax_alpha0_sigma = key_value_pair[1];
				if(determine_rmax_alpha0_sigma == "1PPM:"){determine_rmax_alpha0_sigma = "1PPM: r(max)";}
				if(determine_rmax_alpha0_sigma == "2PPM:"){determine_rmax_alpha0_sigma = "2PPM: r(max), "+alpha0char+", "+sigmachar;}				
			}	


			if(key_value_pair[0] == "starting_polarization_angle" ){starting_pol_angle = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "polarization_angle_increment"){pol_angle_increment = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "F_cutoff_for_LD_display"	  ){F_cutoff_for_LD_display = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "F_cutoff_for_azimuth_display"){F_cutoff_for_azimuth_display = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "azimuth_pixel_bin_size"      ){azimuth_pixel_bin_size = parseFloat(key_value_pair[1]);}
			if(key_value_pair[0] == "show_azimuths"               ){show_azimuths = key_value_pair[1];
				if(show_azimuths == "No"){show_azimuths = "No                   ";}}
			if(key_value_pair[0] == "show_LD"               	  ){show_LD = key_value_pair[1];}
			if(key_value_pair[0] == "azimuth_color"               ){azimuth_color = key_value_pair[1];}

		}
	}
}


function save_params_to_config_file(){
	macro_path = getDirectory("macros");
	config_file_handle = File.open(macro_path+"pol_deconvolution_config.txt");
	print(config_file_handle, "bleedthrough\t"+bleedthrough+"\r\n");
	print(config_file_handle, "polarization_direction\t"+polarization_direction+"\r\n");
	print(config_file_handle, "polarization_purity\t"+polarization_purity+"\r\n");
	print(config_file_handle, "polarization_purity_hor\t"+polarization_purity_hor+"\r\n");
	print(config_file_handle, "polarization_purity_ver\t"+polarization_purity_ver+"\r\n");
	print(config_file_handle, "look_for_pol_angles\t"+look_for_pol_angles+"\r\n");
	print(config_file_handle, "phase_fixed\t"+phase_fixed+"\r\n");
	print(config_file_handle, "phase1_deg\t"+phase1_deg+"\r\n");
	print(config_file_handle, "phase2_deg\t"+phase2_deg+"\r\n");
	print(config_file_handle, "pol_modulation\t"+pol_modulation+"\r\n");
	print(config_file_handle, "first_pixel_pol\t"+first_pixel_pol+"\r\n");
	print(config_file_handle, "pixel_shape\t"+pixel_shape+"\r\n");
	//print(config_file_handle, "color_scheme\t"+color_scheme+"\r\n");
	print(config_file_handle, "two_pol_color_scheme\t"+two_pol_color_scheme+"\r\n");
	print(config_file_handle, "multipol_color_scheme\t"+multipol_color_scheme+"\r\n");
	print(config_file_handle, "make_16_bit\t"+make_16_bit+"\r\n");
	print(config_file_handle, "orig_image_close\t"+orig_image_close+"\r\n");
	print(config_file_handle, "do_segmentation\t"+do_segmentation+"\r\n");
	print(config_file_handle, "ratiometric_fitting\t"+ratiometric_fitting+"\r\n");
	print(config_file_handle, "look_for_vesicles\t"+look_for_vesicles+"\r\n");
	print(config_file_handle, "membranes_or_filaments\t"+membranes_or_filaments+"\r\n");
	print(config_file_handle, "save_files\t"+save_files+"\r\n");
	print(config_file_handle, "determine_rmax_alpha0_sigma\t"+determine_rmax_alpha0_sigma+"\r\n");
	print(config_file_handle, "cutoff_distance_from_snake\t"+cutoff_distance_from_snake+"\r\n");

	print(config_file_handle, "starting_polarization_angle\t"+starting_pol_angle+"\r\n");
	print(config_file_handle, "polarization_angle_increment\t"+pol_angle_increment+"\r\n");
	print(config_file_handle, "F_cutoff_for_LD_display\t"+F_cutoff_for_LD_display+"\r\n");
	print(config_file_handle, "F_cutoff_for_azimuth_display\t"+F_cutoff_for_azimuth_display+"\r\n");
	print(config_file_handle, "azimuth_pixel_bin_size\t"+azimuth_pixel_bin_size+"\r\n");
	print(config_file_handle, "show_LD\t"+show_LD+"\r\n");
	print(config_file_handle, "show_azimuths\t"+show_azimuths+"\r\n");
	print(config_file_handle, "azimuth_color\t"+azimuth_color+"\r\n");

	File.close(config_file_handle);

}


function turn_RPM_image_into_a_composite(raw_image_id, bleedthrough, pol_modulation, first_pixel_pol, pixel_shape){
	selectImage(raw_image_id);
	raw_image_name = getTitle();
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);

	if(endsWith(raw_image_name, ".tif") || endsWith(raw_image_name, ".oib")){
		raw_composite_name = replace(raw_image_name, ".tif", "_rcomp.tif");
		raw_composite_name = replace(raw_image_name, ".oib", "_rcomp.tif");
	}
	else{
		raw_composite_name = raw_image_name+"_rcomp.tif";		
	}
	
	correct_negative_values(raw_image_id);


	run("Duplicate...", "title=processed");
	processed_image_id = getImageID();
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Left");
	}

	correct_bleedthrough(processed_image_id, bleedthrough);  //img_id, bleedthrough in % (0.05 = 5%)

	selectImage(processed_image_id);
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);

	
	first_image_id = make_image_from_columns(processed_image_id, 0, "first");  //img_id, offset, name of new image
	first_image_name = getTitle();

	selectImage(processed_image_id);
	second_image_id = make_image_from_columns(processed_image_id, 1, "second");  //img_id, offset, name of new image
	second_image_name = getTitle();

	selectImage(processed_image_id);
	close();


	//-----------------------------------correct for artifacts from different pixel locations in odd/even images----------------------
	for (col= 0; col<w/2; col++) {
		for (line = 0; line < h; line++){
			setPixel(col, line, (getPixel(col,line)+getPixel((col+1),line))/2 );
		}
	}


	if(first_pixel_pol == "horizontal"){
		hor_image_id = first_image_id;
		hor_image_name = first_image_name;
		ver_image_id = second_image_id;
		ver_image_name = second_image_name;
	}
	else{
		hor_image_id = second_image_id;
		hor_image_name = second_image_name;
		ver_image_id = first_image_id;
		ver_image_name = first_image_name;
	}


	//-----------------------------------correct for polarization purity of each image ----------------------

	updateDisplay();		 
	//--------------------------------------------------------------------------------------------------------------------------------

	run("Merge Channels...", "red="+hor_image_name+" green="+ver_image_name+" blue=*None* gray=*None* create keep");
	rename(raw_composite_name);


	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Right");
	}

	if(pixel_shape == "square"){
		run("Size...", "width="+w+" height="+h+" interpolation=None");
	}

	composite_img_id = getImageID();

	if(pixel_width > pixel_height){pixel_height = pixel_width;}  // if image is 'stretched' because of pol. modulation, use the coarser resolution
	else{pixel_width = pixel_height;}
	run("Properties...", "unit="+pixel_size_unit+" pixel_width="+pixel_width+" pixel_height="+pixel_height);



	
	selectImage(hor_image_id);
	close();
	selectImage(ver_image_id);
	close();
	return(composite_img_id);
}





function guess_params(composite_id){
	w = getWidth;
	h = getHeight;
	

	selectImage(composite_id);
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Left");
	}

	selectImage(composite_id);
	setSlice(1);
	run("Duplicate...", "title=hor_image.tif");
	hor_image_id = getImageID();
	
	selectImage(composite_id);
	setSlice(2);
	run("Duplicate...", "title=ver_image.tif");
	ver_image_id = getImageID();

	imageCalculator("average create", hor_image_id, ver_image_id);
	avg_image_id = getImageID();
	rename("avg_image");

	min_mean = 100000;
	for(x = w/20; x < w - w/10; x = x + w/20){
		for(y = h/20; y < h - h/10; y = y + h/20){
			makeRectangle(x, y, w/10, h/10);
			getStatistics(area, mean);
			if(mean < min_mean){
				min_mean = mean;
				min_x = x;
				min_y = y;
			}
		}
	}
	
	selectImage(hor_image_id);
	makeRectangle(min_x, min_y, w/10, h/10);
	getStatistics(area, hor_mean, hor_min, hor_max, hor_stdev);
	hor_bkgd = round(hor_mean);
	run("Select All");
	run("Subtract...", "value="+hor_bkgd);
	run("Add...", "value=1");

	selectImage(ver_image_id);
	makeRectangle(min_x, min_y, w/10, h/10);
	getStatistics(area, ver_mean, ver_min, ver_max, ver_stdev);
	ver_bkgd = round(ver_mean);
	run("Select All");
	run("Subtract...", "value="+ver_bkgd);
	run("Add...", "value=1");





	//  make a mask ------------------------------------------------------------------------
	imageCalculator("average create", hor_image_id, ver_image_id);
	rename("avg_image");
	setAutoThreshold();
	run("Convert to Mask");
	run("Divide...", "value=255");
	rename("mask");
	mask_image_id = getImageID();

	imageCalculator("multiply",hor_image_id, mask_image_id);
	imageCalculator("multiply",ver_image_id, mask_image_id);


	//---stripe correction ----------------------------------------------------------
	stripe_style = "stripes";
	stripe_adj = 1.00;
	hor_stripe_adj = 1.00;
	ver_stripe_adj = 1.00;
	hor_checker_adj = 1.00;
	ver_checker_adj = 1.00;


	selectImage(ver_image_id);
	w = getWidth();
	h = getHeight();
	even_pixel_values = 0.00;
	odd_pixel_values = 0.00;
	top_left_diag_pixel_values = 0.00;   //diagonal, including pixel 0,0
	top_left_1_diag_pixel_values = 0.00; //diagonal, including pixel 1,0
	for (col= 0; col < w; col++) {
		for (row = 0; row < h; row++){
			pixel_value = getPixel(col,row);
			if(col/2 == floor(col/2)){
				even_pixel_values += pixel_value;
			}
			else{
				odd_pixel_values += pixel_value;
			}
			if( (row+col)/2 == floor((row+col)/2)){
				top_left_diag_pixel_values += pixel_value;
			}
			else{
				top_left_1_diag_pixel_values += pixel_value;
			}
		}
	}
	ver_stripe_adj = sqrt(odd_pixel_values/even_pixel_values);
	ver_checker_adj = sqrt(top_left_1_diag_pixel_values/top_left_diag_pixel_values);


	selectImage(hor_image_id);
	even_pixel_values = 0.00;
	odd_pixel_values = 0.00;
	top_left_diag_pixel_values = 0.00;   //diagonal, including pixel 0,0
	top_left_1_diag_pixel_values = 0.00; //diagonal, including pixel 1,0
	for (col= 0; col < w; col++) {
		for (row = 0; row < h; row++){
			pixel_value = getPixel(col,row);
			if(col/2 == floor(col/2)){
				even_pixel_values += pixel_value;
			}
			else{
				odd_pixel_values += pixel_value;
			}
			if( (row+col)/2 == floor((row+col)/2)){
				top_left_diag_pixel_values += pixel_value;
			}
			else{
				top_left_1_diag_pixel_values += pixel_value;
			}
		}
	}
	hor_stripe_adj = sqrt(odd_pixel_values/even_pixel_values);
	hor_checker_adj = sqrt(top_left_1_diag_pixel_values/top_left_diag_pixel_values);


	if (abs(log(hor_stripe_adj)) > abs(log(ver_stripe_adj))){
		ver_stripe_adj = 1.00;
		stripe_adj = hor_stripe_adj;
	}
	else{
		hor_stripe_adj = 1.00;
		stripe_adj = ver_stripe_adj;

	}

	if (abs(log(hor_checker_adj)) > abs(log(ver_checker_adj))){
		ver_checker_adj = 1.00;
		checker_adj = hor_checker_adj;
	}
	else{
		hor_checker_adj = 1.00;
		checker_adj = ver_checker_adj;
	}


	if (abs(log(stripe_adj)) > abs(log(checker_adj))){
		stripe_style = "stripes";
	}
	else{
		stripe_style = "checkerboard";
		hor_stripe_adj = hor_checker_adj;
		ver_stripe_adj = ver_checker_adj;
	}


	if( (abs(log(checker_adj)) < abs(log(1.01))) && (abs(log(checker_adj)) < abs(log(1.01))) ){
		stripe_style = "none";
		hor_stripe_adj = 1.00;
		ver_stripe_adj = 1.00;
		hor_checker_adj = 1.00;
		ver_checker_adj = 1.00;

	}



	//---end of stripe correction ----------------------------------------------------------


	//---overall brightness adjustment correction ----------------------------------------------------------
	selectImage(avg_image_id);
	getMinAndMax(avg_image_min_value,avg_image_max_value);
	bright_adj = round((avg_image_max_value-avg_image_min_value)/1.2);
	selectImage(avg_image_id);
	close();

	//-------------------------------------------------------
	//---hor/ver brightness adjustment correction ----------------------------------------------------------
	imageCalculator("divide create 32-bit", hor_image_id, ver_image_id);
	run("Log");
	rename("log_ratio_image");
	log_ratio_image_id = getImageID();

	getStatistics(area, mean_rg_logratio, min_rg_logratio, max_rg_logratio, logratio_stdev);
	hor_adj = 1;
	ver_adj = 1;
	if(mean_rg_logratio > 0){ver_adj = exp(mean_rg_logratio);}  // Fh > Fv
	if(mean_rg_logratio < 0){hor_adj = exp(-mean_rg_logratio);}




	//-------------------------------------------------------
	//---R/G range calculation ----------------------------------------------------------

	if(range_lock == true){rg_range = locked_rg_range;}
	else{
		if((mean_rg_logratio - min_rg_logratio) > (max_rg_logratio - mean_rg_logratio)){
			log_halfrange = mean_rg_logratio-min_rg_logratio;
		}
		else{
			log_halfrange = max_rg_logratio - mean_rg_logratio;
		}
	
		getHistogram(hist_values,hist_counts,100);
		
		total_counts = 0;
		for (index = 0; index < 100; index++){
			total_counts = total_counts + hist_counts[index];
		}
		cutoff_counts = total_counts * 0.1;
	
		incr = 0;
		out_of_range_percentage = 0;
		while(out_of_range_percentage < 0.025){
			lower_limit = mean_rg_logratio - log_halfrange * (1 - incr/100);
			upper_limit = mean_rg_logratio + log_halfrange * (1 - incr/100);
			num_of_out_of_range = 0;
			for (index = 0; index < 50; index++){
				
				if ((hist_values[index] < lower_limit) || (hist_values[index] > upper_limit)){
					num_of_out_of_range = num_of_out_of_range + hist_counts[index];
				}
				if ((hist_values[99-index] < lower_limit) || (hist_values[99-index] > upper_limit)){
					num_of_out_of_range = num_of_out_of_range + hist_counts[99-index];
				}
			}
			out_of_range_percentage = num_of_out_of_range/total_counts;
			incr++;
		}
		rg_range = (round((exp(mean_rg_logratio-lower_limit))*100))/100;
	}

	selectImage(composite_id);
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Right");
	}
	selectImage(log_ratio_image_id);
	close();
	selectImage(mask_image_id);
	close();
	selectImage(hor_image_id);
	close();
	selectImage(ver_image_id);
	close();



	function_array = newArray(9);
	function_array[0] = hor_bkgd;
	function_array[1] = ver_bkgd;
	function_array[2] = bright_adj;
	function_array[3] = rg_range;
	function_array[4] = hor_adj;
	function_array[5] = ver_adj;
	function_array[6] = stripe_style;
	function_array[7] = hor_stripe_adj;
	function_array[8] = ver_stripe_adj;

	return function_array;
}





function get_manual_params(params, rg_image_id){
//pol_modulation, first_pixel_pol, bleedthrough, pixel_shape, color_scheme, composite_img_id, hor_min, ver_min, bright_adj, rg_range, hor_adj, ver_adj, stripe_style, hor_stripe_adj, ver_stripe_adj, polarization_purity_hor, polarization_purity_ver;

	pol_modulation  = params[0];
	first_pixel_pol = params[1];
	bleedthrough    = params[2];
	pixel_shape     = params[3];
	color_scheme    = params[4];
	polarization_purity_hor = params[5];
	polarization_purity_ver = params[6];

	hor_bkgd_adj    = params[7];
	ver_bkgd_adj    = params[8];
	bright_adj      = params[9];
	rg_range        = params[10];
	hor_adj         = params[11];
	ver_adj         = params[12];
	stripe_style    = params[13];
	hor_stripe_adj  = params[14];
	ver_stripe_adj  = params[15];

	selectImage(rg_image_id);
	getLocationAndSize(x, y, width, height);

	Dialog.create("Image processing parameters");
	if(screenHeight > 800){Dialog.addMessage("Linear dichroism visualization parameters:");}
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addChoice("Color scheme (Fh, Fv):",newArray("red-green","magenta-green","magenta-yellow","cyan-yellow", "red-blue","yellow-blue","rainbow"), color_scheme);
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Pure color means dichroic ratio of", rg_range, 2,5, "or higher");
	Dialog.setInsets(0, 40, 0);
	Dialog.addCheckbox("Use the above value for the next image", range_lock);
	Dialog.addMessage("Range of intensities shown:");
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Shown as black:", hor_bkgd_adj,0,5,"(hor. polarization)");
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Shown as black:", ver_bkgd_adj,0,5,"(ver. polarization)");
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Shown as 100% bright:", bright_adj+(hor_bkgd_adj+ver_bkgd_adj)/2,0,5,"(both polarizations)");
	Dialog.setInsets(5, 40, 0);
	if(mixed_image == true){
		Dialog.addCheckboxGroup(1,2,newArray("Process whole stack", "Close orig. image"), newArray(process_whole_stack, orig_image_close));	
	}
	else{
		Dialog.addCheckbox("Close orig. image", orig_image_close);
	}

	if(screenHeight > 900){
		Dialog.setInsets(10, 20, 0);
		Dialog.addMessage("_________________________________________________");	
		Dialog.setInsets(0, 20, 0);
		Dialog.addMessage(" \nCorrections to underlying image:");
	}
	else{
		Dialog.setInsets(0, 20, 5);
		Dialog.addMessage("---------------------------------------------");		
	}
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Hor. intensity multiplied by", hor_adj);
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Ver. intensity multiplied by", ver_adj);
	Dialog.setInsets(10, 0, 0);


	if(screenHeight > 900){
		Dialog.setInsets(10, 20, 0);
		Dialog.addMessage("_________________________________________________");	
		Dialog.setInsets(0, 20, 0);
		Dialog.addMessage(" \nMicroscope parameters:");
	}
	else{
		Dialog.setInsets(0, 20, 0);
		Dialog.addMessage("---------------------------------------------");		
	}

	if(mixed_image == true){
		Dialog.setInsets(0, 0, 0);
		Dialog.addChoice("Polarization modulation",newArray("columns","rows"), pol_modulation);
		Dialog.setInsets(0, 0, 0);
		Dialog.addChoice("Polarization in top left pixel",newArray("horizontal","vertical"), first_pixel_pol);
		Dialog.setInsets(0, 0, 0);
		Dialog.addChoice("Pixel shape",newArray("square","rectangle"), pixel_shape);
		if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
		Dialog.addNumber("Bleedthrough:", bleedthrough*100,2,5,"%");
	}
	
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Hor. polarization purity:", polarization_purity_hor,2,5,"%");
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Ver. polarization purity:", polarization_purity_ver,2,5,"%");


	if(screenHeight > 900){
		Dialog.setInsets(10, 20, 0);
		Dialog.addMessage("_________________________________________________");	
		Dialog.setInsets(10, 0, 0);
	}
	else{
		Dialog.setInsets(0, 20, 5);
		Dialog.addMessage("---------------------------------------------");		
	}

	Dialog.setInsets(0, 0, 0);
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addChoice("Perform quantitative analysis:", newArray("none", "1PPM: r(max)", "2PPM: r(max), "+alpha0char+", "+sigmachar), determine_rmax_alpha0_sigma);
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	else{Dialog.setInsets(0, 80, 0);}
	Dialog.addRadioButtonGroup("", newArray("Analyze a membrane", "Analyze a filament"), 1, 2, membranes_or_filaments);
	Dialog.setInsets(0, 80, 0);

	Dialog.addCheckbox("Attempt automated vesicle recognition", look_for_vesicles);
	Dialog.addNumber("Membrane/filament selection width:", cutoff_distance_from_snake,3,5,muchar+"m");

//	Dialog.addCheckbox("Attempt automated vesicle recognition", look_for_vesicles);
	Dialog.setInsets(0, 0, 0);


//	phase_fixed_choices = newArray("found from fit", "set by user");
	phase_choices = newArray("Identical, found from fit", "Identical, set by user", "Distinct, found from fit", "Distinct, set by user");
//	Dialog.addChoice("Deviations from hor./ver. direction:", phase_fixed_choices, phase_fixed_choices[phase_fixed]);
	Dialog.addChoice("Pol. deviations from hor./ver.:", phase_choices, phase_choice);

	/*
	Dialog.setInsets(0, 0, 0);
	look_for_pol_angles_choices = newArray("identical", "distinct");
	if(phase_fixed == true){
		if(phase1_deg == phase2_deg){look_for_pol_angles = 1;}
		else{look_for_pol_angles = 2;}
	}
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addChoice("Deviations from hor./ver. direction:", look_for_pol_angles_choices, look_for_pol_angles_choices[look_for_pol_angles-1]);

	*/

	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Deviation from hor. direction:", phase1_deg,3,5,"degrees");
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addNumber("Deviation from ver. direction:", phase2_deg,3,5,"degrees");
	Dialog.setInsets(0, 40, 0);
	Dialog.addCheckbox("Save result files", save_files);
	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}

	Dialog.setLocation(x + width,y);

	Dialog.show();

	color_scheme = Dialog.getChoice();
	rg_range = Dialog.getNumber();
	range_lock = Dialog.getCheckbox();
	if(range_lock == true){
		locked_rg_range = rg_range;
	}
	hor_bkgd_adj = Dialog.getNumber();
	ver_bkgd_adj = Dialog.getNumber();
	bright_adj = Dialog.getNumber();
	if(mixed_image == true){
		process_whole_stack = Dialog.getCheckbox();	
	}
	else{
		process_whole_stack = false;
	}
	orig_image_close = Dialog.getCheckbox();

	hor_adj = Dialog.getNumber();
	ver_adj = Dialog.getNumber();
	bright_adj = bright_adj - (hor_bkgd_adj + ver_bkgd_adj)/2;

	if(mixed_image == true){
		pol_modulation = Dialog.getChoice();
		first_pixel_pol = Dialog.getChoice();
		pixel_shape = Dialog.getChoice();
		bleedthrough = Dialog.getNumber()/100;	
	}
	else{
		pol_modulation = "columns";
		first_pixel_pol = "horizontal";
		pixel_shape = "square";
		bleedthrough = 0;		
	}
	
	polarization_purity_hor = Dialog.getNumber();
	polarization_purity_ver = Dialog.getNumber();

	determine_rmax_alpha0_sigma = Dialog.getChoice();
	if(determine_rmax_alpha0_sigma == "none"){
		make_16_bit = false;		
		further_processing = false;
		do_segmentation = false;
	}
	else{
		make_16_bit = true;
		further_processing = true;
		do_segmentation = true;
	}
	membranes_or_filaments = Dialog.getRadioButton();
	look_for_vesicles = Dialog.getCheckbox();
	cutoff_distance_from_snake = Dialog.getNumber();
	phase_choice = Dialog.getChoice();
//	look_for_phase = Dialog.getChoice();
//	if(look_for_phase == "found from fit"){phase_fixed = false;}
	if(phase_choice == "Identical, found from fit" || phase_choice == "Distinct, found from fit"){phase_fixed = false;}
	else{phase_fixed = true;}
//	phase_choice = Dialog.getChoice();
//	if(phase_choice == "identical"){look_for_pol_angles = 1;}
//	if(phase_choice == "distinct" ){look_for_pol_angles = 2;}
	if(phase_choice == "Identical, found from fit" || phase_choice == "Identical, set by user"){look_for_pol_angles = 1;}
	else{look_for_pol_angles = 1;}
	
	
	phase1_deg = Dialog.getNumber();
	phase2_deg = Dialog.getNumber();
	save_files = Dialog.getCheckbox();

	nparams = newArray(17);
	nparams[0]  = pol_modulation;
	nparams[1]  = first_pixel_pol;
	nparams[2]  = bleedthrough;
	nparams[3]  = pixel_shape;
	nparams[4]  = color_scheme;
	nparams[5]  = polarization_purity_hor;
	nparams[6]  = polarization_purity_ver;

	nparams[7]  = hor_bkgd_adj;
	nparams[8]  = ver_bkgd_adj;
	nparams[9]  = bright_adj;
	nparams[10] = rg_range;
	nparams[11] = hor_adj;
	nparams[12] = ver_adj;
	nparams[13] = stripe_style;
	nparams[14] = hor_stripe_adj;
	nparams[15] = ver_stripe_adj;

	return nparams;

}


function get_manual_multipol_params(params){


	bkgd 							 = params[0];
	bright_adj      				 = params[1];
	show_LD  						 = params[2];
	color_scheme 					 = params[3];
	rg_range 						 = params[4];
	F_cutoff_for_LD_azimuth_display	 = params[5];
	show_azimuths        			 = params[6];
	azimuth_color					 = params[7];
	starting_pol_angle				 = params[8];
	pol_angle_increment				 = params[9];
	azimuth_pixel_bin_size			 = params[10];
	azimuth_scaling_factor           = params[11];
	F_cutoff_for_LD_azimuth_display	 = params[12];
	use_fitting_for_LD				 = params[13];

	
	Dialog.create("Image processing parameters");

	Dialog.addMessage("Range of fluorescence intensities shown:");
	Dialog.addNumber("Shown as black:", bkgd, 0, 5, "");
	Dialog.addNumber("Shown as 100% bright:", bright_adj,0,5, "");	
	Dialog.addNumber("Starting polarization angle: ", starting_pol_angle, 0, 3, degreechar+" (0"+degreechar+" = horizontal)");
	Dialog.addNumber("Polarization angle increment: ", pol_angle_increment, 0, 3, degreechar+" (clockwise)");
	Dialog.addMessage("__________________________________________________");	
	Dialog.setInsets(0, 10, 0);
	Dialog.addRadioButtonGroup("Linear dichroism visualization:", newArray("None", "Crude/fast", "Accurate/slow"), 1, 3, show_LD);
	Dialog.addNumber("Show linear dichroism for F >", round(F_cutoff_for_LD_display*100),0,2,"% of Fmax");

	Dialog.addChoice("Color scheme",newArray("white-red","white-green","white-blue","white-cyan","white-magenta","rainbow"), color_scheme);
	Dialog.addNumber("Saturated color means dichroic ratio of", rg_range, 2,5, "or higher");
	Dialog.setInsets(0, 40, 0);
	Dialog.addCheckbox("Use the above value for the next image", range_lock);

	Dialog.addMessage("   ");	
	Dialog.setInsets(0, 10, 0);
	Dialog.addRadioButtonGroup("Azimuth visualization:", newArray("No                   ", "Yes"), 1, 3, show_azimuths);
	Dialog.addNumber("Show azimuths for F >", round(F_cutoff_for_azimuth_display*100),0,2,"% of Fmax");
	Dialog.addNumber("Pixel binning: ",azimuth_pixel_bin_size,0,2,"pixels");
	Dialog.addChoice("Azimuth display color", newArray("direction","white","black","gray","red","green","blue","magenta","cyan","yellow"), azimuth_color);
	Dialog.addNumber("Azimuth size scaling factor", azimuth_scaling_factor, 2,5,"");
	//Dialog.setInsets(0, 40, 0);
	//Dialog.addCheckbox("Use the above value for the next image", azimuth_size_lock);
	Dialog.addMessage("__________________________________________________");	
	Dialog.addMessage("Linear dichroism quantitation/analysis:");
	Dialog.addMessage("(Results do not depend on visualization settings)");
	Dialog.setInsets(10, 0, 0);

	if(screenHeight < 801){Dialog.setInsets(0, 0, 0);}
	Dialog.addChoice("Perform quantitative analysis:", newArray("none", "1PPM: r(max)", "2PPM: r(max), "+alpha0char+", "+sigmachar), determine_rmax_alpha0_sigma);
	Dialog.setInsets(0, 80, 0);
	Dialog.addRadioButtonGroup("", newArray("Analyze a membrane", "Analyze a filament"), 1, 2, membranes_or_filaments);
	Dialog.setInsets(0, 80, 0);
	Dialog.addCheckbox("Attempt automated vesicle recognition", look_for_vesicles);
	Dialog.addNumber("Membrane/filament selection width:", cutoff_distance_from_snake,3,5,muchar+"m");
	Dialog.setInsets(0, 80, 0);
	Dialog.addCheckbox("Show additional plots", show_additional_plots);
	Dialog.setInsets(0, 80, 0);
	Dialog.addCheckbox("Save result files", save_files);
	Dialog.setInsets(0, 80, 0);
	Dialog.addCheckbox("Close orig. image when done", orig_image_close);

	//Dialog.setLocation(x + width,y);

	Dialog.show();
	bkgd = Dialog.getNumber();
	bright_adj = Dialog.getNumber();
	starting_pol_angle = Dialog.getNumber();
	pol_angle_increment = Dialog.getNumber();

	show_LD = Dialog.getRadioButton();
	F_cutoff_for_LD_display = Dialog.getNumber();
	F_cutoff_for_LD_display = F_cutoff_for_LD_display/100;
	multipol_color_scheme = Dialog.getChoice();
	color_scheme = multipol_color_scheme;
	rg_range = Dialog.getNumber();
	range_lock = Dialog.getCheckbox();
	if(range_lock == true){
		locked_rg_range = rg_range;
	}

	
	show_azimuths = Dialog.getRadioButton();
	F_cutoff_for_azimuth_display = Dialog.getNumber();
	F_cutoff_for_azimuth_display = F_cutoff_for_azimuth_display/100;
	azimuth_pixel_bin_size = Dialog.getNumber();
	azimuth_color = Dialog.getChoice();
	azimuth_scaling_factor = Dialog.getNumber();
	//azimuth_size_lock = Dialog.getChoice();
	
	determine_rmax_alpha0_sigma = Dialog.getChoice();
	if(determine_rmax_alpha0_sigma == "none"){	
		further_processing = false;
		do_segmentation = false;
	}
	else{
		further_processing = true;
		do_segmentation = true;
	}
	membranes_or_filaments = Dialog.getRadioButton();
	look_for_vesicles = Dialog.getCheckbox();
	cutoff_distance_from_snake = Dialog.getNumber();
	show_additional_plots = Dialog.getCheckbox();
	save_files = Dialog.getCheckbox();
	orig_image_close = Dialog.getCheckbox();


	nparams = newArray(14);
	nparams[0] = bkgd;
	nparams[1] = bright_adj;
	nparams[2] = show_LD;
	nparams[3] = color_scheme;
	nparams[4] = rg_range;
	// range_lock not sent: doesn't affect rendering of the current image
	nparams[5] = F_cutoff_for_LD_display;
	nparams[6] = show_azimuths;
	nparams[7] = azimuth_color;
	nparams[8] = starting_pol_angle;
	nparams[9] = pol_angle_increment;
	nparams[10] = azimuth_pixel_bin_size;
	nparams[11] = azimuth_scaling_factor;
	nparams[12] = F_cutoff_for_azimuth_display;
	nparams[13] = use_fitting_for_LD;

	// orig_image_close not sent: doesn't affect rendering of the current image
	// determine_alpha0_sigma not sent
	// cutoff_distance_from snake not sent
	// look_for_vesicles not sent
	// show_additional_plots not sent
	// save_files not sent
	
	return nparams;

}


function make_processed_composite(raw_composite_id, params){
	hor_bkgd_adj    = params[0];
	ver_bkgd_adj    = params[1];
	bright_adj      = params[2];
	rg_range        = params[3];
	hor_adj         = params[4];
	ver_adj         = params[5];
	stripe_style    = params[6];
	hor_stripe_adj  = params[7];
	ver_stripe_adj  = params[8];


	setBatchMode(true);
	
	selectImage(raw_composite_id);
	h = getHeight();
	w = getWidth();
	raw_composite_image_name = getTitle();

	processed_composite_name = replace(raw_composite_image_name, "_rcomp.tif", "_pcomp.tif");

	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Left");
	}

	setSlice(1);
	run("Duplicate...", "title=hor_image.tif");
	hor_image_id = getImageID();
	hor_image_name = getTitle();
	run("Subtract...", "value="+hor_bkgd_adj);
	run("Multiply...", "value="+hor_adj);
	run("Add...", "value=1");
	
	selectImage(raw_composite_id);
	setSlice(2);
	run("Duplicate...", "title=ver_image.tif");
	ver_image_id = getImageID();
	ver_image_name = getTitle();
	run("Subtract...", "value="+ver_bkgd_adj);
	run("Multiply...", "value="+ver_adj);
	run("Add...", "value=1");

	selectImage(raw_composite_id);
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Right");
	}


	//--adjust for striping artifact in hor image:-------------------------------------------------------------------------------
	// unidirectional scanning: stripes


	if((hor_stripe_adj != 1) && (stripe_style == "stripes")){
		selectImage(ver_image_id);
		for (col= 0; col < w-1; col+=2) {
			for (row = 0; row < h; row++){
				setPixel(col, row, (getPixel(col,row)*hor_stripe_adj ) );
				setPixel(col+1, row, (getPixel(col+1,row)/hor_stripe_adj ) );
			}
		}
	}
	

	if((ver_stripe_adj != 1) && (stripe_style == "stripes")){	
		selectImage(ver_image_id);
		for (col= 0; col < w-1; col+=2) {
			for (row = 0; row < h; row++){
				setPixel(col, row, (getPixel(col,row)*ver_stripe_adj ) );
				setPixel(col+1, row, (getPixel(col+1,row)/ver_stripe_adj ) );
			}
		}
	}
	


	// bidirectional scanning: checkerboard
	if((hor_stripe_adj != 1) && (stripe_style == "checkerboard")){
		selectImage(hor_image_id);
		for (col= 0; col < w; col++) {
			for (row = 0; row < h; row++){
				if( (row+col)/2 == floor((row+col)/2)){
					setPixel(col, row, (getPixel(col,row)*hor_stripe_adj ) );
				}
				else{
					setPixel(col, row, (getPixel(col,row)/hor_stripe_adj ) );
				}
			}
		}
	}

	if((ver_stripe_adj != 1) && (stripe_style == "checkerboard")){
		selectImage(ver_image_id);
		for (col= 0; col < w; col++) {
			for (row = 0; row < h; row++){
				if( (row+col)/2 == floor((row+col)/2)){
					setPixel(col, row, (getPixel(col,row)*ver_stripe_adj ) );
				}
				else{
					setPixel(col, row, (getPixel(col,row)/ver_stripe_adj ) );
				}
			}
		}
	}

// polarization purity correction

	Fh_purity_float = parseFloat(polarization_purity_hor)/100;
	Fv_purity_float = parseFloat(polarization_purity_ver)/100;


	Fh_expression = "("+toString(Fv_purity_float) +" * A + (" + toString(Fh_purity_float) + " - 1) * B) / " + toString((Fh_purity_float + Fv_purity_float - 1));
	Fv_expression = "("+toString((Fv_purity_float - 1)) +" * A + " + toString(Fh_purity_float) + " * B) / " + toString((Fh_purity_float + Fv_purity_float - 1));

	run("Image Expression Parser (Macro)", "expression=["+ Fh_expression +"] a="+hor_image_name+" b="+ver_image_name+" c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	pure_Fh_image_id = getImageID();
	rename(hor_image_name);
	run("Properties...", "channels=1 slices=1 frames=1 unit="+pixel_size_unit+" pixel_width="+pixel_width+" pixel_height="+pixel_height+" voxel_depth=0.005");


	run("Image Expression Parser (Macro)", "expression=["+ Fv_expression +"] a="+hor_image_name+" b="+ver_image_name+" c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	pure_Fv_image_id = getImageID();
	rename(ver_image_name);
	run("Properties...", "channels=1 slices=1 frames=1 unit="+pixel_size_unit+" pixel_width="+pixel_width+" pixel_height="+pixel_height+" voxel_depth=0.005");

	selectImage(hor_image_id);
	close();

	selectImage(ver_image_id);
	close();

	hor_image_id = pure_Fh_image_id;
	ver_image_id = pure_Fv_image_id;
	

//
	

	selectImage(raw_composite_id);
	run("Duplicate...", "duplicate");
	processed_composite_image_id = getImageID();
	rename(processed_composite_name);

	selectImage(hor_image_id);
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Right");
	}
	run("Copy");
	selectImage(processed_composite_image_id);
	setSlice(1);
	run("Paste");
	selectImage(hor_image_id);
	close();


	selectImage(ver_image_id);
	if(pol_modulation == "rows"){
		run("Rotate 90 Degrees Right");
	}
	run("Copy");
	selectImage(processed_composite_image_id);
	setSlice(2);
	run("Paste");
	selectImage(ver_image_id);
	close();

	selectImage(processed_composite_image_id);

	setSlice(2);
	run("Green");
	setSlice(1);
	run("Red");

	return(processed_composite_image_id);
}





function make_rg_image(composite_img_id, bright_adj, rg_range, color_scheme){


	setBatchMode(true);

	
	selectImage(composite_img_id);
	processed_composite_image_name = getTitle();


	setSlice(1);
	run("Duplicate...", "title=hor_image.tif");
	hor_image_id = getImageID();

	
	selectImage(composite_img_id);
	setSlice(2);
	run("Duplicate...", "title=ver_image.tif");
	ver_image_id = getImageID();

	//  make a ratio image, 32-bit   ---------------------------------------------------------------------------------------------
	imageCalculator("divide create 32-bit", hor_image_id, ver_image_id);
	rename("ratio_image");
	ratio_image_id = getImageID();




	//  make a normalized average intensity image, 32-bit ------------------------------------------------------------------------
	imageCalculator("average create", hor_image_id, ver_image_id);
	rename("avg_image");
	avg_image_id = getImageID();
	run("32-bit");
	run("Divide...", "value="+toString(bright_adj-bkgd));
	//now the intensity corresponding to bright_adj - bkgd is equal to 1


	//  modify the ratio image to use a chosen rg range------------------------------------------------------------------------------
	selectImage(ratio_image_id);
	run("Log");
	run("Multiply...", "value="+log(3)/log(rg_range));
	run("Exp");				

	//  make a 'blank' rgb image (RGB) -------------------------------------------------------------------------------------------
	run("Duplicate...", "title=rgb_image");
	run("RGB Color");
	rgb_image_id = getImageID();


		

	scale_factor = 0.5;
	
	
	// fill up the RGB image
	if(color_scheme == "red-green"){
		color_scheme_suffix = "_RG.tif";
		selectImage(ratio_image_id);
		run("Duplicate...", "title=temp_ratio_image");
		temp_ratio_image_id = getImageID();
		run("Macro...", "code=[if(v > 3) v=1000]");  
		run("Macro...", "code=[if(v < 1/3) v=0.001]");  

		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");



		red_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";		
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);
		

		
		green_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";	
		run("Image Expression Parser (Macro)", "expression=["+green_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		green_image_id = getImageID();
		rename("green image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(2);
		imageCalculator("Multiply", rgb_image_id, green_image_id);


		
		blue_expression = "0 * A";	
		run("Image Expression Parser (Macro)", "expression=["+blue_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		blue_image_id = getImageID();
		rename("blue image");

		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, blue_image_id);


		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(green_image_id);
		close();
		selectImage(blue_image_id);
		close();

		selectImage(temp_ratio_image_id);
		close();

	}


	if(color_scheme == "magenta-green"){
		color_scheme_suffix = "_MG.tif";
		selectImage(ratio_image_id);
		run("Duplicate...", "title=temp_ratio_image");
		temp_ratio_image_id = getImageID();
		run("Macro...", "code=[if(v > 3) v=1000]");  
		run("Macro...", "code=[if(v < 1/3) v=0.001]");  

		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");



		red_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";		
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);
		

		
		green_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";	
		run("Image Expression Parser (Macro)", "expression=["+green_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		green_image_id = getImageID();
		rename("green image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(2);
		imageCalculator("Multiply", rgb_image_id, green_image_id);



		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, red_image_id);


		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(green_image_id);
		close();

		selectImage(temp_ratio_image_id);
		close();
	}

	if(color_scheme == "magenta-yellow"){
		color_scheme_suffix = "_MY.tif";
		selectImage(ratio_image_id);
		run("Duplicate...", "title=temp_ratio_image");
		temp_ratio_image_id = getImageID();
		run("Macro...", "code=[if(v > 3) v=1000]");  
		run("Macro...", "code=[if(v < 1/3) v=0.001]");  

		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");



		red_expression = "0 * A + 255";		
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);
		

		
		green_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";	
		run("Image Expression Parser (Macro)", "expression=["+green_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		green_image_id = getImageID();
		rename("green image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(2);
		imageCalculator("Multiply", rgb_image_id, green_image_id);



		blue_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";	
		run("Image Expression Parser (Macro)", "expression=["+blue_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		blue_image_id = getImageID();
		rename("blue image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, blue_image_id);


		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(green_image_id);
		close();
		selectImage(blue_image_id);
		close();

		selectImage(temp_ratio_image_id);
		close();
	
	}


	if(color_scheme == "cyan-yellow"){
		color_scheme_suffix = "_CY.tif";
		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");

		red_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";	
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);


		
		selectImage(rgb_image_id);
		setSlice(2);
		run("Set...", "value=255");

		blue_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";	
		run("Image Expression Parser (Macro)", "expression=["+ blue_expression +"] a=ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		blue_image_id = getImageID();
		rename("blue image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, blue_image_id);


		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(blue_image_id);
		close();

	}


	if(color_scheme == "red-blue"){
		color_scheme_suffix = "_RB.tif";
		selectImage(ratio_image_id);
		run("Duplicate...", "title=temp_ratio_image");
		temp_ratio_image_id = getImageID();
		run("Macro...", "code=[if(v > 3) v=1000]");  
		run("Macro...", "code=[if(v < 1/3) v=0.001]");  

		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");



		red_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";		
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);
		

		
		blue_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";	
		run("Image Expression Parser (Macro)", "expression=["+blue_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		blue_image_id = getImageID();
		rename("blue image");
		run("Macro...", "code=[if(v > 255) v=255]");  
		
		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, blue_image_id);



		imageCalculator("Subtract create 32-bit", red_image_id, blue_image_id);
		green_image_id = getImageID();
		rename("green image");
		run("Abs");
		run("Multiply...", "value=-1");
		run("Add...", "value=255");

		selectImage(rgb_image_id);
		setSlice(2);
		imageCalculator("Multiply", rgb_image_id, green_image_id);


		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(green_image_id);
		close();
		selectImage(blue_image_id);
		close();

		selectImage(temp_ratio_image_id);
		close();

	}


	if(color_scheme == "yellow-blue"){
		color_scheme_suffix = "_YB.tif";
		selectImage(ratio_image_id);
		run("Duplicate...", "title=temp_ratio_image");
		temp_ratio_image_id = getImageID();
		run("Macro...", "code=[if(v > 3) v=1000]");  
		run("Macro...", "code=[if(v < 1/3) v=0.001]");  

		selectImage(rgb_image_id);
		run("RGB Stack");
		run("Set...", "value=1 stack");


		red_expression = "255 * (A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+")";		
		run("Image Expression Parser (Macro)", "expression=["+ red_expression +"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		red_image_id = getImageID();
		rename("red image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(1);
		imageCalculator("Multiply", rgb_image_id, red_image_id);
		

		selectImage(rgb_image_id);
		setSlice(2);
		imageCalculator("Multiply", rgb_image_id, red_image_id);


		blue_expression = "255 / ((A + "+toString(scale_factor)+")/(1 + "+toString(scale_factor)+"))";	
		run("Image Expression Parser (Macro)", "expression=["+blue_expression+"] a=temp_ratio_image b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		blue_image_id = getImageID();
		rename("blue image");
		run("Macro...", "code=[if(v > 255) v=255]");  

		selectImage(rgb_image_id);
		setSlice(3);
		imageCalculator("Multiply", rgb_image_id, blue_image_id);

		imageCalculator("Multiply stack", rgb_image_id, avg_image_id);
		run("RGB Color");

		selectImage(red_image_id);
		close();
		selectImage(blue_image_id);
		close();

		selectImage(temp_ratio_image_id);
		close();
		
	
	}



	if(color_scheme == "rainbow"){
		color_scheme_suffix = "_rainbow.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);
				new_ratio_pixel_value= getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){new_ratio_pixel_value = 1;}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);	
				
				hue = round(log(new_ratio_pixel_value)/log(3)*120 + 120);
				if(hue < 10){hue = 0;}
				if(hue > 230){hue = 240;}

		
				hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
				if(hue >= 0 && hue < 60){
					red = 255;
					green = round(hsv_x * 255);
					blue = 0;
				}
		
				if(hue >= 60 && hue < 120){
					red = round(hsv_x * 255);
					green = 255;
					blue = 0;
				}

				if(hue == 120){
					red = 255;
					green = 255;
					blue = 255;
				}					
		
				if(hue > 120 && hue < 180){
					red = 0;
					green = 255;
					blue = round(hsv_x * 255);
				}
				if(hue >= 180 && hue < 240){
					red = 0;
					green = round(hsv_x * 255);
					blue = 255;
				}
		
				if(hue >= 240 && hue < 300){
					red = round(hsv_x * 255);
					green = 0;
					blue = 255;
				}
				
				if(hue >= 300 && hue < 360){
					red = 255;
					green = 0;
					blue = round(hsv_x * 255);
				}

				
				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}		
	}


//"white-red","white-green","white-blue","white-cyan","white-magenta"



	if(color_scheme == "white-red"){
		color_scheme_suffix = "_WR.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);  //   Fmax/Fmin
				new_ratio_pixel_value= 1/getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){
					new_ratio_pixel_value = 1;
				}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);				
						
				if(new_ratio_pixel_value <= 1){
					red = 255;
					green = 255;
					blue = 255;					
				}

				if(new_ratio_pixel_value > 1 && new_ratio_pixel_value < 3){
					red   = 255;
					green = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					blue  = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
				}
				if(new_ratio_pixel_value >= 3){
					red = 255;
					green = 0;
					blue = 0;
				}

				// rgb_intensity = red + green + blue;
				rgb_intensity = 255;
				red = red * 255/rgb_intensity;
				green = green * 255/rgb_intensity;
				blue = blue * 255/rgb_intensity;

				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}	
	}




	if(color_scheme == "white-green"){
		color_scheme_suffix = "_WG.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);  //   Fmax/Fmin
				new_ratio_pixel_value= 1/getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){
					new_ratio_pixel_value = 1;
				}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);				
						
				if(new_ratio_pixel_value <= 1){
					red = 255;
					green = 255;
					blue = 255;					
				}

				if(new_ratio_pixel_value > 1 && new_ratio_pixel_value < 3){
					red   = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					green = 255;
					blue  = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
				}
				if(new_ratio_pixel_value >= 3){
					red = 0;
					green = 255;
					blue = 0;
				}

				// rgb_intensity = red + green + blue;
				rgb_intensity = 255;
				red = red * 255/rgb_intensity;
				green = green * 255/rgb_intensity;
				blue = blue * 255/rgb_intensity;

				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}	
	}



	if(color_scheme == "white-blue"){
		color_scheme_suffix = "_WB.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);  //   Fmax/Fmin
				new_ratio_pixel_value= 1/getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){
					new_ratio_pixel_value = 1;
				}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);				
						
				if(new_ratio_pixel_value <= 1){
					red = 255;
					green = 255;
					blue = 255;					
				}

				if(new_ratio_pixel_value > 1 && new_ratio_pixel_value < 3){
					red   = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					green = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					blue  = 255;
				}
				if(new_ratio_pixel_value >= 3){
					red = 0;
					green = 0;
					blue = 255;
				}

				// rgb_intensity = red + green + blue;
				rgb_intensity = 255;
				red = red * 255/rgb_intensity;
				green = green * 255/rgb_intensity;
				blue = blue * 255/rgb_intensity;

				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}	
	}



	if(color_scheme == "white-cyan"){
		color_scheme_suffix = "_WC.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);  //   Fmax/Fmin
				new_ratio_pixel_value= 1/getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){
					new_ratio_pixel_value = 1;
				}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);				
						
				if(new_ratio_pixel_value <= 1){
					red = 255;
					green = 255;
					blue = 255;					
				}

				if(new_ratio_pixel_value > 1 && new_ratio_pixel_value < 3){
					red   = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					green = 255;
					blue  = 255;
				}
				if(new_ratio_pixel_value >= 3){
					red = 0;
					green = 255;
					blue = 255;
				}

				// rgb_intensity = red + green + blue;
				rgb_intensity = 255;
				red = red * 255/rgb_intensity;
				green = green * 255/rgb_intensity;
				blue = blue * 255/rgb_intensity;

				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}	
	}




	if(color_scheme == "white-magenta"){
		color_scheme_suffix = "_WM.tif";
		for (row= 0; row<h; row++) {
			for (col = 0; col < w; col++){
				selectImage(ratio_image_id);  //   Fmax/Fmin
				new_ratio_pixel_value= 1/getPixel(col,row);
				if(isNaN(new_ratio_pixel_value)){
					new_ratio_pixel_value = 1;
				}
						
				selectImage(avg_image_id);
				intensity_factor = getPixel(col,row);				
						
				if(new_ratio_pixel_value <= 1){
					red = 255;
					green = 255;
					blue = 255;					
				}

				if(new_ratio_pixel_value > 1 && new_ratio_pixel_value < 3){
					red   = 255;
					green = (3 - new_ratio_pixel_value*0.8)/2.2 * 255;
					blue  = 255;
				}
				if(new_ratio_pixel_value >= 3){
					red = 255;
					green = 0;
					blue = 255;
				}

				// rgb_intensity = red + green + blue;
				rgb_intensity = 255;
				red = red * 255/rgb_intensity;
				green = green * 255/rgb_intensity;
				blue = blue * 255/rgb_intensity;

				red = round(red*intensity_factor);
				green = round(green*intensity_factor);
				blue = round(blue*intensity_factor);
																	
				setColor(red,green,blue);

				selectImage(rgb_image_id);
				drawLine(col,row,col,row);


			}
		}	
	}




	rg_image_name = replace(processed_composite_image_name, "_pcomp.tif", "_range_") + rg_range + color_scheme_suffix;
	selectImage(rgb_image_id);
	rename(rg_image_name);



	selectImage(ratio_image_id);
	close();
	selectImage(avg_image_id);
	close();
	selectImage(hor_image_id);
	close();
	selectImage(ver_image_id);
	close();

	selectImage(rgb_image_id);




	updateDisplay();
	return(rgb_image_id);
}








function correct_negative_values(image_id){
	selectImage(image_id);
	number_of_slices = nSlices;

	for(slice_number = 1; slice_number < number_of_slices; slice_number++){
		setSlice(slice_number);
		for (line= 0; line<h; line++) {
			for (col = 0; col < w; col++){
				if (getPixel(col, line) > 65000){
					setPixel(col, line, 65535-getPixel(col, line));
				}
			}		 
		}
		updateDisplay();
	}
}


function correct_bleedthrough(image_id, bleedthrough){
	selectImage(image_id);
	for (line= 0; line<h; line++) {
		for (col = 1; col < w; col++){
			setPixel(col, line, getPixel(col, line) - bleedthrough * getPixel((col-1),line));
			setPixel(col-1, line, getPixel((col-1), line) * (1 + bleedthrough));
		}		 
	}
	updateDisplay();
}


function make_image_from_columns(image_id, offset, name){
	selectImage(image_id);
	run("Duplicate...", "title=column_image");
	column_image_id = getImageID();
	
	w=getWidth;
	h=getHeight;
	//set values of unwanted pixels to 0
	for (col= 1-offset; col<w; col+=2) {
		for (line = 0; line < h; line++){
				setPixel(col, line, getPixel((col-1),line));
		}		 
	}

	updateDisplay();
	run("Size...", "width="+w/2+ " height="+h+" interpolation=None");
	run("Multiply...", "value=2");
	column_image_id = getImageID();
	
	name = name+column_image_id;
	rename(name);
	return(column_image_id);
}


function get_main_h_dialog_values(h_image_id, t_bkgd_intensity){
		selectImage(h_image_id);
		run("Select All");
		run("Duplicate...", " ");
		run("Subtract...", "value="+t_bkgd_intensity);
		bkgd_subbed_image_id = getImageID();
		setBatchMode("show");

		Dialog.create("Image processing parameters");

		Dialog.setInsets(10, 10, 10);
		Dialog.addMessage("Enter/select parameters for processing an image\nacquired using a single polarization of excitation light:");
		Dialog.addNumber("Background intensity:", t_bkgd_intensity,0,5,"");
	
		Dialog.setInsets(10, 20, 0);
		Dialog.addMessage("___________________________________________");	
		Dialog.setInsets(0, 20, 0);
		Dialog.addMessage(" \nMicroscope parameters:");
	
		Dialog.setInsets(0, 0, 0);
		Dialog.addChoice("Polarization direction:",newArray("horizontal","vertical"), polarization_direction);
		Dialog.setInsets(0, 0, 0);
		Dialog.addNumber("Polarization purity:", polarization_purity,2,6,"%");
		Dialog.setInsets(0, 30, 0);
		Dialog.addCheckbox("Set deviation from hor./ver. directions manually", phase_fixed);
		Dialog.setInsets(0, 0, 0);
		Dialog.addNumber("Deviation ["+degreechar+"]:", phase1_deg,2,6,"");
	
	
		Dialog.setInsets(10, 20, 0);
		Dialog.addMessage("___________________________________________");	
		Dialog.setInsets(0, 20, 0);
		Dialog.addMessage(" \nPerform quantitative analysis:");
		Dialog.setInsets(5, 0, 0)
		Dialog.addChoice("Analysis type:", newArray("1PPM: r(max)", "2PPM: r(max), "+alpha0char+", "+sigmachar), determine_rmax_alpha0_sigma);
		Dialog.setInsets(0, 0, 0);
		if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um" || pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm" || pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
			Dialog.addNumber("Cell/vesicle outline width:", cutoff_distance_from_snake,3,5,muchar+"m");
		}
		else{
			Dialog.addNumber("Cell/vesicle outline width:", cutoff_distance_from_snake,3,5,"pixels");
		}
		Dialog.setInsets(0, 20, 0);
		ratiometric_fitting_choices = newArray("non-ratiometric", "ratiometric");
		Dialog.addRadioButtonGroup("Data fitting:", ratiometric_fitting_choices, 1, 2, ratiometric_fitting_choices[ratiometric_fitting]);
		Dialog.setInsets(0, 30, 0);
		Dialog.addCheckbox("Attempt automated vesicle recognition", look_for_vesicles);

		
		Dialog.setInsets(0, 30, 0);		
		Dialog.addCheckbox("Close original image when done", orig_image_close);
		Dialog.setInsets(0, 30, 0);
		Dialog.addCheckbox("Save result files", save_files);

		
		Dialog.show();
	
		new_bkgd_intensity = Dialog.getNumber();
		polarization_direction = Dialog.getChoice();
		polarization_purity = Dialog.getNumber;
		polarization_purity_hor = polarization_purity;
		polarization_purity_ver = polarization_purity;		
		phase_fixed = Dialog.getCheckbox();
		if(phase_fixed == true){look_for_pol_angles = 0;}
		else{look_for_pol_angles = 1;}
		phase1_deg = Dialog.getNumber();
		phase2_deg = phase1_deg;
		
		determine_rmax_alpha0_sigma = Dialog.getChoice();
		cutoff_distance_from_snake = Dialog.getNumber();
		ratiometric_fitting = Dialog.getRadioButton();
		if(ratiometric_fitting == "ratiometric"){ratiometric_fitting = true;}
		else{ratiometric_fitting = false;}
		look_for_vesicles = Dialog.getCheckbox();
		orig_image_close = Dialog.getCheckbox();
		save_files = Dialog.getCheckbox();

		if(new_bkgd_intensity == t_bkgd_intensity){check_with_user = false;}
		else{
			bkgd_intensity = new_bkgd_intensity;
			check_with_user = true;
		}

		dialog_values_text = "";
		dialog_values_text = dialog_values_text+"bkgd_subbed_image_id:"+bkgd_subbed_image_id+"\n";
		dialog_values_text = dialog_values_text+"bkgd_intensity:"+bkgd_intensity+"\n";
		dialog_values_text = dialog_values_text+"polarization_direction:"+polarization_direction+"\n";
		dialog_values_text = dialog_values_text+"polarization_purity:"+polarization_purity+"\n";
		dialog_values_text = dialog_values_text+"determine_rmax_alpha0_sigma:"+determine_rmax_alpha0_sigma+"\n";
		dialog_values_text = dialog_values_text+"cutoff_distance_from_snake:"+cutoff_distance_from_snake+"\n";
		dialog_values_text = dialog_values_text+"look_for_vesicles:"+look_for_vesicles+"\n";
		dialog_values_text = dialog_values_text+"phase_fixed:"+phase_fixed+"\n";
		dialog_values_text = dialog_values_text+"phase1_deg:"+phase1_deg+"\n";
		dialog_values_text = dialog_values_text+"ratiometric_fitting:"+ratiometric_fitting+"\n";
		dialog_values_text = dialog_values_text+"orig_image_close:"+orig_image_close+"\n";
		dialog_values_text = dialog_values_text+"save_files:"+save_files+"\n";
		dialog_values_text = dialog_values_text+"check_with_user:"+check_with_user+"\n";

		return dialog_values_text;
}



function make_BIN_image(t_segmentation_image_id, t_bkgd_subbed_image_id, t_cutoff_distance_from_snake_in_pixels){
	selectImage(t_segmentation_image_id);
	getSelectionCoordinates(t_snake_points_x, t_snake_points_y);	// x, y coordinates of the spline
	getPixelSize(t_pixel_size_unit, t_pixel_width, t_pixel_height);
	run("Select None");

	newImage("snake_ROI", "8-bit black", w, h, 1);
	t_snake_ROI_image_id = getImageID();
	setVoxelSize(t_pixel_width, t_pixel_height, 1, t_pixel_size_unit);
	makeSelection("polyline", t_snake_points_x, t_snake_points_y);
	setForegroundColor(2);
	Roi.setStrokeWidth(t_cutoff_distance_from_snake_in_pixels);
	run("Draw", "slice");
	run("Select None");

	selectImage(t_segmentation_image_id);
	if(bitDepth() == 24){
		// number_of_polarizations = 2;

		selectImage(t_segmentation_image_id);
		run("Duplicate...", "title=temp_16_bit_segmentation_image");
		run("HSB Stack");
		run("Delete Slice");
		run("Delete Slice");
		run("16-bit");
		temp_16_bit_segmentation_image_id = getImageID();

		selectImage(t_bkgd_subbed_image_id);
		run("Duplicate...", "title=temp_16_bit_bkgd_subbed_image");
		run("HSB Stack");
		run("Delete Slice");
		run("Delete Slice");
		run("16-bit");
		temp_16_bit_bkgd_subbed_image_id = getImageID();
		
		imageCalculator("Difference", temp_16_bit_segmentation_image_id, temp_16_bit_bkgd_subbed_image_id);
		getStatistics(area, mean, min, max);
		run("Macro...", "code=[if(v == 0) v=255; else v=0]");  // non-modified parts = 255, parts to be rejected = 0

		imageCalculator("Multiply", temp_16_bit_segmentation_image_id, t_snake_ROI_image_id);
		imageCalculator("Divide create 32-bit", temp_16_bit_segmentation_image_id, temp_16_bit_segmentation_image_id);
		BIN_image_id = getImageID();
		rename("bin");

//

		run("Line Width...", "line=1"); 

		make_HSV_snake_overlay(BIN_image_id, t_snake_points_x, t_snake_points_y);

		run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 


//
		
		makeSelection("polyline", t_snake_points_x, t_snake_points_y);
//		Overlay.addSelection("orange");	




		selectImage(temp_16_bit_segmentation_image_id);
		close();
		selectImage(temp_16_bit_bkgd_subbed_image_id);
		close();

	}
	else{
		// number_of_polarizations = 1;
		
		temp_16_bit_snake_image_id = 0;

		imageCalculator("Difference", t_segmentation_image_id, t_bkgd_subbed_image_id);
		run("Macro...", "code=[if(v == 0) v=255; else v=0]");  // non-modified parts = 255, parts to be rejected = 0

		imageCalculator("Multiply", t_segmentation_image_id, t_snake_ROI_image_id);
		imageCalculator("Divide create 32-bit", t_segmentation_image_id, t_segmentation_image_id);
		BIN_image_id = getImageID();
		rename("bin");
		makeSelection("polygon", t_snake_points_x, t_snake_points_y);
		Overlay.addSelection("orange");	
	
	
	}
		
	selectImage(t_snake_ROI_image_id);
	close();

	selectImage(BIN_image_id);		

	
	
	return(BIN_image_id);

}


function make_theta_array(t_snake_points_x, t_snake_points_y){
	theta = newArray(t_snake_points_x.length);
	//vertex_angles = Array.getVertexAngles(t_snake_points_x, t_snake_points_y, 1);

	//theta[1] = atan( (t_snake_points_x[2]-t_snake_points_x[0]) / (t_snake_points_y[2]-t_snake_points_y[0]) );
	//for (i=2; i<t_snake_points_x.length-1; i++){
	//	theta[i] = theta[i-1] + vertex_angles[i]/180*PI;
	//}
	
	
	//-------------------------------------
	for (i=1; i<t_snake_points_x.length-1; i++){
		theta[i] = atan( (t_snake_points_x[i+1]-t_snake_points_x[i-1]) / (t_snake_points_y[i+1]-t_snake_points_y[i-1]) );
	}
	//-------------------------------------

	// is the snake a loop? Are the endpoints > 5x farther apart than points #0 and #1?

	if( (t_snake_points_x[0]-t_snake_points_x[t_snake_points_x.length-1])*(t_snake_points_x[0]-t_snake_points_x[t_snake_points_x.length-1]) + (t_snake_points_y[0]-t_snake_points_y[t_snake_points_x.length-1])*(t_snake_points_y[0]-t_snake_points_y[t_snake_points_x.length-1]) > 25*( (t_snake_points_x[0]-t_snake_points_x[1])*(t_snake_points_x[0]-t_snake_points_x[1]) + (t_snake_points_y[0]-t_snake_points_y[1])*(t_snake_points_y[0]-t_snake_points_y[1]))  ) {
		// snake is open-ended
		theta[0] = theta[1] - (theta[2]-theta[1]);
		theta[t_snake_points_x.length-1] = theta[t_snake_points_x.length-2] - (theta[t_snake_points_x.length-3]-theta[t_snake_points_x.length-2]);
	}
	else{
		// snake makes a loop
		theta[0] = atan( (t_snake_points_x[1]-t_snake_points_x[t_snake_points_x.length-1]) / (t_snake_points_y[1]-t_snake_points_y[t_snake_points_x.length-1]) );
		theta[t_snake_points_x.length-1] = atan( (t_snake_points_x[0]-t_snake_points_x[t_snake_points_x.length-2]) / (t_snake_points_y[0]-t_snake_points_y[t_snake_points_x.length-2]) );
	}	

	return theta;
}



function t_make_a_tiled_stack(B_C_b_major_stack_id, num_of_alpha0s_minor, num_of_sigmas_minor){

	selectImage(B_C_b_major_stack_id);
	num_of_slices = nSlices;
	num_of_alpha0s_major = getWidth;
	num_of_sigmas_major  = getHeight;


	newImage("B_C_b_tiled_stack", "32-bit black", num_of_alpha0s_major*num_of_alpha0s_minor, num_of_sigmas_major*num_of_sigmas_minor, num_of_slices);
	tiled_stack_id = getImageID();

	for(i=1; i<=num_of_slices; i++){
		selectImage(B_C_b_major_stack_id);
		setSlice(i);
		run("Select All");
		run("Copy");

		selectImage(tiled_stack_id);
		setSlice(i);
		run("Specify...", "width="+num_of_alpha0s_major+" height="+num_of_sigmas_major+" x=0 y=0 slice="+i);
		run("Paste");

		last_column = num_of_alpha0s_major;

		do{
			run("Specify...", "width="+last_column+" height="+num_of_sigmas_major+" x=0 y=0 slice="+i);
			run("Copy");
						
			run("Specify...", "width="+last_column+" height="+num_of_sigmas_major+" x="+last_column+" y=0 slice="+i);
			run("Paste");
			last_column = last_column + last_column;
		}while(last_column < num_of_alpha0s_major * num_of_alpha0s_minor)

		last_row = num_of_sigmas_major;
		do{
			run("Specify...", "width="+toString(num_of_alpha0s_minor*num_of_alpha0s_major)+" height="+last_row+" x=0 y=0 slice="+i);
			run("Copy");
						
			run("Specify...", "width="+toString(num_of_alpha0s_minor*num_of_alpha0s_major)+" height="+last_row+" x=0 y="+last_row+" slice="+i);
			run("Paste");
			last_row = last_row + last_row;
		}while(last_row < num_of_sigmas_major * num_of_sigmas_minor)

		run("Select None");
		
	}


	return(tiled_stack_id);
}





//////////////////////////////////////////////////////////////////////////////


function t_make_a_stretched_stack(B_C_b_minor_stack_id, num_of_alpha0s_major, num_of_sigmas_major){  //should be correct and fast... work in progress

	selectImage(B_C_b_minor_stack_id);
	run("Duplicate...", "title=B_C_b_stretched_minor_stack duplicate");
	stretched_tiled_image_id = getImageID();
	
	num_of_alpha0s_minor = getWidth;
	num_of_sigmas_minor  = getHeight;

	run("Size...", "width="+toString(num_of_alpha0s_minor*num_of_alpha0s_major)+" height="+toString(num_of_sigmas_minor*num_of_sigmas_major)+" interpolation=None");

	return(stretched_tiled_image_id);
}


//////////////////////////////////////////////////////////////////////////////////////


function t_format_output(t_major_fraction_percentage, t_alpha0_major, t_sigma_major, t_minor_fraction_percentage, t_alpha0_minor, t_sigma_minor, t_min_RMSD, t_entropy, t_predicted_rmax_1PPM, t_predicted_rmax_2PPM, t_predicted_b_1PPM, t_predicted_B_2PPM, t_predicted_C_2PPM){
		
	major_fraction_percentage_string = toString(t_major_fraction_percentage*100)+"%";
	while(lengthOf(major_fraction_percentage_string) < 4){
		major_fraction_percentage_string = "  "+major_fraction_percentage_string;
	}
	alpha0_major_string = toString(t_alpha0_major) + degreechar;
	while(lengthOf(alpha0_major_string) < 3){
		alpha0_major_string = "  "+alpha0_major_string;
	}
	sigma_major_string = toString(t_sigma_major) + degreechar;
	while(lengthOf(sigma_major_string) < 3){
		sigma_major_string = "  "+sigma_major_string;
	}

	minor_fraction_percentage_string = toString(t_minor_fraction_percentage*100)+"%";
	while(lengthOf(minor_fraction_percentage_string) < 5){
		minor_fraction_percentage_string = "  "+minor_fraction_percentage_string;
	}
	alpha0_minor_string = toString(t_alpha0_minor) + degreechar;
	while(lengthOf(alpha0_minor_string) < 3){
		alpha0_minor_string = "  "+alpha0_minor_string;
	}
	sigma_minor_string = toString(t_sigma_minor) + degreechar;
	while(lengthOf(sigma_minor_string) < 3){
		sigma_minor_string = "  "+sigma_minor_string;
	}

	min_RMSD_string = toString(t_min_RMSD);
	while(lengthOf(min_RMSD_string) < 10){
		min_RMSD_string = min_RMSD_string+"0";
	}
	
	entropy_string = toString(t_entropy);
	if(entropy_string == "-Infinity"){
		entropy_string = " "+entropy_string;
	}
	if(t_entropy == t_entropy){
		while(lengthOf(entropy_string) < 7){
			entropy_string = entropy_string+"0";
		}	
	}
	else{
		entropy_string = "NaN     ";
	}

	predicted_rmax_1PPM_string = toString(t_predicted_rmax_1PPM);
	while(lengthOf(predicted_rmax_1PPM_string) < 6){
		predicted_rmax_1PPM_string = predicted_rmax_1PPM_string+"0";
	}
	
	predicted_rmax_2PPM_string = toString(t_predicted_rmax_2PPM);
	while(lengthOf(predicted_rmax_2PPM_string) < 6){
		predicted_rmax_2PPM_string = predicted_rmax_2PPM_string+"0";
	}

	predicted_b_1PPM_string = toString(t_predicted_b_1PPM);
	while(lengthOf(predicted_b_1PPM_string) < 6){
		predicted_b_1PPM_string = predicted_b_1PPM_string+"0";
	}

	predicted_B_2PPM_string = toString(t_predicted_B_2PPM);
	while(lengthOf(predicted_B_2PPM_string) < 6){
		predicted_B_2PPM_string = predicted_B_2PPM_string+"0";
	}
	
	predicted_C_2PPM_string = toString(t_predicted_C_2PPM);
	while(lengthOf(predicted_C_2PPM_string) < 9){
		predicted_C_2PPM_string = predicted_C_2PPM_string+"0";
	}

	output_t_string = major_fraction_percentage_string+"  "+alpha0_major_string+", "+sigma_major_string+"     "+minor_fraction_percentage_string+"  "+alpha0_minor_string+", "+sigma_minor_string+"   "+min_RMSD_string+"   "+entropy_string+"    "+predicted_rmax_1PPM_string+"    "+predicted_rmax_2PPM_string+"    "+predicted_b_1PPM_string+"     "+predicted_B_2PPM_string+"     "+predicted_C_2PPM_string;
	return output_t_string;

}



function t_make_a_sum_of_B_C_b_major_minor_values_stack(t_B_C_b_tiled_stack_id, t_B_C_b_stretched_stack_id, t_minor_fraction_percentage){
	t_major_fraction_percentage = 1 - t_minor_fraction_percentage;
	selectImage(t_B_C_b_tiled_stack_id);
	run("Select None");		
	run("Duplicate...", "title=t_major_B_C_b_tiled_stack duplicate");
	rename("major_B_C_b_tiled_stack"+toString(t_minor_fraction_percentage));
	t_major_B_C_b_tiled_stack_id = getImageID();
	run("Multiply...", "value="+t_major_fraction_percentage+" stack");
	
	selectImage(t_B_C_b_stretched_stack_id);
	run("Select None");
	run("Duplicate...", "title=t_minor_B_C_b_tiled_stack duplicate");
	t_minor_B_C_b_tiled_stack_id = getImageID();
	rename("t_minor_B_C_b_tiled_stack"+toString(t_minor_fraction_percentage));
	run("Multiply...", "value="+t_minor_fraction_percentage+" stack");

	imageCalculator("Add 32-bit stack", t_major_B_C_b_tiled_stack_id, t_minor_B_C_b_tiled_stack_id);	//the result now has the image_id of major_B_C_b....
	rename("t_sum_of_B_C_b_major_minor_values"+toString(t_minor_fraction_percentage));
	t_sum_of_B_C_b_major_minor_values_stack_id = getImageID();
	
	selectImage(t_minor_B_C_b_tiled_stack_id);
	close();

	return t_sum_of_B_C_b_major_minor_values_stack_id;
}


function t_make_rmax_1PPM_image(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b)/(1-b) = 1 + 2b/(1-b)

	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);
		
	setSlice(3);
	run("Duplicate...", "title=[l2rmax1P_pred]");								// b
	t_rmax_1P_image_id = getImageID();


	run("Multiply...", "value=-1");												// -b
	run("Add...", "value=1");													// 1 - b
	run("Reciprocal");															// 1/(1 - b)
																				// b/(1-b)
	imageCalculator("Multiply 32-bit", t_rmax_1P_image_id, t_sum_of_B_C_b_major_minor_values_stack_id);	
	run("Multiply...", "value=2");												// 2b/(1-b)
	run("Add...", "value=1");													// 1 + 2b/(1-b) = Rmax_1P

	return t_rmax_1P_image_id;
	
}

function t_make_rmax_1PPM_image_15(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b*0.866)/(1-b*0.866) = 1 + 2*0.866b/(1-0.866b)

	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);
	setSlice(3);
	run("Duplicate...", "title=[l2rmax1P_pred]");								// b
	t_rmax_1P_image_id = getImageID();

	//run("Macro...", "code=[v = 1 + 2*v/(1-v)]"); 		//~2x slower than the code below

	run("Multiply...", "value=-0.866");												// -0.866b
	run("Add...", "value=1");													// 1 - 0.866b
	run("Reciprocal");															// 1/(1 - 0.866b)
																				// b/(1-b)
	imageCalculator("Multiply 32-bit", t_rmax_1P_image_id, t_sum_of_B_C_b_major_minor_values_stack_id);	
	run("Multiply...", "value=1.732");												// 2*0.866b/(1-0.866b)
	run("Add...", "value=1");													// 1 + 2*0.866b/(1-0.866b) = r @ 15 degrees

	return t_rmax_1P_image_id;
	
}

function t_make_rmax_1PPM_image_30(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b*0.5)/(1-b*0.5) = 1 + b/(1-0.5b)

	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);
	setSlice(3);
	run("Duplicate...", "title=[l2rmax1P_pred]");								// b
	t_rmax_1P_image_id = getImageID();

	//run("Macro...", "code=[v = 1 + 2*v/(1-v)]"); 		//~2x slower than the code below

	run("Multiply...", "value=-0.5");												// -0.5b
	run("Add...", "value=1");													// 1 - 0.5b
	run("Reciprocal");															// 1/(1 - 0.5b)
																				// b/(1-b)
	imageCalculator("Multiply 32-bit", t_rmax_1P_image_id, t_sum_of_B_C_b_major_minor_values_stack_id);	
	//run("Multiply...", "value=1.732");												// 2*0.866b/(1-0.866b)
	run("Add...", "value=1");													// 1 + b/(1-0.5b) = r @ 30 degrees

	return t_rmax_1P_image_id;
	
}


function t_make_rmax_2PPM_image(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b+c)/(1-b+c) = 1 + 2b/(1-b+c)


	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);

	run("Duplicate...", "duplicate");
	temp_stack_image_id = getImageID();

	setSlice(3);
	run("Set...", "value=1 slice");															// set B1P slice to 1, so that it works with the following steps 
	setSlice(1);																			// select image B2P
	run("Duplicate...", "title=[t_predicted B]");
	t_predicted_B_2PPM_image_id = getImageID();
	selectImage(temp_stack_image_id);
	setSlice(1);																			// select image B2P
	run("Multiply...", "value=-1 slice");													// -B2P
	run("Z Project...", "projection=[Sum Slices]");		//makes a new image					// 1 - B2P + C2P
	rename("t_rmax2P_pred");
	t_rmax_2P_image_id = getImageID();
	run("Reciprocal");																		// 1/(1 - B2P + C2P)
	imageCalculator("Multiply 32-bit", t_rmax_2P_image_id, t_predicted_B_2PPM_image_id);	// B2P/(1 - B2P + C2P)
	run("Multiply...", "value=2");															// 2 * B2P/(1 - B2P + C2P)
	run("Add...", "value=1");																// 1 + 2 * B2P/(1 - B2P + C2P) = Rmax

	selectImage(t_predicted_B_2PPM_image_id);
	close();

	selectImage(temp_stack_image_id);
	close();

	
	return(t_rmax_2P_image_id);
	
}


function t_make_rmax_2PPM_image_15(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b+c)/(1-b+c) = 1 + 2b/(1-b+c)
	//	rmax = (1+0.866b+0.5c)/(1-0.866b+0.5c) = 1 + 1.732b/(1-0.866b+0.5c)


	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);

	run("Duplicate...", "duplicate");
	temp_stack_image_id = getImageID();
	setSlice(3);
	run("Set...", "value=1 slice");															// set B1P slice to 1, so that it works with the following steps 
	setSlice(1);																			// select image B2P
	run("Duplicate...", "title=[t_predicted B]");
	t_predicted_B_2PPM_image_id = getImageID();
	selectImage(temp_stack_image_id);
	setSlice(1);																			// select image B2P
	run("Multiply...", "value=-0.866 slice");													// -0.866 B2P
	setSlice(2);																			// select image B2P
	run("Multiply...", "value=0.5 slice");													// 0.5 C2P
	run("Z Project...", "projection=[Sum Slices]");		//makes a new image					// 1 -0.866 B2P + 0.5 C2P
	rename("t_rmax2P_pred");
	t_rmax_2P_image_id = getImageID();
	run("Reciprocal");																		// 1/(1 -0.866 B2P + 0.5 C2P)
	imageCalculator("Multiply 32-bit", t_rmax_2P_image_id, t_predicted_B_2PPM_image_id);	// B2P/(1 -0.866 B2P + 0.5 C2P)
	run("Multiply...", "value=1.732");															// 1.732 * B2P/(1 -0.866 B2P + 0.5 C2P)
	run("Add...", "value=1");																// 1 + 1.732 * B2P/(1 -0.866 B2P + 0.5 C2P) = R @ 15 degrees

	selectImage(t_predicted_B_2PPM_image_id);
	close();

	selectImage(temp_stack_image_id);
	close();

	
	return(t_rmax_2P_image_id);
	
}


function t_make_rmax_2PPM_image_30(t_sum_of_B_C_b_major_minor_values_stack_id){

	//	rmax = (1+b+c)/(1-b+c) = 1 + 2b/(1-b+c)
	//	rmax = (1+0.5b-0.5c)/(1-0.5b-0.5c) = 1 + b/(1-0.5b-0.5c)


	selectImage(t_sum_of_B_C_b_major_minor_values_stack_id);
	run("Duplicate...", "duplicate");
	temp_stack_image_id = getImageID();

	setSlice(3);
	run("Set...", "value=1 slice");															// set B1P slice to 1, so that it works with the following steps 
	setSlice(1);																			// select image B2P
	run("Duplicate...", "title=[t_predicted B]");
	t_predicted_B_2PPM_image_id = getImageID();
	
	selectImage(temp_stack_image_id);
	setSlice(1);																			// select image B2P
	run("Multiply...", "value=-0.5 slice");													// -0.5 B2P
	setSlice(2);																			// select image B2P
	run("Multiply...", "value=-0.5 slice");													// -0.5 C2P
	run("Z Project...", "projection=[Sum Slices]");		//makes a new image					// 1 -0.5 B2P - 0.5 C2P
	rename("t_rmax2P_pred");
	t_rmax_2P_image_id = getImageID();
	run("Reciprocal");																		// 1/(1 -0.5 B2P - 0.5 C2P)
	imageCalculator("Multiply 32-bit", t_rmax_2P_image_id, t_predicted_B_2PPM_image_id);	// B2P/(1 -0.5 B2P - 0.5 C2P)
	run("Add...", "value=1");																// 1 + B2P/(1 -0.5 B2P - 0.5 C2P) = R @ 30 degrees

	selectImage(t_predicted_B_2PPM_image_id);
	close();

	selectImage(temp_stack_image_id);
	close();

	
	return(t_rmax_2P_image_id);
	
}



function find_bkgd_intensity(image_id_f){
	selectImage(image_id_f);
	w=getWidth;
	h=getHeight;
	getStatistics(area, mean, min, max);
	bkgd_intensity = max;
	min_x = 0;
	min_y = 0;
	for(x = w/20; x < w - w/10; x = x + w/20){
		for(y = h/20; y < h - h/10; y = y + h/20){
			makeRectangle(x, y, w/10, h/10);
			getStatistics(area, mean);
			if(mean < bkgd_intensity){
				bkgd_intensity = mean;
				min_x = x;
				min_y = y;
			}
		}
	}
	run("Select None");

	return(bkgd_intensity);
}



function make_membrane_snake(make_membrane_snake_image_id){	// should receive a background subtracted 16-bit image (1 pol) or RGB image (2 pols)
	selectImage(make_membrane_snake_image_id);

	if(bitDepth() == 24){  	// if the image is RGB: two polarizations, create an image that encodes only brightness
		t_number_of_polarizations = 2;
		run("Select None");
		run("Duplicate...", "title=temp_16_bit_image");
		run("HSB Stack");
		run("Delete Slice");
		run("Delete Slice");
		run("16-bit");
		temp_16_bit_snake_image_id = getImageID();
		selectImage(temp_16_bit_snake_image_id);
		// do the shape fitting
	}
	else{
		t_number_of_polarizations = 1;
		temp_16_bit_snake_image_id = 0;
	}
	// 16-bit B/W image

	center_x = w/2;
	center_y = h/2;
	radius = (h+w)/4-1;

	circle_found = false;
	snake_threshold = 0;


	getLine(x1, y1, x2, y2, original_line_width);

	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2));
	}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2) * 1000);
	}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2) / 1000000);
	}
	run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

	if(look_for_vesicles == true){
		List.setCommands;
		
		if(List.indexOf("Hough Circle Transform") == -1){
			showMessage("Please install the Hough Circle Transform plugin for markedly improved automated vesicle recognition.");		
		}
		else{
			run("Duplicate...", "vesicle finding image");
			vesicle_finding_image_ID = getImageID();
			
			getStatistics(thresholding_image_area, thresholding_image_mean, thresholding_image_min, thresholding_image_max);
			snake_threshold = (thresholding_image_max + thresholding_image_min)/2;
			setAutoThreshold("Default dark");
			run("Convert to Mask");

			//resize circle-finding image to max 120x120 pixels
			if(h > w){size_reduction_factor = h/120;}
			else{size_reduction_factor = w/120;}
			run("Size...", "width="+toString(round(w/size_reduction_factor))+" height="+toString(round(h/size_reduction_factor))+" average interpolation=Bilinear");		
			run("Multiply...", "value=255");
			setBatchMode("show");
			print("\nbefore Hough transform");
			run("Clear Results");

			min_radius = 5; //5 microns
			print("min_radius="+min_radius);
			toUnscaled(min_radius);
			print("min_radius="+min_radius);
			min_radius=round(min_radius);

	
			run("Hough Circle Transform","minRadius="+min_radius+", maxRadius=100, inc=1, minCircles=1, maxCircles=5, threshold=0.5, resolution=656, ratio=1.0, bandwidth=10, local_radius=10,  reduce results_table");

			wait(1000);

			selectImage(vesicle_finding_image_ID);
			close();
	
			if(nResults > 0){
				circle_found = true;	
				print("a circle found by Hough Circle Transform", nResults);
				if(	getResult("X (microns)",0) == getResult("X (microns)",0) ){
					center_x = getResult("X (microns)",0);
					center_y = getResult("Y (microns)",0);
					radius   = getResult("Radius (microns)",0);	
		
					print("center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	
					
					toUnscaled(center_x, center_y);
					toUnscaled(radius);
					
					print("unscaled center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	
					center_x = round(center_x * size_reduction_factor);
					center_y = round(center_y * size_reduction_factor);;
					radius   = round(radius * size_reduction_factor);;	
	
					print("expanded center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	
					if(radius > w/2 || radius > h/2){
						center_x = round(center_x / size_reduction_factor);
						center_y = round(center_y / size_reduction_factor);;
						radius   = round(radius / size_reduction_factor);;					
					}				
				}

				if(	getResult("X (pixels)",0) == getResult("X (pixels)",0) ){
					center_x = getResult("X (pixels)",0);
					center_y = getResult("Y (pixels)",0);
					radius   = getResult("Radius (pixels)",0);	
		
					print("center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);

					center_x = round(center_x * size_reduction_factor);
					center_y = round(center_y * size_reduction_factor);;
					radius   = round(radius * size_reduction_factor);;	
	
					print("expanded center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	
				}

				
				if(t_number_of_polarizations == 2){selectImage(temp_16_bit_snake_image_id);}
				else{selectImage(make_membrane_snake_image_id);}
				
				run("Specify...", "width="+toString(2*radius)+" height="+toString(2*radius)+" x="+toString(center_x-radius)+" y="+toString(center_y-radius)+" oval");
			
			}
			else{
				circle_found = false;	
				print("No circle found by Hough Circle Transform plugin", nResults);
			}
		}	
	}		
	print("center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	setBatchMode(true);

	if(t_number_of_polarizations == 2){
		fit_outline_by_spline(temp_16_bit_snake_image_id, center_x, center_y, radius, t_number_of_polarizations, circle_found);
		selectImage(make_membrane_snake_image_id);
		run("Restore Selection");
		selectImage(temp_16_bit_snake_image_id);
		close();
	}
	else{
		fit_outline_by_spline(make_membrane_snake_image_id, center_x, center_y, radius, t_number_of_polarizations, circle_found);
	}

	selectImage(make_membrane_snake_image_id);

}



function segmentation_function(stack_for_segmentation_image_ID){  // [s]
	print("entering segmentation function");
	setBatchMode(true);
	hor_image_ID = 0;
	ver_image_ID = 0;
	w=getWidth;
	h=getHeight;


	// prepare an image for thresholding
	if (nSlices == 1){   // use a single polarization image; 
		qp_image_name = getTitle();
		t_number_of_polarizations = 1;	

		run("Select None");
		run("Duplicate...", " ");
		rename("bkgd subtracted duplicate");
		
			
		setBatchMode("exit and display");

		
		image_for_thresholding_id = getImageID();
		hv_composite_avg_ID = image_for_thresholding_id;


	}

	else{  		// use a hor/ver polarization composite

		hv_composite_ID = getImageID();
		hv_composite_name = getTitle();
		t_number_of_polarizations = 2;		

		file_name_base = replace(hv_composite_name,"_pcomp.tif","");
		file_name_base = replace(file_name_base,".tif","");
	

		setSlice(1);
		run("Select None");
		run("Duplicate...", " ");
		hor_image_name = file_name_base+"_comp_HOR.tif";
		rename(hor_image_name);
		hor_image_ID = getImageID();


		selectImage(hv_composite_ID);
		setSlice(2);
		run("Duplicate...", " ");
		ver_image_name = file_name_base + "_comp_VER.tif";
		rename(ver_image_name);
		ver_image_ID = getImageID();

		setForegroundColor(0,0,0);
		setTool("Paintbrush Tool");

		image_for_thresholding_id = 0;

		//create the average of the 16-bit composite stack.

		imageCalculator("Average create", hor_image_ID, ver_image_ID);
		run("Grays");
		run("Enhance Contrast", "saturated=0.35");
		hv_composite_avg_ID = getImageID();
		image_for_thresholding_id = hv_composite_avg_ID;

		getPixelSize(pixel_size_unit, pixel_width, pixel_height);
		if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um"){
			// nothing
		}
		else{
			run("Properties...", "unit=micron pixel_width="+pixel_size+" pixel_height="+pixel_size);
		}


		//check if the RG version of the 16-bit composite stack is open. If yes, use a copy of it for spline fitting and segmentation.

		image_titles = getList("image.titles");	
		RG_image_of_LD_image_id = 0;
		for (i = 0; i<image_titles.length; i++){
			if(startsWith(image_titles[i],file_name_base) && (endsWith(image_titles[i],"RG.tif"))){
				selectImage(image_titles[i]);
				RG_image_of_LD_image_id  = getImageID();
				run("Duplicate...", " ");
				image_for_thresholding_id = getImageID();
				i = image_titles.length;
				if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometer" || pixel_size_unit == "um"){
					// nothing
				}
				else{
					run("Properties...", "unit=micron pixel_width="+pixel_size+" pixel_height="+pixel_size);
				}
			}
		}
	}

	setTool("rectangle");

	// do thresholding/manual segmentation


	//--------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------- spline fitting -----------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------


	if(image_for_thresholding_id != hv_composite_avg_ID){selectImage(image_for_thresholding_id);} // ... use the copy of the RG image for thresholding
	else{selectImage(hv_composite_avg_ID);}

	center_x = w/2;
	center_y = h/2;
	radius = (h+w)/4-1;

	circle_found = false;
	snake_threshold = 0;


	getLine(x1, y1, x2, y2, original_line_width);

	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2));
	}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2) * 1000);
	}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_height + pixel_width)/2) / 1000000);
	}
	run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

	if(look_for_vesicles == true){
		List.setCommands;
		
		if(List.indexOf("Hough Circle Transform") == -1){
			showMessage("Please install the Hough Circle Transform plugin for markedly improved automated vesicle recognition.");		
		}
		else{
			run("Duplicate...", "vesicle finding image");
			vesicle_finding_image_ID = getImageID();
			image_info = getInfo("image.subtitle");
			if(indexOf(image_info, "RGB") > 0){  //if the image is RGB, create an image that encodes only brightness
				run("HSB Stack");
				run("Delete Slice");
				run("Delete Slice");
				run("16-bit");
			}
			
			getStatistics(thresholding_image_area, thresholding_image_mean, thresholding_image_min, thresholding_image_max);
			snake_threshold = (thresholding_image_max + thresholding_image_min)/2;
			setAutoThreshold("Default dark");
	
			run("Convert to Mask");
			//resize circle-finding image to max 120x120 pixels
			if(h > w){size_reduction_factor = h/120;}
			else{size_reduction_factor = w/120;}
			run("Size...", "width="+toString(round(w/size_reduction_factor))+" height="+toString(round(h/size_reduction_factor))+" average interpolation=Bilinear");		
			run("Multiply...", "value=255");
			setBatchMode("exit and display");
			print("\nbefore Hough transform");
			run("Clear Results");

			min_radius = 5; //5 microns
			toUnscaled(min_radius);
			min_radius=round(min_radius);

	
			run("Hough Circle Transform","minRadius="+min_radius+", maxRadius=100, inc=1, minCircles=1, maxCircles=5, threshold=0.5, resolution=656, ratio=1.0, bandwidth=10, local_radius=10,  reduce results_table");

			wait(2000);

			selectImage(vesicle_finding_image_ID);
			close();
	
			if(nResults > 0){
				circle_found = true;	
				print("a circle found by Hough Circle Transform", nResults);	
				center_x = getResult("X (microns)",0);
				center_y = getResult("Y (microns)",0);
				radius   = getResult("Radius (microns)",0);	
	
				print("center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);

				
				toUnscaled(center_x, center_y);
				toUnscaled(radius);
				
				print("unscaled center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);

				center_x = round(center_x * size_reduction_factor);
				center_y = round(center_y * size_reduction_factor);;
				radius   = round(radius * size_reduction_factor);;	

				print("expanded center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);

				if(radius > w/2 || radius > h/2){
					center_x = round(center_x / size_reduction_factor);
					center_y = round(center_y / size_reduction_factor);;
					radius   = round(radius / size_reduction_factor);;					
				}
				
		
				if(image_for_thresholding_id != hv_composite_avg_ID){selectImage(image_for_thresholding_id);} // ... use the copy of the RG image for thresholding
				else{selectImage(hv_composite_avg_ID);}
				run("Specify...", "width="+toString(2*radius)+" height="+toString(2*radius)+" x="+toString(center_x-radius)+" y="+toString(center_y-radius)+" oval");
			
			}
			else{
				circle_found = false;	
				print("No circle found by Hough Circle Transform plugin", nResults);
			}
		}	
	}		
	print("center_x = "+center_x+", center_y = "+center_y+", radius = "+radius);
	setBatchMode(true);


	// make an initial polygon selection, either by using Hough circle finding parameters or center in the middle of image, radius = (h/2+v/2)/2
	if(image_for_thresholding_id != hv_composite_avg_ID){selectImage(image_for_thresholding_id);} // ... use the copy of the RG image for thresholding
	else{selectImage(hv_composite_avg_ID);}
	fit_outline_by_spline(t_number_of_polarizations, center_x, center_y, radius, circle_found);
	setBatchMode("exit and display");


	// if user modifies selection, look for a new cell (no node fixing in fit_guv_by_spline) or for a new guv (make nodes ~circular in fit_guv_by_spline)
	check_with_user = true;
	while(check_with_user == true){
		check_with_user = false;		
		if(image_for_thresholding_id != hv_composite_avg_ID){selectImage(image_for_thresholding_id);} // ... use the copy of the RG image for thresholding
		else{selectImage(hv_composite_avg_ID);}
		run("Copy");
		setTool("multipoint");
		waitForUser("Modify selection or make a new selection, such as by using the polygon or segmented line tool, \n\nfollowed by a spline fit (pressing 'f'). For precise shape fitting of round vesicles, select \n3 points belonging to the vesicle outline. For fitting the shape of a cell, select 1 point\n approximately in the center of the cell. Erase unwanted pixels by the paintbrush tool.\n \nWhen done, click OK.");
	
		getSelectionCoordinates(xpoints, ypoints);
		//print("xpoints", xpoints.length);

		if(selectionType() == 6){ // (poly)line
			makeSelection("polygon", xpoints, ypoints);
		}
		

		if(selectionType() == 10 && xpoints.length == 1){  // a single point
			check_with_user = true;		
			center_x = xpoints[0];
			center_y = ypoints[0];
			radius = (h+w)/4-1;
			look_for_vesicles = false;
			fit_outline_by_spline(t_number_of_polarizations, center_x, center_y, radius, false); //'false' = do not look for a circle by my algorithm
		}
		
		if(selectionType() == 10 && xpoints.length == 3){ // three points
			check_with_user = true;		
			run("Fit Circle");		
			Roi.getBounds(x, y, width, height)
			center_x = x + width/2;
			center_y = y + height/2;
			radius   = width/2;
			look_for_vesicles = true;
			fit_outline_by_spline(t_number_of_polarizations, center_x, center_y, radius, true); //'true' = look for a circle by my algorithm
		}

		if(selectionType() == 1){  // circle/oval
			run("Interpolate");
			run("Fit Spline");
		}

	}

	run("Fit Spline");

	setColor("orange");
	Overlay.addSelection;	

		
	if(t_number_of_polarizations == 1){		
		snake_image_name = file_name_base + "_SNK.txt";
	}
	else{
		snake_image_name = file_name_base + "_comp_SNK.txt";
	}
	rename(snake_image_name);
	//save(snake_image_name);	
	saveAs("XY Coordinates");

	//if the user saves it under a name that does not end with _SNK.txt, fix it
	snake_image_name = File.name;
	snake_image_directory = File.directory;
	if(endsWith(snake_image_name,"_SNK.txt")){}
	else{
		File.delete(snake_image_directory+snake_image_name);
		if(endsWith(snake_image_name,".txt")){
			if(t_number_of_polarizations == 1){snake_image_name = replace(snake_image_name, "\\.txt", "_SNK\\.txt");}
			else{snake_image_name = replace(snake_image_name, "\\.txt", "_comp_SNK\\.txt");}
		}
		else{
			if(t_number_of_polarizations == 1){snake_image_name = snake_image_name + "_SNK.txt";}
			else{snake_image_name = snake_image_name + "_comp_SNK.txt";}
		}
		save(snake_image_directory+snake_image_name);	
	}


	if(hor_image_ID < 0){	
		selectImage(hor_image_ID);
		run("Restore Selection");
		setColor("orange");
		Overlay.addSelection("orange");	
		hor_image_name = replace(snake_image_name, "SNK\\.txt", "HOR\\.tif");
		save(File.directory+hor_image_name);
		close();
	}
	if(ver_image_ID < 0){
		selectImage(ver_image_ID);
		run("Restore Selection");
		setColor("orange");
		Overlay.addSelection("orange");	
		ver_image_name = replace(snake_image_name, "SNK\\.txt", "VER\\.tif");
		save(File.directory+ver_image_name);
		close();
	}



	if(image_for_thresholding_id != hv_composite_avg_ID){
		selectImage(image_for_thresholding_id); // ... use the copy of the RG image for thresholding
		run("Select None");
		imageCalculator("Difference", image_for_thresholding_id, RG_image_of_LD_image_id);
		getStatistics(area, mean, min, max);
		run("16-bit");
		run("Macro...", "code=[if(v == 0) v=255; else v = 0]");
		run("Convert to Mask");
		if(max > 0){run("Invert");} // if there is no touch-up, use the whole image
	}
	else{selectImage(hv_composite_avg_ID);}


	mask_image_name = replace(snake_image_name, "SNK\\.txt", "BIN\\.tif");
	rename(mask_image_name);	
	save(File.directory+mask_image_name);
	close();

	selectImage(hv_composite_avg_ID);
	close();

	run("Line Width...", "line="+original_line_width); 

	segmented_file_path = File.directory+mask_image_name;
	return segmented_file_path;
}


function fit_outline_by_spline(spline_fitting_image_id, center_x_f, center_y_f, radius_f, t_number_of_polarizations, circle_found_f){ //[s]
	
	selectImage(spline_fitting_image_id);
	w=getWidth;
	h=getHeight;
	Overlay.remove;
	setFont("SansSerif", 10);

	if(bitDepth() == 24){  	// if the image is RGB: two polarizations, create an image that encodes only brightness
		t_number_of_polarizations = 2;
		run("Select None");
		run("Duplicate...", "title=temp_16_bit_image");
		run("HSB Stack");
		run("Delete Slice");
		run("Delete Slice");
		run("16-bit");
		temp_16_bit_snake_image_id = getImageID();
		selectImage(temp_16_bit_snake_image_id);
		// do the shape fitting
	}
	else{
		t_number_of_polarizations = 1;
		temp_16_bit_snake_image_id = 0;
	}


	// shoot lines in theta increments, find corresponding radii
	
	num_of_theta_increments = 36;
	theta_increment_size = 2 * PI/num_of_theta_increments;
	radius_array         = newArray(num_of_theta_increments);
	thickness_array      = newArray(num_of_theta_increments);
	sorted_radii         = newArray(num_of_theta_increments);
	sorted_thicknesses   = newArray(num_of_theta_increments);
	spline_node_x_values = newArray(num_of_theta_increments);
	spline_node_y_values = newArray(num_of_theta_increments);
	nodes_ok_by_radius   = newArray(num_of_theta_increments);		
	nodes_ok_by_angle    = newArray(num_of_theta_increments);		
	nodes_ok_by_thickness= newArray(num_of_theta_increments);		
	corrected_spline_node_x_values = newArray(num_of_theta_increments);
	corrected_spline_node_y_values = newArray(num_of_theta_increments);
	corrected_radii      = newArray(num_of_theta_increments);
		
	for(i = 0; i < num_of_theta_increments; i++){
		theta = theta_increment_size * i;
		end_x = center_x_f + cos(theta)*radius_f*1.2;
		end_y = center_y_f + sin(theta)*radius_f*1.2;
		makeLine(center_x_f, center_y_f, end_x, end_y, 30);

		profile = getProfile();
		Array.getStatistics(profile, min, max, mean, stdDev);

		temp_radius = radius_f;
		if( temp_radius > profile.length-2 ){
			temp_radius = profile.length-2;
		}
		if(profile[temp_radius] == max){   // search for outline boundaries starting from the center of the outline
			temp_radius = radius_f;
			while(profile[temp_radius] > 0.9*max && temp_radius > 0){temp_radius--;}
			inner_radius = temp_radius;
			temp_radius = radius_f;
			while(profile[temp_radius] > 0.7*max && temp_radius < profile.length-1){temp_radius++;}
			outer_radius = temp_radius;			
		}
		else{								// search for outline boundaries starting from the center of the circle
			temp_radius = 0;
			while(profile[temp_radius] < 0.9*max && temp_radius < profile.length-1){temp_radius++;}
			inner_radius = temp_radius;
			while(profile[temp_radius] > 0.7*max && temp_radius < profile.length-1){temp_radius++;}
			outer_radius = temp_radius;			
		}
		
		mid_radius = (inner_radius + outer_radius)/2;
		radius_array[i] = mid_radius;
		thickness_array[i] = outer_radius - inner_radius;

		//print(i, mid_radius);

		spline_node_x_values[i] = center_x_f + cos(theta)*mid_radius;
		spline_node_y_values[i] = center_y_f + sin(theta)*mid_radius;
		makePoint(spline_node_x_values[i], spline_node_y_values[i]);
				
		//wait(200);
		
	}
	vertex_angles = Array.getVertexAngles(spline_node_x_values, spline_node_y_values, 1); 

	setColor("orange");
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));
	}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);
	}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){
		cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);
	}

	run("Line Width...", "line="+toString(cutoff_distance_from_snake_in_pixels)); 

	makeSelection("polyline",Array.concat(spline_node_x_values,spline_node_x_values[0]), Array.concat(spline_node_y_values,spline_node_y_values[0]));
	

	if(look_for_vesicles == true || circle_found_f == true){	// make more circular
	
		// identify nodes that need to be adjusted
		
		for(i=0; i < num_of_theta_increments; i++){
			sorted_radii[i] = radius_array[i];
			sorted_thicknesses[i] = thickness_array[i];
		}
		sorted_radii = Array.sort(sorted_radii);
		sorted_thicknesses = Array.sort(sorted_thicknesses);
		
		if(radius_array.length/2 == round(radius_array.length/2)){	
			median_radius = (sorted_radii[radius_array.length/2-1] + sorted_radii[radius_array.length/2])/2;
			median_thickness = (sorted_thicknesses[radius_array.length/2-1] + sorted_thicknesses[radius_array.length/2])/2;
		}
		else{
			median_radius = sorted_radii[radius_array.length/2-1];
			median_thickness = sorted_thicknesses[radius_array.length/2-1];
		}

		radius_cutoff = 0.05;      // < 5% deviation is allowed
		thickness_cutoff = 0.20;   // < 20% deviation is allowed
		vertex_angle_cutoff = 5;   // < 5 degree deviation is allowed
		setFont("SansSerif", 10);
	
		for(i=0; i < num_of_theta_increments; i++){
			if( abs(radius_array[i] - median_radius)/median_radius <= radius_cutoff){nodes_ok_by_radius[i] = 1;}
			else{nodes_ok_by_radius[i] = 0;}
			
			if( abs(thickness_array[i] - median_thickness)/median_thickness > thickness_cutoff && t_number_of_polarizations == 2){nodes_ok_by_thickness[i] = 0;}
			else{nodes_ok_by_thickness[i] = 1;}

			if( abs(vertex_angles[i] - 360/num_of_theta_increments) <= vertex_angle_cutoff){nodes_ok_by_angle[i] = 1;}
			else{nodes_ok_by_angle[i] = 0;}
		

		}
	
			

		completion_needed = false;		// are there any nodes that need to be fixed?
		for(i=0; i < num_of_theta_increments; i++){
			if(nodes_ok_by_radius[i] * nodes_ok_by_angle[i] * nodes_ok_by_thickness[i] == 0){completion_needed = true;}
		}

		completion_possible = false;	// is there any consecutive pair of nodes that is good, and therefore can be used for fixing others?
		for(i=0; i < num_of_theta_increments-1; i++){
			if(nodes_ok_by_radius[i] * nodes_ok_by_angle[i] * nodes_ok_by_thickness[i] == 1 && nodes_ok_by_radius[i+1] * nodes_ok_by_angle[i+1] * nodes_ok_by_thickness[i+1] == 1){completion_possible = true;}
		}
	

		
		if(completion_needed == true && completion_possible == true){
			// find the longest stretch of good nodes - initialization
			longest_good_stretch_start = 0;
			longest_good_stretch_end = 0;
			longest_good_stretch_length = 0;
			current_stretch_start = 0;
			current_stretch_end = 0;
			current_stretch_length = 0;

			nodes_ok_by_radius_twice = Array.concat(nodes_ok_by_radius, nodes_ok_by_radius);
			nodes_ok_by_angle_twice = Array.concat(nodes_ok_by_angle, nodes_ok_by_angle);
			nodes_ok_by_thickness_twice = Array.concat(nodes_ok_by_thickness, nodes_ok_by_thickness);
			radius_array_twice = Array.concat(radius_array,radius_array);

			// find the longest stretch of good nodes - actual searching			
			for(i=0; i < 2*num_of_theta_increments-1; i++){
				// if current node is bad and the next one is good, note the start of a stretch, look for its end
				if(nodes_ok_by_radius_twice[i] * nodes_ok_by_angle_twice[i] * nodes_ok_by_thickness_twice[i] == 0 && nodes_ok_by_radius_twice[i+1] * nodes_ok_by_angle_twice[i+1] * nodes_ok_by_thickness_twice[i+1] == 1){
					current_stretch_start = i+1;
					while(nodes_ok_by_radius_twice[i+1] * nodes_ok_by_angle_twice[i+1] * nodes_ok_by_thickness_twice[i+1] == 1 && i < (2*num_of_theta_increments-2)){
						current_stretch_end = i+1;
						current_stretch_length = current_stretch_end - current_stretch_start + 1;
						i++;
					}
					if (current_stretch_length > longest_good_stretch_length){
						longest_good_stretch_start = current_stretch_start;
						longest_good_stretch_end = current_stretch_end;   // this is different from how it was implemented earlier, where the end was < num_of_theta_increments
						if(longest_good_stretch_end > num_of_theta_increments - 1){longest_good_stretch_end = longest_good_stretch_end - num_of_theta_increments - 1;}
						longest_good_stretch_length = current_stretch_length;
					}
				}
			}


			// starting from the end of the longest stretch of good nodes, correct the bad nodes by interpolating between
			// values of the node before the stretch of bad nodes and the node after the stretch of bad nodes
			corrected_spline_node_x_values_twice = Array.concat(spline_node_x_values, spline_node_x_values);
			corrected_spline_node_y_values_twice = Array.concat(spline_node_y_values, spline_node_y_values);
			corrected_radii_twice = Array.concat(radius_array, radius_array);
	
			for(i = longest_good_stretch_start; i < longest_good_stretch_start + num_of_theta_increments - 1; i++){
				// find the end of the current good stretch: if the current node is ok and the next one is ok as well, increase i by 1
				while(nodes_ok_by_radius_twice[i] * nodes_ok_by_angle_twice[i] * nodes_ok_by_thickness_twice[i] == 1 && nodes_ok_by_radius_twice[i+1] * nodes_ok_by_angle_twice[i+1] * nodes_ok_by_thickness_twice[i+1] == 1 && i < longest_good_stretch_start + num_of_theta_increments - 1){
					i = i + 1;
				}
				// now i is the index of a node that is good, node i + 1 is bad
				// find the next node that is good, judging by radius onlyhou
				j = i + 1;
				while (nodes_ok_by_radius_twice[j] == 0  && j < num_of_theta_increments * 2 - 1){
					j = j + 1;
				}
				// now I have the index (i) of the last good node of the current good stretch
				// and the index (j) of the first good node of the next good stretch 

				// Interpolate node (i + 1)
				corrected_radii_twice[i+1] = corrected_radii_twice[i] + (corrected_radii_twice[j] - corrected_radii_twice[i])/(j - i);
				corrected_spline_node_x_values_twice[i+1] = center_x_f + cos((i+1) * 360/num_of_theta_increments/180*PI)*corrected_radii_twice[i+1];
				corrected_spline_node_y_values_twice[i+1] = center_y_f + sin((i+1) * 360/num_of_theta_increments/180*PI)*corrected_radii_twice[i+1];
			}

			corrected_radii                = Array.slice(corrected_radii_twice,                longest_good_stretch_start, longest_good_stretch_start + num_of_theta_increments);
			corrected_spline_node_x_values = Array.slice(corrected_spline_node_x_values_twice, longest_good_stretch_start, longest_good_stretch_start + num_of_theta_increments);
			corrected_spline_node_y_values = Array.slice(corrected_spline_node_y_values_twice, longest_good_stretch_start, longest_good_stretch_start + num_of_theta_increments);

		 

			setColor("magenta");
			makeSelection("polyline",Array.concat(corrected_spline_node_x_values,corrected_spline_node_x_values[0]), Array.concat(corrected_spline_node_y_values,corrected_spline_node_y_values[0]));

			setColor("white");

		}
	}
	run("Fit Spline");
	if(t_number_of_polarizations == 2){
		selectImage(spline_fitting_image_id);
		run("Restore Selection");
		selectImage(temp_16_bit_snake_image_id);
		close();
	}
	
}


function show_LD_azimuths_along_the_snake(t_BIN_image_id, t_bkgd_subbed_multipol_stack_id){
	selectImage(t_bkgd_subbed_multipol_stack_id);
	number_of_polarizations = nSlices;
	w=getWidth;
	h=getHeight;
	getPixelSize(pixel_size_unit, pixel_width, pixel_height);

	polarization_list = newArray(number_of_polarizations);
	for (i = 0; i < number_of_polarizations; i++) {
		polarization_list[i] = starting_pol_angle + i * pol_angle_increment;
	}


	selectImage(t_BIN_image_id);
	BIN_file_name = getTitle();
	getSelectionCoordinates(snake_points_x, snake_points_y);
	run("Select None");

	file_name_base_f = replace(BIN_file_name,"_BIN.tif","");

	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}

	imageCalculator("Multiply create 32-bit stack", t_bkgd_subbed_multipol_stack_id, t_BIN_image_id);
	BINned_multipol_stack_image_id = getImageID();	

	makeSelection("polyline", snake_points_x, snake_points_y);
	run("Reslice [/]...", "output=0.1 start=Top avoid");
	snake_intensity_image_id = getImageID();
	snake_intensity_list_length = getWidth();

	snake_point_intensities = newArray(number_of_polarizations);
	snake_log2rmax_values = newArray(snake_intensity_list_length);
	snake_azimuths = newArray(snake_intensity_list_length);
	snake_fit_r_squared_values = newArray(snake_intensity_list_length);
	snake_point_distances = newArray(snake_intensity_list_length);
	List.setMeasurements;
	snake_length_in_microns = List.getValue("Width");
	for (i = 0; i < snake_intensity_list_length; i++) {
		snake_point_distances[i] = i * snake_length_in_microns/snake_intensity_list_length;
	}
	//print("snake length = "+snake_length_in_microns);
	//snake_point_index = newArray(snake_intensity_list_length);

	if(determine_rmax_alpha0_sigma == "1PPM: r(max)"){		//, "2PPM: r(max), "+alpha0char+", "+sigmachar

		fitting_equation = "y = a*a + b*b * cos(x/180*3.14159 - c*c)*cos(x/180*3.14159 - c*c)";		
	
		for (i = 0; i < snake_intensity_list_length; i++) {
			showStatus("Fitting...");
			showProgress(i/snake_intensity_list_length);
			makeLine(i, 0, i, number_of_polarizations, 1);
			snake_point_intensity_profile = getProfile();

			snake_intensities_real = true;
	
			for(j = 0; j < number_of_polarizations; j++){
				snake_point_intensities[j] = snake_point_intensity_profile[j];
				if(isNaN(snake_point_intensities[j])){snake_intensities_real = false;}
			}
	
			// print("snake point #: "+i);
			// Array.print(snake_point_intensities);

			Fmin = NaN;
			Fmax = NaN;
			azimuth = NaN;
			snake_azimuths[i] = NaN;

			log2rmax = NaN;
			snake_log2rmax_values[i] = NaN;
			snake_fit_r_squared_values[i] = NaN;


			if(snake_intensities_real == true){
				Fit.doFit(fitting_equation, polarization_list, snake_point_intensities);
				//Fit.logResults;
				if(Fit.rSquared > 0){
					Fmin = Fit.p(0) * Fit.p(0);
					Fmax = Fmin + Fit.p(1) * Fit.p(1);
					azimuth = (Fit.p(2) * Fit.p(2))/PI*180 - floor((Fit.p(2) * Fit.p(2))/PI*180/180)*180;
					// snake_azimuths[i] = azimuth;
					snake_azimuths[i] = (-azimuth + 270) - floor((-azimuth + 270)/180)*180;
			
					log2rmax = log(Fmax/Fmin)/log(2);
					snake_log2rmax_values[i] = log2rmax;
					snake_fit_r_squared_values[i] = Fit.rSquared;
					
				}
			}
	
	
			//print("log2rmax = "+log2rmax+", azimuth = "+azimuth);
			//print("--------------------");
	
			//snake_point_index[i] = i;
			
		}
	
		Array.getStatistics(Array.deleteValue(snake_log2rmax_values, NaN), log2rmax_min, log2rmax_max, log2rmax_mean, log2rmax_stdev); 
	
		
		Plot.create("Log"+subtwochar+"(rmax) values along the selected outline", "Position along the outline, in microns", "Log"+subtwochar+"(rmax)", snake_point_distances, snake_log2rmax_values);	
		Plot.setColor("black");
		Plot.add("point", snake_point_distances, snake_fit_r_squared_values);
		Plot.setColor("black");
		Plot.setLegend("Log2(rmax); Log2(rmax) = "+toString(round(log2rmax_mean*1000)/1000)+" ± "+toString(round((log2rmax_stdev/sqrt(snake_intensity_list_length))*1000)/1000)+" (mean ± 95%CI)\tGoodness of fit (r-squared values)");

	}
	else{    // 2PPM fitting
		fitting_equation = "y = a * (1 + b * cos(x/180*3.14159 - d*d)*cos(x/180*3.14159 - d*d) + c * cos(x/180*3.14159 - d*d)*cos(x/180*3.14159 - d*d)*cos(x/180*3.14159 - d*d)*cos(x/180*3.14159 - d*d)  )";		
	
		for (i = 0; i < snake_intensity_list_length; i++) {
			showStatus("Fitting...");
			showProgress(i/snake_intensity_list_length);
			makeLine(i, 0, i, number_of_polarizations, 1);
			snake_point_intensity_profile = getProfile();

			snake_intensities_real = true;

			for(j = 0; j < number_of_polarizations; j++){
				snake_point_intensities[j] = snake_point_intensity_profile[j];
				if(isNaN(snake_point_intensities[j])){snake_intensities_real = false;}
			}
	
			//print("snake point #: "+i);
			//Array.print(snake_point_intensities);

			Fmin = NaN;
			Fmax = NaN;
			azimuth = NaN;
			snake_azimuths[i] = NaN;

			log2rmax = NaN;
			snake_log2rmax_values[i] = NaN;
			snake_fit_r_squared_values[i] = NaN;

	
			if(snake_intensities_real == true){
				Fit.doFit(fitting_equation, polarization_list, snake_point_intensities);
				//Fit.logResults;
				if(Fit.rSquared > 0){
					Fmin = 1;
					Fmax = 1 + Fit.p(1) + Fit.p(2);
					azimuth = (Fit.p(3) * Fit.p(3))/PI*180 - floor((Fit.p(3) * Fit.p(3))/PI*180/180)*180;
					// snake_azimuths[i] = azimuth;
					snake_azimuths[i] = (-azimuth + 270) - floor((-azimuth + 270)/180)*180;
			
					log2rmax = log(Fmax/Fmin)/log(2);
					snake_log2rmax_values[i] = log2rmax;
					snake_fit_r_squared_values[i] = Fit.rSquared;
					
				}
			}
	
	
			//print("log2rmax = "+log2rmax+", azimuth = "+azimuth);
			//print("--------------------");
	
			//snake_point_index[i] = i;
			
		}
	
		Array.getStatistics(Array.deleteValue(snake_log2rmax_values, NaN), log2rmax_min, log2rmax_max, log2rmax_mean, log2rmax_stdev); 
	
		
		Plot.create("Log"+subtwochar+"(rmax) values along the selected outline", "Position along the outline, in microns", "Log"+subtwochar+"(rmax)", snake_point_distances, snake_log2rmax_values);	
		Plot.setColor("black");
		Plot.add("point", snake_point_distances, snake_fit_r_squared_values);
		Plot.setColor("black");
		Plot.setLegend("Log2(rmax); Log2(rmax) = "+toString(round(log2rmax_mean*1000)/1000)+" ± "+toString(round((log2rmax_stdev/sqrt(snake_intensity_list_length))*1000)/1000)+" (mean ± 95%CI)\tGoodness of fit (r-squared values)");
		
	}
	
	for (i = 0; i < snake_intensity_list_length; i++) {
		hue = round(360/snake_intensity_list_length * i);
		hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
		if(hue >= 0 && hue < 60){
			red = 255;
			green = round(hsv_x * 255);
			blue = 0;
		}

		if(hue >= 60 && hue < 120){
			red = round(hsv_x * 255);
			green = 255;
			blue = 0;
		}

		if(hue >= 120 && hue < 180){
			red = 0;
			green = 255;
			blue = round(hsv_x * 255);
		}
		if(hue >= 180 && hue < 240){
			red = 0;
			green = round(hsv_x * 255);
			blue = 255;
		}

		if(hue >= 240 && hue < 300){
			red = round(hsv_x * 255);
			green = 0;
			blue = 255;
		}
		
		if(hue >= 300 && hue < 360){
			red = 255;
			green = 0;
			blue = round(hsv_x * 255);
		}


		
		// for some reason the line below throws an error that I cannot solve
		// I think the IJ.pad() function cannot accept the output of the toHex function, no matter what I try
		// plot_color = "#"+IJ.pad(toHex(red),2)+IJ.pad(toHex(green),2)+IJ.pad(toHex(blue),2);

		// This works:
		plot_color = "#"+convert_to_hex(red)+convert_to_hex(green)+convert_to_hex(blue);
		
		Plot.setColor(plot_color,plot_color);
		snake_point_distance = newArray(snake_point_distances[i],snake_point_distances[i]);
		value = newArray(snake_point_distances[i]-snake_point_distances[i],snake_point_distances[i]-snake_point_distances[i]);
		Plot.add("box", snake_point_distance, value);

	}
	if(log2rmax_max > 1){
		Plot.setLimits(0, snake_length_in_microns, 0, log2rmax_max * 1.1);		
	}
	else{
		Plot.setLimits(0, snake_length_in_microns, 0, 1.1);	
	}
	Plot.show();
	log2rmax_r_squared_plot_image_id = getImageID();

	theta_array = make_theta_array(snake_points_x, snake_points_y);
	//print("theta_array:");
	//Array.print(theta_array);
	theta_deg_array = newArray(theta_array.length);
	theta_distances_array = newArray(theta_array.length);
	
	for (i = 0; i < theta_array.length; i++) {
		theta_deg_array[i] = theta_array[i]/PI*180 + 90;
	}
	//print("theta_deg_array:");
	//Array.print(theta_deg_array);

	theta_distances_array[0] = 0;
	for (i = 1; i < theta_array.length; i++) {
		theta_distances_array[i] = theta_distances_array[i-1] + sqrt((snake_points_x[i]-snake_points_x[i-1])*(snake_points_x[i]-snake_points_x[i-1])*pixel_width*pixel_width + (snake_points_y[i]-snake_points_y[i-1])*(snake_points_y[i]-snake_points_y[i-1])*pixel_height*pixel_height); 
	}
	//print("theta_distances_array:");
	//Array.print(theta_distances_array);


	Plot.setColor("black");
	Plot.create("Azimuth values along the selected outline", "Position along the outline, in microns", "Azimuth (°)");
	Plot.setLineWidth(2);
	Plot.add("dot", snake_point_distances, snake_azimuths);	
	Plot.setColor("red");
	Plot.setLineWidth(1);
	Plot.add("circle", theta_distances_array, theta_deg_array);
	Plot.setLegend("Azimuth angles calculated from fluorescence intensity profiles\tOutline directions (values of the angle "+thetachar+")");
	for (i = 0; i < snake_intensity_list_length; i++) {
		hue = round(360/snake_intensity_list_length * i);
		hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
		if(hue >= 0 && hue < 60){
			red = 255;
			green = round(hsv_x * 255);
			blue = 0;
		}

		if(hue >= 60 && hue < 120){
			red = round(hsv_x * 255);
			green = 255;
			blue = 0;
		}

		if(hue >= 120 && hue < 180){
			red = 0;
			green = 255;
			blue = round(hsv_x * 255);
		}
		if(hue >= 180 && hue < 240){
			red = 0;
			green = round(hsv_x * 255);
			blue = 255;
		}

		if(hue >= 240 && hue < 300){
			red = round(hsv_x * 255);
			green = 0;
			blue = 255;
		}
		
		if(hue >= 300 && hue < 360){
			red = 255;
			green = 0;
			blue = round(hsv_x * 255);
		}


		
		// for some reason the line below throws an error that I cannot solve
		// I think the IJ.pad() function cannot accept the output of the toHex function, no matter what I try
		// plot_color = "#"+IJ.pad(toHex(red),2)+IJ.pad(toHex(green),2)+IJ.pad(toHex(blue),2);

		// This works:
		plot_color = "#"+convert_to_hex(red)+convert_to_hex(green)+convert_to_hex(blue);
			
		Plot.setColor(plot_color,plot_color);
		snake_point_distance = newArray(snake_point_distances[i],snake_point_distances[i]);
		value = newArray(snake_point_distances[i]-snake_point_distances[i],snake_point_distances[i]-snake_point_distances[i]);
		Plot.add("box", snake_point_distance, value);

	}
	Plot.setLimits(0, snake_length_in_microns, 0, 180);
	Plot.setColor("black");
	Plot.show();
	azimuths_thetas_plot_image_id = getImageID();

	selectImage(BINned_multipol_stack_image_id);
	close();
	selectImage(snake_intensity_image_id);
	close();

	//setBatchMode("exit and display");


	/*
		plot_image_id =       parseInt(LD_quantitation_results[0]);
		hyperstack_image_id = parseInt(LD_quantitation_results[1]);
		theta_image_id =      parseInt(LD_quantitation_results[2]);
		fitting_data_text =  	       LD_quantitation_results[3];
		fitting_results_text =         LD_quantitation_results[4];

	*/
	setTool(0);

	plot_image_id_array = newArray(log2rmax_r_squared_plot_image_id, azimuths_thetas_plot_image_id);
	
	return plot_image_id_array;
}

function convert_to_hex(dec_number){
	hex_number = toHex(dec_number);
	if(dec_number < 16){hex_number = "0"+hex_number;}
	return hex_number;
}


function perform_LD_quantitation_1PPM(t_BIN_image_id, t_HOR_VER_image_id){  // the BIN image contains values of 1 and NaN, and carries a selection (= snake points)
	selectImage(t_BIN_image_id);
	BIN_file_name = getTitle();
	getSelectionCoordinates(snake_points_x, snake_points_y);
	run("Select None");

	file_name_base_f = replace(BIN_file_name,"_BIN.tif","");

	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}

	
	selectImage(t_HOR_VER_image_id);
	t_number_of_polarizations = nSlices;
	HOR_VER_image_name = getTitle();
	if(t_number_of_polarizations == 1){
		if(endsWith(HOR_VER_image_name,"_HOR.tif")){hor_exists = true;  ver_exists = false;}
		if(endsWith(HOR_VER_image_name,"_VER.tif")){hor_exists = false; ver_exists = true; }	
		imageCalculator("Multiply create 32-bit", t_HOR_VER_image_id, t_BIN_image_id);	
		HOR_VER_snk_image_id = getImageID();
		rename(HOR_VER_image_name);
	}
	if(t_number_of_polarizations == 2){
		hor_exists = true;
		ver_exists = true;
		imageCalculator("Multiply create 32-bit stack", t_HOR_VER_image_id, t_BIN_image_id);	
		HOR_VER_snk_image_id = getImageID();
		rename(HOR_VER_image_name);
	}
	if(t_number_of_polarizations > 2){
		hor_exists = true;
		ver_exists = true;
		imageCalculator("Multiply create 32-bit stack", t_HOR_VER_image_id, t_BIN_image_id);	
		MULTIPOL_snk_image_id = getImageID();
		rename("MULTIPOL_snk_image_name");
	}


	
	w=getWidth;
	h=getHeight;


	//for each point on the snake calculate theta
	theta = make_theta_array(snake_points_x, snake_points_y);

	if(membranes_or_filaments == "Analyze a filament"){
		for(i = 0; i < theta.length; i++){
			print("before: "+theta[i]);
			theta[i] = theta[i] + PI/2;
			if(theta[i] > PI/2){theta[i] = theta[i] - PI;}		
			print("after: "+theta[i]);		
		}
	}



	// For each pixel of interest:
	// 1) find a closest point on the snake to the pixel of interest, 
	// 2) add hor intensity to the snake point's hor count, 
	// 3) add ver intensity to the snake point's ver count, 
	// 4) increment the snake point's pixel count
	// 5) put the snake point's theta value in the pixel of interest of a 32-bit theta value image 


	hor_sum = newArray(snake_points_x.length);
	ver_sum = newArray(snake_points_x.length);
	pixel_count = newArray(snake_points_x.length);
	r = newArray(snake_points_x.length);
	log2r = newArray(snake_points_x.length);


	showStatus("Calculating membrane orientations ("+thetachar+" values");

	

	Array.getStatistics(snake_points_x, min_snake_x, max_snake_x);
	Array.getStatistics(snake_points_y, min_snake_y, max_snake_y);

	col_min = floor(min_snake_x - cutoff_distance_from_snake_in_pixels - 1);
	col_max = floor(max_snake_x + cutoff_distance_from_snake_in_pixels + 1);
	row_min = floor(min_snake_y - cutoff_distance_from_snake_in_pixels - 1);
	row_max = floor(max_snake_y + cutoff_distance_from_snake_in_pixels + 1);


	if(t_number_of_polarizations <= 2){
	
		showStatus("Associating pixels with spline points");
		selectImage(HOR_VER_snk_image_id);
		
		for (col = col_min; col <= col_max; col++) {
			showProgress(col, col_max - col_min);		
			for (row = row_min; row <= row_max; row++){
				pixel_value = getPixel(col, row);
				if(pixel_value == pixel_value){		// false if pixel_value = NaN
					min_distance = 100000;
									
					for (i=0; i < snake_points_x.length; i++){
						distance_from_snake_point = sqrt((col-snake_points_x[i])*(col-snake_points_x[i]) + (row-snake_points_y[i])*(row-snake_points_y[i]));
						if(distance_from_snake_point < min_distance){
							min_distance = distance_from_snake_point;
							closest_snake_point_index = i;
						}
					}
	
					if(t_number_of_polarizations == 2){
						setSlice(1);
						hor_sum[closest_snake_point_index] = hor_sum[closest_snake_point_index] + getPixel(col, row);				
						setPixel(col, row, theta[closest_snake_point_index]);
						setSlice(2);
						ver_sum[closest_snake_point_index] = ver_sum[closest_snake_point_index] + getPixel(col, row);				
					}
					if(t_number_of_polarizations == 1){
						if(hor_exists == true){hor_sum[closest_snake_point_index] = hor_sum[closest_snake_point_index] + getPixel(col, row);}
						else{ver_sum[closest_snake_point_index] = ver_sum[closest_snake_point_index] + getPixel(col, row);}
					}
					
					pixel_count[closest_snake_point_index] = pixel_count[closest_snake_point_index]+1;
	
					setPixel(col, row, theta[closest_snake_point_index]);
				
				}
			}
		}
	
	
		theta_image_file_name = replace(getTitle(),"BIN.tif", "THT.tif");
		rename(theta_image_file_name);
		theta_image_id = getImageID();			// the HOR_VER_snk image is now THT image
		if(t_number_of_polarizations == 2){
			setSlice(2);
			run("Delete Slice", "delete=channel");
			run("phase");
			setMinAndMax(-PI/2,PI/2);
		}

		
		//for each point on the snake, find the theta value, find the pixels that have it, sum up their Fh and Fv
	
	
		for(i=0; i<snake_points_x.length; i++){
	
			if (hor_exists == true && ver_exists == true){
				r[i] = hor_sum[i]/ver_sum[i];
				log2r[i] = log(r[i])/log(2);
			}
			else{
				if (hor_exists == true){
					r[i]=hor_sum[i]/pixel_count[i];
					log2r[i]=hor_sum[i]/pixel_count[i];
				}
				if (ver_exists == true){
					r[i]=ver_sum[i]/pixel_count[i];
					log2r[i]=ver_sum[i]/pixel_count[i];
				}
			}
		}
	
	

		// for VER polarization only shift the theta angles by 90o (so that we can treat the image as acquired with HOR polarization
		if (hor_exists == false && ver_exists == true && ratiometric_fitting == false){	
			for (i=0; i<snake_points_x.length-1; i++){
				theta[i] = theta[i]-PI/2;
				if(theta[i] < -PI/2){
					theta[i] = theta[i]+PI;
				}
			}
		}


		// find number of snake points with valid theta values
		// put the 'good' snake points and associated log2r values into arrays x, y
	
		x = newArray();
		y = newArray();
		n = 0;
		if(t_number_of_polarizations == 2 || (t_number_of_polarizations == 1 && ratiometric_fitting == false) ){
			for(i=0; i<snake_points_x.length; i++){
				if(theta[i] == theta[i] && log2r[i] == log2r[i] && log2r[i] < 100000  && log2r[i] > -100000){
					x = Array.concat(x, theta[i]);
					y = Array.concat(y, log2r[i]); 
					n = n + 1;
				}
			}
		}
		Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 
	}
	


	if(t_number_of_polarizations > 2){
		// make a 'BIN' image, in which the pixels that are not NaN hold as their value the index of the nearest snake point
		showStatus("Associating pixels with spline points...");
		selectImage(t_BIN_image_id);
		run("Duplicate...", "title=snake_point_indices");
		snake_point_index_image_id = getImageID();
		for (col = col_min; col <= col_max; col++) {
			showProgress(col, col_max - col_min);		
			for (row = row_min; row <= row_max; row++){
				pixel_value = getPixel(col, row);
				if(pixel_value == pixel_value){		// false if pixel_value = NaN
					min_distance = 100000;
					for (i=0; i < snake_points_x.length; i++){
						distance_from_snake_point = sqrt((col-snake_points_x[i])*(col-snake_points_x[i]) + (row-snake_points_y[i])*(row-snake_points_y[i]));
						if(distance_from_snake_point < min_distance){
							min_distance = distance_from_snake_point;
							closest_snake_point_index = i;
						}
					}
					setPixel(col, row, closest_snake_point_index);				
				}
			}
		}
	
		newImage("snake point intensities pol 1", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_intensity_pol_1_image_id = getImageID();  // this image will hold intensities for each snake point (col) and polarization (row)
	
	
	
		for (i = 0; i < snake_points_x.length; i++) {
			// for each snake point, make a 'BIN' image in which NaN pixels are only the ones closest to that snake point
			showStatus("Analyzing data for individual spline points...");
			showProgress(i, snake_points_x.length);
			selectImage(snake_point_index_image_id);
			run("Duplicate...", "title=snake_point_"+toString(i)+"_BIN");
			particular_snake_point_BIN_image_id = getImageID();
			run("Macro...", "code=[if(v == "+toString(i)+") v=1; else v=NaN]");
			showStatus("Analyzing data for individual spline points...");
			pixel_count[i] = getValue("RawIntDen");
			if(matches(pixel_count[i],"NaN")){pixel_count[i] = 0;}

			imageCalculator("Multiply create 32-bit stack", MULTIPOL_snk_image_id, particular_snake_point_BIN_image_id);
			rename("snake_point_#"+toString(i)+"_intensities");
			particular_snake_point_intensity_image_id = getImageID();
	
			for(j = 1; j <= t_number_of_polarizations; j++){
				selectImage(particular_snake_point_intensity_image_id);
				setSlice(j);
				getStatistics(area, mean, min, max, std, histogram);
				selectImage(snake_point_intensity_pol_1_image_id);
				setPixel(i, j-1, mean/1000*area);
			}
	
			selectImage(particular_snake_point_intensity_image_id);
			close();

			selectImage(particular_snake_point_BIN_image_id);
			close();
			
		}



	
		newImage("snake point intensities pol 2", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_intensity_pol_2_image_id = getImageID();
	
		list_of_polarizations = newArray(t_number_of_polarizations);
		for(i = 0; i < t_number_of_polarizations; i++){
			list_of_polarizations[i] = starting_pol_angle + i * pol_angle_increment;
		}
	
		for(i = 0; i < t_number_of_polarizations; i++){  // for each pol1, find pol2 that is either 90° larger or 90° smaller than pol1
			pol2_index = NaN;
			pol1 = list_of_polarizations[i];
			for(j = 0; j < t_number_of_polarizations; j++){   // look for pol2 values that are 90° larger than pol1
				if(list_of_polarizations[j] == pol1 + 90){
					pol2_index = j;
				}
			}
			if(isNaN(pol2_index)){                                 // if no pol2 has been found...
				for(j = 0; j < t_number_of_polarizations; j++){    // ... look for pol2 that is 90° smaller than pol1
					if(list_of_polarizations[j] == pol1 - 90){
						pol2_index = j;
					}
				}				
			}
			if(isNaN(pol2_index)){
				showMessage("Unable to find a perpendicular pair of polarizations! Polarization #1 = "+toString(pol1)+degreechar);
			}
			else{                                                    
				selectImage(snake_point_intensity_pol_1_image_id);
				makeLine(0, pol2_index, snake_points_x.length-1, pol2_index, 1);
				row_of_F_values = getProfile();
		
				selectImage(snake_point_intensity_pol_2_image_id);
				for(j = 0; j < snake_points_x.length; j++){
					setPixel(j, i, row_of_F_values[j]);
				}
				//makeLine(0, i, snake_points_x.length-1, i, 1);
				//run("Paste");		
			}
				
		}


	
		imageCalculator("Divide 32-bit create", snake_point_intensity_pol_1_image_id, snake_point_intensity_pol_2_image_id);
		run("Log");
		run("Divide...", "value=0.693");
		rename("Log2(r) image");
		log2r_image_id = getImageID();	
	
		x = newArray();
		y = newArray();
	
	
	
		newImage("snake point theta values", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_theta_image_id = getImageID();
	
		// calculate theta and log2r for each snake point and 90° polarization pair
		for (col = 0; col < snake_points_x.length; col++) {

			if(pixel_count[col] > 0){
				
				for (row = 0; row < t_number_of_polarizations; row++){
	//				theta_value = abs((theta[col]-list_of_polarizations[row]/180*PI) - floor((theta[col]-list_of_polarizations[row]/180*PI)/PI)*PI);			
					theta_value = abs((theta[col]+list_of_polarizations[row]/180*PI) - floor((theta[col]+list_of_polarizations[row]/180*PI)/PI)*PI);			
					if(theta_value > PI/2){theta_value = theta_value - PI;}
					selectImage(snake_point_theta_image_id);
					setPixel(col, row, theta_value);				
	
					selectImage(log2r_image_id);
					log2r = getPixel(col, row);
					if(isNaN(log2r) || toString(log2r) == "Infinity" || toString(log2r) == "-Infinity"){}
					else{
						x = Array.concat(x, theta_value); 
						y = Array.concat(y, log2r);
					}
				}
			}			
		}

		//print("x array:");
		//Array.print(x);
		//print("y array:");
		//Array.print(y);
	
		n = x.length;
		Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 

		selectImage(MULTIPOL_snk_image_id);
		close();

		if(save_files == false){
			selectImage(log2r_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_1_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_2_image_id);
			close();
	
			selectImage(snake_point_theta_image_id);
			close();				
		}



		selectImage(snake_point_index_image_id);
		for (i = 0; i < snake_points_x.length; i++) {
			run("Macro...", "code=[if(v == "+toString(i)+") v="+toString(theta[i])+";]");	
			showStatus("Associating pixels with spline points...");
			showProgress(i, snake_points_x.length);		
		}
		rename("theta image");
		theta_image_id = getImageID();
		//close();


	//setBatchMode("exit and display");
	}

	
	selectImage(t_BIN_image_id);
	makeSelection("polyline", snake_points_x, snake_points_y);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1P fitting and plotting //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// ratiometric fitting of data obtained with a single polarization
	if(t_number_of_polarizations <= 2 && ratiometric_fitting == true){
		binned_theta_deg_array = newArray(181);

		binned_hor = newArray(181);
		binned_hor_counts = newArray(181);

		binned_ver = newArray(181);
		binned_ver_counts = newArray(181);

		binned_log2r = newArray(181);

//
		number_of_nonzero_theta_bins = 0;

		for(i=0; i<snake_points_x.length; i++){
			binned_theta_deg = floor(theta[i]/PI*180);   // binned_theta_deg runs from -90 to +90 
			if(hor_exists == true && log2r[i]==log2r[i]){
				binned_hor[binned_theta_deg + 90] = binned_hor[binned_theta_deg + 90] + log2r[i];  // I am increasing the indices by 90 just so that they would not be negative
				binned_hor_counts[binned_theta_deg + 90] = binned_hor_counts[binned_theta_deg + 90] + 1;
			}
			if(ver_exists == true && log2r[i]==log2r[i]){
				binned_ver[binned_theta_deg + 90] = binned_ver[binned_theta_deg + 90] + log2r[i];
				binned_ver_counts[binned_theta_deg + 90] = binned_ver_counts[binned_theta_deg +90] + 1;
			}
		}
		for(i=0; i<=180; i++){
			binned_theta_deg_array[i] = i-90;
			if(binned_hor_counts[i] > 0){ binned_hor[i] = binned_hor[i]/binned_hor_counts[i]; }
			if(binned_ver_counts[i] > 0){ binned_ver[i] = binned_ver[i]/binned_ver_counts[i]; }
			if(binned_hor_counts[i] > 0 || binned_ver_counts[i] > 0){ number_of_nonzero_theta_bins++; }
		}
	}

//


	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		if(hor_exists == true){
			for(i=0; i<=90; i++){binned_ver[i+90] = binned_hor[i]; }
			for(i=90; i<=180; i++){binned_ver[i-90] = binned_hor[i]; }
		}
		if(ver_exists == true){
			for(i=0; i<=90; i++){binned_hor[i+90] = binned_ver[i]; }
			for(i=90; i<=180; i++){binned_hor[i-90] = binned_ver[i]; }
		}
		for(i=0; i<=180; i++){
			if(binned_hor[i] > 0 && binned_ver[i] > 0){
				binned_log2r[i] = log(  binned_hor[i]/binned_ver[i]  )/log(2);
				x = Array.concat(x,i/180*PI - PI/2);
				y = Array.concat(y,binned_log2r[i]);
				n = n + 1;			
			}
		}
		hor_exists = true;
		ver_exists = true;
	}

	Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 


	n_of_snake_points_with_positive_theta = 0;
	n_of_snake_points_with_negative_theta = 0;
	fraction_of_snake_points_with_data = n/snake_points_x.length;
	
	rmax = 0;

	// Make a header if there are two or more polarizations...
	if ((hor_exists == true && ver_exists == true) || (t_number_of_polarizations == 1 && ratiometric_fitting == true)){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tr \tpixel_count";}
	// Make a header if there is only one polarization...
	else{
		if(hor_exists == true && ratiometric_fitting == false){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tFh_avg \tpixel_count";}
		if(ver_exists == true && ratiometric_fitting == false){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tFv_avg \tpixel_count";}
	}


	// List data used for fitting, if there is only one polarization, but ratiometric fitting (which relies on theta binning) is being used 
	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		for(i=0; i<=180; i++){
			if(binned_hor[i] > 0 && binned_ver[i] > 0){
				data_for_fitting_text = data_for_fitting_text + "\r\n" + toString((binned_theta_deg_array[i]/180*PI))+" \t"+toString(binned_hor[i])+" \t"+toString(binned_ver[i])+" \t"+"0"+" \t"+"0"+" \t"+toString((binned_hor[i]/binned_ver[i]))+" \t"+"1";
			}
		}
	}
	if(t_number_of_polarizations == 1 && ratiometric_fitting == false){
		for(i=0; i<snake_points_x.length; i++){
			data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta[i])+" \t"+toString(hor_sum[i])+" \t"+toString(ver_sum[i])+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r[i])+" \t"+toString(pixel_count[i]);	 		
		}		
	}
	if(t_number_of_polarizations == 2){
		for(i=0; i<snake_points_x.length; i++){
			data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta[i])+" \t"+toString(hor_sum[i])+" \t"+toString(ver_sum[i])+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r[i])+" \t"+toString(pixel_count[i]);	 		
		}		
	}
	if(t_number_of_polarizations > 2){
		if(save_files == true){
			for(i=0; i<snake_points_x.length; i++){
				for(j=0; j<number_of_polarizations; j++){
					selectImage(log2r_image_id);
					log2r = getPixel(i, j);
					if(isNaN(log2r)){}
					else{
						r = pow(2, log2r);
						selectImage(snake_point_theta_image_id);
						theta = getPixel(i, j);
						selectImage(snake_point_intensity_pol_1_image_id);
						hor_sum = getPixel(i, j);
						selectImage(snake_point_intensity_pol_2_image_id);
						ver_sum = getPixel(i, j);
						data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta)+" \t"+toString(hor_sum)+" \t"+toString(ver_sum)+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r)+" \t"+toString(pixel_count[i]);	 			
					}
				}				
			}	
			selectImage(log2r_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_1_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_2_image_id);
			close();
	
			selectImage(snake_point_theta_image_id);
			close();
		}					
	}


	
	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		if(phase_fixed == false){look_for_pol_angles = 1;}
	}


	
	// look for Rmax and phase
	if(t_number_of_polarizations > 2){

		/*
		number_of_fitting_params = 3;
		initialGuesses = newArray(1, 0, 0);
		offset_string = "c + ";
		offset = 0;

		fitting_equation = "y ="+offset_string+"log((1 + a*cos(2*x-2*b)  )/(1 - a*cos(2*x-2*b))  )/log(2)";				
		index_of_param_cos2theta = 0;
		index_of_param_phase = 1;
		index_of_param_phase2 = -1;
		index_of_param_offset = 2;				
		index_of_param_normalization_factor = -1;	
		*/
		
		number_of_fitting_params = 2;
		initialGuesses = newArray(1, 0);
		offset_string = "b + ";
		offset = 0;
		phase = 0;
		phase1_deg = 0;

		fitting_equation = "y ="+offset_string+"log((1 + a*cos(2*x)  )/(1 - a*cos(2*x))  )/log(2)";				
		index_of_param_cos2theta = 0;
		index_of_param_phase = -1;
		index_of_param_phase2 = -1;
		index_of_param_offset = 1;				
		index_of_param_normalization_factor = -1;	
		
	}
	if(t_number_of_polarizations <= 2){
		number_of_fitting_params = 0;
		if (hor_exists == true && ver_exists == true){  // 2 polarizations
			if(phase_fixed == true){         // look for 0 phases
				if(look_for_pol_angles == 1){
					phase2_deg = phase1_deg;
				}
				if(number_of_nonzero_theta_bins < 90){
					number_of_fitting_params = 1;
					initialGuesses = newArray(1);
					offset_string = "";
					index_of_param_offset = -1;
					offset = 0;
				}
				else{
					number_of_fitting_params = 2;
					initialGuesses = newArray(1, 0);
					offset_string = "b + ";
					index_of_param_offset = 1;				
				}			
				fitting_equation = "y ="+offset_string+"log((1 + a*cos(2*x-2*("+toString(phase1_deg/180*PI)+"))    )/(1 -       a*cos(2*x-2*("+toString(phase2_deg/180*PI)+"))    ))/log(2)";				
				index_of_param_cos2theta = 0;
				index_of_param_phase = -1;
				index_of_param_phase2 = -1;
				index_of_param_normalization_factor = -1;
			}
			if(phase_fixed == false && look_for_pol_angles == 1){         // look for 1 phase  
				if(number_of_nonzero_theta_bins < 90){
					number_of_fitting_params = 2;
					initialGuesses = newArray(1, 0);
					offset_string = "";
					index_of_param_offset = -1;
					offset = 0;
				}
				else{
					number_of_fitting_params = 3;
					initialGuesses = newArray(1, 0, 0);
					offset_string = "c + ";
					index_of_param_offset = 2;				
				}			
				fitting_equation = "y ="+offset_string+"log((1 + a*cos(2*x-2*b)  )/(1 - a*cos(2*x-2*b))  )/log(2)";				
				index_of_param_cos2theta = 0;
				index_of_param_phase = 1;
				index_of_param_phase2 = -1;
				index_of_param_normalization_factor = -1;
			}
			if(phase_fixed == false && look_for_pol_angles == 2){         // look for 2 phase
				if(number_of_nonzero_theta_bins < 90){
					number_of_fitting_params = 3;
					initialGuesses = newArray(1, 0, 0);
					offset_string = "";
					index_of_param_offset = -1;
					offset = 0;
				}
				else{
					number_of_fitting_params = 4;
					initialGuesses = newArray(1, 0, 0, 0);
					offset_string = "d + ";
					index_of_param_offset = 3;				
				}			
				fitting_equation = "y ="+offset_string+"log(  (1 +  a*cos(2*x-2*b) )/(1 - a*cos(2*x-2*c)  )  )/log(2)";				
				index_of_param_cos2theta = 0;
				index_of_param_phase = 1;
				index_of_param_phase2 = 2;
				index_of_param_normalization_factor = -1;
			}
		}
		else{                                    // 1 polarization
			if(phase_fixed == true){        // look for 0 phases			
				initialGuesses = newArray(1, 0);
				fitting_equation = "y = b * (1 + ("+toString(polarization_purity/100)+")*a*cos(2*x-2*("+toString(phase1_deg/180*PI)+")) + ("+toString(1-polarization_purity/100)+")*a*cos(2*x-2*(("+toString(phase1_deg/180*PI)+")-"+toString(PI)+"/2)))";		
				index_of_param_cos2theta = 0;
				index_of_param_offset = -1;
				index_of_param_phase = -1;
				index_of_param_phase2 = -1;
				index_of_param_normalization_factor = 1;
				number_of_fitting_params = 2;
			}
			if(phase_fixed == false){        // look for 1 phase1
				initialGuesses = newArray(1, 0, 1);
				fitting_equation = "y = c * (1 + ("+toString(polarization_purity/100)+")*a*cos(2*x-2*b) + ("+toString(1-polarization_purity/100)+")*a*cos(2*x-2*(b-"+toString(PI)+"/2)))";		
				index_of_param_cos2theta = 0;
				index_of_param_offset = -1;
				index_of_param_phase = 1;
				index_of_param_phase2 = -1;
				index_of_param_normalization_factor = 2;
				number_of_fitting_params = 3;
			}
		}
	}

	showStatus("Fitting 1P data...");
	Fit.logResults;
	Fit.doFit(fitting_equation, x, y, initialGuesses);

	r_square = Fit.rSquared;
	
	B_from_fit = Fit.p(index_of_param_cos2theta);
	
	if(index_of_param_phase > -1){phase = Fit.p(index_of_param_phase);}
	else{phase = parseFloat(phase1_deg)/180*PI;}
	
	if(index_of_param_phase2 > -1){phase2 = Fit.p(index_of_param_phase2);}
	else{phase2 = phase;}
	//else{phase2 = parseFloat(phase2_deg)/180*PI;}
	
	if(index_of_param_offset > -1){offset = Fit.p(index_of_param_offset);}
	else{offset = 0;}
	
	if(index_of_param_normalization_factor > -1){normalization_factor = Fit.p(index_of_param_normalization_factor);}
	else{normalization_factor = 1;}


	phase  = phase  % PI;
	phase2 = phase2 % PI;

	if (phase  < 0){phase  = phase  + PI;}		// from here on, phase should be between 0 and Pi
	if (phase2 < 0){phase2 = phase2 + PI;}		// from here on, phase2 should be between 0 and Pi	 
	if (phase  > PI/2){phase  = phase  - PI;}	// from here on, phase should be between -Pi/2 and Pi/2 
	if (phase2 > PI/2){phase2 = phase2 - PI;}	// from here on, phase2 should be between -Pi/2 and Pi/2 
	
	if (phase > PI/4){
		phase = phase - PI/2; 	
		phase2 = phase2 - PI/2;			
		B_from_fit = - B_from_fit;
	}
	
	if (phase < -PI/4){
		phase = phase + PI/2;
		phase2 = phase2 + PI/2;
		B_from_fit = - B_from_fit;
	}
	
	if (phase2 >  PI/2){phase2 = phase2 - PI;}	 
	if (phase2 < -PI/2){phase2 = phase2 + PI;}


	rmax = (1+B_from_fit)/(1-B_from_fit);		
	log2_rmax = log(rmax)/log(2);
	
	phase1_deg  = round((phase/PI*180)*10)/10;
	phase2_deg = round((phase2/PI*180)*10)/10;

	if(index_of_param_phase > -1){
		if(look_for_pol_angles == 2){phase_text = toString(phase1_deg)+degreechar+"/"+toString(phase2_deg)+degreechar;}
		if(look_for_pol_angles == 1){phase_text = toString(phase1_deg)+degreechar;}	
	}
	else{
		if(look_for_pol_angles == 2){phase_text = toString(phase1_deg)+degreechar+"/"+toString(phase2_deg)+degreechar+" (fixed)";}
		if(look_for_pol_angles == 1){phase_text = toString(phase1_deg)+degreechar+" (fixed)";}	
		if(look_for_pol_angles == 0){
			if(phase1_deg != phase2_deg){phase_text = toString(phase1_deg)+degreechar+"/"+toString(phase2_deg)+degreechar+" (fixed)";}
			else{phase_text = toString(phase1_deg)+degreechar+" (fixed)";}		
		}		
		phase = phase1_deg/180*PI;
		phase2 = phase2_deg/180*PI;
	}





	results_of_fitting_text = "Results of 1PPM fitting: \r\nFitting equation: "+fitting_equation+"\r\n";
	if(	number_of_fitting_params == 2 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\n\r\n";}
	if(	number_of_fitting_params == 3 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\nc = "+Fit.p(2)+"\r\n\r\n";}
	if(	number_of_fitting_params == 4 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\nc = "+Fit.p(2)+"\r\nd = "+Fit.p(3)+"\r\n\r\n";}

	if(index_of_param_offset > -1){results_of_fitting_text = results_of_fitting_text + "offset = "+offset+"\r\n";}
	if(index_of_param_normalization_factor > -1){results_of_fitting_text = results_of_fitting_text + "normalization_factor = "+normalization_factor+"\r\n";}
	if(index_of_param_cos2theta > -1){results_of_fitting_text = results_of_fitting_text + "B1P = "+B_from_fit+"\r\n";}
	if(index_of_param_phase > -1){results_of_fitting_text = results_of_fitting_text + "phase = "+phase+" (from fit)\r\n";}
	else{results_of_fitting_text = results_of_fitting_text + "phase = "+phase+" (from user)\r\n";}
	if(index_of_param_phase2 > -1){results_of_fitting_text = results_of_fitting_text + "phase2 = "+phase2+" (from fit)\r\n";}
	else{
		if(phase1_deg != phase2_deg){
			results_of_fitting_text = results_of_fitting_text + "phase2 = "+phase2+" (from user)\r\n";		
		}
	}
	print("Here:"+results_of_fitting_text);
	if(hor_exists == false || ver_exists == false){results_of_fitting_text = results_of_fitting_text + "polarization_direction = "+polarization_direction+"\r\n";}
	results_of_fitting_text = results_of_fitting_text+"\r\nnumber_of_polarizations = "+t_number_of_polarizations+"\r\n\r\nr(max) \t1/r(max) \tlog2(r(max)) \tphase \tr"+suptwochar+" \r\n"+rmax+" \t"+1/rmax+" \t"+log2_rmax+" \t"+phase_text+" \t"+r_square+" \r\n--------------------------------------------------\r\n";




	print(results_of_fitting_text);


	// generate x, y values of the fit by using the fitting function

	num_of_samples=200;
	fit_x_values = newArray(num_of_samples); 
	fit_y_values = newArray(num_of_samples); 
	xmin=-3.14159/2; 
	xmax= 3.14159/2;
	step=(xmax-xmin)/num_of_samples;
	for (i=0; i< num_of_samples; i++) {
		fit_x_values[i]=xmin+(i+1)*step;
		fit_y_values[i]=Fit.f(fit_x_values[i]);
	}

	Array.getStatistics(y, y_min, y_max);
	Array.getStatistics(fit_y_values, fit_y_min, fit_y_max);

	if (hor_exists == true && ver_exists == true){
		if(toString(fit_y_min) == "Infinity" || toString(fit_y_min) == "-Infinity" || toString(fit_y_max) == "Infinity" || toString(fit_y_max) == "-Infinity"){
			lower_plot_limit = y_min * 1.1;
			upper_plot_limit = y_max * 1.1;
		}
		else{
			upper_plot_limit = y_max * 1.1;
			if(abs(log2_rmax) > y_max){ upper_plot_limit = abs(log2_rmax) * 1.1; }
			lower_plot_limit = y_min * 1.6;
			if(-abs(log2_rmax) < y_min){ lower_plot_limit = -abs(log2_rmax) * 1.6;}		
		}

	}
	else{
		lower_plot_limit = 0;
		upper_plot_limit = y_max * 1.1;
	}

	Plot.create(file_name_base_f+"_1P_FIT", "Membrane orientation, "+thetachar+" ["+degreechar+"]", "Log"+subtwochar+"(r)");

	Plot.setLimits(-90, 90, lower_plot_limit, upper_plot_limit);
	Plot.setLineWidth(1.5);
	Plot.setColor("Blue");
	fit_x_values_deg = newArray(fit_x_values.length);
	for(i = 0; i < fit_x_values.length; i++){
		fit_x_values_deg[i] = fit_x_values[i]/3.14159*180;
	}
	Plot.add("line", fit_x_values_deg, fit_y_values);
	Plot.setLineWidth(1);
	Plot.setColor("gray");
	x_deg = newArray(x.length);
	for(i = 0; i < x.length; i++){
		x_deg[i] = x[i]/3.14159*180;
	}
	Plot.add("circle", x_deg, y);
	Plot.setLineWidth(1.5);
	Plot.setColor("Blue");
	Plot.add("line", fit_x_values_deg, fit_y_values);
//	text = "r(max) = "+rmax+" ("+1/rmax+"), phase = "+phase_text+", B = "+B_from_fit+", r"+suptwochar+" = "+r_square;
	text = "log2(r(max)) = "+log2_rmax+", phase = "+phase_text+", B = "+B_from_fit+", r"+suptwochar+" = "+r_square;
	Plot.addLegend(text, "Auto");
	Plot.setFormatFlags("11000100110011");
	Plot.show;
	plot_window_image_id = getImageID();

	theta0_hor_sum = 0;
	theta0_ver_sum = 0;
	theta0_pixel_count = 0;
	theta0_r = 0;
	theta0_log2r = 0;
	thetaPi2_hor_sum = 0;
	thetaPi2_ver_sum = 0;
	thetaPi2_pixel_count = 0;
	thetaPi2_r = 0;
	thetaPi2_log2r = 0;
	thetaPi4_hor_sum = 0;
	thetaPi4_ver_sum = 0;
	thetaPi4_pixel_count = 0;
	thetaPi4_r = 0;
	thetaPi4_log2r = 0;




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1P alpha0, sigma /////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();

	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();

	num_of_alpha0s = getWidth;
	num_of_sigmas  = getHeight;

	showStatus("Calculating "+alpha0char+", "+sigmachar);

	run("Image Expression Parser (Macro)", "expression=[(4*C - 7*B)/(-10 + B + 2*C)] a=None b=2PPM_B.tif c=2PPM_C.tif d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
	rename("1PPM_B.tif");
	B_1PPM_image_id = getImageID();


	newImage(file_name_base_f+"_1P_RMSD_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	RMSD_image_id = getImageID();

	newImage(file_name_base_f+"_1P_chi_sq_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	chi_sq_image_id = getImageID();

	for(i = 0; i < x.length; i++){
		showProgress(i, x.length);
		showStatus("Calculating "+alpha0char+", "+sigmachar);

		if (hor_exists == true && ver_exists == true){
			expectation_expression = toString(offset) + " + log( (1 + A * cos(" + toString(2*x[i]-2*phase) + ") ) / (1 - A * cos(" + toString(2*x[i]-2*phase2) + ") ) ) / log(2)";				
		}
		else{
			expectation_expression = toString(normalization_factor) + " *  (1 + A * "+toString(polarization_purity/100)+" * cos(" + toString(2*x[i]-2*phase) + ") + A * ("+toString(1-polarization_purity/100)+") * cos(" + toString(2*x[i]-2*(phase-3.14159/2)) + "))";
		}
		
		run("Image Expression Parser (Macro)", "expression=["+ expectation_expression +"] a=1PPM_B.tif b=None c=None d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		expectation_image_id = getImageID();
		rename("expectation image");

		//RMSD
		run("Duplicate...", "title=sq_dev");
		sq_dev_image_id = getImageID();
		run("Subtract...", "value="+toString(y[i]));
		run("Square");
		imageCalculator("Add 32-bit", RMSD_image_id, sq_dev_image_id); 
		showStatus("Calculating "+alpha0char+", "+sigmachar);

		//Chi-squared
		selectImage(expectation_image_id);
		run("Abs");
		run("Add...", "value=0.01");
		imageCalculator("Divide 32-bit", sq_dev_image_id, expectation_image_id);
		imageCalculator("Add 32-bit", chi_sq_image_id, sq_dev_image_id);
		showStatus("Calculating "+alpha0char+", "+sigmachar);


		selectImage(sq_dev_image_id);
		close();
		
		selectImage(expectation_image_id);
		close();
	}

	selectImage(chi_sq_image_id);
	run("Divide...", "value="+toString(n*1.2));
	run("phase");
	setMinAndMax(0, 1);


	selectImage(RMSD_image_id);
	run("Divide...", "value="+n);
	run("Square Root");
	getRawStatistics(nPixels, mean, min, max);


	run("Find Maxima...", "noise="+max+" output=[Point Selection] light");
	getSelectionBounds(alpha0, sigma_row_number, w, h);
	sigma = num_of_sigmas - 1 - sigma_row_number;
	selectImage(B_1PPM_image_id);
	B1P = getPixel(alpha0, sigma_row_number);
	log2rmax = abs(log((1+B1P)/(1-B1P))/log(2));
	selectImage(RMSD_image_id);
	run("Select None");

	run("phase");
	setMinAndMax(min, 5*min);


////////////////////////////////////////////////////////////////////////////////////
//////////////////////// make hyperstack of 1P results /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



	hyperstack_name = file_name_base_f+"_1P_goodness_of_fit";
	
	newImage(hyperstack_name, "32-bit color-mode", num_of_alpha0s, num_of_sigmas, 4, 1, 1);
	hyperstack_image_id = getImageID();
	setColor(140,215,255);  // light blue

	selectImage(RMSD_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(1);
	Overlay.drawString("1P r"+suptwochar, 5, 20);
	Overlay.setPosition(1,1,1);
	Overlay.show;	
	run("Set Label...", "label=[1P r"+suptwochar+"]");
	run("Paste");
	run("Divide...", "value="+y_stdev+" slice");
	run("Square", "slice");
	run("Multiply...", "value=-1 slice");
	run("Add...", "value=1 slice");
	run("Cyan Hot");
	setMinAndMax(0, 1);


	selectImage(RMSD_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(2);
	Overlay.drawString("1P RMSD", 5, 20);
	Overlay.setPosition(2,1,1);
	Overlay.show;	
	run("Set Label...", "label=[1P RMSD]");
	run("Paste");
	run("phase");
	setMinAndMax(min,5*min);


	
	selectImage(chi_sq_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(3);
	Overlay.drawString("1P Chi"+suptwochar, 5, 20);
	Overlay.setPosition(3,1,1);
	Overlay.show;
	run("Set Label...", "label=[1P Chi"+suptwochar+"]");
	run("Paste");
	run("phase");
	setMinAndMax(0, 1);


	selectImage(RMSD_image_id);
	run("Duplicate...", "title=RMSD");
	low_1P_RMSD_image_id = getImageID();
	getStatistics(temp_area, temp_mean, temp_min, temp_max, temp_stdev);
	setThreshold(temp_min, 1.01*temp_min);
	run("Convert to Mask");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(4);
	Overlay.drawString("1P "+alpha0char+", "+sigmachar, 5, 20);
	Overlay.setPosition(4,1,1);
	Overlay.show;
	run("Paste");
	run("Grays");
	run("Select None");
	
	selectImage(low_1P_RMSD_image_id);
	close();

	selectImage(hyperstack_image_id);
	setSlice(1);

	run("Properties...", "channels=4 slices=1 frames=1 unit="+degreechar+" pixel_width=1.0000 pixel_height=-1 voxel_depth=1.0000 origin=0,90");


	
	selectImage(RMSD_image_id);
	close();
	selectImage(chi_sq_image_id);
	close();


////////////////////////////////////////////////////


	selectImage(B_1PPM_image_id);
	close();

	selectImage(B_2PPM_image_id);
	close();

	selectImage(C_2PPM_image_id);
	close();



	results_of_1P_fitting_array = newArray(5);
	results_of_1P_fitting_array[0] = toString(plot_window_image_id);
	results_of_1P_fitting_array[1] = toString(hyperstack_image_id);
	results_of_1P_fitting_array[2] = toString(theta_image_id);
	results_of_1P_fitting_array[3] = data_for_fitting_text;
	results_of_1P_fitting_array[4] = results_of_fitting_text;



	return results_of_1P_fitting_array;
}















function perform_LD_quantitation_2PPM(t_BIN_image_id, t_HOR_VER_image_id){  // the BIN image contains values of 1 and NaN, and carries a selection (= snake points)
	selectImage(t_BIN_image_id);
	BIN_file_name = getTitle();
	getSelectionCoordinates(snake_points_x, snake_points_y);
	run("Select None");

	file_name_base_f = replace(BIN_file_name,"_BIN.tif","");

	getPixelSize(pixel_size_unit, pixel_width, pixel_height);
	if(pixel_size_unit == "micron" || pixel_size_unit == "microns" || pixel_size_unit == "micrometer" || pixel_size_unit == "micrometers" || pixel_size_unit == "um"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2));}
	if(pixel_size_unit == "nanometer" || pixel_size_unit == "nanometers" || pixel_size_unit == "nm"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) * 1000);}
	if(pixel_size_unit == "meter" || pixel_size_unit == "meters" || pixel_size_unit == "m"){cutoff_distance_from_snake_in_pixels = round(cutoff_distance_from_snake/((pixel_width + pixel_height)/2) / 1000000);}

	
	selectImage(t_HOR_VER_image_id);
	t_number_of_polarizations = nSlices;
	HOR_VER_image_name = getTitle();
	if(t_number_of_polarizations == 1){
		if(endsWith(HOR_VER_image_name,"_HOR.tif")){hor_exists = true;  ver_exists = false;}
		if(endsWith(HOR_VER_image_name,"_VER.tif")){hor_exists = false; ver_exists = true; }	
		imageCalculator("Multiply create 32-bit", t_HOR_VER_image_id, t_BIN_image_id);	
		HOR_VER_snk_image_id = getImageID();
		rename(HOR_VER_image_name);
	}
	if(t_number_of_polarizations == 2){
		hor_exists = true;
		ver_exists = true;
		imageCalculator("Multiply create 32-bit stack", t_HOR_VER_image_id, t_BIN_image_id);	
		HOR_VER_snk_image_id = getImageID();
		rename(HOR_VER_image_name);
	}
	if(t_number_of_polarizations > 2){
		hor_exists = true;
		ver_exists = true;
		imageCalculator("Multiply create 32-bit stack", t_HOR_VER_image_id, t_BIN_image_id);	
		MULTIPOL_snk_image_id = getImageID();
		rename("MULTIPOL_snk_image_name");
	}


	
	w=getWidth;
	h=getHeight;


	//for each point on the snake calculate theta
	theta = make_theta_array(snake_points_x, snake_points_y);

	if(membranes_or_filaments == "Analyze a filament"){
		for(i = 0; i < theta.length; i++){
			theta[i] = theta[i] + PI/2;
			if(theta[i] > PI/2){theta[i] = theta[i] - PI;}
		}
	}


	// For each pixel of interest:
	// 1) find a closest point on the snake to the pixel of interest, 
	// 2) add hor intensity to the snake point's hor count, 
	// 3) add ver intensity to the snake point's ver count, 
	// 4) increment the snake point's pixel count
	// 5) put the snake point's theta value in the pixel of interest of a 32-bit theta value image 


	hor_sum = newArray(snake_points_x.length);
	ver_sum = newArray(snake_points_x.length);
	pixel_count = newArray(snake_points_x.length);
	r = newArray(snake_points_x.length);
	log2r = newArray(snake_points_x.length);


	showStatus("Calculating membrane orientations ("+thetachar+" values");

	

	Array.getStatistics(snake_points_x, min_snake_x, max_snake_x);
	Array.getStatistics(snake_points_y, min_snake_y, max_snake_y);

	col_min = floor(min_snake_x - cutoff_distance_from_snake_in_pixels - 1);
	col_max = floor(max_snake_x + cutoff_distance_from_snake_in_pixels + 1);
	row_min = floor(min_snake_y - cutoff_distance_from_snake_in_pixels - 1);
	row_max = floor(max_snake_y + cutoff_distance_from_snake_in_pixels + 1);
	

	if(t_number_of_polarizations <= 2){

		showStatus("Associating pixels with spline points");
		selectImage(HOR_VER_snk_image_id);
	
		for (col = col_min; col <= col_max; col++) {
			showProgress(col, col_max - col_min);		
			for (row = row_min; row <= row_max; row++){
				pixel_value = getPixel(col, row);
				if(pixel_value == pixel_value){		// false if pixel_value = NaN
					min_distance = 100000;
									
					for (i=0; i < snake_points_x.length; i++){
						distance_from_snake_point = sqrt((col-snake_points_x[i])*(col-snake_points_x[i]) + (row-snake_points_y[i])*(row-snake_points_y[i]));
						if(distance_from_snake_point < min_distance){
							min_distance = distance_from_snake_point;
							closest_snake_point_index = i;
						}
					}
	
					if(t_number_of_polarizations == 2){
						setSlice(1);
						hor_sum[closest_snake_point_index] = hor_sum[closest_snake_point_index] + getPixel(col, row);				
						setSlice(2);
						ver_sum[closest_snake_point_index] = ver_sum[closest_snake_point_index] + getPixel(col, row);				
					}
					if(t_number_of_polarizations == 1){
						if(hor_exists == true){hor_sum[closest_snake_point_index] = hor_sum[closest_snake_point_index] + getPixel(col, row);}
						else{ver_sum[closest_snake_point_index] = ver_sum[closest_snake_point_index] + getPixel(col, row);}
					}
					
					pixel_count[closest_snake_point_index] = pixel_count[closest_snake_point_index]+1;
	
					setPixel(col, row, theta[closest_snake_point_index]);
				
				}
			}
		}
	
	
		theta_image_file_name = replace(getTitle(),"BIN.tif", "THT.tif");
		rename(theta_image_file_name);
		theta_image_id = getImageID();			// the HOR_VER_snk image is now THT image
		if(t_number_of_polarizations == 2){
			setSlice(2);
			run("Delete Slice", "delete=channel");
			run("phase");
			setMinAndMax(-PI/2,PI/2);
		}
	
	
	
	
	
		//for each point on the snake, find the theta value, find the pixels that have it, sum up their Fh and Fv
	
	
		for(i=0; i<snake_points_x.length; i++){
	
			if (hor_exists == 1 && ver_exists == 1){
				r[i] = hor_sum[i]/ver_sum[i];
				log2r[i] = log(r[i])/log(2);
			}
			else{
				if (hor_exists == 1){
					r[i]=hor_sum[i]/pixel_count[i];
					log2r[i]=hor_sum[i]/pixel_count[i];
				}
				if (ver_exists == 1){
					r[i]=ver_sum[i]/pixel_count[i];
					log2r[i]=ver_sum[i]/pixel_count[i];
				}
			}
		}
	
	
	
		// for VER polarization only shift the theta angles by 90o (so that we can treat the image as acquired with HOR polarization
		if (hor_exists == false && ver_exists == true && ratiometric_fitting == false){	
			for (i=0; i<snake_points_x.length-1; i++){
				theta[i] = theta[i]-PI/2;
				if(theta[i] < -PI/2){
					theta[i] = theta[i]+PI;
				}
			}
		}
	
	
		// find number of snake points with valid theta values
		// put the 'good' snake points and associated log2r values into arrays x, y
	
		x = newArray();
		y = newArray();
		n = 0;
		if(t_number_of_polarizations == 2 || (t_number_of_polarizations == 1 && ratiometric_fitting == false) ){
			for(i=0; i<snake_points_x.length; i++){
				if(theta[i] == theta[i] && log2r[i] == log2r[i] && log2r[i] < 100000  && log2r[i] > -100000){
					x = Array.concat(x, theta[i]);
					y = Array.concat(y, log2r[i]); 
					n = n + 1;
				}
			}
		}
		Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 
	}



	if(t_number_of_polarizations > 2){
		showStatus("Associating pixels with spline points");
		selectImage(t_BIN_image_id);
		run("Duplicate...", "title=snake_point_indices");
		snake_point_index_image_id = getImageID();
		for (col = col_min; col <= col_max; col++) {
			showProgress(col, col_max - col_min);		
			for (row = row_min; row <= row_max; row++){
				pixel_value = getPixel(col, row);
				if(pixel_value == pixel_value){		// false if pixel_value = NaN
					min_distance = 100000;
					for (i=0; i < snake_points_x.length; i++){
						distance_from_snake_point = sqrt((col-snake_points_x[i])*(col-snake_points_x[i]) + (row-snake_points_y[i])*(row-snake_points_y[i]));
						if(distance_from_snake_point < min_distance){
							min_distance = distance_from_snake_point;
							closest_snake_point_index = i;
						}
					}
					setPixel(col, row, closest_snake_point_index);				
				}
			}
		}
	
		newImage("snake point intensities pol 1", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_intensity_pol_1_image_id = getImageID();
	
	
	
		for (i = 0; i < snake_points_x.length; i++) {
			// for each snake point, make a 'BIN' image in which NaN pixels are only the ones closest to that snake point
			showStatus("Analyzing data for individual spline points...");
			showProgress(i, snake_points_x.length);
			selectImage(snake_point_index_image_id);
			run("Duplicate...", "title=snake_point_"+toString(i)+"_BIN");
			particular_snake_point_BIN_image_id = getImageID();
			run("Macro...", "code=[if(v == "+toString(i)+") v=1; else v=NaN]");
			showStatus("Analyzing data for individual spline points...");
			pixel_count[i] = getValue("RawIntDen");
			if(matches(pixel_count[i],"NaN")){pixel_count[i] = 0;}
	
			imageCalculator("Multiply create 32-bit stack", MULTIPOL_snk_image_id, particular_snake_point_BIN_image_id);
			rename("snake_point_#"+toString(i)+"_intensities");
			particular_snake_point_intensity_image_id = getImageID();
	
			for(j = 1; j <= t_number_of_polarizations; j++){
				selectImage(particular_snake_point_intensity_image_id);
				setSlice(j);
				getStatistics(area, mean, min, max, std, histogram);
				selectImage(snake_point_intensity_pol_1_image_id);
				setPixel(i, j-1, mean/1000*area);
			}
	
			selectImage(particular_snake_point_intensity_image_id);
			close();
			selectImage(particular_snake_point_BIN_image_id);
			close();
			
		}
	
		newImage("snake point intensities pol 2", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_intensity_pol_2_image_id = getImageID();
	
		list_of_polarizations = newArray(t_number_of_polarizations);
		for(i = 0; i < t_number_of_polarizations; i++){
			list_of_polarizations[i] = starting_pol_angle + i * pol_angle_increment;
		}
	
		for(i = 0; i < t_number_of_polarizations; i++){
			pol2_index = NaN;
			pol1 = list_of_polarizations[i];
			for(j = 0; j < t_number_of_polarizations; j++){
				if(list_of_polarizations[j] == pol1 + 90){
					pol2_index = j;
				}
			}
			if(isNaN(pol2_index)){
				for(j = 0; j < t_number_of_polarizations; j++){
					if(list_of_polarizations[j] == pol1 - 90){
						pol2_index = j;
					}
				}				
			}
			if(isNaN(pol2_index)){
				showMessage("Unable to find a perpendicular pair of polarizations! Polarization #1 = "+toString(pol1)+degreechar);
			}
			else{
				selectImage(snake_point_intensity_pol_1_image_id);
				makeLine(0, pol2_index, snake_points_x.length-1, pol2_index, 1);
				row_of_F_values = getProfile();
		
				selectImage(snake_point_intensity_pol_2_image_id);
				for(j = 0; j < snake_points_x.length; j++){
					setPixel(j, i, row_of_F_values[j]);
				}
				//makeLine(0, i, snake_points_x.length-1, i, 1);
				//run("Paste");		
			}
				
		}
	
		imageCalculator("Divide 32-bit create", snake_point_intensity_pol_1_image_id, snake_point_intensity_pol_2_image_id);
		run("Log");
		run("Divide...", "value=0.693");
		rename("Log2(r) image");
		log2r_image_id = getImageID();	
	
		x = newArray();
		y = newArray();
	
	
	
		newImage("snake point theta values", "32-bit black", snake_points_x.length, t_number_of_polarizations, 1);
		snake_point_theta_image_id = getImageID();
	
		// calculate theta and log2r for each snake point and 90° polarization pair		
		for (col = 0; col < snake_points_x.length; col++) {
			if(pixel_count[col] > 0){
				for (row = 0; row < t_number_of_polarizations; row++){
					//theta_value = abs((theta[col]-list_of_polarizations[row]/180*PI) - floor((theta[col]-list_of_polarizations[row]/180*PI)/PI)*PI);			
					theta_value = abs((theta[col]+list_of_polarizations[row]/180*PI) - floor((theta[col]+list_of_polarizations[row]/180*PI)/PI)*PI);			
					if(theta_value > PI/2){theta_value = theta_value - PI;}
					selectImage(snake_point_theta_image_id);
					setPixel(col, row, theta_value);				
	
					selectImage(log2r_image_id);
					log2r = getPixel(col, row);
					if(isNaN(log2r) || toString(log2r) == "Infinity" || toString(log2r) == "-Infinity"){}
					else{
						x = Array.concat(x, theta_value); 
						y = Array.concat(y, log2r);
					}
				}
			}
		}
	
		n = x.length;
		Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 

		selectImage(MULTIPOL_snk_image_id);
		close();

		if(save_files == false){
			selectImage(log2r_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_1_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_2_image_id);
			close();
	
			selectImage(snake_point_theta_image_id);
			close();				
		}



		selectImage(snake_point_index_image_id);
		for (i = 0; i < snake_points_x.length; i++) {
			run("Macro...", "code=[if(v == "+toString(i)+") v="+toString(theta[i])+";]");			
		}
		rename("theta image");
		theta_image_id = getImageID();
		//close();

	//setBatchMode("exit and display");
	}


	
	selectImage(t_BIN_image_id);
	makeSelection("polyline", snake_points_x, snake_points_y);



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2P fitting and plotting //////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// ratiometric fitting of data obtained with a single polarization
	//if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		binned_theta_deg_array = newArray(181);

		binned_hor = newArray(181);
		binned_hor_counts = newArray(181);

		binned_ver = newArray(181);
		binned_ver_counts = newArray(181);

		binned_log2r = newArray(181);
		number_of_nonzero_theta_bins = 0;

/*

		for(i=0; i<snake_points_x.length; i++){
			binned_theta_deg = floor(theta[i]/PI*180);   // binned_theta_deg runs from -90 to +90 
			if(hor_exists == true && log2r[i]==log2r[i]){
				binned_hor[binned_theta_deg + 90] = binned_hor[binned_theta_deg + 90] + log2r[i];  // I am increasing the indices by 90 just so that they would not be negative
				binned_hor_counts[binned_theta_deg + 90] = binned_hor_counts[binned_theta_deg + 90] + 1;
			}
			if(ver_exists == true && log2r[i]==log2r[i]){
				binned_ver[binned_theta_deg + 90] = binned_ver[binned_theta_deg + 90] + log2r[i];
				binned_ver_counts[binned_theta_deg + 90] = binned_ver_counts[binned_theta_deg +90] + 1;
			}
		}
		for(i=0; i<=180; i++){
			binned_theta_deg_array[i] = i-90;
			if(binned_hor_counts[i] > 0){ binned_hor[i] = binned_hor[i]/binned_hor_counts[i]; }
			if(binned_ver_counts[i] > 0){ binned_ver[i] = binned_ver[i]/binned_ver_counts[i]; }
			if(binned_hor_counts[i] > 0 || binned_ver_counts[i] > 0){ number_of_nonzero_theta_bins++; }
		}


*/
	
	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		if(hor_exists == true){
			for(i=0; i<=90; i++){binned_ver[i+90] = binned_hor[i]; }
			for(i=90; i<=180; i++){binned_ver[i-90] = binned_hor[i]; }
		}
		if(ver_exists == true){
			for(i=0; i<=90; i++){binned_hor[i+90] = binned_ver[i]; }
			for(i=90; i<=180; i++){binned_hor[i-90] = binned_ver[i]; }
		}
		for(i=0; i<=180; i++){
			if(binned_hor[i] > 0 && binned_ver[i] > 0){
				binned_log2r[i] = log(  binned_hor[i]/binned_ver[i]  )/log(2);
				x = Array.concat(x,i/180*PI - PI/2);
				y = Array.concat(y,binned_log2r[i]);
				n = n + 1;			
			}
		}
		hor_exists = true;
		ver_exists = true;
	}

	Array.getStatistics(y, y_min, y_max, y_mean, y_stdev); 

	index_of_param_cos2theta = -1;
	index_of_param_cos4theta = -1;
	index_of_param_phase = -1;
	index_of_param_phase2 = -1;
	index_of_param_normalization_factor = -1;
	index_of_param_offset = -1;

	//decide whether the fit should look for phase shift or just Rmax

	n_of_snake_points_with_positive_theta = 0;
	n_of_snake_points_with_negative_theta = 0;
	fraction_of_snake_points_with_data = n/snake_points_x.length;

	rmax = 0;

	// Make a header if there are two or more polarizations...
	if ((hor_exists == true && ver_exists == true) || (t_number_of_polarizations == 1 && ratiometric_fitting == true)){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tr \tpixel_count";}
	// Make a header if there is only one polarization...
	else{
		if(hor_exists == true && ratiometric_fitting == false){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tFh_avg \tpixel_count";}
		if(ver_exists == true && ratiometric_fitting == false){data_for_fitting_text = "theta \thorizontal_polarization_intensity \tvertical_polarization_intensity \tx_coordinate \ty_coordinate \tFv_avg \tpixel_count";}
	}


	// List data used for fitting, if there is only one polarization, but ratiometric fitting (which relies on theta binning) is being used 
	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		for(i=0; i<=180; i++){
			if(binned_hor[i] > 0 && binned_ver[i] > 0){
				data_for_fitting_text = data_for_fitting_text + "\r\n" + toString((binned_theta_deg_array[i]/180*PI))+" \t"+toString(binned_hor[i])+" \t"+toString(binned_ver[i])+" \t"+" "+" \t"+" "+" \t"+toString((binned_hor[i]/binned_ver[i]))+" \t"+" ";
			}
		}				
	}
	if(t_number_of_polarizations == 1 && ratiometric_fitting == false){
		for(i=0; i<snake_points_x.length; i++){
			data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta[i])+" \t"+toString(hor_sum[i])+" \t"+toString(ver_sum[i])+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r[i])+" \t"+toString(pixel_count[i]);	 		
		}		
	}
	if(t_number_of_polarizations == 2){
		for(i=0; i<snake_points_x.length; i++){
			data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta[i])+" \t"+toString(hor_sum[i])+" \t"+toString(ver_sum[i])+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r[i])+" \t"+toString(pixel_count[i]);	 		
		}		
	}
	if(t_number_of_polarizations > 2){
		if(save_files == true){
			for(i=0; i<snake_points_x.length; i++){
				for(j=0; j<number_of_polarizations; j++){
					selectImage(log2r_image_id);
					log2r = getPixel(i, j);
					if(isNaN(log2r)){}
					else{
						r = pow(2, log2r);
						selectImage(snake_point_theta_image_id);
						theta = getPixel(i, j);
						selectImage(snake_point_intensity_pol_1_image_id);
						hor_sum = getPixel(i, j);
						selectImage(snake_point_intensity_pol_2_image_id);
						ver_sum = getPixel(i, j);
						data_for_fitting_text = data_for_fitting_text + "\r\n" + toString(theta)+" \t"+toString(hor_sum)+" \t"+toString(ver_sum)+" \t"+toString(snake_points_x[i])+" \t"+toString(snake_points_y[i])+" \t"+toString(r)+" \t"+toString(pixel_count[i]);	 			
					}
				}				
			}	
			selectImage(log2r_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_1_image_id);
			close();
	
			selectImage(snake_point_intensity_pol_2_image_id);
			close();
	
			selectImage(snake_point_theta_image_id);
			close();
		}					
	}
	

	if(t_number_of_polarizations == 1 && ratiometric_fitting == true){
		if(phase_fixed == false){look_for_pol_angles = 1;}
	}

	r_square = NaN;





	while(r_square != r_square){
		//print("t_number_of_polarizations = "+toString(t_number_of_polarizations));
		if(hor_exists == true && ver_exists == true){		//2P with 2 polarizations

			if(t_number_of_polarizations > 2){
				number_of_fitting_params = 4;
				initialGuesses = newArray(1, 0.1, 0, 0);
				offset_string = "d + ";
				index_of_param_offset = 3;				
		
				fitting_equation = "y = "+offset_string+"log(     (1 + a*cos(2*x-2*c) + b*cos(4*x-4*c))         /           (1 - a*cos(2*x-2*c) + b*cos(4*x-4*c))       )/log(2)";
				index_of_param_cos2theta = 0;
				index_of_param_cos4theta = 1;
				index_of_param_phase = 2;
				index_of_param_phase2 = -1;
				index_of_param_normalization_factor = -1;
				phase_suffix = " (from fit)";
			}

			if(t_number_of_polarizations <= 2){
				if(phase_fixed == true){         // look for 0 phases
					if(look_for_pol_angles == 1){
						phase2_deg = phase1_deg;
						polarization_purity_ver = polarization_purity_hor;
					}
					sign = "-";
					if(phase1_deg < 0){sign = "+";}
					b_phase_string = sign+toString(abs(2*phase1_deg/180*PI));
					c_phase_string = sign+toString(abs(4*phase1_deg/180*PI));			
					if(phase1_deg == 0){
						b_phase_string = "";
						c_phase_string = "";
					}
	
					if(number_of_nonzero_theta_bins < 90){
						number_of_fitting_params = 2;
						initialGuesses = newArray(0, 0);
						offset_string = "";
						index_of_param_offset = -1;
						offset = 0;
					}
					else{
						number_of_fitting_params = 3;
						initialGuesses = newArray(0, 0, 0);
						offset_string = "c + ";
						index_of_param_offset = 2;				
					}			
					fitting_equation = "y = "+offset_string+"log(   (1 + a*cos(2*x"+b_phase_string+") + b*cos(4*x"+c_phase_string+") )    /    (1 - a*cos(2*x"+b_phase_string+") + b*cos(4*x"+c_phase_string+")  )     )/log(2)";
					index_of_param_cos2theta = 0;
					index_of_param_cos4theta = 1;
					index_of_param_phase = -1;
					index_of_param_phase2 = -1;
					index_of_param_normalization_factor = -1;
					phase_suffix = " (fixed)";
				}
				if(phase_fixed == false && look_for_pol_angles == 1){         			// look for 1 phase
					if(number_of_nonzero_theta_bins < 90){
						number_of_fitting_params = 3;
						initialGuesses = newArray(1, 0.1, 0);
						offset_string = "";
						index_of_param_offset = -1;
						offset = 0;
					}
					else{
						number_of_fitting_params = 4;
						initialGuesses = newArray(1, 0.1, 0, 0);
						offset_string = "d + ";
						index_of_param_offset = 3;				
					}
					fitting_equation = "y = "+offset_string+"log(     (1 + a*cos(2*x-2*c) + b*cos(4*x-4*c))         /           (1 - a*cos(2*x-2*c) + b*cos(4*x-4*c))       )/log(2)";
					index_of_param_cos2theta = 0;
					index_of_param_cos4theta = 1;
					index_of_param_phase = 2;
					index_of_param_phase2 = -1;
					index_of_param_normalization_factor = -1;
					phase_suffix = " (from fit)";
				}
				if(phase_fixed == false && look_for_pol_angles == 2){         			// look for 2 phases
					if(number_of_nonzero_theta_bins < 90){
						number_of_fitting_params = 4;
						initialGuesses = newArray(1, 1, 0, 0);
						offset_string = "";
						index_of_param_offset = -1;
						offset = 0;
					}
					else{
						number_of_fitting_params = 5;
						initialGuesses = newArray(1, 1, 0, 0, 0);
						offset_string = "e + ";
						index_of_param_offset = 4;				
					}
					fitting_equation = "y = "+offset_string+"log(     (1 + a*cos(2*x-2*c) + b*cos(4*x-4*c))         /           (1 - a*cos(2*x-2*d) + b*cos(4*x-4*d))   )/log(2)";
					index_of_param_cos2theta = 0;
					index_of_param_cos4theta = 1;
					index_of_param_phase = 2;
					index_of_param_phase2 = 3;
					index_of_param_normalization_factor = -1;
					phase_suffix = " (from fit)";
				}
			}

		}
		if(hor_exists == false || ver_exists == false){		//2P with 1 polarizations
			if(look_for_pol_angles > 1){look_for_pol_angles = 1;}
			if(look_for_pol_angles == 0){         			// look for 0 phases
				fitting_equation = "y = c * (1 + a*(2*("+toString(polarization_purity_hor/100)+")  -1)*cos(2*x-2*("+phase+")) + b*cos(4*x-4*("+phase+"))       )";
				number_of_fitting_params = 3;
				initialGuesses = newArray(1, 0.1, 1);
				index_of_param_cos2theta = 0;
				index_of_param_cos4theta = 1;
				index_of_param_normalization_factor = 2;				
				index_of_param_phase = -1;
				index_of_param_phase2 = -1;
				index_of_param_offset = -1;
				phase_suffix = " (fixed)";
			}
			if(look_for_pol_angles == 1){         			// look for 1 phase
				fitting_equation = "y = d * (1 + a*(2*("+toString(polarization_purity_hor/100)+")  -1)*cos(2*x-2*c) + b*cos(4*x-4*c)       )";
				number_of_fitting_params = 4;
				initialGuesses = newArray(1, 1, 0, 1);
				index_of_param_cos2theta = 0;
				index_of_param_cos4theta = 1;
				index_of_param_phase = 2;
				index_of_param_normalization_factor = 3;				
				index_of_param_phase2 = -1;
				index_of_param_offset = -1;
				phase_suffix = " (from fit)";
			}
		}
		
		showStatus("Fitting 2P data...");
	
		//Fit.logResults;
		Fit.doFit(fitting_equation, x, y);
	
		r_square = Fit.rSquared;

		if(r_square != r_square){
			Dialog.create("Fitting failure");
			Dialog.addMessage("The data fitting algorithm failed to converge. In order \nto aid the fitting procedure, please enter the values of deviation \nof the used polarizations from the horizontal/vertical directions.\n "); 
			Dialog.addNumber("Deviation from horizontal direction (in degrees):", phase1_deg);
			Dialog.addNumber("Deviation from vertical direction (in degrees):", phase2_deg);
			Dialog.show();
			phase1_deg = Dialog.getNumber();
			phase2_deg = Dialog.getNumber();
			phase_fixed = true;
			
		}
		if(r_square == r_square && r_square < 0.8){
			showMessage("Warning: The identified fit is not very good (r"+suptwochar+" = "+r_square+"). If you believe \na better fit exists, please run this script again, with values of deviation \nof the used polarizations from the horizontal/vertical directions set manually.\n ");
			
		}

	
	}


	B2P_from_fit = Fit.p(index_of_param_cos2theta);
	C2P_from_fit = Fit.p(index_of_param_cos4theta);


	if(index_of_param_phase  > -1){phase  = Fit.p(index_of_param_phase);}	
	if(index_of_param_phase2 > -1){phase2 = Fit.p(index_of_param_phase2);}
	else{phase2 = phase;}
	if(index_of_param_offset > -1){offset = Fit.p(index_of_param_offset);}
	if(index_of_param_normalization_factor > -1){normalization_factor = Fit.p(index_of_param_normalization_factor);}


	if(index_of_param_phase > -1){
		phase  = phase  % PI;
		phase2 = phase2 % PI;

		if (phase  < 0){phase  = phase  + PI;}		// from here on, phase should be between 0 and Pi
		if (phase2 < 0){phase2 = phase2 + PI;}		// from here on, phase2 should be between 0 and Pi	 
		if (phase  > PI/2){phase  = phase  - PI;}	// from here on, phase should be between -Pi/2 and Pi/2 
		if (phase2 > PI/2){phase2 = phase2 - PI;}	// from here on, phase2 should be between -Pi/2 and Pi/2 
		
		if (phase > PI/4){
			phase = phase - PI/2; 	
			phase2 = phase2 - PI/2;			
			B2P_from_fit = - B2P_from_fit;
		}
		
		if (phase < -PI/4){
			phase = phase + PI/2;
			phase2 = phase2 + PI/2;
			B2P_from_fit = - B2P_from_fit;
		}
		
		if (phase2 >  PI/2){phase2 = phase2 - PI;}	 
		if (phase2 < -PI/2){phase2 = phase2 + PI;}
		
		phase1_deg = round((phase/PI*180)*10)/10;
		phase_text = toString(phase1_deg)+degreechar;

		if(index_of_param_phase2 > -1){
			phase2_deg = round((phase2/PI*180)*10)/10;
			phase_text = phase_text+"/"+toString(phase2_deg)+degreechar;
		}

	}
	else{
		if(look_for_pol_angles == 2){
			phase_text = toString(phase1_deg)+degreechar+"/"+toString(phase2_deg)+degreechar+" (fixed)";	
		}
		if(look_for_pol_angles == 1){
			phase_text = toString(phase1_deg)+degreechar+" (fixed)";	
		}
		if(look_for_pol_angles == 0){
			if(phase1_deg != phase2_deg){phase_text = toString(phase1_deg)+degreechar+"/"+toString(phase2_deg)+degreechar+" (fixed)";}
			else{phase_text = toString(phase1_deg)+degreechar+" (fixed)";}		
		}
		phase = phase1_deg/180*PI;
		print("phase2_deg="+phase2_deg);
		phase2 = phase2_deg/180*PI;
	}

	rmax = (1 + B2P_from_fit + C2P_from_fit)/(1 - B2P_from_fit + C2P_from_fit);			
	log2_rmax = log(rmax)/log(2);


	results_of_fitting_text = "Results of 2PPM fitting: \r\nFitting equation: "+fitting_equation+"\r\n";
	if(	number_of_fitting_params == 3 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\nc = "+Fit.p(2)+"\r\n\r\n";}
	if(	number_of_fitting_params == 4 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\nc = "+Fit.p(2)+"\r\nd = "+Fit.p(3)+"\r\n\r\n";}
	if(	number_of_fitting_params == 5 ){results_of_fitting_text = results_of_fitting_text + "a = "+Fit.p(0)+"\r\nb = "+Fit.p(1)+"\r\nc = "+Fit.p(2)+"\r\nd = "+Fit.p(3)+"\r\ne = "+Fit.p(4)+"\r\n\r\n";}

	

	if(hor_exists == true && ver_exists == true){
		if(index_of_param_phase2 == -1){
			results_of_fitting_text = results_of_fitting_text + "\r\noffset = "+offset+"\r\nphase = "+phase+" "+phase_suffix+"\r\nB2P = "+B2P_from_fit+"\r\nC2P = "+C2P_from_fit;	
			results_of_fitting_text = results_of_fitting_text + "\r\nnumber_of_polarizations = 2";
		}
		else{
			results_of_fitting_text = results_of_fitting_text + "\r\noffset = "+offset+"\r\nphase = "+phase+" "+phase_suffix+"\r\nphase2 = "+phase2+"\r\nB2P = "+B2P_from_fit+"\r\nC2P = "+C2P_from_fit;	
			results_of_fitting_text = results_of_fitting_text + "\r\nnumber_of_polarizations = 2";			
		}
	}
	else{
		results_of_fitting_text = results_of_fitting_text+"\r\nnormalization_factor = "+normalization_factor+"\r\nphase = "+phase+" "+phase_suffix+"\r\nB2P = "+B2P_from_fit+"\r\nC2P = "+C2P_from_fit+"\r\npolarization_purity = "+toString(polarization_purity/100)+"\r\n";						
		if(hor_exists == 1){
			results_of_fitting_text = results_of_fitting_text + "\r\npolarization_direction = horizontal";
			results_of_fitting_text = results_of_fitting_text + "\r\nnumber_of_polarizations = 1";
		}
		else{
			results_of_fitting_text = results_of_fitting_text + "\r\npolarization_direction = vertical";				
			results_of_fitting_text = results_of_fitting_text + "\r\nnumber_of_polarizations = 1";
		}
	}

	results_of_fitting_text = results_of_fitting_text+"\r\n\r\nr(max) \t1/r(max) \tlog"+subtwochar+"(r(max)) \tphase \tr"+suptwochar+" \r\n"+rmax+" \t"+1/rmax+" \t"+log2_rmax+" \t"+phase+" \t"+r_square+" \r\n";
	results_of_fitting_text = results_of_fitting_text+"--------------------------------------------------\r\n";



	//make a plot

	num_of_samples=100;
	fit_x_values = newArray(num_of_samples); 
	fit_y_values = newArray(num_of_samples); 
	xmin=-3.14159/2; 
	xmax= 3.14159/2;
	step=(xmax-xmin)/(num_of_samples-1);
	for (i=0; i< num_of_samples; i++) {
		fit_x_values[i]=xmin + i*step;
		fit_y_values[i]=Fit.f(fit_x_values[i]);
	}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2P alpha0, sigma /////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	setBatchMode(true);


	macro_path = getDirectory("macros");

	open(macro_path+"2PPM_B.tif");
	B_2PPM_image_id = getImageID();

	open(macro_path+"2PPM_C.tif");
	C_2PPM_image_id = getImageID();

	num_of_alpha0s = getWidth;
	num_of_sigmas  = getHeight;

	newImage(file_name_base_f+"_2P_RMSD_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	RMSD_image_id = getImageID();

	newImage(file_name_base_f+"_2P_chi_sq_a0_s.tif", "32-bit black", num_of_alpha0s, num_of_sigmas, 1);
	chi_sq_image_id = getImageID();


	showStatus("Calculating "+alpha0char+", "+sigmachar);


	for(i = 0; i < n; i++){
		showProgress(i, n);
		if(hor_exists == 1 && ver_exists == 1){expectation_expression = toString(offset) + " + log( (1 + B * cos(" + toString(2*x[i]-2*phase) + ") + C * cos(" + toString(4*x[i]-4*phase) + ")) / (1 - B * cos(" + toString(2*x[i]-2*phase2) + ") + C * cos(" + toString(4*x[i]-4*phase2) + ") ) ) / log(2)";}
		else{expectation_expression = toString(normalization_factor) + " *  (1 + B * "+toString(polarization_purity/100)+" * cos(" + toString(2*x[i]-2*phase) + ") + C * "+toString(polarization_purity/100)+" * cos(" + toString(4*x[i]-4*phase) + ")           + B * (1-"+toString(polarization_purity/100)+") * cos(" + toString(2*x[i]-2*(phase-3.14159/2)) + ") + C * (1-"+toString(polarization_purity/100)+") * cos(" + toString(4*x[i]-4*(phase-3.14159/2)) + ")               ) ";}


		run("Image Expression Parser (Macro)", "expression=["+ expectation_expression +"] a=None b=2PPM_B.tif c=2PPM_C.tif d=None e=None f=None g=None h=None i=None j=None k=None l=None m=None n=None o=None p=None q=None r=None s=None t=None u=None v=None w=None x=None y=None z=None");
		expectation_image_id = getImageID();
		rename("expectation image");
			
		//RMSD
		run("Duplicate...", "title=sq_dev");
		sq_dev_image_id = getImageID();
		run("Subtract...", "value="+toString(y[i]));
		run("Square");
		imageCalculator("Add 32-bit", RMSD_image_id, sq_dev_image_id); 
		showStatus("Calculating "+alpha0char+", "+sigmachar);

		//Chi-squared
		selectImage(expectation_image_id);
		run("Abs");
		run("Add...", "value=0.01");
		imageCalculator("Divide 32-bit", sq_dev_image_id, expectation_image_id);
		imageCalculator("Add 32-bit", chi_sq_image_id, sq_dev_image_id);
		showStatus("Calculating "+alpha0char+", "+sigmachar);
		
		selectImage(sq_dev_image_id);
		close();
			
		selectImage(expectation_image_id);
		close();
	}

	selectImage(chi_sq_image_id);
	run("Divide...", "value="+toString(n*1.2));
	getRawStatistics(nPixels, mean, min_chi_sq, max_chi_sq);
	for (col = 0; col <= num_of_alpha0s; col++) {
		showProgress(col, num_of_alpha0s);
		for (row = 0; row <= num_of_sigmas; row++){
			if(getPixel(col,row) != getPixel(col,row)){
				setPixel(col,row,max_chi_sq);
			}
		}
	}

	run("phase");
	setMinAndMax(0, 1);



	selectImage(RMSD_image_id);
	run("Divide...", "value="+n);
	run("Square Root");
	getRawStatistics(nPixels, mean, min_RMSD_from_alpha0_sigma_search, max_RMSD_from_alpha0_sigma_search);

	for (col = 0; col <= num_of_alpha0s; col++) {
		showProgress(col, num_of_alpha0s);
		for (row = 0; row <= num_of_sigmas; row++){
			if(getPixel(col,row) != getPixel(col,row)){
				setPixel(col,row,max_RMSD_from_alpha0_sigma_search);
			}
		}
	}

	run("phase");
	setMinAndMax(sqrt(1-r_square*r_square)*y_stdev, 5*sqrt(1-r_square*r_square)*y_stdev);

	run("Find Maxima...", "noise="+max_RMSD_from_alpha0_sigma_search+" output=[Point Selection] light");
	getSelectionBounds(alpha0, sigma_row_number, w, h);
	sigma = num_of_sigmas - 1 - sigma_row_number;


	run("Select None");

	selectImage(B_2PPM_image_id);
	alpha0_B = getPixel(alpha0,sigma_row_number);

	selectImage(C_2PPM_image_id);
	alpha0_C = getPixel(alpha0,sigma_row_number);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////   add alpha0, sigma curve to the plot   //////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	alpha0_fit_x = newArray(100);
	alpha0_fit_y = newArray(100);

	for(i = 0; i < 100; i++){
		alpha0_fit_x[i] = -3.14159/2 + i/99 * 3.14159;  //makes alpha0_fit_x values from -Pi/2 to Pi/2
		if(hor_exists == 1 && ver_exists == 1){
			alpha0_fit_y[i] = (offset + log((1 +alpha0_B*cos(2*alpha0_fit_x[i]-2*phase) + alpha0_C*cos(4*alpha0_fit_x[i]-4*phase))/(1 + alpha0_B*cos(2*(alpha0_fit_x[i]-3.14159/2)-2*phase2) + alpha0_C*cos(4*(alpha0_fit_x[i]-3.14159/2)-4*phase2)))/log(2));
		}
		else{
			alpha0_fit_y[i] = normalization_factor*(1 + alpha0_B*cos(2*alpha0_fit_x[i]-2*phase) + alpha0_C*cos(4*alpha0_fit_x[i]-4*phase));
		}
	}


	results_of_fitting_text = results_of_fitting_text + "alpha0 = "+alpha0+degreechar+"\r\n";
	results_of_fitting_text = results_of_fitting_text + "sigma = "+sigma+degreechar+"\r\n";
	print(results_of_fitting_text);




////////////////////////////////////////////////////////////////////////////////////
//////////////////////// make hyperstack of results ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



	hyperstack_name = file_name_base_f+"_2P_goodness_of_fit";
	
	newImage(hyperstack_name, "32-bit color-mode", num_of_alpha0s, num_of_sigmas, 4, 1, 1);
	hyperstack_image_id = getImageID();
	setColor("orange");

	selectImage(RMSD_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(2);
	Overlay.drawString("2P RMSD", 5, 20);
	Overlay.setPosition(2,1,1);
	Overlay.show;	
	run("Set Label...", "label=[2P RMSD]");
	run("phase");
	run("Paste");
	setMinAndMax(sqrt(1-r_square*r_square)*y_stdev, 5*sqrt(1-r_square*r_square)*y_stdev);

	selectImage(RMSD_image_id);
	run("Select All");
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(1);
	Overlay.drawString("2P r"+suptwochar, 5, 20);
	Overlay.setPosition(1,1,1);
	Overlay.show;	
	run("Set Label...", "label=[2P r"+suptwochar+"]");
	run("Paste");
	run("Divide...", "value="+y_stdev+" slice");
	run("Square", "slice");
	run("Multiply...", "value=-1 slice");
	run("Add...", "value=1 slice");
	run("Cyan Hot");
	setMinAndMax(0, 1);
	getRawStatistics(nPixels, mean_r_square_from_alpha0_sigma_search, min_r_square_from_alpha0_sigma_search, max_r_square_from_alpha0_sigma_search);

	
	selectImage(chi_sq_image_id);
	run("Copy");
	selectImage(hyperstack_image_id);
	setSlice(3);
	Overlay.drawString("2P Chi"+suptwochar, 5, 20);
	Overlay.setPosition(3,1,1);
	Overlay.show;
	run("Set Label...", "label=[2P Chi"+suptwochar+"]");
	run("Paste");
	run("phase");
	setMinAndMax(0, 1);

	selectImage(hyperstack_image_id);
	setSlice(4);	
	Overlay.drawString("2P "+alpha0char+", "+sigmachar, 5, 20);
	Overlay.setPosition(4,1,1);
	Overlay.show;
	run("Set Label...", "label=[2P alpha0, sigma]");
	setPixel(alpha0, sigma_row_number, 255);
	run("Grays");
	
	setSlice(1);
	setColor(139,0,21);  // = #8b0015, used for the alphs0, sigma prediction in the plot
	Overlay.setStrokeWidth(1);
	getLine(x1, y1, x2, y2, orig_line_width);
	setLineWidth(1);
	Overlay.drawRect(alpha0, sigma_row_number, 1, 1);
	Overlay.show;
	setLineWidth(orig_line_width);

	run("Properties...", "channels=4 slices=1 frames=1 unit="+degreechar+" pixel_width=1.0000 pixel_height=-1 voxel_depth=1.0000 origin=0,90");
	run("Select None");

	selectImage(RMSD_image_id);
	close();
	selectImage(chi_sq_image_id);
	close();	

//////////////////////////////////////////////////////////////////////////////////////



//

	Array.getStatistics(fit_y_values, fit_y_min, fit_y_max);




	if (hor_exists == true && ver_exists == true){
		if(toString(fit_y_min) == "Infinity" || toString(fit_y_min) == "-Infinity" || toString(fit_y_max) == "Infinity" || toString(fit_y_max) == "-Infinity"){
			lower_plot_limit = y_min * 1.1;
			upper_plot_limit = y_max * 1.1;
		}
		else{
			upper_plot_limit = y_max * 1.1;
			if(abs(log2_rmax) > y_max){ upper_plot_limit = abs(log2_rmax) * 1.1; }
			lower_plot_limit = y_min * 2;
			if(-abs(log2_rmax) < y_min){ lower_plot_limit = -abs(log2_rmax) * 1.6;}
		}
		Plot.create(file_name_base_f+"_2P_FIT", "Membrane orientation, "+thetachar+" ["+degreechar+"]", "Log"+subtwochar+"(r)");	
	}
	else{
		upper_plot_limit = y_max * 1.5;
		lower_plot_limit = 0;

		Plot.create(file_name_base_f+"_2P_FIT", "Membrane orientation ("+thetachar+")", "F");
	}

	// print("Lower plot limit: "+lower_plot_limit+", upper plot limit: "+upper_plot_limit);


	Plot.setLimits(-90, 90, lower_plot_limit, upper_plot_limit);

	Plot.setLineWidth(1);
	Plot.setColor("gray");
	x_deg = newArray(x.length);
	for(i = 0; i < x.length; i++){
		x_deg[i] = x[i]/3.14159*180;
	}

	Plot.add("circle", x_deg, y);
	plot_legend = "\n";

	Plot.setLineWidth(1.5);
	Plot.setColor("#8b0015");
//dark red/marroon
	fit_x_values_deg = newArray(fit_x_values.length);
	for(i = 0; i < fit_x_values.length; i++){
		fit_x_values_deg[i] = fit_x_values[i]/3.14159*180;
	}
	Plot.add("line", fit_x_values_deg, fit_y_values);

	if(hor_exists == 1 && ver_exists == 1){
		//plot_legend = plot_legend + "r(max) = "+rmax+" ("+1/rmax+"), phase = "+phase_text+", B = "+B2P_from_fit+", C = "+B2P_from_fit+", r"+suptwochar+" = "+r_square+"\n";	
		plot_legend = plot_legend + "log2(r(max)) = "+log2_rmax+", phase = "+phase_text+", B = "+B2P_from_fit+", C = "+B2P_from_fit+", r"+suptwochar+" = "+r_square+"\n";	
	}
	else{
		//plot_legend = plot_legend + "r(max) = "+rmax+" ("+1/rmax+"), phase = "+toString(round(10*phase/3.14159*180)/10)+degreechar+", B = "+B2P_from_fit+", C = "+C2P_from_fit+", r"+suptwochar+" = "+r_square+"\n";			
		plot_legend = plot_legend + "log2(r(max)) = "+log2_rmax+", phase = "+toString(round(10*phase/3.14159*180)/10)+degreechar+", B = "+B2P_from_fit+", C = "+C2P_from_fit+", r"+suptwochar+" = "+r_square+"\n";			
	}



	Plot.setLineWidth(1.5);
	Plot.setColor("Red");
	alpha0_fit_x_values_deg = newArray(alpha0_fit_x.length);
	for(i = 0; i < alpha0_fit_x.length; i++){
		alpha0_fit_x_values_deg[i] = alpha0_fit_x[i]/3.14159*180;
	}

	Plot.add("line", alpha0_fit_x_values_deg, alpha0_fit_y);

	plot_legend = plot_legend + alpha0char+" = "+alpha0+degreechar+", "+sigmachar+" = "+sigma+degreechar+", B = "+alpha0_B+", C = "+alpha0_C+", r"+suptwochar+" = "+max_r_square_from_alpha0_sigma_search;

	Plot.addLegend(plot_legend, "Auto");

	Plot.setFormatFlags("11000100110011");
	Plot.show;
	plot_window_image_id = getImageID();



	selectImage(B_2PPM_image_id);
	close();

	selectImage(C_2PPM_image_id);
	close();


	print("Finding "+alpha0char+", "+sigmachar+":");
	print("("+alpha0char+", "+sigmachar+") = ("+alpha0+degreechar+", "+sigma+degreechar+")\r\nr"+suptwochar+" = "+max_r_square_from_alpha0_sigma_search); 
	print("B = "+alpha0_B+"\r\nC = "+alpha0_C);

	results_of_2P_fitting_array = newArray(5);
	results_of_2P_fitting_array[0] = toString(plot_window_image_id);
	results_of_2P_fitting_array[1] = toString(hyperstack_image_id);
	results_of_2P_fitting_array[2] = toString(theta_image_id);
	results_of_2P_fitting_array[3] = data_for_fitting_text;
	results_of_2P_fitting_array[4] = results_of_fitting_text;

	return results_of_2P_fitting_array;	

}





function guess_multipol_params(multipol_stack_id){
	w = getWidth;
	h = getHeight;
	number_of_polarizations = nSlices;
	

	selectImage(multipol_stack_id);
	run("Duplicate...", "duplicate");
	temp_multipol_stack_id = getImageID();
	
	run("Z Project...", "projection=[Average Intensity]");
	rename("Avg_F");
	temp_avg_F_image_id = getImageID();

	selectImage(multipol_stack_id);
	run("Z Project...", "projection=[Min Intensity]");
	rename("Min_F");
	temp_min_F_image_id = getImageID();

	selectImage(multipol_stack_id);
	run("Z Project...", "projection=[Max Intensity]");
	rename("Max_F");
	temp_max_F_image_id = getImageID();



	

	min_mean = 100000;
	for(x = w/20; x < w - w/10; x = x + w/20){
		for(y = h/20; y < h - h/10; y = y + h/20){
			makeRectangle(x, y, w/10, h/10);
			getStatistics(area, mean);
			if(mean < min_mean){
				min_mean = mean;
				min_x = x;
				min_y = y;
			}
		}
	}

	bkgd = round(min_mean);
	
	selectImage(temp_multipol_stack_id);
	run("Subtract...", "value="+bkgd+" stack");
	run("Add...", "value=1 stack");



	//---overall brightness adjustment correction ----------------------------------------------------------
	selectImage(temp_avg_F_image_id);
	run("Subtract...", "value="+min_mean);

	getMinAndMax(avg_image_min_value,avg_image_max_value);
	bright_adj = round((avg_image_max_value-avg_image_min_value)/1.2);

	//  make a mask ------------------------------------------------------------------------
	selectImage(temp_avg_F_image_id);
	setAutoThreshold();
	run("Convert to Mask");
	run("Divide...", "value=255");
	rename("mask");
	temp_mask_image_id = getImageID();

	imageCalculator("multiply stack",temp_multipol_stack_id, temp_mask_image_id);
	imageCalculator("multiply", temp_min_F_image_id, temp_mask_image_id);
	imageCalculator("multiply", temp_max_F_image_id, temp_mask_image_id);

	imageCalculator("Divide create 32-bit", temp_max_F_image_id, temp_min_F_image_id);
	run("Log");
	run("Divide...", "value=0.693");
	temp_log2r_image_id = getImageID();
	rename("log2r");


	//-------------------------------------------------------
	//---R/G range calculation ----------------------------------------------------------


	getStatistics(area, mean_rg_logratio, min_rg_logratio, max_rg_logratio, logratio_stdev);

	if(range_lock == true){rg_range = locked_rg_range;}
	else{
	
		log_halfrange = max_rg_logratio - mean_rg_logratio;

		total_counts = 0;
		
		getHistogram(hist_values,hist_counts,100);
		for (index = 0; index < 100; index++){
			total_counts = total_counts + hist_counts[index];
		}
		cutoff_counts = total_counts * 0.1;
	
		incr = 0;
		out_of_range_percentage = 0;
		while(out_of_range_percentage < 0.025){
			lower_limit = mean_rg_logratio - log_halfrange * (1 - incr/100);
			upper_limit = mean_rg_logratio + log_halfrange * (1 - incr/100);
			num_of_out_of_range = 0;
			for (index = 0; index < 50; index++){
				
				if ((hist_values[index] < lower_limit) || (hist_values[index] > upper_limit)){
					num_of_out_of_range = num_of_out_of_range + hist_counts[index];
				}
				if ((hist_values[99-index] < lower_limit) || (hist_values[99-index] > upper_limit)){
					num_of_out_of_range = num_of_out_of_range + hist_counts[99-index];
				}
			}
			out_of_range_percentage = num_of_out_of_range/total_counts;
			incr++;
		}

		rg_range = (round((exp(mean_rg_logratio-lower_limit))*100))/100;
	}

	selectImage(temp_multipol_stack_id);
	close();

	selectImage(temp_log2r_image_id);
	close();
	selectImage(temp_mask_image_id);
	close();
	selectImage(temp_min_F_image_id);
	close();
	selectImage(temp_max_F_image_id);
	close();


	function_array = newArray(3);
	function_array[0] = bkgd;
	function_array[1] = bright_adj;
	function_array[2] = rg_range;


	return function_array;
}


function make_processed_multipol_img(temp_multipol_stack_image_id, temp_multipol_params, log2r_azimuth_r_squared_image_id, make_new_log2r_azimuth_image) {
	setBatchMode(true);

/*
		bkgd 							 = new_processing_params[0];
		bright_adj      				 = new_processing_params[1];
		show_LD  						 = new_processing_params[2];
		color_scheme 					 = new_processing_params[3];
		rg_range 						 = new_processing_params[4];
		F_cutoff_for_LD_display			 = new_processing_params[5];
		show_azimuths        			 = new_processing_params[6];
		azimuth_color					 = new_processing_params[7];
		starting_pol_angle				 = new_processing_params[8];
		pol_angle_increment				 = new_processing_params[9];
		azimuth_pixel_bin_size			 = new_processing_params[10];
		F_cutoff_for_azimuth_display	 = new_processing_params[11];
		use_fitting_for_LD				 = new_processing_params[12];

 */


	temp_bkgd = temp_multipol_params[0];
	temp_bright_adj = temp_multipol_params[1];
	temp_show_LD = temp_multipol_params[2];
	temp_rg_range = temp_multipol_params[4];
	temp_LD_cutoff = temp_multipol_params[5];
	temp_show_azimuths = temp_multipol_params[6];
	temp_azimuth_color_scheme = temp_multipol_params[7];
	temp_azimuth_bin_size = temp_multipol_params[10];
	temp_azimuth_scaling_factor = temp_multipol_params[11];
	temp_azimuth_cutoff = temp_multipol_params[12];
	temp_use_fitting_for_LD = temp_multipol_params[13];

	polarization_list = newArray(number_of_polarizations);
	for (i = 0; i < polarization_list.length; i++){
		polarization_list[i] = starting_pol_angle + i * pol_angle_increment;
	}

	if(temp_azimuth_cutoff > temp_LD_cutoff){lower_of_LD_azimuth_cutoffs = temp_LD_cutoff;}
	else{lower_of_LD_azimuth_cutoffs = temp_azimuth_cutoff;}

	selectImage(temp_multipol_stack_image_id);
	temp_multipol_stack_name = getTitle();
	temp_file_name_base = substring(temp_multipol_stack_name,0,lengthOf(temp_multipol_stack_name)-4);


	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////    Crude/fast LD visualization     //////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	if(temp_show_LD == "Crude/fast"){
		// make temporal average image, use it to obtain max_F
		selectImage(temp_multipol_stack_image_id);
		run("Z Project...", "projection=[Average Intensity]");
		run("Subtract...", "value="+temp_bkgd);
		rename("Avg_F");
		temp_avg_F_image_id = getImageID();
		getStatistics(area, mean_F, min_F, max_F, std, histogram);
		setMinAndMax(0, temp_bright_adj - temp_bkgd);	

		// make a shrunk mask for LD
		selectImage(temp_avg_F_image_id);
		run("Duplicate...", "title=scaled_down_avg_F_image2");
		run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=Bilinear");
		run("Macro...", "code=[if(v < "+toString(temp_LD_cutoff * max_F)+") v=0]");
		temp_scaled_down_avg_F_image_id = getImageID();
		imageCalculator("Divide create 32-bit", temp_scaled_down_avg_F_image_id, temp_scaled_down_avg_F_image_id);
		rename("scaled_down_LD_mask_image");
		temp_scaled_down_LD_mask_image_id = getImageID();
		selectImage(temp_scaled_down_avg_F_image_id);
		close();

		// make a shrunk mask for azimuths
		if(temp_show_azimuths == "Yes"){		
			selectImage(temp_avg_F_image_id);
			run("Duplicate...", "title=scaled_down_avg_F_image3");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=Bilinear");
			run("Macro...", "code=[if(v < "+toString(temp_azimuth_cutoff * max_F)+") v=0]");
			temp_scaled_down_avg_F_image_id = getImageID();
			imageCalculator("Divide create 32-bit", temp_scaled_down_avg_F_image_id, temp_scaled_down_avg_F_image_id);
			rename("scaled_down_azimuth_mask_image");
			temp_scaled_down_azimuth_mask_image_id = getImageID();
			selectImage(temp_scaled_down_avg_F_image_id);
			close();
		}
		
		// make a shrunk mask for LD/azimuths
		if(temp_show_azimuths == "Yes"){		
			selectImage(temp_avg_F_image_id);
			run("Duplicate...", "title=scaled_down_avg_F_image");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=Bilinear");
			run("Macro...", "code=[if(v < "+toString(lower_of_LD_azimuth_cutoffs * max_F)+") v=0]");
			temp_scaled_down_avg_F_image_id = getImageID();
			imageCalculator("Divide create 32-bit", temp_scaled_down_avg_F_image_id, temp_scaled_down_avg_F_image_id);
			rename("scaled_down_LD_azimuth_mask_image");
			temp_scaled_down_LD_azimuth_mask_image_id = getImageID();
			selectImage(temp_scaled_down_avg_F_image_id);
			close();
		}
		else{
			selectImage(temp_scaled_down_LD_mask_image_id);
			run("Duplicate...", "title=scaled_down_LD_azimuth_mask_image");
			temp_scaled_down_LD_azimuth_mask_image_id = getImageID();			
		}
		
		// make a shrunk t-stack, apply the shared LD/azimuth mask, reslice the stack			
		selectImage(temp_multipol_stack_image_id);				
		run("Duplicate...", "duplicate title=scaled_down_bkgd_subbed_stack");
		run("Subtract...", "value="+toString(temp_bkgd)+" stack");
		run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth="+toString(number_of_polarizations)+" interpolation=Bilinear");
		scaled_down_stack_image_id = getImageID();
		imageCalculator("Multiply create 32-bit stack", scaled_down_stack_image_id, temp_scaled_down_LD_azimuth_mask_image_id);
		temp_scaled_down_masked_stack_image_id = getImageID();
		rename("scaled_down_masked_stack");

		scaled_down_width = getWidth();
		scaled_down_height = getHeight();
		
		// reslice the shrunk t-stack so that intensity profile of a line selections = F(pol)		
		run("Reslice [/]...", "output=1.000 start=Top avoid");
		resliced_scaled_down_masked_stack_image_id = getImageID();

		makeLine(1, 1, 2, 2);
		getLine(x1, y1, x2, y2, orig_line_width);
		run("Select None");

		selectImage(scaled_down_stack_image_id);
		run("Duplicate...", "title=Fh_Fv_multiplier_composite");
		run("Add Slice");
		run("32-bit");
		run("Set...", "value=1 stack");
		Fh_Fv_multiplier_image_id = getImageID();

		// create a shrunk log2r_azimuth_r_squared 3-slice stack to hold the fitting results

		if(make_new_log2r_azimuth_image == true){   // make a new log2r_azimuth_r_squared image
			if(log2r_azimuth_r_squared_image_id < 0){
				selectImage(log2r_azimuth_r_squared_image_id);
				close();
			}
			selectImage(temp_scaled_down_LD_azimuth_mask_image_id);
			run("Duplicate...", "title=log2r_azimuth_r_squared_composite");
			run("Add Slice");
			run("Add Slice");
			log2r_azimuth_r_squared_image_id = getImageID();
		}

		selectImage(temp_scaled_down_masked_stack_image_id);

		// for all pixels in the shared LD/azimuth mask, do the cos^2 fitting
		// - store the fitting parameters in the shrunk log2r_azimuth_r_squared stack
		// - use the fitting parameters to make the Fh_Fv_multiplier composite
		fitting_equation = "y = a*a + b*b * cos(x/180*3.14159 - c*c)*cos(x/180*3.14159 - c*c)";	
			
		for (row = 0; row < scaled_down_height; row++){
			showStatus("Fitting row "+toString(row)+" of "+toString(scaled_down_height)+" ("+toString(round(row/scaled_down_height*100))+"%)");		
			showProgress(row/scaled_down_height);
			for (col = 0; col < scaled_down_width; col++) {
				selectImage(resliced_scaled_down_masked_stack_image_id);
				setSlice(row + 1);

				if(isNaN(getPixel(col, 0))){}
				else{
					if(make_new_log2r_azimuth_image == true){ 					
						makeLine(col, 0, col, number_of_polarizations-1);	
						run("Line Width...", "line=1"); 
					
						//run("Properties... ", " width=1");
						F_over_time = getProfile();
						Fit.doFit(fitting_equation, polarization_list, F_over_time);
						Fmin = Fit.p(0) * Fit.p(0);
						Fmax = Fmin + Fit.p(1) * Fit.p(1);
	
						log2r = log(Fmax/Fmin)/log(2);
						azimuth = (Fit.p(2) * Fit.p(2))/PI*180 - floor((Fit.p(2) * Fit.p(2))/PI*180/180)*180;
						rsquared = Fit.rSquared;
			
						selectImage(log2r_azimuth_r_squared_image_id);
						setSlice(1);
						setPixel(col, row, log2r);
						setSlice(2);
						setPixel(col, row, azimuth);
						setSlice(3);
						setPixel(col, row, rsquared);
	
						selectImage(Fh_Fv_multiplier_image_id);
						setSlice(1);
						setPixel(col, row, 1/sqrt(pow(2, log2r)));
						setSlice(2);
						setPixel(col, row, sqrt(pow(2, log2r)));
					}
					else{
						selectImage(log2r_azimuth_r_squared_image_id);
						setSlice(1);
						log2r = getPixel(col, row);

						selectImage(Fh_Fv_multiplier_image_id);
						setSlice(1);
						setPixel(col, row, 1/sqrt(pow(2, log2r)));
						setSlice(2);
						setPixel(col, row, sqrt(pow(2, log2r)));

					}					
				}
			}
		}
		selectImage(resliced_scaled_down_masked_stack_image_id);
		close();	
		run("Line Width...", "line="+orig_line_width); 
		//setLineWidth(orig_line_width);
	
		// modify the Fh_Fv_multiplier composite to have non-1 values only in pixels matchng the LD mask
		// (so far it has non-1 values in pixels matching the shared LD/azimuth mask)
		selectImage(Fh_Fv_multiplier_image_id);

		imageCalculator("Multiply stack", Fh_Fv_multiplier_image_id, temp_scaled_down_LD_mask_image_id);  // applies a mask to show LD only in the selected parts of the image
		run("Macro...", "code=[if(isNaN(v)) v=1] stack");	
		run("Size...", "width="+w+" height="+h+" depth=2 constrain interpolation=Bilinear");
		Fh_Fv_multiplier_image_id = getImageID();

	
		// Finally, use the Fh_Fv_multiplier composite to make a Fh, Fv composite that can be used to make an RG image
		imageCalculator("Multiply stack",Fh_Fv_multiplier_image_id, temp_avg_F_image_id);  // the Fh_Fv_multiplier composite should now hold Fh, Fv... but for a threshold that is the lower threshold of azimuth, LD threshold
		temp_composite_image_id = getImageID();
		rename(temp_file_name_base+"_");
	
		//... and turn the Fh, Fv composite into an RG image
		rg_image_id = make_rg_image(temp_composite_image_id, temp_bright_adj, temp_rg_range, color_scheme);		

		selectImage(temp_composite_image_id);
		//rename("composite_for_making_multipol_RG_image");
		close();
	
		// that should take care of the crude LD visualization... now the azimuths

		if(temp_show_azimuths == "Yes"){
			for (col = 0; col < scaled_down_width; col++) {
				showProgress(col/scaled_down_width);
				for (row = 0; row < scaled_down_height; row++) {
					
					pixel_x_coord = floor(w/scaled_down_width * col);// + floor(temp_azimuth_bin_size/2) + 1;
					pixel_y_coord = floor(h/scaled_down_height * row);// + floor(temp_azimuth_bin_size/2) + 1;
					/* I only need these next few lines if I want the azimuth sizes to correspond to fluorescence intensities
					// coordinates in the unscaled images:
					selectImage(temp_avg_F_image_id);
					rel_F = getPixel(pixel_x_coord, pixel_y_coord)/max_F;  
					*/
		
					// check that the investigated pixel matches the azimuth mask
					selectImage(temp_scaled_down_azimuth_mask_image_id);

					if(isNaN(getPixel(col, row))){}
					else{
						selectImage(log2r_azimuth_r_squared_image_id);
						setSlice(1);
						log2r = getPixel(col, row);
						if(isNaN(log2r)){log2r = 0;}
						
						rel_log2r = log2r/(log(temp_rg_range)/log(2))/2;
						if(rel_log2r > log(temp_rg_range)/log(2) ){rel_log2r = 1;}
						setSlice(2);
						azimuth = getPixel(col, row);
						if(isNaN(azimuth)){azimuth = 0;}
						
						//print("col: "+col+", row: "+row+", azimuth: "+azimuth);
						selectImage(rg_image_id);

						if(temp_azimuth_color_scheme == "direction"){

							//hue = 360 - azimuth * 2;
							hue = 300 - azimuth * 2;
							if(hue < 0){hue = hue + 360;}
							hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
							if(hue >= 0 && hue < 60){
								red = 255;
								green = round(hsv_x * 255);
								blue = 0;
							}
					
							if(hue >= 60 && hue < 120){
								red = round(hsv_x * 255);
								green = 255;
								blue = 0;
							}
					
							if(hue >= 120 && hue < 180){
								red = 0;
								green = 255;
								blue = round(hsv_x * 255);
							}
							if(hue >= 180 && hue < 240){
								red = 0;
								green = round(hsv_x * 255);
								blue = 255;
							}
					
							if(hue >= 240 && hue < 300){
								red = round(hsv_x * 255);
								green = 0;
								blue = 255;
							}
							if(hue >= 300 && hue < 360){
								red = 255;
								green = 0;
								blue = round(hsv_x * 255);
							}
							setColor(red, green, blue);

							
						}
						else{
							temp_azimuth_color = temp_azimuth_color_scheme;
							setColor(temp_azimuth_color);  //  overlay color
						}
						
						setLineWidth(2);
		
						// azimuth length = F intensity
						// Overlay.drawLine(pixel_x_coord - cos(azimuth/180*PI) * temp_azimuth_bin_size/2 * rel_F, pixel_y_coord - sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F, pixel_x_coord + cos(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F, pixel_y_coord + sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F);			
	
						// azimuth length = LD
						Overlay.drawLine(pixel_x_coord - cos(azimuth/180*PI) * temp_azimuth_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_y_coord - sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_x_coord + cos(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_y_coord + sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor);			
						Overlay.show;		

							
					}
				}
			}
			selectImage(temp_scaled_down_azimuth_mask_image_id);
			close();
		}		
	
	//	selectImage(temp_azimuth_image_id);
	//	close();
	
		selectImage(temp_avg_F_image_id);
		close();
		selectImage(temp_scaled_down_LD_mask_image_id);
		close();
		selectImage(temp_scaled_down_LD_azimuth_mask_image_id);
		close();
		selectImage(scaled_down_stack_image_id);
		close();
		selectImage(temp_scaled_down_masked_stack_image_id);
		close();
		//selectImage(log2r_azimuth_r_squared_image_id);
		//close();
	}

	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////   Accurate/slow LD visualization   //////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	if(temp_show_LD == "Accurate/slow"){

		// slow, accurate ----------------------------------------------------------------------

		//print("Slow and accurate LD visualization...");

		// make avgF image
		selectImage(temp_multipol_stack_image_id);
		run("Z Project...", "projection=[Average Intensity]");
		run("Subtract...", "value="+temp_bkgd);
		rename("Avg_F");
		temp_avg_F_image_id = getImageID();
		getStatistics(area, mean_F, min_F, max_F, std, histogram);


		if(make_new_log2r_azimuth_image == true){   // make a new log2r_azimuth_r_squared image
			if(log2r_azimuth_r_squared_image_id < 0){
				selectImage(log2r_azimuth_r_squared_image_id);
				close();
			}
			newImage("log2r_azimuth_r_squared_composite", "32-bit black", w, h, 3);
			log2r_azimuth_r_squared_image_id = getImageID();
		}



		if(temp_show_azimuths == "Yes"){		
			selectImage(temp_avg_F_image_id);
			run("Duplicate...", "title=scaled_down_avg_F_image3");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=None");
			run("Macro...", "code=[if(v < "+toString(temp_azimuth_cutoff * max_F)+") v=0]");
			temp_scaled_down_avg_F_image_id = getImageID();
			imageCalculator("Divide create 32-bit", temp_scaled_down_avg_F_image_id, temp_scaled_down_avg_F_image_id);
			rename("scaled_down_azimuth_mask_image");
			temp_scaled_down_azimuth_mask_image_id = getImageID();
			selectImage(temp_scaled_down_avg_F_image_id);
			close();
		}




		// initialize an Fh, Fv composite image with temporal averages of fluorescence intensities
		selectImage(temp_avg_F_image_id);
		run("Duplicate...", "composite image");
		rename("composite image");
		run("Select All");
		run("Copy");	
		run("Add Slice");
		run("Paste");
		temp_composite_image_id = getImageID();
		rename(temp_file_name_base+"_");

		makeLine(1, 1, 2, 2);
		getLine(x1, y1, x2, y2, orig_line_width);
		run("Select None");

	
		// for each pixel over the cutoff intensity, do cos^2 fitting of F(t)
	
		fitting_equation = "y = a*a + b*b * cos(x/180*3.14159 - c*c)*cos(x/180*3.14159 - c*c)";	
		if(show_azimuths == "Yes"){
			if(temp_azimuth_cutoff > temp_LD_cutoff){lower_of_LD_azimuth_cutoffs = temp_LD_cutoff;}
			else{lower_of_LD_azimuth_cutoffs = temp_azimuth_cutoff;}
		}
		else{
			lower_of_LD_azimuth_cutoffs = temp_LD_cutoff;
		}
		
		selectImage(temp_multipol_stack_image_id);
		run("Reslice [/]...", "output=1.000 start=Top avoid");
		resliced_multipol_image_id = getImageID();
		run("Subtract...", "value="+toString(temp_bkgd)+" stack");


		for (row = 1; row <= h; row++){
			showProgress(-row/h);
			showStatus("Fitting row "+toString(row)+" of "+toString(h)+" ("+toString(round(row/h*100))+"%)");		
			for (col = 0; col <= w; col++) {
				//showProgress(row/h);

				selectImage(temp_composite_image_id);
				if(getPixel(col, row-1)/max_F > lower_of_LD_azimuth_cutoffs){
					if(make_new_log2r_azimuth_image == true){ 					
	
			
						selectImage(resliced_multipol_image_id);
						setSlice(row);
						makeLine(col, 0, col, number_of_polarizations-1);
						run("Line Width...", "line=1"); 
						//run("Properties... ", " width=1");
						F_over_time = getProfile();
						Fit.doFit(fitting_equation, polarization_list, F_over_time);
						Fmin = Fit.p(0) * Fit.p(0);
						Fmax = Fmin + Fit.p(1) * Fit.p(1);
	
						log2r = log(Fmax/Fmin)/log(2);
						azimuth = (Fit.p(2) * Fit.p(2))/PI*180 - floor((Fit.p(2) * Fit.p(2))/PI*180/180)*180;
						rsquared = Fit.rSquared;
			
						selectImage(temp_composite_image_id);
						setSlice(1);
						setPixel(col, row-1, Fmin);
						setSlice(2);
						setPixel(col, row-1, Fmax);
	
						selectImage(log2r_azimuth_r_squared_image_id);
						setSlice(1);
						setPixel(col, row-1, log2r);
						setSlice(2);
						setPixel(col, row-1, azimuth);
						setSlice(3);
						setPixel(col, row-1, rsquared);
					}
					else{
						selectImage(log2r_azimuth_r_squared_image_id);
						setSlice(1);
						log2r = getPixel(col, row);

						selectImage(temp_avg_F_image_id);
						F_intensity = getPixel(col, row);

						selectImage(temp_composite_image_id);
						setSlice(1);
						setPixel(col, row, F_intensity / sqrt(pow(2, log2r)));
						setSlice(2);
						setPixel(col, row, F_intensity * sqrt(pow(2, log2r)));
							
					}
				}
			}
		}
		selectImage(resliced_multipol_image_id);
		close();
		run("Line Width...", "line="+orig_line_width); 



		// now I have the Fh, Fv composite, and also the log2r_azimuth_rsquared composite
		// This turns the Fh, Fv composite into an RG image
		rg_image_id = make_rg_image(temp_composite_image_id, temp_bright_adj, temp_rg_range, color_scheme);		
		selectImage(temp_composite_image_id);
		close();

		// Now I need to add the azimuths to the RG image
		if(temp_show_azimuths == "Yes"){
			
			selectImage(log2r_azimuth_r_squared_image_id);
			setSlice(2);
			run("Duplicate...", "title=azimuths");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=None");
			num_of_azimuth_rows = getHeight();
			num_of_azimuth_cols = getWidth();
			//run("Size...", "width="+toString(w)+" height="+toString(h)+" depth=1 constrain interpolation=None");
			temp_azimuth_image_id = getImageID();

			imageCalculator("Multiply", temp_azimuth_image_id, temp_scaled_down_azimuth_mask_image_id);
			
	
			selectImage(log2r_azimuth_r_squared_image_id);
			setSlice(1);
			run("Duplicate...", "title=log2rmax");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=None");
			temp_log2rmax_image_id = getImageID();
	
	
			//selectImage(rg_image_id);
			//Overlay.drawString("Hide azimuths: Image"+rightarrowchar+"Overlay"+rightarrowchar+"Hide Overlay", 10, 30);
		
			
			for (col = 0; col < num_of_azimuth_cols; col++) {
				for (row = 0; row < num_of_azimuth_rows; row++) {
					selectImage(temp_azimuth_image_id);
					azimuth = getPixel(col, row);
					if(azimuth > 0){
						pixel_x_coord = floor(w/num_of_azimuth_cols * col); // + floor(temp_azimuth_bin_size/2) + 1;
						pixel_y_coord = floor(h/num_of_azimuth_rows * row); // + floor(temp_azimuth_bin_size/2) + 1;
	
						selectImage(temp_log2rmax_image_id);
						//relative_log2rmax = getPixel(col, row)/(log(temp_rg_range)/log(2));
						//if(relative_log2rmax > 1){relative_log2rmax = 1;}

						relative_log2rmax = getPixel(col, row)/(log(temp_rg_range)/log(2))/2;
						if(relative_log2rmax > log(temp_rg_range)/log(2) ){relative_log2rmax = 1;}

	
						selectImage(rg_image_id);
						if(temp_azimuth_color_scheme == "direction"){

							//hue = 360 - azimuth * 2;
							hue = 300 - azimuth * 2;
							if(hue < 0){hue = hue + 360;}
							hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
							if(hue >= 0 && hue < 60){
								red = 255;
								green = round(hsv_x * 255);
								blue = 0;
							}
					
							if(hue >= 60 && hue < 120){
								red = round(hsv_x * 255);
								green = 255;
								blue = 0;
							}
					
							if(hue >= 120 && hue < 180){
								red = 0;
								green = 255;
								blue = round(hsv_x * 255);
							}
							if(hue >= 180 && hue < 240){
								red = 0;
								green = round(hsv_x * 255);
								blue = 255;
							}
					
							if(hue >= 240 && hue < 300){
								red = round(hsv_x * 255);
								green = 0;
								blue = 255;
							}
							if(hue >= 300 && hue < 360){
								red = 255;
								green = 0;
								blue = round(hsv_x * 255);
							}
							setColor(red, green, blue);
							setLineWidth(2);
							
						}
						else{
							temp_azimuth_color = temp_azimuth_color_scheme;
							setColor(temp_azimuth_color);  //  overlay color
							setLineWidth(2);
						}
						
						Overlay.drawLine(pixel_x_coord - cos(azimuth/180*PI) * temp_azimuth_bin_size/2 * relative_log2rmax * azimuth_scaling_factor, pixel_y_coord - sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * relative_log2rmax * azimuth_scaling_factor, pixel_x_coord + cos(azimuth/180*PI) * azimuth_pixel_bin_size/2 * relative_log2rmax * azimuth_scaling_factor, pixel_y_coord + sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * relative_log2rmax * azimuth_scaling_factor);			
						Overlay.show;			
						
					}
				}
			}
	
			selectImage(temp_azimuth_image_id);
			close();
			selectImage(temp_log2rmax_image_id);
			close();
			selectImage(temp_scaled_down_azimuth_mask_image_id);
			close();
		}
		//selectImage(log2r_azimuth_r_squared_image_id);
		//close();
		selectImage(temp_avg_F_image_id);
		close();
		
		selectImage(rg_image_id);
			
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////         No LD visualization        //////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	if(temp_show_LD == "None"){
		selectImage(temp_multipol_stack_image_id);
		run("Z Project...", "projection=[Average Intensity]");
		run("Subtract...", "value="+temp_bkgd);
		rename("Avg_F");
		temp_avg_F_image_id = getImageID();
		getStatistics(area, mean_F, min_F, max_F, std, histogram);
		setMinAndMax(0, temp_bright_adj - temp_bkgd);	

		run("Duplicate...", "rg image");
		rename(temp_file_name_base);
		rg_image_id = getImageID();

		if(temp_show_azimuths == "Yes"){
			// make an azimuth mask image
			selectImage(temp_avg_F_image_id);
			run("Duplicate...", "title=scaled_down_avg_F_image");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth=1 constrain interpolation=None");
			run("Macro...", "code=[if(v < "+toString(temp_azimuth_cutoff * max_F)+") v=0]");
			temp_scaled_down_avg_F_image_id = getImageID();
			imageCalculator("Divide create 32-bit", temp_scaled_down_avg_F_image_id, temp_scaled_down_avg_F_image_id);
			rename("scaled_down_azimuth_mask_image");
			temp_scaled_down_azimuth_mask_image_id = getImageID();
			selectImage(temp_scaled_down_avg_F_image_id);
			close();

			// make a shrunk t-stack, apply the shared LD/azimuth mask, reslice the stack			
			selectImage(temp_multipol_stack_image_id);				
			run("Duplicate...", "duplicate title=scaled_down_bkgd_subbed_stack");
			run("Subtract...", "value="+toString(temp_bkgd)+" stack");
			run("Size...", "width="+toString(round(w/temp_azimuth_bin_size))+" height="+toString(round(h/temp_azimuth_bin_size))+" depth="+toString(number_of_polarizations)+" interpolation=None");
			scaled_down_stack_image_id = getImageID();
			imageCalculator("Multiply create 32-bit stack", scaled_down_stack_image_id, temp_scaled_down_azimuth_mask_image_id);
			temp_scaled_down_masked_stack_image_id = getImageID();
			rename("scaled_down_masked_stack");
	
			scaled_down_width = getWidth();
			scaled_down_height = getHeight();
			
			// reslice the shrunk t-stack so that intensity profile of a line selections = F(pol)		
			run("Reslice [/]...", "output=1.000 start=Top avoid");
			resliced_scaled_down_masked_stack_image_id = getImageID();
	
			makeLine(1, 1, 2, 2);
			getLine(x1, y1, x2, y2, orig_line_width);
			run("Select None");

			// create a shrunk log2r_azimuth_r_squared 3-slice stack to hold the fitting results

			if(make_new_log2r_azimuth_image == true){   // make a new log2r_azimuth_r_squared image
				if(log2r_azimuth_r_squared_image_id < 0){
					selectImage(log2r_azimuth_r_squared_image_id);
					close();
				}

				selectImage(temp_scaled_down_azimuth_mask_image_id);
				run("Duplicate...", "title=log2r_azimuth_r_squared_composite");
				run("Add Slice");
				run("Add Slice");
				log2r_azimuth_r_squared_image_id = getImageID();
			}
	
			selectImage(temp_scaled_down_masked_stack_image_id);
	
			// for all pixels in the azimuth mask, do the cos^2 fitting
			// - store the fitting parameters in the shrunk log2r_azimuth_r_squared stack
			// - use the fitting parameters to make the Fh_Fv_multiplier composite
			fitting_equation = "y = a*a + b*b * cos(x/180*3.14159 - c*c)*cos(x/180*3.14159 - c*c)";	
				
			for (row = 0; row < scaled_down_height; row++){
				showStatus("Fitting row "+toString(row)+" of "+toString(scaled_down_height)+" ("+toString(round(row/scaled_down_height*100))+"%)");		
				showProgress(row/scaled_down_height);
				for (col = 0; col < scaled_down_width; col++) {
					selectImage(resliced_scaled_down_masked_stack_image_id);
					setSlice(row + 1);
	
					if(isNaN(getPixel(col, 0))){}
					else{
						makeLine(col, 0, col, number_of_polarizations-1);				
						run("Properties... ", " width=1");
						F_over_time = getProfile();
						Fit.doFit(fitting_equation, polarization_list, F_over_time);
						if(isNaN(Fit.p(0))){}
						else{						
							Fmin = Fit.p(0) * Fit.p(0);
							Fmax = Fmin + Fit.p(1) * Fit.p(1);
		
							log2r = log(Fmax/Fmin)/log(2);
							rel_log2r = log2r/(log(temp_rg_range)/log(2))/2;
							if(rel_log2r > log(temp_rg_range)/log(2) ){rel_log2r = 1;}
	
							
							azimuth = (Fit.p(2) * Fit.p(2))/PI*180 - floor((Fit.p(2) * Fit.p(2))/PI*180/180)*180;
							rsquared = Fit.rSquared;
	
							pixel_x_coord = floor(w/scaled_down_width * col) + floor(temp_azimuth_bin_size/2) + 1;
							pixel_y_coord = floor(h/scaled_down_height * row) + floor(temp_azimuth_bin_size/2) + 1;
	
	
							selectImage(rg_image_id);
							if(temp_azimuth_color_scheme == "direction"){
	
								//hue = 360 - azimuth * 2;
								hue = 300 - azimuth * 2;
								if(hue < 0){hue = hue + 360;}
								hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
								if(hue >= 0 && hue < 60){
									red = 255;
									green = round(hsv_x * 255);
									blue = 0;
								}
						
								if(hue >= 60 && hue < 120){
									red = round(hsv_x * 255);
									green = 255;
									blue = 0;
								}
						
								if(hue >= 120 && hue < 180){
									red = 0;
									green = 255;
									blue = round(hsv_x * 255);
								}
								if(hue >= 180 && hue < 240){
									red = 0;
									green = round(hsv_x * 255);
									blue = 255;
								}
						
								if(hue >= 240 && hue < 300){
									red = round(hsv_x * 255);
									green = 0;
									blue = 255;
								}
								if(hue >= 300 && hue < 360){
									red = 255;
									green = 0;
									blue = round(hsv_x * 255);
								}
								setColor(red, green, blue);
								setLineWidth(2);
								
							}
							else{
								temp_azimuth_color = temp_azimuth_color_scheme;
								setColor(temp_azimuth_color);  //  overlay color
								setLineWidth(2);
							}
			
							// azimuth length = F intensity
							// Overlay.drawLine(pixel_x_coord - cos(azimuth/180*PI) * temp_azimuth_bin_size/2 * rel_F, pixel_y_coord - sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F, pixel_x_coord + cos(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F, pixel_y_coord + sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_F);			
			
							// azimuth length = LD
							Overlay.drawLine(pixel_x_coord - cos(azimuth/180*PI) * temp_azimuth_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_y_coord - sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_x_coord + cos(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor, pixel_y_coord + sin(azimuth/180*PI) * azimuth_pixel_bin_size/2 * rel_log2r * azimuth_scaling_factor);			
							Overlay.show;		

						}							
					}
				}
			}
			selectImage(resliced_scaled_down_masked_stack_image_id);
			close();	
		
	
	
		//	selectImage(temp_azimuth_image_id);
		//	close();
		
			selectImage(temp_avg_F_image_id);
			close();
			selectImage(scaled_down_stack_image_id);
			close();
			selectImage(temp_scaled_down_azimuth_mask_image_id);
			close();
			selectImage(temp_scaled_down_masked_stack_image_id);
			close();
		//	selectImage(log2r_azimuth_r_squared_image_id);
		//	close();

		}	
	}

	return_array = newArray(rg_image_id, log2r_azimuth_r_squared_image_id);

	return return_array;

}





function trace_a_filament(t_segmentation_image_id){
	selectImage(t_segmentation_image_id);
	getSelectionCoordinates(user_selected_points_x, user_selected_points_y);	// x, y coordinates of the spline

	makeSelection("polyline", user_selected_points_x, user_selected_points_y);
	setColor(255,255,255);	
	Overlay.addSelection;
	Overlay.show;
	run("Select None");

	traced_points_x = newArray("Null");
	traced_points_y = newArray("Null");


	for(segment_number = 0; segment_number < user_selected_points_x.length - 1; segment_number ++){
		traced_points_x = Array.deleteIndex(traced_points_x, traced_points_x.length-1);
		traced_points_y = Array.deleteIndex(traced_points_y, traced_points_y.length-1);

		
		segment_length = sqrt( (user_selected_points_x[segment_number+1] - user_selected_points_x[segment_number])*(user_selected_points_x[segment_number+1] - user_selected_points_x[segment_number]) +  (user_selected_points_y[segment_number+1] - user_selected_points_y[segment_number])*(user_selected_points_y[segment_number+1] - user_selected_points_y[segment_number]));
		segment_direction_x = (user_selected_points_x[segment_number+1] - user_selected_points_x[segment_number])/segment_length;
		segment_direction_y = (user_selected_points_y[segment_number+1] - user_selected_points_y[segment_number])/segment_length;
	
		tangent_direction_x = segment_direction_y;
		tangent_direction_y = -segment_direction_x;
		
	
		// let me start with a 2-point selection first
		// split the selection into 10 segments
	
	
		number_of_subsegments = 3;
		subsegment_points_x = newArray(number_of_subsegments + 1);
		subsegment_points_y = newArray(number_of_subsegments + 1);
	
		for(i = 0; i <= number_of_subsegments; i++){
			subsegment_points_x[i] = user_selected_points_x[segment_number] + i/number_of_subsegments*(user_selected_points_x[segment_number+1] - user_selected_points_x[segment_number]);
			subsegment_points_y[i] = user_selected_points_y[segment_number] + i/number_of_subsegments*(user_selected_points_y[segment_number+1] - user_selected_points_y[segment_number]);		
		}
	
		Array.print(user_selected_points_x);
		Array.print(user_selected_points_y);
		print(" ");
		Array.print(subsegment_points_x);
		Array.print(subsegment_points_y);
	
		run("Select None");
		makeSelection("point",subsegment_points_x,subsegment_points_y);
	
	
		for(i = 0; i <= number_of_subsegments; i++){
			processed_point_x = subsegment_points_x[i];
			processed_point_y = subsegment_points_y[i];

			run("Line Width...", "line=5"); 
	
			makeLine(processed_point_x - tangent_direction_x * segment_length/10, processed_point_y - tangent_direction_y * segment_length/10, processed_point_x + tangent_direction_x * segment_length/10, processed_point_y + tangent_direction_y * segment_length/10);	
	
			setColor(255,0,0);
			Overlay.addSelection;
			Overlay.show;
	
			profile = getProfile();
			Array.print(profile);
			maximum_intensity_positions = Array.findMaxima(profile, 1);
			if(maximum_intensity_positions.length > 0){
				print("maximum position = "+maximum_intensity_positions[0]+"/"+profile.length);
				maximum_intensity_coordinate_x = (processed_point_x - tangent_direction_x * segment_length/10) + maximum_intensity_positions[0]/profile.length * 2 * tangent_direction_x * segment_length/10;
				maximum_intensity_coordinate_y = (processed_point_y - tangent_direction_y * segment_length/10) + maximum_intensity_positions[0]/profile.length * 2 * tangent_direction_y * segment_length/10;
		
				setColor(0,255,0);
				makePoint(maximum_intensity_coordinate_x, maximum_intensity_coordinate_y);
				Overlay.addSelection;
				Overlay.show;

				Array.getStatistics(profile, min, max, mean, stdDev);


				traced_points_x = Array.concat(traced_points_x, maximum_intensity_coordinate_x);
				traced_points_y = Array.concat(traced_points_y, maximum_intensity_coordinate_y);


			
			}
			else{
				setColor(0,255,0);
				makePoint(processed_point_x, processed_point_y);
				Overlay.addSelection;
				Overlay.show;
				traced_points_x = Array.concat(traced_points_x, processed_point_x);
				traced_points_y = Array.concat(traced_points_y, processed_point_y);
	
			}
		}
	}

	selectImage(t_segmentation_image_id);
	Overlay.hide;
	run("Select None");
	Array.print(traced_points_x);
	Array.print(traced_points_y);
	makeSelection("polyline", traced_points_x, traced_points_y);
	run("Fit Spline");
}





function make_HSV_snake_overlay(t_image_id, t_snake_points_x, t_snake_points_y){

	selectImage(t_image_id);

	for (i = 0; i < t_snake_points_x.length; i++) {
		hue = i * 360/t_snake_points_x.length;
		hsv_x = 1 - abs(2 * ((hue/60 / 2) - floor(hue/60 /2)) - 1);
		if(hue >= 0 && hue < 60){
			red = 255;
			green = round(hsv_x * 255);
			blue = 0;
		}

		if(hue >= 60 && hue < 120){
			red = round(hsv_x * 255);
			green = 255;
			blue = 0;
		}

		if(hue >= 120 && hue < 180){
			red = 0;
			green = 255;
			blue = round(hsv_x * 255);
		}
		if(hue >= 180 && hue < 240){
			red = 0;
			green = round(hsv_x * 255);
			blue = 255;
		}

		if(hue >= 240 && hue < 300){
			red = round(hsv_x * 255);
			green = 0;
			blue = 255;
		}
		if(hue >= 300 && hue < 360){
			red = 255;
			green = 0;
			blue = round(hsv_x * 255);
		}
		setColor(red, green, blue);
		setLineWidth(1);
		if(i < t_snake_points_x.length - 1){
			Overlay.drawLine(t_snake_points_x[i], t_snake_points_y[i], t_snake_points_x[i+1], t_snake_points_y[i+1]);		
			Overlay.show;
		}
	}
}
