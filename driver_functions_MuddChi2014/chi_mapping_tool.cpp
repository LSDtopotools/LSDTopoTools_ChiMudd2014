//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chi_mapping_tool
//
// This program takes two arguments, the path name and the driver name
// The driver file has a number of options that allow the user to calculate 
// different kinds of chi analysis
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Copyright (C) 2016 Simon M. Mudd 2016
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../LSDStatsTools.hpp"
#include "../LSDChiNetwork.hpp"
#include "../LSDRaster.hpp"
#include "../LSDRasterInfo.hpp"
#include "../LSDIndexRaster.hpp"
#include "../LSDFlowInfo.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDIndexChannelTree.hpp"
#include "../LSDBasin.hpp"
#include "../LSDChiTools.hpp"
#include "../LSDShapeTools.hpp"

int main (int nNumberofArgs,char *argv[])
{
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the chi mapping tool!                    ||" << endl;
    cout << "|| This program has a number of options to make chi    ||" << endl;
    cout << "|| plots and to map out slopes in chi space.           ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd         ||" << endl;
    cout << "||  at the University of Edinburgh                     ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./chi_mapping_tool.exe home/fieldwork/Chile/CRN/ My_analysis.param" << endl;
    cout << "In windows (the slash directions will change and there is no leading ./)" << endl;
    cout << "chi_mapping_tool.exe c:\\fieldwork\\Chile\\CRN\\ My_analysis.param" << endl;
    cout << "=========================================================" << endl;
    cout << "For more documentation on the parameter file, " << endl;
    cout << " see readme and online documentation." << endl;
    cout << " http://lsdtopotools.github.io/LSDTT_book/#_chi_analysis_part_3_getting_chi_gradients_for_the_entire_landscape" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  cout << "The path is: " << path_name << " and the filename is: " << f_name << endl;

  string full_name = path_name+f_name;
    
  ifstream file_info_in;
  file_info_in.open(full_name.c_str());
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: The parameter \"" << full_name
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  // strings for filenames. Some set to defaults
  string DATA_DIR,OUTPUT_DIR,DEM_ID;
  string CHeads_file = "Null";
  string raster_ext = "bil";

  // initialise variables to be assigned from .driver file
  // These will all be assigned default values
  float A_0 = 1000;
  float movern = 0.5;
  float Minimum_Slope = 0.001;
  int n_iterations = 20;
  int minimum_segment_length = 10;
  int n_nodes_to_visit = 10;             // when constructing channel network, this 
  float sigma = 20;
  int target_nodes = 80;
  int skip = 2;
  int threshold_pixels_for_chi = 1000;
  int threshold_contributing_pixels = 1000;
  int minimum_basin_size_pixels = 1000;
  int basic_Mchi_regression_nodes = 11;
  bool test_drainage_boundaries = true;
  bool only_take_largest_basin = true;

  // bools for what the program will print
  bool only_check_parameters = false;
  bool print_stream_order_raster = false; 
  bool print_junction_index_raster = false; 
  bool print_fill_raster = false; 
  bool print_DrainageArea_raster = false; 
  bool print_chi_coordinate_raster = false;
  bool print_simple_chi_map_to_csv = false; 
  bool print_segmented_M_chi_map_to_csv = false;
  bool print_basic_M_chi_map_to_csv = false;  

  //============================================================================
  // Ingest parameters
  //============================================================================
  string parameter, value, lower, lower_val;
  string bc;
  while (file_info_in.good())
  {
    parse_line(file_info_in, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = tolower(parameter[i]);
    }

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);

    if (lower == "read path")
    {
      DATA_DIR = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      DATA_DIR = RemoveControlCharactersFromEndOfString(DATA_DIR);
    }
    else if (lower == "write path")
    {
      OUTPUT_DIR = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      OUTPUT_DIR = RemoveControlCharactersFromEndOfString(OUTPUT_DIR);
    }
    else if (lower == "read fname")
    {
      DEM_ID = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      DEM_ID = RemoveControlCharactersFromEndOfString(DEM_ID);
      //cout << "Got the read name, it is: " << read_fname << endl;
    }
    else if (lower == "channel heads fname")
    {
      CHeads_file = value;
      // get rid of any control characters from the end (if param file was made in DOS)
      CHeads_file = RemoveControlCharactersFromEndOfString(CHeads_file);
      //cout << "Got the channel heads name, it is: " << CHeads_file << endl;
    }
    
    //============================================================
    // Parameters for fill
    //============================================================
    else if (lower == "min_slope_for_fill")
    {
      Minimum_Slope = atof(value.c_str());
    }
    
    //============================================================
    // Parameters for chi
    //============================================================
    else if (lower == "a_0")
    {
      A_0 = atof(value.c_str());
    }
    else if (lower == "m_over_n")
    {
      movern = atof(value.c_str());
    }
    else if (lower == "threshold_pixels_for_chi")
    {
      threshold_pixels_for_chi = atoi(value.c_str());
    }
    
    //============================================================
    // Parameters for segment fitting
    //============================================================
    else if (lower == "n_iterations")
    {
      n_iterations = atoi(value.c_str());
    }
    else if (lower == "minimum_segment_length")
    {
      minimum_segment_length = atoi(value.c_str());
    }
    else if (lower == "skip")
    {
      skip = atoi(value.c_str());
    }
    else if (lower == "sigma")
    {
      sigma = atof(value.c_str());
    }
    else if (lower == "target_nodes")
    {
      target_nodes = atoi(value.c_str());
    }
    
    //============================================================
    // Parameters for selecting basins
    //============================================================
    else if (lower == "threshold_contributing_pixels")
    {
      threshold_contributing_pixels = atoi(value.c_str());
    }
    else if (lower == "minimum_basin_size_pixels")
    {
      minimum_basin_size_pixels = atoi(value.c_str());
    }
    else if (lower == "test_drainage_boundaries")
    {
      test_drainage_boundaries = atobool(value.c_str());
    }
    else if (lower == "only_take_largest_basin")
    {
      only_take_largest_basin = atobool(value.c_str());
    }
    else if (lower == "n_nodes_to_visit")
    {
      n_nodes_to_visit = atoi(value.c_str());
    }


    //============================================================
    // Parameters for deciding what you want to print
    //============================================================
    else if (lower == "only_check_parameters")
    {
      only_check_parameters = atobool(value.c_str());
    }
    else if (lower == "print_stream_order_raster")
    {
      print_stream_order_raster = atobool(value.c_str());
    }
    else if (lower == "print_junction_index_raster")
    {
      print_junction_index_raster = atobool(value.c_str());
    }
    else if (lower == "print_fill_raster")
    {
      print_fill_raster = atobool(value.c_str());
    }
    else if (lower == "print_drainagearea_raster")
    {
      print_DrainageArea_raster = atobool(value.c_str());
    }
    else if (lower == "print_chi_coordinate_raster")
    {
      print_chi_coordinate_raster = atobool(value.c_str());
    }
    else if (lower == "print_simple_chi_map_to_csv")
    {
      print_simple_chi_map_to_csv = atobool(value.c_str());
    }
    else if (lower == "print_segmented_m_chi_map_to_csv")
    {
      print_segmented_M_chi_map_to_csv = atobool(value.c_str());
    }
    else if (lower == "print_basic_m_chi_map_to_csv")
    {
      print_basic_M_chi_map_to_csv = atobool(value.c_str());
    }
  }


  file_info_in.close();

  
  // set no flux boundary conditions
  cout << "Setting boundary conditions" << endl;
  vector<string> boundary_conditions(4);
  boundary_conditions[0] = "No";
  boundary_conditions[1] = "no flux";
  boundary_conditions[2] = "no flux";
  boundary_conditions[3] = "No flux";

  // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);  
        
  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;
  
  float Resolution = topography_raster.get_DataResolution();
  map<string,string> GRS = topography_raster.get_GeoReferencingStrings();

  float thresh_area_for_chi = float(threshold_pixels_for_chi)*Resolution*Resolution;
  float thresh_area_for_channel = float(threshold_contributing_pixels)*Resolution*Resolution;
  float min_basin_area =  float(minimum_basin_size_pixels)*Resolution*Resolution;

  //============================================================================
  // Print the parameters to an input file that you can read later
  //============================================================================
  string parameter_values_name = OUTPUT_DIR+DEM_ID+"_Input.param";
  ofstream param_input;
  param_input.open(parameter_values_name.c_str());
  
  param_input << "# You have been using the ChiToolsDriver.exe program." << endl
              << "# This was developed by the Land Surface Dyanmics Group" << endl
              << "#  at the University of Edinburgh. " << endl
              << "# Algorithms used to create this data can be found in" << endl
              << "# Mudd et al., 2014, JGR-ES: http://onlinelibrary.wiley.com/doi/10.1002/2013JF002981/full" << endl
              << endl
              << "# ================================================" << endl
              << "# The file information is: "<< endl
              << "read path: " << DATA_DIR<< endl
              << "write path: " << OUTPUT_DIR<< endl
              << "read fname: " << DEM_ID << endl 
              << "channel heads fname: " << CHeads_file<< endl
              << "raster_ext: " << raster_ext<< endl
              << "#Note that if the CHeads file is Null then a pixel threshold will be used to determine channel sources." << endl
              << endl
              << "# ================================================" << endl
              << "# The parameters for the analysis are: " << endl
              << "# Note: A_0 in m^2, sigma in m"
              << "min_slope_for_fill: " << Minimum_Slope << endl
              << "A_0: " <<  A_0 << endl
              << "m_over_n: " <<  movern << endl
              << "#And these are parameters for segmentation. Only used if you have selected the segmentation algorithm. " << endl
              << "n_iterations: " << n_iterations << endl
              << "minimum_segment_length: " <<  minimum_segment_length << endl
              << "target nodes: " <<  target_nodes << endl
              << "sigma: " <<  sigma << endl
              << "skip: " <<  skip << endl
              << "# This parameter is for the basic M_chi calculations that use a fixed window length." << endl
              << "basic_Mchi_regression_nodes: " << basic_Mchi_regression_nodes << endl
              << endl
              << "# ================================================" << endl
              << "# These are parameters for controlling thresholds for chi and basins:" << endl 
              << "# Areas are in m^2"
              << "threshold_pixels_for_chi: " << threshold_pixels_for_chi << endl
              << "thresh_area_for_chi: " << thresh_area_for_chi << endl
              << "threshold_contributing_pixels: " << threshold_contributing_pixels << endl
              << "thresh_area_for_channel: " << thresh_area_for_channel << endl
              << "test_drainage_boundaries: " << test_drainage_boundaries << endl 
              << "minimum_basin_size_pixels: " << minimum_basin_size_pixels << endl
              << "min_basin_area: " << min_basin_area << endl
              << "n_nodes_to_visit: " << n_nodes_to_visit  << endl
              << endl
              << "# ================================================" << endl
              << "#  These are all the things that the program will print: " << endl
              << "#  Note: " << endl
              << "#    1 == true" << endl
              << "#    0 == false" << endl
              << "#  I am only going to check parameters and not actually print anything: " << only_check_parameters << endl
              << "print_fill_raster: " << print_fill_raster << endl
              << "print_DrainageArea_raster: " << print_DrainageArea_raster  << endl
              << "print_stream_order_raster: " << print_stream_order_raster  << endl
              << "print_junction_index_raster: " << print_junction_index_raster << endl
              << "print_chi_coordinate_raster: " << print_chi_coordinate_raster << endl
              << "print_simple_chi_map_to_csv: " << print_simple_chi_map_to_csv << endl
              << "print_segmented_M_chi_map_to_csv: " << print_segmented_M_chi_map_to_csv << endl
              << "print_basic_M_chi_map_to_csv: " << print_basic_M_chi_map_to_csv << endl
              << endl << endl;

  param_input.close();

  if(only_check_parameters)
  {
    cout << "You set the only_check_parameters flag to true; I have only printed" << endl;
    cout << "the parameters to file and am now exiting." << endl;
    exit(0);
  }

  //============================================================================
  // Start gathering necessary rasters
  //============================================================================
  cout << "Filling topography." << endl;
  LSDRaster filled_topography = topography_raster.fill(Minimum_Slope);
  
  if (print_fill_raster)
  {
    string filled_raster_name = OUTPUT_DIR+DEM_ID+"_Fill";
    filled_topography.write_raster(filled_raster_name,raster_ext);
  }
  
  
  cout << "\t Flow routing..." << endl;
  // get a flow info object
  LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

  // calculate the flow accumulation
  cout << "\t Calculating flow accumulation (in pixels)..." << endl;
  LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
  
  cout << "\t Converting to flow area..." << endl;
  LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

  if (print_DrainageArea_raster)
  {
    string DA_raster_name = OUTPUT_DIR+DEM_ID+"_DArea";
    DrainageArea.write_raster(DA_raster_name,raster_ext);
  }

  // calcualte the distance from outlet
  cout << "\t Calculating flow distance..." << endl;
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file != "null")
  {
    cout << endl << endl << endl << "==================================" << endl;
    cout << "The channel head file is null. " << endl;
    cout << "Getting sources from a threshold of "<< threshold_contributing_pixels << " pixels." <<endl;
    sources = FlowInfo.get_sources_index_threshold(FlowAcc, threshold_contributing_pixels);
    
    cout << "The number of sources is: " << sources.size() << endl;
    
  }
  else
  {
    cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
    sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
    cout << "\t Got sources!" << endl;
  }

  // now get the junction network
  LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);
  
  if (print_stream_order_raster)
  { 
    LSDIndexRaster SOArray = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SO_raster_name = OUTPUT_DIR+DEM_ID+"_SO";
    SOArray.write_raster(SO_raster_name,raster_ext);
  }
  if (print_junction_index_raster)
  { 
    LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
    string JI_raster_name = OUTPUT_DIR+DEM_ID+"_JI";
    JIArray.write_raster(JI_raster_name,raster_ext);
  }

  //Array2D<float> ChiGradientArray(NRows,NCols,NoData);

  // need to get base-level nodes , otherwise these catchments will be missed!
  vector< int > BaseLevelJunctions_Initial = JunctionNetwork.get_BaseLevelJunctions();
  vector< int > BaseLevelJunctions;
  int N_BaseLevelJuncs = BaseLevelJunctions_Initial.size();
  cout << "The number of base level junctions is: " << N_BaseLevelJuncs << endl; 

  // remove basins drainage from edge if that is what the user wants
  
  cout << "Right, let me check the drainage basins. " << endl;
  if (test_drainage_boundaries)
  {
    cout << "Test_drainage_boundaries: " << test_drainage_boundaries << endl;
  
    cout << endl << endl << "I am going to remove any basins draining to the edge." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Edge(BaseLevelJunctions_Initial,FlowInfo); 
  }
  else
  {
    cout << "I'm not going to remove any drainage basins drainage to the edge." << endl;
    BaseLevelJunctions = BaseLevelJunctions_Initial;
  }
  
  // remove base level junctions for which catchment is too small
  BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions, 
                                              FlowInfo, FlowAcc, minimum_basin_size_pixels);

  if (only_take_largest_basin)
  {
    cout << "I am only going to take the largest basin." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Largest(BaseLevelJunctions, FlowInfo, FlowAcc);
  }

  // Correct number of base level junctions
  N_BaseLevelJuncs = BaseLevelJunctions.size();
  cout << "The number of basins I will analyse is: " << N_BaseLevelJuncs << endl;

  // calculate chi for the entire DEM
  cout << "Calculating the chi coordinate for A_0: " << A_0 << " and m/n: " << movern << endl; 
  LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);
                                                                            
  if(print_chi_coordinate_raster)                                           
  {                                                                         
    string chi_coord_string = OUTPUT_DIR+DEM_ID+"_chi_coord";               
    chi_coordinate.write_raster(chi_coord_string,raster_ext);               
  }                                                                         
                                                                            
  // now use a ChiTool object to print the chi tree to csv                  
  LSDChiTools ChiTool(FlowInfo);
  
  if (print_simple_chi_map_to_csv)
  {
    cout <<"I am printing a simple chi map for you to csv." << endl;
    string chi_csv_fname = OUTPUT_DIR+DEM_ID+"_chi_coord.csv";
    ChiTool.chi_map_to_csv(FlowInfo, chi_csv_fname, chi_coordinate);
  }


  // now calculate segmentation
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  if (print_segmented_M_chi_map_to_csv || print_basic_M_chi_map_to_csv)
  // now get the overlapping channels from the junction network file
  {
    cout << "I am getting the source and outlet nodes for the overlapping channels" << endl;
    cout << "The n_nodes to visit are: " << n_nodes_to_visit << endl;
    JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, DistanceFromOutlet, 
                                    source_nodes,outlet_nodes,n_nodes_to_visit);
  }

  if (print_segmented_M_chi_map_to_csv)
  {
    cout << "I am calculating the segmented channels" << endl;
    if (source_nodes.size() == 0)
    {
      cout << "I don't seem to have any source nodes!" << endl;
    }
    ChiTool.chi_map_automator(FlowInfo, source_nodes, outlet_nodes,
                            filled_topography, DistanceFromOutlet, 
                            DrainageArea, chi_coordinate, target_nodes, 
                            n_iterations, skip, minimum_segment_length, sigma);
  
    string csv_full_fname = OUTPUT_DIR+DEM_ID+"_MChiSegmented.csv";
    ChiTool.print_data_maps_to_file_full(FlowInfo, csv_full_fname);
  }

  if (print_basic_M_chi_map_to_csv)
  {
    LSDChiTools ChiTool2(FlowInfo);
    ChiTool2.chi_map_automator_rudimentary(FlowInfo, source_nodes,outlet_nodes,
                                    filled_topography, DistanceFromOutlet, DrainageArea, 
                                    chi_coordinate, basic_Mchi_regression_nodes);
    string csv_full_fname = OUTPUT_DIR+DEM_ID+"_MChiBasic.csv";
    ChiTool2.print_data_maps_to_file_full(FlowInfo, csv_full_fname);
  }

}
