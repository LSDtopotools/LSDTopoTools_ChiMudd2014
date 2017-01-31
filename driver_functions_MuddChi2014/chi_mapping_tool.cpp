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
#include "../LSDParameterParser.hpp"
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

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);
  
  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // set default float parameters
  int_default_map["n_iterations"] = 20;
  int_default_map["minimum_segment_length"] = 10;
  int_default_map["n_nodes_to_visit"] = 10;
  int_default_map["target_nodes"] = 80;
  int_default_map["skip"] = 2;
  int_default_map["threshold_pixels_for_chi"] = 1000;
  int_default_map["minimum_basin_size_pixels"] = 1000;
  int_default_map["threshold_contributing_pixels"] = 1000;
  
  // set default in parameter
  float_default_map["min_slope_for_fill"] = 0.0001;
  float_default_map["A_0"] = 1;
  float_default_map["m_over_n"] = 0.5;
  float_default_map["sigma"] = 20;
  
  // set default methods
  bool_default_map["only_check_parameters"] = false;
  bool_default_map["print_stream_order_raster"] = false;
  bool_default_map["print_junction_index_raster"] = false;
  bool_default_map["print_fill_raster"] = false;
  bool_default_map["print_DrainageArea_raster"] = false;
  bool_default_map["print_chi_coordinate_raster"] = false;
  bool_default_map["print_simple_chi_map_to_csv"] = false;
  bool_default_map["print_segmented_M_chi_map_to_csv"] = false;
  bool_default_map["print_basic_M_chi_map_to_csv"] = false;  
  bool_default_map["test_drainage_boundaries"] = false;  
  bool_default_map["only_take_largest_basin"] = false;  
  bool_default_map["print_source_keys"] = false;
  bool_default_map["print_baselevel_keys"] = false;
  bool_default_map["print_basin_raster"] = false;
  bool_default_map["write hillshade"] = false;
  bool_default_map["print_simple_chi_map_with_basins_to_csv"] = false;
  
  // set default string method
  string_default_map["CHeads_file"] = "NULL";
  string_default_map["averaging_raster_vector"] = "NULL";
  
  // Use the parameter parser to get the maps of the parameters required for the 
  // analysis

  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();
  
  cout << "Read filename is:" <<  DATA_DIR+DEM_ID << endl;
  
    // check to see if the raster exists
  LSDRasterInfo RI((DATA_DIR+DEM_ID), raster_ext);  

  // initialise variables to be assigned from .driver file
  // These will all be assigned default values
  float A_0 = this_float_map["A_0"];
  float movern = this_float_map["m_over_n"];
  float Minimum_Slope = this_float_map["min_slope_for_fill"];
  int n_iterations = this_int_map["n_iterations"];
  int minimum_segment_length = this_int_map["minimum_segment_length"];
  int n_nodes_to_visit = this_int_map["n_nodes_to_visit"];             // when constructing channel network, this 
  float sigma = this_float_map["sigma"];
  int target_nodes = this_int_map["target_nodes"];
  int skip = this_int_map["skip"];
  int threshold_contributing_pixels = this_int_map["threshold_contributing_pixels"];
  int minimum_basin_size_pixels = this_int_map["minimum_basin_size_pixels"];
  int basic_Mchi_regression_nodes = this_int_map["basic_Mchi_regression_nodes"];
  bool test_drainage_boundaries = this_bool_map["test_drainage_boundaries"];

  // load the  DEM
  LSDRaster topography_raster((DATA_DIR+DEM_ID), raster_ext);
  cout << "Got the dem: " <<  DATA_DIR+DEM_ID << endl;
  
  float Resolution = topography_raster.get_DataResolution();
  map<string,string> GRS = topography_raster.get_GeoReferencingStrings();

  float thresh_area_for_chi = float(this_int_map["threshold_pixels_for_chi"])*Resolution*Resolution;

  if(this_bool_map["only_check_parameters"])
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
  
  if (this_bool_map["print_fill_raster"])
  {
    string filled_raster_name = OUT_DIR+DEM_ID+"_Fill";
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

  if (this_bool_map["print_DrainageArea_raster"])
  {
    string DA_raster_name = OUT_DIR+DEM_ID+"_DArea";
    DrainageArea.write_raster(DA_raster_name,raster_ext);
  }

  // calcualte the distance from outlet
  cout << "\t Calculating flow distance..." << endl;
  LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

  cout << "\t Loading Sources..." << endl;
  // load the sources
  vector<int> sources;
  if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
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
  
  if (this_bool_map["print_stream_order_raster"])
  { 
    LSDIndexRaster SOArray = JunctionNetwork.StreamOrderArray_to_LSDIndexRaster();
    string SO_raster_name = OUT_DIR+DEM_ID+"_SO";
    SOArray.write_raster(SO_raster_name,raster_ext);
  }
  if (this_bool_map["print_junction_index_raster"])
  { 
    LSDIndexRaster JIArray = JunctionNetwork.JunctionIndexArray_to_LSDIndexRaster();
    string JI_raster_name = OUT_DIR+DEM_ID+"_JI";
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
  if (this_bool_map["test_drainage_boundaries"])
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
  cout << "Removing basins with fewer than " << minimum_basin_size_pixels << " pixels" << endl;
  BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Area(BaseLevelJunctions, 
                                              FlowInfo, FlowAcc, minimum_basin_size_pixels);

  if (this_bool_map["only_take_largest_basin"])
  {
    cout << "I am only going to take the largest basin." << endl;
    BaseLevelJunctions = JunctionNetwork.Prune_Junctions_Largest(BaseLevelJunctions, FlowInfo, FlowAcc);
  }

  // Correct number of base level junctions
  N_BaseLevelJuncs = BaseLevelJunctions.size();
  cout << "The number of basins I will analyse is: " << N_BaseLevelJuncs << endl;
  
  // This is for debugging
  //for (int BN = 0; BN< N_BaseLevelJuncs; BN++)
  //{
  //  cout << "BL junc is: " << BaseLevelJunctions[BN] << " node is: " << JunctionNetwork.get_Node_of_Junction(BaseLevelJunctions[BN]) << endl;
  //  vector<int> UPSN = FlowInfo.get_upslope_nodes(JunctionNetwork.get_Node_of_Junction(BaseLevelJunctions[BN]));
  //  cout << "Pixels for that node are: " << UPSN.size() << endl;
  //}

  // now use a ChiTool object to print the chi tree to csv                  
  LSDChiTools ChiTool(FlowInfo);

  //============================================================================
  // Print a basin raster if you want it.
  if(this_bool_map["print_basin_raster"])
  {
    string basin_raster_prefix = OUT_DIR+DEM_ID;
    ChiTool.print_basins(FlowInfo, JunctionNetwork, BaseLevelJunctions, basin_raster_prefix);
  }

  // calculate chi for the entire DEM
  cout << "Calculating the chi coordinate for A_0: " << A_0 << " and m/n: " << movern << endl; 
  LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

  if(this_bool_map["print_chi_coordinate_raster"])
  {                                                                         
    string chi_coord_string = OUT_DIR+DEM_ID+"_chi_coord"; 
    chi_coordinate.write_raster(chi_coord_string,raster_ext);
  }                                                                         
  
  if (this_bool_map["print_simple_chi_map_to_csv"])
  {
    cout <<"I am printing a simple chi map for you to csv." << endl;
    string chi_csv_fname = OUT_DIR+DEM_ID+"_chi_coord.csv";
    ChiTool.chi_map_to_csv(FlowInfo, chi_csv_fname, chi_coordinate);
  }

  if (this_bool_map["print_simple_chi_map_with_basins_to_csv"])
  {
    cout <<"I am printing a simple chi map with basins for you to csv." << endl;
    LSDIndexRaster basin_raster = ChiTool.get_basin_raster(FlowInfo, JunctionNetwork, BaseLevelJunctions);
    string chi_csv_fname = OUT_DIR+DEM_ID+"_chi_coord_basins.csv";
    ChiTool.chi_map_to_csv(FlowInfo, chi_csv_fname, chi_coordinate,basin_raster);
  }


  // now calculate segmentation
  vector<int> source_nodes;
  vector<int> outlet_nodes;
  if (this_bool_map["print_segmented_M_chi_map_to_csv"] || this_bool_map["print_basic_M_chi_map_to_csv"])
  // now get the overlapping channels from the junction network file
  {
    cout << "I am getting the source and outlet nodes for the overlapping channels" << endl;
    cout << "The n_nodes to visit are: " << n_nodes_to_visit << endl;
    JunctionNetwork.get_overlapping_channels(FlowInfo, BaseLevelJunctions, DistanceFromOutlet, 
                                    source_nodes,outlet_nodes,n_nodes_to_visit);
  }

  if (this_bool_map["print_segmented_M_chi_map_to_csv"])
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
  
    string csv_full_fname = OUT_DIR+DEM_ID+"_MChiSegmented.csv";
    cout << "Let me print all the data for you into a csv file called " << csv_full_fname << endl;
    ChiTool.print_data_maps_to_file_full(FlowInfo, csv_full_fname);
    cout << "That is your file printed!" << endl;
    
    // These print the source and baselelvel keys if wanted
    if (this_bool_map["print_source_keys"])
    {
      string sources_keys_name = OUT_DIR+DEM_ID+"_SourceKeys.csv";
      ChiTool.print_source_keys(FlowInfo, sources_keys_name);
    }
    if (this_bool_map["print_baselevel_keys"])
    {
      string baselevel_keys_name = OUT_DIR+DEM_ID+"_BaselevelKeys.csv";
      ChiTool.print_baselevel_keys(FlowInfo, JunctionNetwork, baselevel_keys_name);
    }
    
  }

  if (this_bool_map["print_basic_M_chi_map_to_csv"])
  {
    
    LSDChiTools ChiTool2(FlowInfo);
    ChiTool2.chi_map_automator_rudimentary(FlowInfo, source_nodes,outlet_nodes,
                                    filled_topography, DistanceFromOutlet, DrainageArea, 
                                    chi_coordinate, basic_Mchi_regression_nodes);
    string csv_full_fname = OUT_DIR+DEM_ID+"_MChiBasic.csv";
    ChiTool2.print_data_maps_to_file_full(FlowInfo, csv_full_fname);
  }

  if (this_bool_map["write hillshade"])
  {
    float hs_azimuth = 315;
    float hs_altitude = 45;
    float hs_z_factor = 1;
    LSDRaster hs_raster = topography_raster.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

    string hs_fname = OUT_DIR+DEM_ID+"_hs";
    hs_raster.write_raster(hs_fname,raster_ext);
  }

}
