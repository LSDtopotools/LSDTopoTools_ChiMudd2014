//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDChiTools
// Land Surface Dynamics ChiTools object
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for performing various analyses in chi space
//
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
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
// either version 2 of the License, or (at your option) any later version.
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
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// LSDChiTools.cpp
// LSDChiTools object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDChiTools_CPP
#define LSDChiTools_CPP

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDChiTools.hpp"
#include "LSDChiNetwork.hpp"
using namespace std;
using namespace TNT;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDRaster 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDIndexRaster& ThisRaster)
{
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDFlowInfo& ThisFI)
{
  NRows = ThisFI.get_NRows();
  NCols = ThisFI.get_NCols();
  XMinimum = ThisFI.get_XMinimum();
  YMinimum = ThisFI.get_YMinimum();
  DataResolution = ThisFI.get_DataResolution();
  NoDataValue = ThisFI.get_NoDataValue();
  GeoReferencingStrings = ThisFI.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Creates an LSDChiTools from an LSDFlowInfo
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::create(LSDJunctionNetwork& ThisJN)
{
  NRows = ThisJN.get_NRows();
  NCols = ThisJN.get_NCols();
  XMinimum = ThisJN.get_XMinimum();
  YMinimum = ThisJN.get_YMinimum();
  DataResolution = ThisJN.get_DataResolution();
  NoDataValue = ThisJN.get_NoDataValue();
  GeoReferencingStrings = ThisJN.get_GeoReferencingStrings();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This resets all the data maps
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::reset_data_maps()
{
  map<int,float> empty_map;
  vector<int> empty_vec;
  
  M_chi_data_map = empty_map;
  b_chi_data_map = empty_map;
  elev_data_map = empty_map;
  chi_data_map = empty_map;
  flow_distance_data_map = empty_map;
  drainage_area_data_map = empty_map;
  node_sequence = empty_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, double& x_loc, double& y_loc)
{
  
  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
    
  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns the x and y location of a row and column
// Same as above but with floats
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_x_and_y_locations(int row, int col, float& x_loc, float& y_loc)
{
  
  x_loc = XMinimum + float(col)*DataResolution + 0.5*DataResolution;
    
  // Slightly different logic for y because the DEM starts from the top corner
  y_loc = YMinimum + float(NRows-row)*DataResolution - 0.5*DataResolution;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to convert a node position with a row and column to a lat
// and long coordinate
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_lat_and_long_locations(int row, int col, double& lat, 
                   double& longitude, LSDCoordinateConverterLLandUTM Converter)
{
  // get the x and y locations of the node
  double x_loc,y_loc;
  get_x_and_y_locations(row, col, x_loc, y_loc);
  
  // get the UTM zone of the node
  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);
  //cout << endl << endl << "Line 1034, UTM zone is: " << UTM_zone << endl;
  
  
  if(UTM_zone == NoDataValue)
  {
    lat = NoDataValue;
    longitude = NoDataValue;
  }
  else
  {
    // set the default ellipsoid to WGS84
    int eId = 22;
  
    double xld = double(x_loc);
    double yld = double(y_loc);
  
    // use the converter to convert to lat and long
    double Lat,Long;
    Converter.UTMtoLL(eId, yld, xld, UTM_zone, is_North, Lat, Long);
          
  
    lat = Lat;
    longitude = Long;
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);
  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key] ;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;
    
    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;
        
  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname, 
                                 float A_0, float m_over_n, float area_threshold)
{
  
  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());
  
  chi_map_csv_out.precision(9);
  
  float chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  
  chi_map_csv_out << "latitude,longitude,chi" << endl;
  
  LSDRaster Chi = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(m_over_n, A_0, area_threshold);
  
  float NDV = Chi.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      chi_coord =  Chi.get_data_element(row,col);
      
      if (chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out << latitude << "," << longitude  << "," << chi_coord << endl;
      }
    }
  }
  
  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This prints a chi map to csv with an area threshold in m^2. You feed it the chi map
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_to_csv(LSDFlowInfo& FlowInfo, string chi_map_fname, 
                                 LSDRaster& chi_coord)
{
  
  ofstream chi_map_csv_out;
  chi_map_csv_out.open(chi_map_fname.c_str());
  
  
  
  float this_chi_coord;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  
  chi_map_csv_out << "latitude,longitude,chi" << endl;

  float NDV = chi_coord.get_NoDataValue();

  for(int row = 0; row<NRows; row++)
  {
    for(int col = 0; col<NCols; col++)
    {
      this_chi_coord = chi_coord.get_data_element(row,col);
      
      if (this_chi_coord != NDV)
      {
        get_lat_and_long_locations(row, col, latitude, longitude, Converter);
        chi_map_csv_out.precision(9);
        chi_map_csv_out << latitude << "," << longitude  << ",";
        chi_map_csv_out.precision(5);
        chi_map_csv_out << this_chi_coord << endl;
      }
    }
  }
  
  chi_map_csv_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is for calculating segments from all sources in a DEM
// The sources and their outlets are supplied by the source and outlet nodes
// vectors. These are generated from the LSDJunctionNetwork function
// get_overlapping_channels
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator(LSDFlowInfo& FlowInfo, 
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance, 
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate, 
                                    int target_nodes, 
                                    int n_iterations, int skip,
                                    int minimum_segment_length, float sigma)
{
  
  // IMPORTANT THESE PARAMETERS ARE NOT USED BECAUSE CHI IS CALUCALTED SEPARATELY
  // However we need to give something to pass to the Monte carlo functions
  // even through they are not used (they are inherited)
  float A_0 = 1;
  float m_over_n = 0.5;

  // These elements access the chi data
  vector< vector<float> > chi_m_means;
  vector< vector<float> > chi_b_means;
  vector< vector<float> > chi_coordinates;
  vector< vector<int> > chi_node_indices;
  
  // these are for the individual channels
  vector<float> these_chi_m_means;
  vector<float> these_chi_b_means;
  vector<float> these_chi_coordinates;
  vector<int> these_chi_node_indices;
  
  // these are maps that will store the data
  map<int,float> m_means_map;
  map<int,float> b_means_map;
  map<int,float> chi_coord_map;
  map<int,float> elev_map;
  map<int,float> area_map;
  map<int,float> flow_distance_map;
  vector<int> node_sequence_vec;

  // these are for working with the FlowInfo object
  int this_node,row,col;

  // get the number of channels
  int n_channels = int(source_nodes.size());
  
  for(int chan = 0; chan<n_channels; chan++)
  {
    cout << "Sampling channel " << chan+1 << " of " << n_channels << endl;
    
    // get this particualr channel (it is a chi network with only one channel)
    LSDChiNetwork ThisChiChannel(FlowInfo, source_nodes[chan], outlet_nodes[chan], 
                                Elevation, FlowDistance, DrainageArea,chi_coordinate);
    
    // split the channel
    //cout << "Splitting channels" << endl;
    ThisChiChannel.split_all_channels(A_0, m_over_n, n_iterations, skip, target_nodes, minimum_segment_length, sigma);
    
    // monte carlo sample all channels
    //cout << "Entering the monte carlo sampling" << endl;
    ThisChiChannel.monte_carlo_sample_river_network_for_best_fit_after_breaks(A_0, m_over_n, n_iterations, skip, minimum_segment_length, sigma);
  
    // okay the ChiNetwork has all the data about the m vales at this stage. 
    // Get these vales and print them to a raster
    chi_m_means = ThisChiChannel.get_m_means();
    chi_b_means = ThisChiChannel.get_b_means();
    chi_coordinates = ThisChiChannel.get_chis();
    chi_node_indices = ThisChiChannel.get_node_indices();
    
    // now get the number of channels. This should be 1!
    int n_channels = int(chi_m_means.size());
    if (n_channels != 1)
    {
      cout << "Whoa there, I am trying to make a chi map but something seems to have gone wrong with the channel extraction."  << endl;
      cout << "I should only have one channel per look but I have " << n_channels << " channels." << endl;
    }
    
    // now get the m_means out
    these_chi_m_means = chi_m_means[0];
    these_chi_b_means = chi_b_means[0];
    these_chi_coordinates = chi_coordinates[0];
    these_chi_node_indices = chi_node_indices[0];
    
    //cout << "I have " << these_chi_m_means.size() << " nodes." << endl;
    
    int n_nodes_in_channel = int(these_chi_m_means.size());
    for (int node = 0; node< n_nodes_in_channel; node++)
    {
      
      this_node =  these_chi_node_indices[node];
      //cout << "This node is " << this_node << endl;
      
      // only take the nodes that have not been found
      if (m_means_map.find(this_node) == m_means_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      
        //cout << "This is a new node; " << this_node << endl;
        m_means_map[this_node] = these_chi_m_means[node];
        b_means_map[this_node] = these_chi_b_means[node];
        chi_coord_map[this_node] = these_chi_coordinates[node];
        elev_map[this_node] = Elevation.get_data_element(row,col);
        area_map[this_node] = DrainageArea.get_data_element(row,col);
        flow_distance_map[this_node] = FlowDistance.get_data_element(row,col);
        node_sequence_vec.push_back(this_node);
      }
      else
      {
        //cout << "I already have node: " << this_node << endl;
      }
    }
  }
  
  // set the opject data members
  M_chi_data_map =m_means_map; 
  b_chi_data_map = b_means_map;
  elev_data_map = elev_map;
  chi_data_map = chi_coord_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_sequence_vec;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is a much more rudimentary version that mimics the
// channel steepness caluclations.
// chi needs tobe calculated outside of the function 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::chi_map_automator_rudimentary(LSDFlowInfo& FlowInfo, 
                                    vector<int> source_nodes,
                                    vector<int> outlet_nodes,
                                    LSDRaster& Elevation, LSDRaster& FlowDistance, 
                                    LSDRaster& DrainageArea, LSDRaster& chi_coordinate, 
                                    int regression_nodes)
{

  // the data is stored in maps, for easier testing if a node has been
  // visited. 
  // You might consider having these as data elements in the object so you don't 
  // have to pass them
  map<int,float> gradient_data_map;
  map<int,float> intercept_data_map;
  map<int,float> R2_data_map;
  map<int,float> chi_coordinate_data_map;
  map<int,float> elevation_data_map;
  map<int,float> flow_distance_map;
  map<int,float> area_map;
  vector<int> node_order;
  
  // check if the number of nodes are odd .If not add 1
  if (regression_nodes % 2 == 0)
  {
    cout << "Hello user. You need an odd number of regression nodes." << endl;
    regression_nodes = regression_nodes+1;
    cout << " Changing your regression nodes to " << regression_nodes << endl;
  }
  
  // now get the midpoint
  int mp_nodes = (regression_nodes-1)/2;
  
  //cout << "The number of mp nodes is: " << mp_nodes << endl;
  
  // these keep track of the beginning and ending nodes of a given channel
  int channel_start_node;
  int channel_end_node;
  float channel_end_elevation;
  
  // vectors for holding the chi elevation data
  vector<float> chi_vec;
  vector<float> elev_vec;
  vector<float> empty_vec;
  
  // these are extracted from the channel segments using linear regression
  float intercept,gradient,R_squared;
  
  //float this_chi;
  //float this_elev;
  int this_mp_node;
  int this_end_node;
  int this_start_node;
  
  // these are for getting information out of the FlowInfo object
  int row,col, this_node;
  int r_node, r_row,r_col;          // reciever row and column. 

  // The way this works is that it starts at the top of a channel. It then works 
  // its way down and find the node that is the midpoint and the node that is the
  // end point. The midpoint node is where the data will be recorded.
  // It then puts the data from the start node to the end node into a vector 
  // and performs a linear regression of this vector. The regression data from these
  // vectors are recorded at the nodes.
  // We then want to cover all the nodes with data so what happens if some nodes
  // do not become midpoints? 
  // We could start at the top and get the first midpoint. 
  // From there we can work our way down checking if the top of the regression segment
  // is more than one node down from the end point...

  // get the number of channels
  int n_channels = int(source_nodes.size());
  // now loop through the channels
  for(int chan = 0; chan<n_channels; chan++)
  {
    channel_start_node = source_nodes[chan];
    channel_end_node = outlet_nodes[chan];
    
    // Get the elevation of the end node as a secondary check of the ending of the channel
    // segment
    FlowInfo.retrieve_current_row_and_col(channel_end_node,row,col);
    channel_end_elevation = Elevation.get_data_element(row,col);
    
    // reset the flag for ending the channel
    bool is_end_of_channel = false;

    // set the segment start node to the channel start node
    this_start_node = channel_start_node;

    // now retrieve the midpoint node
    this_node = channel_start_node;
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_mp_node = this_node;
    this_node = r_node;
    
    // now go down one step
    FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    this_node = r_node;
    
    // now get the end node
    for(int n = 0; n<mp_nodes; n++)
    {
      FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
      this_node = r_node;
    }
    this_end_node = this_node;
    
    //================================================
    // This loop is for bug checking
    //this_node = this_start_node;
    //do
    //{
    //  // get the elevation and chi vectors by following the flow
    //  cout << "This node is: " << this_node << endl;
    //  FlowInfo.retrieve_current_row_and_col(this_node,row,col);
    //  FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
    //  this_node = r_node;
    //}
    //while(this_node != this_end_node);
    //
    //cout << "And the midpoint node was: " << this_mp_node << endl;
    //================================================  
      
    // we search down the channel, collecting linear regressions at the 
    // midpoint of the intervals
    while (not is_end_of_channel)
    {
      // get a vector of chi and elevation from the start node to the end node
      chi_vec = empty_vec;
      elev_vec = empty_vec;
      
      // copy the data elements into the vecotrs. This is a little stupid
      // because one might just use a deque to pop the first element
      // and push the last, but the linear regression takes vectors, 
      // not deques so you would have to copy the deque element-wise anyway
      // If you wanted, you could speed this up by implementing a linear regression
      // of deques, but that will need to wait for another day. 
      this_node = this_start_node;
      do
      {
        // get the elevation and chi vectors by following the flow
        FlowInfo.retrieve_current_row_and_col(this_node,row,col);
        chi_vec.push_back(chi_coordinate.get_data_element(row,col));
        elev_vec.push_back(Elevation.get_data_element(row,col));
        
        FlowInfo.retrieve_receiver_information(this_node,r_node,r_row,r_col);
        this_node = r_node;
      } while(this_node != this_end_node);
      
      // do a linear regression on the segment
      least_squares_linear_regression(chi_vec,elev_vec, intercept, gradient, R_squared);
      
      // now add the intercept and gradient data to the correct node
      // only take data that has not been calculated before
      // The channels are in order of descending length so data from
      // longer channels take precidence. 
      if (gradient_data_map.find(this_mp_node) == gradient_data_map.end() )
      {
        FlowInfo.retrieve_current_row_and_col(this_mp_node,row,col);
        gradient_data_map[this_mp_node] = gradient; 
        intercept_data_map[this_mp_node] = intercept;
        R2_data_map[this_mp_node] = R_squared;
        chi_coordinate_data_map[this_mp_node] = chi_coordinate.get_data_element(row,col);
        elevation_data_map[this_mp_node] = Elevation.get_data_element(row,col);
        flow_distance_map[this_mp_node] = FlowDistance.get_data_element(row,col);
        area_map[this_mp_node] = DrainageArea.get_data_element(row,col);
        node_order.push_back(this_mp_node);
      }
      else
      {
        is_end_of_channel = true;
      }
      
      // now move all the nodes down one
      FlowInfo.retrieve_receiver_information(this_start_node,r_node,r_row,r_col);
      this_start_node = r_node;
      
      FlowInfo.retrieve_receiver_information(this_mp_node,r_node,r_row,r_col);
      this_mp_node = r_node;
      
      FlowInfo.retrieve_receiver_information(this_end_node,r_node,r_row,r_col);
      this_end_node = r_node;
      
      // check if we are at the end of the channel
      if (this_end_node == channel_end_node)
      {
        is_end_of_channel = true;
      }
      // also check if the end node is lower elevation than the end node,
      // just to try and stop the channel passing the end node
      FlowInfo.retrieve_current_row_and_col(this_end_node,row,col);
      if (channel_end_elevation > Elevation.get_data_element(row,col))
      {
        is_end_of_channel = true;
      }
    }          // This finishes the regression segment loop
  }            // This finishes the channel and resets channel start and end nodes
  
  // set the data objects
  M_chi_data_map = gradient_data_map; 
  b_chi_data_map = intercept_data_map;
  elev_data_map = elevation_data_map;
  chi_data_map = chi_coordinate_data_map;
  flow_distance_data_map = flow_distance_map;
  drainage_area_data_map = area_map;
  node_sequence = node_order;

  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_full(LSDFlowInfo& FlowInfo, string filename)
{
  
  // these are for extracting element-wise data from the channel profiles. 
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  
  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,chi,elevation,flow distance,drainage area,m_chi,b_chi" << endl;
  
  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter); 

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(5);
      chi_data_out << chi_data_map[this_node] << ","
                   << elev_data_map[this_node] << ","
                   << flow_distance_data_map[this_node] << ","
                   << drainage_area_data_map[this_node] << ","
                   << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Print data maps to file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDChiTools::print_data_maps_to_file_basic(LSDFlowInfo& FlowInfo, string filename)
{
  
  // these are for extracting element-wise data from the channel profiles. 
  int this_node, row,col;
  double latitude,longitude;
  LSDCoordinateConverterLLandUTM Converter;
  
  // open the data file
  ofstream  chi_data_out;
  chi_data_out.open(filename.c_str());
  chi_data_out << "latitude,longitude,m_chi,b_chi" << endl;
  
  // find the number of nodes
  int n_nodes = (node_sequence.size());
  if (n_nodes <= 0)
  {
    cout << "Cannot print since you have not calculated channel properties yet." << endl;
  }
  else
  {
    for (int n = 0; n< n_nodes; n++)
    {
      this_node = node_sequence[n];
      FlowInfo.retrieve_current_row_and_col(this_node,row,col);
      get_lat_and_long_locations(row, col, latitude, longitude, Converter); 

      chi_data_out.precision(9);
      chi_data_out << latitude << ","
                   << longitude << ",";
      chi_data_out.precision(6);
      chi_data_out << M_chi_data_map[this_node] << ","
                   << b_chi_data_map[this_node] << "," << endl;
    }
  }

  chi_data_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




#endif