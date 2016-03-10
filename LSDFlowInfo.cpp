//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo
// Land Surface Dynamics FlowInfo
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for organizing flow routing under the Fastscape algorithm
//  (see Braun and Willett, Geomorphology 2013, v180, p 170-179)
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
// Copyright (C) 2013 Simon M. Mudd 2013
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

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDFlowInfo.cpp
// cpp file for the LSDFlowInfo object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// Fiona Clubb, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.1.0		21/10/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------



#ifndef LSDFlowInfo_CPP
#define LSDFlowInfo_CPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <math.h>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDStatsTools.hpp"
//#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create function, this is empty, you need to include a filename
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::create()
{
  cout << endl << "-------------------------------" << endl;
  cout << "I am an empty flow info object. " << endl;
  cout << "Did you forget to give me a DEM?" << endl;
  cout << endl << "-------------------------------" << endl;
  //exit(EXIT_FAILURE);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Create function, this creates from a pickled file
// fname is the name of the pickled flow info file
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(string fname)
{
  unpickle(fname);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This defaults to no flux boundary conditions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(LSDRaster& TopoRaster)
{
  vector<string> BoundaryConditions(4, "No Flux");
  create(BoundaryConditions, TopoRaster);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function calcualtes the receiver nodes
// it returns the receiver vector r_i
// it also returns a flow direction array in this ordering:
//
// 7  0 1
// 6 -1 2
// 5  4 3
//
// note this is different from ArcMap flowdirection
//	int Arc_flowdir;			// flow direction in arcmap format
//								// 32  64  128
//								// 16  --  1
//								// 8    4  2
// one can convert nthese indices using the LSDIndexRaster object
// note in arc the row index increases down (to the south)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::create(vector<string>& temp_BoundaryConditions,
                         LSDRaster& TopoRaster)
{

  // initialize several data members
  BoundaryConditions = temp_BoundaryConditions;
  //cout << "TBC" << endl;
  NRows = TopoRaster.get_NRows();
  //cout << "Rows: " << NRows << endl;
  NCols = TopoRaster.get_NCols();
  //cout << "Cols: " << NCols << endl;
  XMinimum = TopoRaster.get_XMinimum();
  //cout << "Xmin: " << XMinimum << endl;
  YMinimum = TopoRaster.get_YMinimum();
  //cout << "Ymin: " << YMinimum << endl;
  NoDataValue = int(TopoRaster.get_NoDataValue());
  //cout << "NDV: " << NoDataValue << endl;
  DataResolution = TopoRaster.get_DataResolution();
  //cout << "Data resolution: " <<DataResolution << endl; 	
  
  GeoReferencingStrings = TopoRaster.get_GeoReferencingStrings();
  //cout << "GRS" << endl;

  //cout << "1" << endl;

  // Declare matrices for calculating flow routing
  float one_ov_root2 = 0.707106781;
  float target_elev;				// a placeholder for the elevation of the potential receiver
  float slope;
  float max_slope;				// the maximum slope away from a node
  int max_slope_index;			// index into the maximum slope

  int row, col;						// index for the rows and column
  int receive_row,receive_col;
  string::iterator string_iterator;	// used to get characters from string

  // we need logic for all of the boundaries.
  // there are 3 kinds of edge boundaries:
  // no flux
  // base level
  // periodic
  // These are denoted in a vector of strings.
  // the vector has four elements
  // North boundary, East boundary, South bondary and West boundary
  // the strings can be any length, as long as the first letter corresponds to the
  // first letter of the boundary condition. It is not case sensitive.

  // go through the boundaries
  // A NOTE ON CARDINAL DIRECTIONS
  // If one looks at the raster data, the top row in the data corresponds to the NORTH boundary
  // This is the row that is first read into the code
  // so row 0 is the NORTH boundary
  // row NRows-1 is the SOUTH boundary
  // column 0 is the WEST boundary
  // column NCols-1 is the EAST boundary
  vector<float> slopes(8,NoDataValue);
  vector<int> row_kernal(8);
  vector<int> col_kernal(8);
  int ndv = NoDataValue;
  NDataNodes = 0; 			// the number of nodes in the raster that have data
  int one_if_a_baselevel_node;	// this is a switch used to tag baseleve nodes

  // the first thing you need to do is construct a topoglogy matrix
  // the donor, receiver, etc lists are as long as the number of nodes.
  // these are made of vectors that are exactly dimension n, which is the number of nodes with
  // data. Each of these nodes has several index vectors, that point the program to where the node is
  // we construct these index vectors first
  // we need to loop through all the data before we calcualte slopes because the
  // receiver node indices must be known before the slope calculations are run
  vector<int> empty_vec;
  RowIndex = empty_vec;
  ColIndex = empty_vec;
  BaseLevelNodeList = empty_vec;
  ReceiverVector = empty_vec;
  Array2D<int> ndv_raster(NRows,NCols,ndv);

  NodeIndex = ndv_raster.copy();
  FlowDirection = ndv_raster.copy();
  FlowLengthCode = ndv_raster.copy();


  //cout << "2" << endl;

  // loop through the topo data finding places where there is actually data
  for (row = 0; row<NRows; row++)
    {
      for (col = 0; col<NCols; col++)
	{
	  // only do calcualtions if there is data
	  if(TopoRaster.RasterData[row][col] != NoDataValue)
	    {
	      RowIndex.push_back(row);
	      ColIndex.push_back(col);
	      NodeIndex[row][col] = NDataNodes;
	      NDataNodes++;
	    }
	}
    }

  //cout << "3" << endl;
  
  // now the row and col index are populated by the row and col of the node in row i
  // and the node index has the indeces into the row and col vectors
  // next up, make d, delta, and D vectors
  vector<int> ndn_vec(NDataNodes,0);
  vector<int> ndn_nodata_vec(NDataNodes,ndv);
  vector<int> ndn_plusone_vec(NDataNodes+1,0);
  vector<int> w_vector(NDataNodes,0);

  NDonorsVector = ndn_vec;
  DonorStackVector = ndn_vec;
  DeltaVector = ndn_plusone_vec;

  SVector = ndn_nodata_vec;
  BLBasinVector = ndn_nodata_vec;

  // this vector starts out empty and then base level nodes are added to it
  for (row = 0; row<NRows; row++)
    {
      for (col = 0; col<NCols; col++)
	{
	  // only do calcualtions if there is data
	  if(TopoRaster.RasterData[row][col] != NoDataValue)
	    {
	      // calcualte 8 slopes
	      // no slopes mean get NoDataValue entries
	      // the algorithm loops through the neighbors to the cells, collecting
	      // receiver indices. The order is
	      // 7 0 1
	      // 6 - 2
	      // 5 4 3
	      // where the above directions are cardinal directions
	      // do slope 0
	      row_kernal[0] = row-1;
	      row_kernal[1] = row-1;
	      row_kernal[2] = row;
	      row_kernal[3] = row+1;
	      row_kernal[4] = row+1;
	      row_kernal[5] = row+1;
	      row_kernal[6] = row;
	      row_kernal[7] = row-1;

	      col_kernal[0] = col;
	      col_kernal[1] = col+1;
	      col_kernal[2] = col+1;
	      col_kernal[3] = col+1;
	      col_kernal[4] = col;
	      col_kernal[5] = col-1;
	      col_kernal[6] = col-1;
	      col_kernal[7] = col-1;

	      // check for periodic boundary conditions
	      if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
		{
		  if( BoundaryConditions[2].find("P") != 0 && BoundaryConditions[2].find("p") != 0 )
		    {
		      cout << "WARNING!!! North boundary is periodic! Changing South boundary to periodic" << endl;
		      BoundaryConditions[2] = "P";
		    }
		}
	      if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0 )
		{
		  if( BoundaryConditions[3].find("P") != 0 && BoundaryConditions[3].find("p") != 0 )
		    {
		      cout << "WARNING!!! East boundary is periodic! Changing West boundary to periodic" << endl;
		      BoundaryConditions[3] = "P";
		    }
		}
	      if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0 )
		{
		  if( BoundaryConditions[0].find("P") != 0 && BoundaryConditions[0].find("p") != 0 )
		    {
		      cout << "WARNING!!! South boundary is periodic! Changing North boundary to periodic" << endl;
		      BoundaryConditions[0] = "P";
		    }
		}
	      if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0 )
		{
		  if( BoundaryConditions[1].find("P") != 0 && BoundaryConditions[1].find("p") != 0 )
		    {
		      cout << "WARNING!!! West boundary is periodic! Changing East boundary to periodic" << endl;
		      BoundaryConditions[1] = "P";
		    }
		}

	      // reset baselevel switch for boundaries
	      one_if_a_baselevel_node = 0;

	      // NORTH BOUNDARY
	      if (row == 0)
		{
		  if( BoundaryConditions[0].find("B") == 0 || BoundaryConditions[0].find("b") == 0 )
		    {
		      one_if_a_baselevel_node = 1;
		    }
		  else
		    {
		      // if periodic, reflect across to south boundary
		      if( BoundaryConditions[0].find("P") == 0 || BoundaryConditions[0].find("p") == 0 )
			{
			  row_kernal[0] = NRows-1;
			  row_kernal[1] = NRows-1;
			  row_kernal[7] = NRows-1;
			}
		      else
			{
			  row_kernal[0] = ndv;
			  row_kernal[1] = ndv;
			  row_kernal[7] = ndv;
			}
		    }
		}
	      // EAST BOUNDAY
	      if (col == NCols-1)
		{
		  if( BoundaryConditions[1].find("B") == 0 || BoundaryConditions[1].find("b") == 0 )
		    {
		      one_if_a_baselevel_node = 1;
		    }
		  else
		    {
		      if( BoundaryConditions[1].find("P") == 0 || BoundaryConditions[1].find("p") == 0)
			{
			  col_kernal[1] = 0;
			  col_kernal[2] = 0;
			  col_kernal[3] = 0;
			}
		      else
			{
			  col_kernal[1] = ndv;
			  col_kernal[2] = ndv;
			  col_kernal[3] = ndv;
			}
		    }
		}
	      // SOUTH BOUNDARY
	      if (row == NRows-1)
		{
		  if( BoundaryConditions[2].find("B") == 0 || BoundaryConditions[2].find("b") == 0 )
		    {
		      one_if_a_baselevel_node = 1;
		    }
		  else
		    {
		      if( BoundaryConditions[2].find("P") == 0 || BoundaryConditions[2].find("p") == 0)
			{
			  row_kernal[3] = 0;
			  row_kernal[4] = 0;
			  row_kernal[5] = 0;
			}
		      else
			{
			  row_kernal[3] = ndv;
			  row_kernal[4] = ndv;
			  row_kernal[5] = ndv;
			}
		    }
		}
	      // WEST BOUNDARY
	      if (col == 0)
		{
		  if( BoundaryConditions[3].find("B") == 0 || BoundaryConditions[3].find("b") == 0 )
		    {
		      one_if_a_baselevel_node = 1;
		    }
		  else
		    {
		      if( BoundaryConditions[3].find("P") == 0 || BoundaryConditions[3].find("p") == 0)
			{
			  col_kernal[5] = NCols-1;
			  col_kernal[6] = NCols-1;
			  col_kernal[7] = NCols-1;
			}
		      else
			{
			  col_kernal[5] = ndv;
			  col_kernal[6] = ndv;
			  col_kernal[7] = ndv;
			}
		    }
		}

	      // now loop through the surrounding nodes, calculating the slopes
	      // slopes with NoData get NoData slopes
	      // reminder of ordering:
	      // 7 0 1
	      // 6 - 2
	      // 5 4 3
	      // first logic for baselevel node
	      if (one_if_a_baselevel_node == 1)
		{
		  // get reciever index
		  FlowDirection[row][col] = -1;
		  ReceiverVector.push_back(NodeIndex[row][col]);
		  FlowLengthCode[row][col] = 0;
		}
	      // now the rest of the nodes
	      else
		{
		  FlowLengthCode[row][col] = 0;		// set flow length code to 0, this gets reset
		  // if there is a maximum slope
		  max_slope = 0;
		  max_slope_index = -1;
		  receive_row = row;
		  receive_col = col;
		  for (int slope_iter = 0; slope_iter<8; slope_iter++)
		    {
		      if (row_kernal[slope_iter] == ndv || col_kernal[slope_iter] == ndv)
			{
			  slopes[slope_iter] = NoDataValue;
			}
		      else
			{
			  target_elev = TopoRaster.RasterData[ row_kernal[slope_iter] ][ col_kernal[slope_iter] ];
			  if(target_elev == NoDataValue)
			    {
			      slopes[slope_iter] = NoDataValue;
			    }
			  else
			    {
			      if(slope_iter%2 == 0)
				{
				  //cout << "LINE 988, cardinal direction, slope iter = " << slope_iter << endl;
				  slope = TopoRaster.RasterData[row][col]-target_elev;
				}
			      else
				{
				  slope = one_ov_root2*(TopoRaster.RasterData[row][col]-target_elev);
				}

			      if (slope > max_slope)
				{
				  max_slope_index = slope_iter;
				  receive_row = row_kernal[slope_iter];
				  receive_col = col_kernal[slope_iter];
				  max_slope = slope;
				  if(slope_iter%2 == 0)
				    {
				      FlowLengthCode[row][col] = 1;
				    }
				  else
				    {
				      FlowLengthCode[row][col] = 2;
				    }
				}
			    }
			}
		    }
		  // get reciever index
		  FlowDirection[row][col] = max_slope_index;
		  ReceiverVector.push_back(NodeIndex[receive_row][receive_col]);
		}		// end if baselevel boundary  conditional

				// if the node is a base level node, add it to the base level node list
	      if (FlowLengthCode[row][col] == 0)
		{
		  BaseLevelNodeList.push_back(NodeIndex[row][col]);
		}
	    }			// end if there is data conditional
	}				// end col loop
    }					// end row loop



  // first create the number of donors vector
  // from braun and willett eq. 5
  for(int i = 0; i<NDataNodes; i++)
    {
      NDonorsVector[ ReceiverVector[i] ]++;
    }


  // now create the delta vector
  // this starts on the last element and works its way backwards
  // from Braun and Willett eq 7 and 8
  DeltaVector[NDataNodes] = NDataNodes;
  for(int i = NDataNodes; i>0; i--)
    {
      DeltaVector[i-1] = DeltaVector[i] -  NDonorsVector[i-1];
    }


  // now the DonorStack and the r vectors. These come from Braun and Willett
  // equation 9.
  // Note that in the manscript I have there is a typo in eqaution 9
  // (Jean Braun's code is correct)
  // it should be w_{r_i} = w_{r_i}+1
  int r_index;
  int w_index;
  int delta_index;
  for (int i = 0; i<NDataNodes; i++)
    {
      r_index = ReceiverVector[i];
      delta_index = DeltaVector[ r_index ];
      w_index = w_vector[ r_index ];
      DonorStackVector[  delta_index+w_index ] = i;
      w_vector[r_index] += 1;
      //cout << "i: " << i << " r_i: " << r_index << " delta_i: " << delta_index << " w_index: " << w_index << endl;
    }


  // now go through the base level node list, building the drainage tree for each of these nodes as one goes along
  int n_base_level_nodes;
  n_base_level_nodes = BaseLevelNodeList.size();

  int k;
  int j_index;
  int begin_delta_index, end_delta_index;
  int l_index;

  j_index = 0;
  for (int i = 0; i<n_base_level_nodes; i++)
    {
      k = BaseLevelNodeList[i];			// set k to the base level node

      // This doesn't seem to be in Braun and Willet but to get the ordering correct you
      // need to make sure that the base level node appears first in the donorstack
      // of nodes contributing to the baselevel node.
      // For example, if base level node is 4, with 4 donors
      // and the donor stack has 3 4 8 9
      // the code has to put the 4 first.
      if (DonorStackVector[ DeltaVector[k] ] != k)
	{
	  int this_index = DonorStackVector[ DeltaVector[k] ];
	  int bs_node = k;

	  for(int ds_node = 1; ds_node < NDonorsVector[k]; ds_node++)
	    {
	      if( DonorStackVector[ DeltaVector[k] + ds_node ] == bs_node )
		{
		  DonorStackVector[ DeltaVector[k] ] = k;
		  DonorStackVector[ DeltaVector[k] + ds_node ] = this_index;
		}
	    }
	}

      // now run recursive algorithm
      begin_delta_index = DeltaVector[k];
      end_delta_index = DeltaVector[k+1];

      //cout << "base_level_node is: " << k << " begin_index: " << begin_delta_index << " end: " << end_delta_index << endl;

      for (int delta_index = begin_delta_index; delta_index<end_delta_index; delta_index++)
	{
	  l_index = DonorStackVector[delta_index];
	  add_to_stack(l_index, j_index, k);
	}
    }

  // now calcualte the indices
  calculate_upslope_reference_indices();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Checks to see is a point is in the raster
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool LSDFlowInfo::check_if_point_is_in_raster(float X_coordinate, float Y_coordinate)
{
  bool is_in_raster = true;
  
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
    {
      is_in_raster = false;
    }

  return is_in_raster;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the reciever of current_node (its node, row, and column)
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_receiver_information(int current_node,
						int& receiver_node, int& receiver_row,
						int& receiver_col)
{
  int rn, rr, rc;
  rn = ReceiverVector[current_node];
  rr = RowIndex[rn];
  rc = ColIndex[rn];
  receiver_node = rn;
  receiver_row = rr;
  receiver_col = rc;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the row and column of the current node
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::retrieve_current_row_and_col(int current_node,int& curr_row,
					       int& curr_col)
{
  int cr, cc;
  cr = RowIndex[current_node];
  cc = ColIndex[current_node];
  curr_row = cr;
  curr_col = cc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// algorithms for searching the vectors
// This gets the row and column of the current node
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file(vector<int>& nodeindex_vec, string outfilename)
{

  // fid the last '.' in the filename to use in the scv filename
  unsigned dot = outfilename.find_last_of(".");

  string prefix = outfilename.substr(0,dot);
  //string suffix = str.substr(dot);
  string insert = "_nodeindices_for_Arc.csv";
  string outfname = prefix+insert;

  cout << "the Arc filename is: " << outfname << endl;

  int n_nodes = nodeindex_vec.size();
  int n_nodeindeces = RowIndex.size();

  // open the outfile
  ofstream csv_out;
  csv_out.open(outfname.c_str());
  csv_out.precision(8);

  csv_out << "x,y,node,row,col" << endl;

  int current_row, current_col;
  float x,y;

  // loop through node indices in vector
  for (int i = 0; i<n_nodes; i++)
    {
      int current_node = nodeindex_vec[i];

      // make sure the nodeindex isn't out of bounds
      if (current_node < n_nodeindeces)
	{
	  // get the row and column
	  retrieve_current_row_and_col(current_node,current_row,
				       current_col);

	  // get the x and y location of the node
	  // the last 0.0001*DataResolution is to make sure there are no integer data points
	  x = XMinimum + float(current_col)*DataResolution + 0.5*DataResolution + 0.0001*DataResolution;

	  // the last 0.0001*DataResolution is to make sure there are no integer data points
	  // y coord a bit different since the DEM starts from the top corner
	  y = YMinimum + float(NRows-current_row)*DataResolution - 0.5*DataResolution + 0.0001*DataResolution;;
	  csv_out << x << "," << y << "," << current_node << "," << current_row << "," << current_col << endl;
	}
    }

  csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Write nodeindex vector to csv file, and give each row a unique ID
//
// SWDG after SMM 2/2/2016
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_vector_of_nodeindices_to_csv_file_Unique(vector<int>& nodeindex_vec, string outfilename)
{

  // fid the last '.' in the filename to use in the scv filename
  unsigned dot = outfilename.find_last_of(".");

  string prefix = outfilename.substr(0,dot);
  //string suffix = str.substr(dot);
  string insert = "_nodeindices_for_Arc.csv";
  string outfname = prefix+insert;

  cout << "the Arc filename is: " << outfname << endl;

  int n_nodes = nodeindex_vec.size();
  int n_nodeindeces = RowIndex.size();

  // open the outfile
  ofstream csv_out;
  csv_out.open(outfname.c_str());
  csv_out.precision(8);

  csv_out << "x,y,node,row,col,unique_ID" << endl;

  int current_row, current_col;
  float x,y;

  // loop through node indices in vector
  for (int i = 0; i<n_nodes; i++)
    {
      int current_node = nodeindex_vec[i];

      // make sure the nodeindex isn't out of bounds
      if (current_node < n_nodeindeces)
	{
	  // get the row and column
	  retrieve_current_row_and_col(current_node,current_row,
				       current_col);

	  // get the x and y location of the node
	  // the last 0.0001*DataResolution is to make sure there are no integer data points
	  x = XMinimum + float(current_col)*DataResolution + 0.5*DataResolution + 0.0001*DataResolution;

	  // the last 0.0001*DataResolution is to make sure there are no integer data points
	  // y coord a bit different since the DEM starts from the top corner
	  y = YMinimum + float(NRows-current_row)*DataResolution - 0.5*DataResolution + 0.0001*DataResolution;;
	  csv_out << x << "," << y << "," << current_node << "," << current_row << "," << current_col << "," << i << endl;
	}
    }

  csv_out.close();
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function returns the base level node with the greatest drainage area
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::retrieve_largest_base_level()
{
  int n_bl = BaseLevelNodeList.size();		// get the number of baselevel nodes
  int max_bl = 0;
  for (int i = 0; i<n_bl; i++)
    {
      if(NContributingNodes[ BaseLevelNodeList[i] ] > max_bl)
	{
	  max_bl = NContributingNodes[ BaseLevelNodeList[i] ];
	}
    }
  return max_bl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Get the node for a cell at a given row and column
//@author DTM
//@date 08/11/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDFlowInfo::retrieve_node_from_row_and_column(int row, int column)
{
  int Node = NodeIndex[row][column];
  return Node;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// gets a vector of all the donors to a given node
// @author SMM
// @date 19/09/2014
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::retrieve_donors_to_node(int current_node)
{
  // get the numver of donors
  int NDonors = NDonorsVector[current_node];
  
  // create the vecotr of donating nodes
  vector<int> donorvec(NDonors,NoDataValue);
  
  // loop over the donating nodes, getting their nodeindicies
  for(int dnode = 0; dnode<NDonors; dnode++)
    {
      donorvec[dnode] = DonorStackVector[ DeltaVector[current_node]+dnode];  
    }
  return donorvec;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// recursive add_to_stack routine, from Braun and Willett eq. 12 and 13
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDFlowInfo::add_to_stack(int lm_index, int& j_index, int bl_node)
{
  //cout << "j_index: " << j_index << " and s_vec: " << lm_index << endl;

  SVector[j_index] = lm_index;
  BLBasinVector[j_index] = bl_node;
  j_index++;


  int begin_m,end_m;
  int l_index;
  // if donating to itself, need escape hatch
  if ( lm_index == bl_node)
    {
      begin_m = 0;
      end_m = 0;
    }
  else
    {
      begin_m = DeltaVector[lm_index];
      end_m =  DeltaVector[ lm_index+1];
    }
  //cout << "lm_index: " << lm_index << " begin_m: " << begin_m << " end m: " << end_m << endl;
  for( int m_index = begin_m; m_index<end_m; m_index++)
    {
      //cout << "recursion, begin_m: " << begin_m << " and end_m: " << end_m << endl;
      l_index = DonorStackVector[m_index];
      add_to_stack(l_index, j_index, bl_node);
    }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function pickles the data from the flowInfo object into a binary format
// which can be read by the unpickle function later
// the filename DOES NOT include and extension: this is added by the
// function
//
// WARNING: These files are HUGE and testing indicates they don't save much time
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::pickle(string filename)
{
  string ext = ".FIpickle";
  string hdr_ext = ".FIpickle.hdr";

  string hdr_fname = filename+hdr_ext;
  string data_fname = filename+ext;

  ofstream header_out;
  header_out.open(hdr_fname.c_str());

  int contributing_nodes = int(NContributingNodes.size());
  int BLNodes = int(BaseLevelNodeList.size());

  // print the header file
  header_out <<  "ncols         		" << NCols
	     << "\nnrows         		" << NRows
	     << "\nxllcorner     		" << setprecision(14) << XMinimum
	     << "\nyllcorner     		" << setprecision(14) << YMinimum
	     << "\ncellsize      		" << DataResolution
	     << "\nNODATA_value  		" << NoDataValue
	     << "\nNDataNodes    		" << NDataNodes
	     << "\nNBaseLevelNodes    " << BLNodes
	     << "\nNContributingNodes " << contributing_nodes
	     << "\nBoundaryConditions ";
  for(int i = 0; i<4; i++)
    {
      header_out << " " << BoundaryConditions[i];
    }
  header_out << endl;
  header_out.close();


  cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
  cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
  cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
  cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
  cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
  cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;


  // now do the main data
  ofstream data_ofs(data_fname.c_str(), ios::out | ios::binary);
  int temp;
  for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
	{
	  temp = NodeIndex[i][j];
	  data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
    }
  for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
	{
	  temp = FlowDirection[i][j];
	  data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
    }
  for (int i=0; i<NRows; ++i)
    {
      for (int j=0; j<NCols; ++j)
	{
	  temp = FlowLengthCode[i][j];
	  data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
	}
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = RowIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = ColIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<BLNodes; i++)
    {
      temp = BaseLevelNodeList[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = NDonorsVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = ReceiverVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes+1; i++)
    {
      temp = DeltaVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = DonorStackVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = SVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = BLBasinVector[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<NDataNodes; i++)
    {
      temp = SVectorIndex[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }
  for (int i = 0; i<contributing_nodes; i++)
    {
      temp = NContributingNodes[i];
      data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
    }

  data_ofs.close();


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this unpickles a pickled flow info object. It is folded into a create function
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::unpickle(string filename)
{

  string ext = ".FIpickle";
  string hdr_ext = ".FIpickle.hdr";

  string hdr_fname = filename+hdr_ext;
  string data_fname = filename+ext;

  ifstream header_in;
  header_in.open(hdr_fname.c_str());

  string temp_str;
  int contributing_nodes;
  vector<string> bc(4);
  int BLNodes;

  header_in >> temp_str >> NCols >> temp_str >> NRows >> temp_str >> XMinimum
	    >> temp_str >> YMinimum >> temp_str >> DataResolution
	    >> temp_str >> NoDataValue >> temp_str >> NDataNodes
	    >> temp_str >> BLNodes
	    >> temp_str >> contributing_nodes
	    >> temp_str >> bc[0] >> bc[1] >> bc[2] >> bc[3];
  header_in.close();
  BoundaryConditions = bc;


  // now read the data, using the binary stream option
  ifstream ifs_data(data_fname.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
    {
      cout << "\nFATAL ERROR: the data file \"" << data_fname
	   << "\" doesn't exist" << endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      // initialze the arrays
      Array2D<int> data_array(NRows,NCols,NoDataValue);
      NodeIndex = data_array.copy();
      FlowDirection = data_array.copy();
      FlowLengthCode = data_array.copy();

      vector<int> data_vector(NDataNodes,NoDataValue);
      vector<int> BLvector(BLNodes,NoDataValue);
      vector<int> deltaV(NDataNodes+1,NoDataValue);
      vector<int> CNvec(contributing_nodes,NoDataValue);

      int temp;
      for (int i=0; i<NRows; ++i)
	{
	  for (int j=0; j<NCols; ++j)
	    {
	      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	      NodeIndex[i][j] =temp;
	    }
	}
      for (int i=0; i<NRows; ++i)
	{
	  for (int j=0; j<NCols; ++j)
	    {
	      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	      FlowDirection[i][j] =temp;
	    }
	}
      for (int i=0; i<NRows; ++i)
	{
	  for (int j=0; j<NCols; ++j)
	    {
	      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	      FlowLengthCode[i][j] =temp;
	    }
	}
      RowIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  RowIndex[i] =temp;

	}
      ColIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  ColIndex[i] =temp;

	}
      BaseLevelNodeList = BLvector;
      for (int i=0; i<BLNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  BaseLevelNodeList[i] =temp;

	}
      NDonorsVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  NDonorsVector[i] =temp;

	}
      ReceiverVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  ReceiverVector[i] =temp;

	}
      DeltaVector = deltaV;
      for (int i=0; i<NDataNodes+1; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  DeltaVector[i] =temp;

	}
      DonorStackVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  DonorStackVector[i] =temp;

	}
      SVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  SVector[i] =temp;

	}
      BLBasinVector = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  BLBasinVector[i] =temp;

	}
      SVectorIndex = data_vector;
      for (int i=0; i<NDataNodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  SVectorIndex[i] =temp;

	}
      NContributingNodes = CNvec;
      for (int i=0; i<contributing_nodes; ++i)
	{
	  ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	  NContributingNodes[i] =temp;

	}


    }
  ifs_data.close();

  cout << "sizes RC indices: " << RowIndex.size() << " " << ColIndex.size() << endl;
  cout << "BLNL size: " << BaseLevelNodeList.size() << endl;
  cout << "donors: " << NDonorsVector.size() << " Reciev: " << ReceiverVector.size() << endl;
  cout << "delta: " << DeltaVector.size() << " S: " << SVector.size() << endl;
  cout << "donorstack: " << DonorStackVector.size() << " BBasin: " << BLBasinVector.size() << endl;
  cout << "SVectorIndex " << SVectorIndex.size() << " NContrib: " << NContributingNodes.size() << endl;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to ingest the channel heads raster generated using channel_heads_driver.cpp
// into a vector of source nodes so that an LSDJunctionNetwork can be created easily
// from them. Assumes the FlowInfo object has the same dimensions as the channel
// heads raster.
//
// Takes the filename and extension of the channel heads raster.
//
// SWDG 05/12/12
//
// Update: 6/6/14 Happy 3rd birthday Skye!!!!
// SMM
// Now if the file extension is "csv" then the script reads a csv channel heads
// file
//
// Update 30/09/14 Altered structure of function, but key difference is that it
// is now much better in how it goes about reading in channel heads using
// the coordinates, so that channel heads for a region determined usiong one DEM
// can be loaded in to another covering a subsample of the area, or a different
// resolution, which was impossible before.
// DTM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::Ingest_Channel_Heads(string filename, string extension, int input_switch){

  vector<int> Sources;
  int CH_node;

  // if this is a csv file, read its contents directly into the node index vector
  if(extension == "csv")
    {
      if(input_switch != 0 && input_switch != 1 && input_switch != 2)
  {
    cout << "\t Note, you have specified an unsupported value for the input switch.  Note: \n\t\t 0=take node index\n\t\t 1=take row and column indices\n\t\t 2=take x and y coordinates"  << endl;
    cout << "\t ...taking node index by default" << endl;
  }
      ifstream ch_csv_in;
      string fname = filename +"."+extension;
      ch_csv_in.open(fname.c_str());
    
      if(not ch_csv_in.good())
  {
    cout << "Hey DUDE, you are trying to ingest a sources file that doesn't exist!!" << endl;
    cout << "Check your filename" << endl;
    exit(EXIT_FAILURE);
  }
    
      cout << "fname is: " << fname << endl;
    
      string sline = "";
      getline(ch_csv_in,sline);
  
      vector<int> nodeindex,rowindex,colindex;
      vector<float> x_coord,y_coord;    
      while(!ch_csv_in.eof())
  {   
    char name[256];
    ch_csv_in.getline(name,256);
    sline = name;
      
    // a very tedious way to get the right bit of data. There is probably a 
    // better way to do this but this way works
    if (sline.size() > 0)
      {
        // column index
        string prefix = sline.substr(0,sline.size());
        unsigned comma = sline.find_last_of(",");
        string suffix = prefix.substr(comma+1,prefix.size()); 
        colindex.push_back(atoi(suffix.c_str()));         
        // row index
        prefix = sline.substr(0,comma);
        comma = prefix.find_last_of(",");
        suffix = prefix.substr(comma+1,prefix.size());     
        rowindex.push_back(atoi(suffix.c_str()));        
        // node index
        prefix = sline.substr(0,comma);
        comma = prefix.find_last_of(",");
        suffix = prefix.substr(comma+1,prefix.size());
        nodeindex.push_back(atoi(suffix.c_str()));        
        // y coordinate
        prefix = sline.substr(0,comma);
        comma = prefix.find_last_of(",");
        suffix = prefix.substr(comma+1,prefix.size());     
        y_coord.push_back(atof(suffix.c_str()));      
        // x coordinate
        prefix = sline.substr(0,comma);
        comma = prefix.find_last_of(",");
        suffix = prefix.substr(comma+1,prefix.size());
        x_coord.push_back(atof(suffix.c_str()));
      }
  }
      int node;
      // use row and column indices to locate source nodes.
      if(input_switch == 1)
  {
    for(int i = 0; i < int(rowindex.size()); ++i)
      {
        if(rowindex[i]<NRows && rowindex[i]>=0 && colindex[i]<NCols && colindex[i] >=0 && NodeIndex[rowindex[i]][colindex[i]]!=NoDataValue)
    {
      node = retrieve_node_from_row_and_column(rowindex[i],colindex[i]);
      Sources.push_back(node);
    }
      }
  }
      // Use coordinates to locate source nodes. Note that this enables the use 
      // of LiDAR derived channel heads in coarser DEMs of the same area or
      // subsets of the original DEM for more efficient processing.
      else if(input_switch == 2)
  {
    vector<int> Sources_temp;
    int N_coords = x_coord.size();
    int N_sources_1 = 0;
    for(int i = 0; i < N_coords; ++i)
      {
        node = get_node_index_of_coordinate_point(x_coord[i], y_coord[i]);
        if (node != NoDataValue) 
    {
      // Test 1 - Check for channel heads that fall in same pixel
      int test1 = 0;
      N_sources_1 = Sources_temp.size();
      for(int i_test=0; i_test<N_sources_1;++i_test)
        {
          if(node==Sources_temp[i_test]) test1 = 1;
        }
      if(test1==0) Sources_temp.push_back(node);
      else cout << "\t\t ! removed node from sources list - coincident with another source node" << endl; 
    }
      }
    // Test 2 - Need to do some extra checks to load sources correctly. 
    int N_sources_2 = Sources_temp.size();
    for(int i = 0; i<N_sources_2; ++i)
      {
        int test2 = 0;
        for(int i_test = 0; i_test<int(Sources_temp.size()); ++i_test)
    {
      if(i!=i_test)
        {
          if(is_node_upstream(Sources_temp[i],Sources_temp[i_test])==true) test2 = 1;
        }
    }
        if(test2 ==0) Sources.push_back(Sources_temp[i]);
        else cout << "\t\t ! removed node from sources list - other sources upstream" << endl; 
      }
  }
      // Using Node Index directly (default)
      else Sources = nodeindex;

    }
  
  // if not the code assums a sources raster. 
  else
    {
      LSDIndexRaster CHeads(filename, extension);
    
      for (int i = 0; i < NRows; ++i){
  for (int j = 0; j < NCols; ++j){
    if (CHeads.get_data_element(i,j) != NoDataValue){
      CH_node = retrieve_node_from_row_and_column(i,j);
      if (CH_node != NoDataValue){
        Sources.push_back(CH_node);
      } 
    }
  }
      }
    }
  return Sources;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the flownet information
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::print_flow_info_vectors(string filename)
{
  string string_filename;
  string dot = ".";
  string extension = "txt";
  string_filename = filename+dot+extension;
  cout << "The filename is " << string_filename << endl;

  // print out all the donor, reciever and stack info
  ofstream donor_info_out;
  donor_info_out.open(string_filename.c_str());
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << i << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << ReceiverVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << NDonorsVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes+1; i++)
    {
      donor_info_out << DeltaVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << DonorStackVector[i] << " ";
    }
  donor_info_out << endl;
  for(int i = 0; i<NDataNodes; i++)
    {
      donor_info_out << SVector[i] << " ";
    }
  donor_info_out << endl;

  if( int(SVectorIndex.size()) == NDataNodes)
    {
      for(int i = 0; i<NDataNodes; i++)
	{
	  donor_info_out << SVectorIndex[i] << " ";
	}
      donor_info_out << endl;
      for(int i = 0; i<NDataNodes; i++)
	{
	  donor_info_out << NContributingNodes[i] << " ";
	}
      donor_info_out << endl;
    }

  donor_info_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// these functions write the index arrays to index rasters
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndex_to_LSDIndexRaster()
{
  cout << "NRows: " << NRows << " and NCols: " << NCols << endl; 
  LSDIndexRaster temp_nodeindex(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,NodeIndex,GeoReferencingStrings);
  return temp_nodeindex;
}

LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster()
{
  LSDIndexRaster temp_flowdir(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirection,GeoReferencingStrings);
  return temp_flowdir;
}

LSDIndexRaster LSDFlowInfo::write_FlowLengthCode_to_LSDIndexRaster()
{
  LSDIndexRaster temp_flc(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowLengthCode,GeoReferencingStrings);
  return temp_flc;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes an LSDIndesxRaster given a list of node indices
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster(vector<int>& nodeindexvec)
{
  int n_node_indices = nodeindexvec.size();
  //cout << "The number of nodeindices is: " << n_node_indices << endl;
  Array2D<int> chan(NRows,NCols,NoDataValue);
  
  //cout << "Raster nr: " << chan.dim1() << " nc: " << chan.dim2() << endl;

  int curr_row, curr_col;

  for(int i = 0; i<n_node_indices; i++)
    {
      // make sure there is no segmentation fault for bad data
      // Note: bad data is ignored
      if(nodeindexvec[i] <= NDataNodes)
	{
	  retrieve_current_row_and_col(nodeindexvec[i],curr_row,
				       curr_col);

	  if(chan[curr_row][curr_col] == NoDataValue)
	    {
	      chan[curr_row][curr_col] = 1;
	    }
	  else
	    {
	      chan[curr_row][curr_col]++;
	    }
	}
      else
	{
	  cout << "WARNING: LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster"
	       << " node index does not exist!"<< endl;
	}
    }


  LSDIndexRaster temp_chan(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,chan,GeoReferencingStrings);
  return temp_chan;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes an LSDIndesxRaster given a list of node indices, and give every
// pixel its nodeindex value, which is unique.
//
// SWDG after SMM 2/2/16
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster_Unique(vector<int>& nodeindexvec)
{
  int n_node_indices = nodeindexvec.size();
  Array2D<int> chan(NRows,NCols,NoDataValue);

  int curr_row, curr_col;

  for(int i = 0; i<n_node_indices; i++){
    // make sure there is no segmentation fault for bad data
    // Note: bad data is ignored
    if(nodeindexvec[i] <= NDataNodes){
	    retrieve_current_row_and_col(nodeindexvec[i], curr_row, curr_col);

	    if(chan[curr_row][curr_col] == NoDataValue){
	      chan[curr_row][curr_col] = i;
	    }
	    else{
	      //revisted points will be overwritten, most recent id will be pres
	      chan[curr_row][curr_col] = i;
	    }
	  }
    else
	  {
	    cout << "WARNING: LSDFlowInfo::write_NodeIndexVector_to_LSDIndexRaster"
	       << " node index does not exist!"<< endl;
	  }
  }

  LSDIndexRaster temp_chan(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,chan,GeoReferencingStrings);
  return temp_chan;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// This function used the S vector index and NContributing nodes to calculate the area
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_NContributingNodes_to_LSDIndexRaster()
{
  Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
  int row,col;

  // loop through the node vector, adding pixels to receiver nodes
  for(int node = 0; node<NDataNodes; node++)
    {
      row = RowIndex[node];
      col = ColIndex[node];
      contributing_pixels[row][col] = NContributingNodes[node];
    }

  LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels,GeoReferencingStrings);
  return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// writes a index raster in arc format
// LSD format:
// 7  0 1
// 6 -1 2
// 5  4 3
//	int Arc_flowdir;			// flow direction in arcmap format
//								// 32  64  128
//								// 16   0  1
//								// 8    4  2
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::write_FlowDirection_to_LSDIndexRaster_Arcformat()
{

  Array2D<int> FlowDirectionArc(NRows,NCols,NoDataValue);
  for(int row = 0; row<NRows; row++)
    {
      for (int col = 0; col<NCols; col++)
	{
	  if ( FlowDirection[row][col] == -1)
	    {
	      FlowDirectionArc[row][col] = 0;
	    }
	  else if ( FlowDirection[row][col] == 0)
	    {
	      FlowDirectionArc[row][col] = 64;
	    }
	  else if ( FlowDirection[row][col] == 1)
	    {
	      FlowDirectionArc[row][col] = 128;
	    }
	  else if ( FlowDirection[row][col] == 2)
	    {
	      FlowDirectionArc[row][col] = 1;
	    }
	  else if ( FlowDirection[row][col] == 3)
	    {
	      FlowDirectionArc[row][col] = 2;
	    }
	  else if ( FlowDirection[row][col] == 4)
	    {
	      FlowDirectionArc[row][col] = 4;
	    }
	  else if ( FlowDirection[row][col] == 5)
	    {
	      FlowDirectionArc[row][col] = 8;
	    }
	  else if ( FlowDirection[row][col] == 6)
	    {
	      FlowDirectionArc[row][col] = 16;
	    }
	  else if ( FlowDirection[row][col] == 7)
	    {
	      FlowDirectionArc[row][col] = 32;
	    }
	}
    }

  LSDIndexRaster temp_fd(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,FlowDirectionArc,GeoReferencingStrings);
  return temp_fd;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function writes the drainage area (number of contributing nodes
// * DataResolution^2) to an LSDRaster object
// Added by FC 15/11/12
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::write_DrainageArea_to_LSDRaster()
{
  // initialise the 2D array
  int n_i;								// node index
  float ndv = float(NoDataValue);
  Array2D<float> DrainageArea_local(NRows,NCols,ndv);

  //get the contributing nodes
  for (int row = 0; row < NRows; row++)
    {
      for (int col = 0; col < NCols; col++)
	{
	  n_i = NodeIndex[row][col];
	  DrainageArea_local[row][col] = float(NContributingNodes[n_i])*DataResolution*DataResolution;
	}
    }
  // create the LSDRaster object
  LSDRaster DrainageArea(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,DrainageArea_local,GeoReferencingStrings);
  return DrainageArea;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
// In this function a pixel that has no donors has a contributing pixel value of 0
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDFlowInfo::calculate_n_pixels_contributing_from_upslope()
{
  Array2D<int> contributing_pixels(NRows,NCols,NoDataValue);
  int row,col;
  int receive_row, receive_col;
  int receiver_node;

  // loop through the s vector, adding pixels to receiver nodes
  for(int node = NDataNodes-1; node>=0; node--)
    {

      row = RowIndex[SVector[node]];
      col = ColIndex[SVector[node]];
      // if the pixel exists and has no contributing pixels,
      // change from nodata to zero

      if(contributing_pixels[row][col] == NoDataValue)
	{
	  contributing_pixels[row][col] = 0;
	}

      receiver_node = ReceiverVector[ SVector[node] ] ;
      receive_row = RowIndex[ receiver_node ];
      receive_col = ColIndex[ receiver_node ];

      cout << "node " << node << " pixel: " << SVector[node] << " receiver: " << receiver_node << endl;
      cout << "contributing: " << contributing_pixels[row][col] << endl;

      if ( receiver_node  == SVector[node])
	{
	  // do nothing
	}
      else if ( contributing_pixels[receive_row][receive_col] == NoDataValue)
	{
	  contributing_pixels[receive_row][receive_col] =
	    contributing_pixels[row][col]+1;
	}
      else
	{
	  contributing_pixels[receive_row][receive_col] +=
	    contributing_pixels[row][col]+1;
	}

      cout << "recieving: " << contributing_pixels[receive_row][receive_col] << endl;
    }

  LSDIndexRaster temp_cp(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,contributing_pixels,GeoReferencingStrings);
  return temp_cp;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calcualtes the contributing pixels
// it can be converted to contributing area by multiplying by the
// DataResolution^2
//
// In this function a pixel that has no donors contributes its own flow
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::calculate_upslope_reference_indices()
{
  vector<int> vectorized_area(NDataNodes,1);
  SVectorIndex = vectorized_area;

  int receiver_node;
  int donor_node;

  // loop through the s vector, adding pixels to receiver nodes
  for(int node = NDataNodes-1; node>=0; node--)
    {
      donor_node = SVector[node];
      receiver_node = ReceiverVector[ donor_node ];

      // every node is visited once and only once so we can map the
      // unique positions of the nodes to the SVector
      SVectorIndex[donor_node] = node;

      // add the upslope area (note no action is taken
      // for base level nodes since they donate to themselves and
      // we must avoid float counting
      if (donor_node != receiver_node)
	{
	  vectorized_area[ receiver_node ] +=  vectorized_area[ donor_node ];
	}
    }

  NContributingNodes = vectorized_area;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function returns a integer vector containing all the node numbers upslope
// of of the node with number node_number_outlet
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<int> LSDFlowInfo::get_upslope_nodes(int node_number_outlet)
{
  vector<int> us_nodes;

  if(node_number_outlet < 0 || node_number_outlet > NDataNodes-1)
    {
      cout << "the node index does not exist" << endl;
      exit(EXIT_FAILURE);
    }

  int start_SVector_node = SVectorIndex[node_number_outlet];
  int end_SVector_node = start_SVector_node+NContributingNodes[node_number_outlet];

  for(int node = start_SVector_node; node < end_SVector_node; node++)
    {
      us_nodes.push_back(SVector[node]);
    }

  return us_nodes;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of source nodes and creates a raster where 
// the pixels have a value of 1 where there are upslope nodes and nodata otherwise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_node_mask(vector<int> source_nodes)
{
  // initiate the data array
  Array2D<float> this_raster(NRows,NCols,NoDataValue);
  
  // loop through the nodes, collecting upslope nodes
  int n_nodes = int(source_nodes.size());
  int this_node;
  int us_node;
  int curr_row,curr_col;
  float is_us = 1.0;
  
  // go through all the source nodes, find their upslope nodes
  // and set the value of these nodes to 1.0 on the data array
  for (int n = 0; n<n_nodes; n++)
  {
    this_node = source_nodes[n];
    
    // check if it is in the DEM
    if (this_node < NDataNodes)
    {
      vector<int> upslope_nodes = get_upslope_nodes(this_node);
      
      int n_us_nodes =  int(upslope_nodes.size());
      for(int us = 0; us<n_us_nodes; us++)
      {
        retrieve_current_row_and_col(upslope_nodes[us],curr_row,curr_col);
        this_raster[curr_row][curr_col]=is_us;
      }
    }
  }
  
  // now create the raster
  LSDRaster temp_us(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,this_raster,GeoReferencingStrings);
  return temp_us;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of source nodes and creates a raster where 
// the pixels have a value of upslope_value
// where there are upslope nodes and NoData otherwise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_node_mask(vector<int> source_nodes, vector<float> upslope_values)
{
  // initiate the data array
  Array2D<float> this_raster(NRows,NCols,NoDataValue);
  
  if (source_nodes.size() == upslope_values.size())
  {
    // loop through the nodes, collecting upslope nodes
    int n_nodes = int(source_nodes.size());
    int this_node;
    int us_node;
    int curr_row,curr_col;
    
    // go through all the source nodes, find their upslope nodes
    // and set the value of these nodes to 1.0 on the data array
    for (int n = 0; n<n_nodes; n++)
    {
      this_node = source_nodes[n];
      
      // check if it is in the DEM
      if (this_node < NDataNodes)
      {
        vector<int> upslope_nodes = get_upslope_nodes(this_node);
        
        int n_us_nodes =  int(upslope_nodes.size());
        for(int us = 0; us<n_us_nodes; us++)
        {
          retrieve_current_row_and_col(upslope_nodes[us],curr_row,curr_col);
          this_raster[curr_row][curr_col]= upslope_values[n];
        }
      }
    }
  }
  else
  {
    cout << "The uplsope vlaues vector needs to be the same lengths as the sources vector!" << endl;
    cout << "Returning an nodata raster" << endl;
  }
    
  // now create the raster
  LSDRaster temp_us(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,this_raster,GeoReferencingStrings);
  return temp_us;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//  Accumulate some variable (such a precipitation) from an accumulation raster
//
//  This requires summing all upslope nodes for every node. It seems a bit inefficient
//  but the other simple alternative is to do a sort() operation initially and then
//  move from upslope node down. There is probably a more efficient way to do this
//  and this algorithm should be revisited later to see if we can speed it up. 
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::upslope_variable_accumulator(LSDRaster& accum_raster)
{
  int raster_NRows, raster_NCols;
  float raster_XMin, raster_YMin, raster_DataRes;

  // first check to make sure the raster dimensions match that of the 
  // raster upon which LSDFlowInfo is based
  raster_NRows =  accum_raster.get_NRows();
  raster_NCols =  accum_raster.get_NCols();
  raster_XMin  =  accum_raster.get_XMinimum();
  raster_YMin  =  accum_raster.get_YMinimum();
  raster_DataRes  =  accum_raster.get_DataResolution();
  
  if (raster_NRows != NRows || raster_NCols != NCols ||
      raster_XMin != XMinimum || raster_YMin != YMinimum ||
      raster_DataRes != DataResolution)
    {
      cout << "Warning!!, LSDFlowInfo::upslope_area_accumulator\n"
	   << "Accumulation raster does not match dimensions of original raster" << endl;
      return accum_raster; 
    }
  else
    {
      // create the data array
      Array2D<float> accumulated_data_array(NRows,NCols,NoDataValue);
        
      // loop through all the nodes, accumulating the areas
      for(int this_node = 0; this_node <NDataNodes; this_node++)
	{
	  // get the upslope nodes
	  vector<int> node_vec = get_upslope_nodes(this_node);
      
	  // loop through these nodes, adding them to the accumulator
	  float this_node_accumulated = 0;
	  int this_row, this_col;
	  for (int ni = 0; ni<int(node_vec.size()); ni++)
	    {
	      retrieve_current_row_and_col(node_vec[ni],this_row,this_col);
	      this_node_accumulated += accum_raster.get_data_element(this_row, this_col);
	    }
      
	  // write the accumulated variable to the array
	  int curr_row, curr_col; 
	  retrieve_current_row_and_col(this_node,curr_row,curr_col);
	  accumulated_data_array[curr_row][curr_col] = this_node_accumulated; 
	}
      // create the raster
      LSDRaster accumulated_flow(NRows, NCols, XMinimum, YMinimum, 
				 DataResolution, NoDataValue, accumulated_data_array,GeoReferencingStrings);
      return accumulated_flow;      
    }  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function tests whether one node is upstream of another node
//
// FC 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::is_node_upstream(int current_node, int test_node)
{
  int i = 0;

  int start_SVector_node = SVectorIndex[current_node];
  int end_SVector_node = start_SVector_node+NContributingNodes[current_node];

  int SVector_test_node = SVectorIndex[test_node];

  for(int node = start_SVector_node; node < end_SVector_node; node++)
    {
      if (node == SVector_test_node)
	{
	  i = 1;
	}
    }

  return i;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function tests whether a node is a base level node
//
// FC 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::is_node_base_level(int node)
{
  int i = 0;

  for (int j = 0; j < int(BaseLevelNodeList.size()); j++)
  {
    if (node == BaseLevelNodeList[j])
    {
      i = 1;
    }
  }

  return i;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function redurns a vector of node indices to all the donor
// nodes of a particular node
//
// SMM 21/10/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_donor_nodes(int current_node)
{
  int start_D = DeltaVector[current_node];
  int end_D = DeltaVector[current_node+1];

  vector<int> donor_nodes;
  for(int this_node = start_D; this_node<end_D; this_node++)
    {
      //cout << "node " << current_node << " and donor: " << DonorStackVector[ this_node ] << endl;
      donor_nodes.push_back( DonorStackVector[ this_node ] );
    }

  return donor_nodes;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the chi function for all the nodes upslope a given node
// it takes a node list
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(int starting_node, float m_over_n, float A_0)
{
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);
  vector<float> chi_vec = get_upslope_chi(upslope_pixel_list, m_over_n, A_0);
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calculates the chi function for all the nodes upslope a given node
// it takes a node list
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDFlowInfo::get_upslope_chi(int starting_node, float m_over_n, float A_0,
                                           LSDRaster& Discharge)
{
  vector<int> upslope_pixel_list = get_upslope_nodes(starting_node);
  vector<float> chi_vec = get_upslope_chi(upslope_pixel_list, m_over_n, A_0, Discharge);
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the acutal chi values in a vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
vector<float> LSDFlowInfo::get_upslope_chi(vector<int>& upslope_pixel_list, 
                                           float m_over_n, float A_0)
{

  int receiver_node;
  int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  float pixel_area = DataResolution*DataResolution;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];

  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    node = upslope_pixel_list[n_index];
    receiver_node = ReceiverVector[ node ];
    IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
    row = RowIndex[node];
    col = ColIndex[node];

    if (FlowLengthCode[row][col] == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }


    chi_vec[n_index] = dx*(pow( (A_0/ (float(NContributingNodes[node])*pixel_area) ),m_over_n))
                          + chi_vec[IndexOfReceiverInUplsopePList];

    //	cout << "node: " << upslope_pixel_list[n_index] << " receiver: " << receiver_node
    //	     << " SIndexReciever: " << IndexOfReceiverInUplsopePList 
    //       << " and checked: " << upslope_pixel_list[IndexOfReceiverInUplsopePList]
    //	     << " and chi: " << chi_vec[n_index] << endl;

  }
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function is called from the the get_upslope_chi that only has an integer
// it returns the acutal chi values in a vector
// same as above but uses a discharge raster
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
vector<float> LSDFlowInfo::get_upslope_chi(vector<int>& upslope_pixel_list, 
                                           float m_over_n, float A_0,
                                           LSDRaster& Discharge)
{

  int receiver_node;
  int IndexOfReceiverInUplsopePList;
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;
  float dx;
  int node,row,col;
  // get the number of nodes upslope
  int n_nodes_upslope = upslope_pixel_list.size();
  vector<float> chi_vec(n_nodes_upslope,0.0);

  if(n_nodes_upslope != NContributingNodes[ upslope_pixel_list[0] ])
  {
    cout << "LSDFlowInfo::get_upslope_chi, the contributing pixels don't agree" << endl;
    exit(EXIT_FAILURE);
  }

  int start_SVector_node = SVectorIndex[ upslope_pixel_list[0] ];

  for (int n_index = 1; n_index<n_nodes_upslope; n_index++)
  {
    node = upslope_pixel_list[n_index];
    receiver_node = ReceiverVector[ node ];
    IndexOfReceiverInUplsopePList = SVectorIndex[receiver_node]-start_SVector_node;
    row = RowIndex[node];
    col = ColIndex[node];

    if (FlowLengthCode[row][col] == 2)
    {
      dx = diag_length;
    }
    else
    {
      dx = DataResolution;
    }


    chi_vec[n_index] = dx*(pow( (A_0/ ( Discharge.get_data_element(row, col) ) ),m_over_n))
                          + chi_vec[IndexOfReceiverInUplsopePList];

    //	cout << "node: " << upslope_pixel_list[n_index] << " receiver: " << receiver_node
    //	     << " SIndexReciever: " << IndexOfReceiverInUplsopePList 
    //       << " and checked: " << upslope_pixel_list[IndexOfReceiverInUplsopePList]
    //	     << " and chi: " << chi_vec[n_index] << endl;

  }
  return chi_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of starting nodes and calucaltes chi
// it assumes each chi value has the same base level. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes, 
                              float m_over_n, float A_0, float area_threshold)
{
  // some variables for writing to the raster
  int curr_row;
  int curr_col;
  int curr_node;
  
  float PixelArea = DataResolution*DataResolution;
  float DrainArea;
  
  // an array to hold the chi values
  Array2D<float> new_chi(NRows,NCols,NoDataValue);

  int n_starting_nodes = int(starting_nodes.size());
  
  for(int sn = 0; sn<n_starting_nodes; sn++)
  {
    // first check to see if this node has already been visited. If it has 
    // all upslope nodes have also been visited so there is no point continuing 
    // with this node 
    retrieve_current_row_and_col(starting_nodes[sn],curr_row,curr_col);
    if(new_chi[curr_row][curr_col] == NoDataValue)
    {
      vector<float> us_chi = get_upslope_chi(starting_nodes[sn], m_over_n, A_0);
      vector<int> upslope_pixel_list = get_upslope_nodes(starting_nodes[sn]);
       
      int n_chi_nodes = int(us_chi.size());
      for (int cn = 0; cn<n_chi_nodes; cn++)
      {
        // get the current row and column
        curr_node =  upslope_pixel_list[cn];
        retrieve_current_row_and_col(curr_node,curr_row,curr_col);
          
        // check to see if the drainage area is greater than the threshold
        // if so, calcualte chi
        DrainArea = PixelArea*NContributingNodes[curr_node];
        if(DrainArea > area_threshold)
        {
          new_chi[curr_row][curr_col]= us_chi[cn];
        }      
      }
   }       
  }

  LSDRaster chi_map(NRows, NCols, XMinimum, YMinimum, 
                    DataResolution, NoDataValue, new_chi,GeoReferencingStrings);
  return chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a list of starting nodes and calucaltes chi
// it assumes each chi value has the same base level. 
// Same as above but calculates using discharge
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_multiple_starting_nodes(vector<int>& starting_nodes, 
                              float m_over_n, float A_0, float area_threshold,
                              LSDRaster& Discharge )
{
  // some variables for writing to the raster
  int curr_row;
  int curr_col;
  int curr_node;
  
  float PixelArea = DataResolution*DataResolution;
  float DrainArea;
  
  // an array to hold the chi values
  Array2D<float> new_chi(NRows,NCols,NoDataValue);

  int n_starting_nodes = int(starting_nodes.size());
  
  for(int sn = 0; sn<n_starting_nodes; sn++)
  {
    // first check to see if this node has already been visited. If it has 
    // all upslope nodes have also been visited so there is no point continuing 
    // with this node 
    retrieve_current_row_and_col(starting_nodes[sn],curr_row,curr_col);
    if(new_chi[curr_row][curr_col] == NoDataValue)
    {
      vector<float> us_chi = get_upslope_chi(starting_nodes[sn], m_over_n, A_0,Discharge);
      vector<int> upslope_pixel_list = get_upslope_nodes(starting_nodes[sn]);
       
      int n_chi_nodes = int(us_chi.size());
      for (int cn = 0; cn<n_chi_nodes; cn++)
      {
        // get the current row and column
        curr_node =  upslope_pixel_list[cn];
        retrieve_current_row_and_col(curr_node,curr_row,curr_col);
          
        // check to see if the drainage area is greater than the threshold
        // if so, calcualte chi
        DrainArea = PixelArea*NContributingNodes[curr_node];
        if(DrainArea > area_threshold)
        {
          new_chi[curr_row][curr_col]= us_chi[cn];
        }      
      }
   }       
  }

  LSDRaster chi_map(NRows, NCols, XMinimum, YMinimum, 
                    DataResolution, NoDataValue, new_chi,GeoReferencingStrings);
  return chi_map;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assumes all base level nodes are at the same base level
// and calculates chi for them. Essentially it covers the entire map in
// chi values. 
// This function is probably most appropriate for looking at numberical
// model results
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float A_0, 
                  float area_threshold)
{
  LSDRaster all_chi = get_upslope_chi_from_multiple_starting_nodes(BaseLevelNodeList, 
                                      m_over_n, A_0, area_threshold);
  return all_chi;                       
}                                                
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assumes all base level nodes are at the same base level
// and calculates chi for them. Essentially it covers the entire map in
// chi values. 
// This function is probably most appropriate for looking at numberical
// model results
// same as above but calculates with a discharge
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::get_upslope_chi_from_all_baselevel_nodes(float m_over_n, float A_0, 
                  float area_threshold, LSDRaster& Discharge)
{
  LSDRaster all_chi = get_upslope_chi_from_multiple_starting_nodes(BaseLevelNodeList, 
                                m_over_n, A_0, area_threshold, Discharge);
  return all_chi;                       
}                                                
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// distance from outlet function
// this is overloaded.
// if it isn't provided any argument, it calculates the distance from outlet
// of all the base level nodes
// if it is given a node index number or a row and column, then
// the distance from outlet includes all the distances upstream of that node
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDRaster LSDFlowInfo::distance_from_outlet()
{
  // initialize the array2d that will become the LSDRaster
  float ndv = float(NoDataValue);
  Array2D<float> flow_distance(NRows,NCols,ndv);

  // initialize the 1/root(2)
  float root2 = 1.41421356;
  float diag_length = root2*DataResolution;

  int row,col,bl_row,bl_col,receive_row,receive_col;

  //cout << "FlowLengthCode: " << FlowLengthCode << endl;

  int start_node = 0;
  int end_node;
  int nodes_in_bl_tree;
  int baselevel_node;
  // loop through the base level node list
  int n_base_level_nodes = BaseLevelNodeList.size();
  for(int bl = 0; bl<n_base_level_nodes; bl++)
    {

      baselevel_node = BaseLevelNodeList[bl];

      bl_row = RowIndex[baselevel_node];
      bl_col = ColIndex[baselevel_node];
      // get the number of nodes upslope and including this node
      nodes_in_bl_tree = NContributingNodes[baselevel_node];
      //cout << "LINE 938, FlowInfo, base level: " << bl << " with " << nodes_in_bl_tree << " nodes upstream" << endl;

      end_node = start_node+nodes_in_bl_tree;

      // set the distance of the outlet to zero
      flow_distance[bl_row][bl_col] = 0;

      // now loop through stack
      for(int s_node = start_node; s_node < end_node; s_node++)
	{
	  //cout << "Line 953 flow info, s_node is: " << s_node << endl;

	  //cout << SVector.size() << " " << ReceiverVector.size() << " " << RowIndex.size() << " " << ColIndex.size() << endl;
	  row = RowIndex[ SVector[ s_node]  ];
	  col = ColIndex[ SVector[ s_node]  ];
	  //cout << "got rows and columns " << row << " " << col << endl;
	  receive_row = RowIndex[ ReceiverVector[SVector[s_node] ]];
	  receive_col = ColIndex[ ReceiverVector[SVector[s_node] ]];
	  //cout <<  "get receive " << receive_row << " " << receive_col << endl;

	  if ( FlowLengthCode[row][col] == 1)
	    {
	      flow_distance[row][col] = flow_distance[receive_row][receive_col]+DataResolution;
	    }
	  else if ( FlowLengthCode[row][col] == 2 )
	    {
	      flow_distance[row][col] = flow_distance[receive_row][receive_col]
		+ diag_length;
	    }
	  //cout << "Flow distance: " << flow_distance << endl;
	}
      start_node = end_node;
    }
  //cout << "LINE 971 FlowInfo Flow distance complete, flow_distance is: " << endl;
  //cout << flow_distance << endl;
  LSDRaster FlowLength(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,
                       flow_distance,GeoReferencingStrings);
  return FlowLength;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function get the d8 slope. It points downslope from each node
//
// SMM 22/09/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDFlowInfo::calculate_d8_slope(LSDRaster& Elevation)
{
  float ndv = float(NoDataValue);
  Array2D<float> d8_slope(NRows,NCols,ndv);
	
  // these save a bit of computational expense. 
  float root_2 = pow(2, 0.5);
  float dx_root2 = root_2*DataResolution;
  
  int this_row;
  int this_col;
  int r_row;
  int r_col;
  int r_node;
  float dx;
	
  for (int node = 0; node<NDataNodes; node++)
    {
      // get the row and column
      retrieve_current_row_and_col(node,this_row,this_col);
  
      // get the distance between nodes. Depends on flow direction		
      switch (retrieve_flow_length_code_of_node(node))
	{
	case 0:
	  dx = -99;
	  break;
	case 1:
	  dx = DataResolution;
	  break;
	case 2:
	  dx = dx_root2;
	  break;
	default:
	  dx = -99;
	  break;
	} 
    
      // get the reciever information
      retrieve_receiver_information(node,r_node, r_row, r_col);
    
      // now calculate the slope
      if (r_node == node)
	{
	  d8_slope[this_row][this_col] = 0;  
	}
      else
	{
	  d8_slope[this_row][this_col] = (1/dx)*
	    (Elevation.get_data_element(this_row,this_col)
	     -Elevation.get_data_element(r_row,r_col));
	}
       
    }
  
  LSDRaster d8_slope_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,ndv,
			    d8_slope,GeoReferencingStrings); 
  return d8_slope_raster;                      

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//
// this finds the node that is farthest upstream from a given node
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::find_farthest_upslope_node(int node, LSDRaster& DistFromOutlet)
{
  // set the farthest node to the current node; if the node has no contributing pixels
  // the function will just return itself
  int farthest_upslope_node = node;

  // first get the nodes that are upslope
  vector<int> upslope_node_list = get_upslope_nodes(node);

  int row, col;
  float this_flow_distance;

  // now loop through these, looking for the farthest upstream node
  float farthest = 0.0;
  int n_upslope_nodes = upslope_node_list.size();
  for (int i = 0; i<n_upslope_nodes; i++)
    {
      // get the row and col of upslope nodes
      row = RowIndex[ upslope_node_list[i] ];
      col = ColIndex[ upslope_node_list[i] ];

      // get the flow distance
      this_flow_distance = DistFromOutlet.get_data_element(row, col);

      // test if it is the farthest
      if (this_flow_distance > farthest)
	{
	  farthest = this_flow_distance;
	  farthest_upslope_node =  upslope_node_list[i];
	}
    }

  return farthest_upslope_node;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// get node index of point from X and Y coordinates
// this is different from the above function in that it does not snap to the nearest channel
//
// FJC 11/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDFlowInfo::get_node_index_of_coordinate_point(float X_coordinate, float Y_coordinate)
{
  // Shift origin to that of dataset
  float X_coordinate_shifted_origin = X_coordinate - XMinimum - DataResolution*0.5;
  float Y_coordinate_shifted_origin = Y_coordinate - YMinimum - DataResolution*0.5;;

  // Get row and column of point
  int col_point = int(X_coordinate_shifted_origin/DataResolution);
  int row_point = (NRows-1) - int(round(Y_coordinate_shifted_origin/DataResolution));

  // Get node of point
  int CurrentNode;
  if(col_point>=0 && col_point<NCols && row_point>=0 && row_point<NRows) CurrentNode = retrieve_node_from_row_and_column(row_point, col_point);
  else CurrentNode = NoDataValue;
  return CurrentNode;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses the flow accumulation pixels
//
//
// SMM 01/06/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_index_threshold(LSDIndexRaster& FlowPixels, int threshold)
{
  vector<int> sources;
  int row,col;
  //int n_donors;
  int donor_row,donor_col;
  int thresh_switch;
  int donor_node;

  // drop down through the stack
  // if the node is greater than or equal to the threshold
  // check all the donors
  // if there are no donors, then the node is a source
  // if none of the donors are greater than the threshold, then it also is a source
  for (int node = 0; node<NDataNodes; node++)
    {
      row = RowIndex[node];
      col = ColIndex[node];

      // see if node is greater than threshold
      if(FlowPixels.get_data_element(row,col)>=threshold)
  {
    //cout << "node " << node << " is a potential source, it has a value of "
    //     << FlowPixels.get_data_element(row,col)
    //     << "and it has " << NDonorsVector[node] <<" donors " << endl;

    // if it doesn't have donors, it is a source
    if(NDonorsVector[node] == 0)
      {
        sources.push_back(node);
      }
    else
      {
        thresh_switch = 1;
        // figure out where the donor nodes are, and if
        // the donor node is greater than the threshold
        for(int dnode = 0; dnode<NDonorsVector[node]; dnode++)
    {
      donor_node = DonorStackVector[ DeltaVector[node]+dnode];
      donor_row = RowIndex[ donor_node ];
      donor_col = ColIndex[ donor_node ];

      // we don't float count base level nodes, which donate to themselves
      if (donor_node != node)
      {
        // if the donor node is greater than the threshold,
        // then this node is not a threhold
        if(FlowPixels.get_data_element(donor_row,donor_col)>=threshold)
        {
          thresh_switch = 0;
        }
      }

      //cout << "thresh_switch is: " << thresh_switch << endl;
    }
        // if all of the donors are below the threhold, this is a source
        if (thresh_switch == 1)
    {
      sources.push_back(node);
    }
      }
  }
    }
  return sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses a threshold of drainage area * slope^2
//
//
// FJC 11/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_slope_area(LSDIndexRaster& FlowPixels, LSDRaster& Slope, int threshold)
{
  vector<int> sources;
  int row,col;
  //int n_donors;
  int donor_row,donor_col;
  int thresh_switch;
  int donor_node;

  // drop down through the stack
  // if the node is greater than or equal to the threshold
  // check all the donors
  // if there are no donors, then the node is a source
  // if none of the donors are greater than the threshold, then it also is a source
  for (int node = 0; node<NDataNodes; node++)
    {
      row = RowIndex[node];
      col = ColIndex[node];

      float area = FlowPixels.get_data_element(row,col);
      float slope = Slope.get_data_element(row,col);
      float SA_product = area * (slope*slope);
      // see if node is greater than threshold
      if(SA_product >= threshold)
	{
	  //cout << "node " << node << " is a potential source, it has a value of "
	  //     << SA_product
	  //     << "and it has " << NDonorsVector[node] <<" donors " << endl;

	  // if it doesn't have donors, it is a source
	  if(NDonorsVector[node] == 0)
	    {
	      sources.push_back(node);
	    }
	  else
	    {
	      thresh_switch = 1;
	      // figure out where the donor nodes are, and if
	      // the donor node is greater than the threshold
	      for(int dnode = 0; dnode<NDonorsVector[node]; dnode++)
		{
		  donor_node = DonorStackVector[ DeltaVector[node]+dnode];
		  donor_row = RowIndex[ donor_node ];
		  donor_col = ColIndex[ donor_node ];

		  // we don't float count base level nodes, which donate to themselves
		  if (donor_node != node)
		    {
		      // if the donor node is greater than the threshold,
		      // then this node is not a threhold
		      float area_donor = FlowPixels.get_data_element(donor_row,donor_col);
		      float slope_donor = Slope.get_data_element(donor_row,donor_col);
		      float SA_product_donor = area_donor * (slope_donor*slope_donor);
		      if(SA_product_donor >= threshold)
			{
			  thresh_switch = 0;
			}
		    }

		  //cout << "thresh_switch is: " << thresh_switch << endl;
		}
	      // if all of the donors are below the threhold, this is a source
	      if (thresh_switch == 1)
		{
		  sources.push_back(node);
		}
	    }
	}
    }
  return sources;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// a get sources version that uses the X and Y coordinates of mapped channel heads
//
//
// FJC 17/02/14
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::get_sources_from_mapped_channel_heads(vector<float>& X_coords, vector<float>& Y_coords)
{
  vector<int> SourceNodes;
  cout << "N of channel heads: " << X_coords.size() << endl;
  for (unsigned int i = 0; i < X_coords.size(); i++)
    {
      int NI = get_node_index_of_coordinate_point(X_coords[i], Y_coords[i]);
      if (NI != NoDataValue)
	{
	  SourceNodes.push_back(NI);
	}
    }

  return SourceNodes;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Perform a downslope trace using D8 from a given point source (i,j).
//Overwrites input parameters to return a raster of the path, the length of the
//trace and the final pixel coordinates of the trace.
// SWDG 20/1/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::D8_Trace(int i, int j, LSDIndexRaster StreamNetwork, float& length, int& receiver_row, int& receiver_col, Array2D<int>& Path){

  float root_2 = 1.4142135623;

  Array2D<int> stnet = StreamNetwork.get_RasterData();

  length = 0;

  int node;

  int reciever_node = retrieve_node_from_row_and_column(i, j);
  receiver_row = i;
  receiver_col = j;

  Path[receiver_row][receiver_col] = 1;

  while (StreamNetwork.get_data_element(receiver_row, receiver_col) == NoDataValue){  // need to do edge checking

    retrieve_receiver_information(reciever_node, node, receiver_row, receiver_col);

    Path[receiver_row][receiver_col] = 1;

    //update length
    if (retrieve_flow_length_code_of_node(reciever_node) == 1){ length += DataResolution; }
    else if (retrieve_flow_length_code_of_node(reciever_node) == 2){ length += (DataResolution * root_2); }

    reciever_node = node;

  }

}

// Move the location of the channel head downslope by a user defined distance.
// Returns A vector of node indexes pointing to the moved heads.
// SWDG 27/11/15

void LSDFlowInfo::MoveChannelHeadDown(vector<int> Sources, float MoveDist, vector<int>& DownslopeSources, vector<int>& FinalHeads){

  float root_2 = 1.4142135623;

  float length;

  int receiver_row;
  int receiver_col;
  int node;

  //for loop goes here to iterate over the Sources vector
  for(int q=0; q<int(Sources.size());++q){
    length = 0;

    int reciever_node = Sources[q];

    while (length < MoveDist){

      retrieve_receiver_information(reciever_node, node, receiver_row, receiver_col);

      //update length
      if (retrieve_flow_length_code_of_node(reciever_node) == 1){ length += DataResolution; }
      else if (retrieve_flow_length_code_of_node(reciever_node) == 2){ length += (DataResolution * root_2); }
      else if (retrieve_flow_length_code_of_node(reciever_node) == 0){break;}

      reciever_node = node;


    }

    DownslopeSources.push_back(reciever_node);
    FinalHeads.push_back(Sources[q]);

  }

  //end of for loop

}

// Move the location of the channel head upslope by a user defined distance.
// Returns A vector of node indexes pointing to the moved heads.
// SWDG 27/11/15
void LSDFlowInfo::MoveChannelHeadUp(vector<int> Sources, float MoveDist, LSDRaster DEM, vector<int>& UpslopeSources, vector<int>& FinalHeads){

  float root_2 = 1.4142135623;

  float length;

  Array2D<float> Elevation = DEM.get_RasterData();

  int new_node;

  int new_i;
  int new_j;

  //for loop goes here to iterate over the Sources vector
  for(int q=0; q<int(Sources.size());++q){
    length = 0;

    int i;
    int j;

    retrieve_current_row_and_col(Sources[q], i, j);

    //test for channel heads at edges
    if (i == 0 || i == NRows - 1 || j == 0 || j == NCols - 1){
      cout << "Hit an edge, skipping" << endl;        
    }

    else{

    while (length < MoveDist){

      float currentElevation;
      int Direction; //1 is cardinal 2 is diagonal

      //find the neighbour with the maximum Elevation

      currentElevation = Elevation[i][j];

      //top left
      if(currentElevation < Elevation[i-1][j-1]){
        currentElevation = Elevation[i-1][j-1];
        new_i = i-1;
        new_j = j-1;
        Direction = 2;
      }
      //top
      if(currentElevation < Elevation[i][j-1]){
        currentElevation = Elevation[i][j-1];
        new_i = i;
        new_j = j-1;
        Direction = 1;
      }
      //top right
      if(currentElevation < Elevation[i+1][j-1]){
        currentElevation = Elevation[i+1][j-1];
        new_i = i+1;
        new_j = j-1;
        Direction = 2;
      }
      //right
      if(currentElevation < Elevation[i+1][j]){
        currentElevation = Elevation[i+1][j];
        new_i = i+1;
        new_j = j;
        Direction = 1;
      }
      //botom right
      if(currentElevation < Elevation[i+1][j+1]){
        currentElevation = Elevation[i+1][j+1];
        new_i = i+1;
        new_j = j+1;
        Direction = 2;
      }
      //bottom
      if(currentElevation < Elevation[i][j+1]){
        currentElevation = Elevation[i][j+1];
        new_i = i;
        new_j = j+1;
        Direction = 1;
      }
      //bottom left
      if(currentElevation < Elevation[i-1][j+1]){
        currentElevation = Elevation[i-1][j+1];
        new_i = i-1;
        new_j = j+1;
        Direction = 2;
      }
      //left
      if(currentElevation < Elevation[i-1][j]){
        currentElevation = Elevation[i-1][j];
        new_i = i-1;
        new_j = j;
        Direction = 1;
      }
      //test that we have not hit the ridge
      //this will exit the loop and add the final visited
      //node to the upper soureces vector
      if (currentElevation == Elevation[i][j]){
        cout << "Warning, unable to move channel head " << Sources[q] << " up by user defined distance" << endl;
        break;
      }

      //update length
      if (Direction == 1){ length += DataResolution; }
      else if (Direction == 2){ length += (DataResolution * root_2); }

      i = new_i;
      j = new_j;

    }
}
    new_node = retrieve_node_from_row_and_column(i,j);

    UpslopeSources.push_back(new_node);
    FinalHeads.push_back(Sources[q]);

  } //end of for loop

}

void LSDFlowInfo::HilltopFlowRoutingOriginal(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope, LSDRaster Aspect, LSDIndexRaster StreamNetwork)
{                                                         
  //Declare parameters
  int i,j,a,b;
  //double X,Y;
  float dem_res = DataResolution;
  //float mean_slope, relief;
  int ndv = NoDataValue;
  double slope_total, length, d;
  int flag, count;
  double PI = 3.14159265;
  double degs, degs_old, degs_new, theta;
  double s_local, s_edge;
  double xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;
  double xmin = XMinimum;
  double ymin = YMinimum;
  double ymax = ymin + NRows*dem_res;
  //double xmax = xmin + NCols*dem_res;

  //Declare Arrays
  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols);
  Array2D<float> blank(NRows,NCols,NoDataValue);

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  //for creating file names for hillslope profiles
  string file_part_1 = "prof_";
  string file_part_2;
  string file_part_3 = ".txt";
  string filename;
  char filename_c[256];

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*dem_res - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = xmin + j*dem_res + 0.5;
  }

  int ht_count = 0;

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    for (j=1; j<NCols-1; ++j) {
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != ndv) {
				
	//reset slope, length, hillslope flag and pixel counter
	slope_total = 0;
	length = 0;
	flag = true;
	count = 1;
				
	//copt blank raster to map hillslope trace
	path = blank.copy();

	//update hilltop counter
	++ht_count;
				
	//get aspect in radians
	degs = aspect[i][j];
	theta = (M_PI/180.)*((-1*degs)+90.);
				
	//setup indices
	a = i;
	b = j;
	path[a][b] = 1;

	//add first pixel to easting northing vectors
	east_vec[0] = easting[b];
	north_vec[0] = northing[a];

	//get local slope
	s_local = slope[a][b];
				
	//test direction, calculate outlet coordinates and update indicies
	// easterly, dir == 1
	if (degs >= 45 && degs < 135) {
	  //find location where trace exits the cell and distance
	  xo = 1., yo = (1.+tan(theta))/2.;
	  d = abs(1./(2.*cos(theta)));
	  //transmit to next cell over
	  xi = 0., yi = yo;
	  dir = 1;
	  //add to vector
	  east_vec[count] = easting[b] + 0.5*dem_res;
	  north_vec[count] = northing[a] + yo - 0.5*dem_res;
	  //increment indices
	  ++b;
	  //check we're not right in the corner!
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//southerly
	else if (degs >= 135 && degs < 225) {
	  //find location where trace exits the cell and distance
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  //transmit to next cell over
	  xi = xo, yi = 1;
	  dir = 2;
	  //add to vector
	  east_vec[count] = easting[b] + xo - 0.5*dem_res;
	  north_vec[count] = northing[a] - 0.5*dem_res;
	  //increment indices
	  ++a;
	  //check we're not right in the corner!
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	// westerly
	else if (degs >= 225 && degs < 315) {
	  //find location where trace exits the cell and distance
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  //transmit to next cell over
	  xi = 1,	yi = yo;
	  dir = 3;
	  //add to vector
	  east_vec[count] = easting[b] -0.5*dem_res;
	  north_vec[count] = northing[a] + yo - 0.5*dem_res;
	  //increment indices
	  --b;
	  //check we're not right in the corner!
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//northerly
	else if (degs >= 315 || degs < 45) {
	  //find location where trace exits the cell and distance
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  //transmit to next cell over
	  xi = xo, yi = 0;
	  dir = 4;
	  //add to vector
	  east_vec[count] = easting[b] + xo - 0.5*dem_res;
	  north_vec[count] = northing[a] + 0.5*dem_res;
	  //increment indices
	  --a;
	  //check we're not right in the corner!
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	else {
	  cout << "FATAL ERROR, Kinematic routing algorithm encountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

	//collect slopes and totals weighted by path length
	//slope_total += s_local*d;
	//length += d;

	s_local = slope[a][b];


	//continue trace until a stream node is encountered
	while (flag == true) {

	  path[a][b] = 1;

	  degs_old = degs;
	  degs_new = aspect[a][b];
	  theta = (M_PI/180.)*((-1*degs_new)+90.);
	  ++ count;

	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      ||	(dir == 2 && degs_new > 90 && degs_new < 270)
	      ||	(dir == 3 && degs_new > 180 && degs_new < 360)
	      ||	((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

	    //DO NORMAL FLOW PATH
	    //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	    temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
	    temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
	    temp_yo2 = yi - xi*tan(theta);			// xo = 0
	    temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0

	    // can't outlet at same point as inlet
	    if (dir == 1) temp_yo2 = -1;
	    else if (dir == 2) temp_xo1 = -1;
	    else if (dir == 3) temp_yo1 = -1;
	    else if (dir == 4) temp_xo2 = -1;

	    s_local = slope[a][b];

	    if (temp_yo1 <= 1 && temp_yo1 > 0) {
	      xo = 1, yo = temp_yo1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 0, yi = yo,
		dir = 1;
	      east_vec[count] = easting[b] + 0.5*dem_res;
	      north_vec[count] = northing[a] + yo - 0.5*dem_res;
	      ++b;
	      if (xi== 0 && yi == 0) yi = 0.00001;
	      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
	    }
	    else if (temp_xo2 <= 1 && temp_xo2 > 0) {
	      xo = temp_xo2, yo = 0;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 1,
		dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*dem_res;
	      north_vec[count] = northing[a] - 0.5*dem_res;
	      ++a;
	      if (xi== 0 && yi == 1) xi = 0.00001;
	      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
	    }
	    else if (temp_yo2 <= 1 && temp_yo2 > 0) {
	      xo = 0, yo = temp_yo2;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 1, yi = yo,
		dir = 3;
	      east_vec[count] = easting[b] -0.5*dem_res;
	      north_vec[count] = northing[a] + yo - 0.5*dem_res;
	      --b;
	      if (xi== 1 && yi == 0) yi = 0.00001;
	      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
	    }

	    else if (temp_xo1 <= 1 && temp_xo1 > 0) {
	      xo = temp_xo1, yo = 1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 0,
		dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*dem_res;
	      north_vec[count] = northing[a] + 0.5*dem_res;
	      --a;
							
	      if (xi == 0 && yi == 0) xi = 0.00001;
	      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
	    }
	    slope_total += s_local*d;
	  }

	  else {

	    // ROUTE ALONG EDGES
	    if (dir	== 1) {
	      if 	(degs_old <= 90 || degs_new >= 270) {
		xo = 0.00001, yo = 1;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*dem_res;
		north_vec[count] = northing[a] + 0.5*dem_res;
		--a;
	      }
	      else if (degs_old > 90 && degs_new < 270) {
		xo = 0.00001, yo = 0;
		s_edge = abs(s_local*sin((PI/2)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*dem_res;
		north_vec[count] = northing[a] - 0.5*dem_res;
		++a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;
		exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 2) {
	      if 	(degs_old <= 180 && degs_new >= 0) {
		xo = 1, yo = 1-0.00001;
		s_edge = abs(s_local*sin((2/PI)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*dem_res;
		north_vec[count] = northing[a] + yo - 0.5*dem_res;
		++b;
	      }
	      else if (degs_old > 180 && degs_new < 360) {
		xo = 0, yo = 1-0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*dem_res;
		north_vec[count] = northing[a] + yo - 0.5*dem_res;
		--b;

	      }
	      else {
		cout << "Flow unable to route E or W" << endl;
		exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 3) {
	      if 	(degs_old <= 270 && degs_new >= 90) {
		xo = 1-0.00001, yo = 0;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*dem_res;
		north_vec[count] = northing[a] - 0.5*dem_res;
		++a;
	      }
	      else if (degs_old > 270 || degs_new < 90) {
		xo = 1-0.00001, yo = 1;
		s_edge = abs(s_local*sin((2/PI) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1- yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*dem_res;
		north_vec[count] = northing[a] + 0.5*dem_res;
		--a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;
		exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 4) {
	      if 	(degs_old <= 360 && degs_new >= 180) {
		xo = 0, yo = 0.00001;
		s_edge = abs(s_local*sin((PI/2) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*dem_res;
		north_vec[count] = northing[a] + yo - 0.5*dem_res;
		--b;
	      }
	      else if (degs_old > 0 && degs_new < 180) {
		xo = 1, yo = 0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*dem_res;
		north_vec[count] = northing[a] + yo - 0.5*dem_res;
		++b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl;
		exit(EXIT_FAILURE);
	      }
	    }
	    slope_total += s_edge*d;

	  }
	  length += d;
	  degs = degs_new;

	  //cout << "[a][b]: " << a << " " << b << endl;

	  if (a <= 0 || b <= 0 ||	a >= NRows-1 || b >= NCols-1) flag = false;
	  else if (stnet[a][b] != ndv || path[a][b] == 1) flag = false;
	}

	//if trace finished at a stream, print hillslope info.
	if (a <= 0 || b <= 0 ||	a >= NRows-1 || b >= NCols-1) continue;
	else
	  {
	    if (path[a][b] == 1) 
	      {
		cout << "Didn't make it to a channel!" << endl;
	      }

	    //					// PRINT TO FILE Cht Sbar Relief Lh
	    //					X = xmin + j*dem_res;
	    //					Y = ymin + (NRows-i)*dem_res;
	    //					relief = zeta[i][j] - zeta[a][b];
	    //					length = length*dem_res;
	    //					mean_slope = slope_total/(length/dem_res);

	    //					ofs << X << " " << Y << " " << seg[i][j] << " "
	    //						<< cht[i][j] << " " << mean_slope << " "
	    //						<< relief << " " << length << " " << "/n"; //area[a][b] << "\n";

	    //PRINT FILE OF PATH NODES FOR EACH HILLSLOPE VECTOR
	    stringstream s;
	    s << ht_count;
	    file_part_2 = s.str();
	    filename = file_part_1;
	    filename.append(file_part_2), filename.append(file_part_3);
	    strcpy(filename_c,filename.c_str());

	    ofstream prof_out;
	    prof_out.open(filename_c);
	    prof_out << "Easting " << "Northing" << endl;
	    prof_out.precision(10);

	    for (int c=0;c<count;++c) {
	      prof_out 	<< east_vec[c] << " "
			<< north_vec[c] << endl;
	    }
	    prof_out.close();
	  }
      }
      //return condition for debugging purposes only
      //if (ht_count > 50) return;
    }
  }
}  

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
                                                         LSDIndexRaster StreamNetwork, LSDRaster Aspect, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
                                                         bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
                                                         vector<int> Target_Basin_Vector){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);
  
  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope,DivergentCount\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {
      
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {
        
	length = 0;
	flag = true;
	count = 1;
	path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        PlanarCountFlag = 0;
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.
        
        E_Star = 0;
        R_Star = 0;
        EucDist = 0;
        
	++ht_count;

	degs = aspect[i][j];
	theta = rads[i][j];
	a = i;
	b = j;
	path[a][b] += 1;
	east_vec[0] = easting[b];
	north_vec[0] = northing[a];
	s_local = slope[a][b];

	//test direction, calculate outlet coordinates and update indicies
	// easterly
	if (degs >= 45 && degs < 135) {
	  //cout << "\neasterly" << endl;
	  xo = 1, yo = (1+tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 0, yi = yo;
	  dir = 1;
	  east_vec[count] = easting[b] + 0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  ++b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//southerly
	else if (degs >= 135 && degs < 225) {
	  //cout << "\nsoutherly" << endl;
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  xi = xo, yi = 1;
	  dir = 2;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] - 0.5*DataResolution;
	  ++a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	// westerly
	else if (degs >= 225 && degs < 315) {
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 1,	yi = yo;
	  dir = 3;
	  east_vec[count] = easting[b] -0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  --b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//northerly
	else if (degs >= 315 || degs < 45) {
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  xi = xo, yi = 0;
	  dir = 4;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] + 0.5*DataResolution;
	  --a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	else {
	  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

	//collect slopes and totals weighted by path length
	length += d;
	s_local = slope[a][b];

	//continue trace until a stream node is encountered
	while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the	edge of the DEM tile.

	  path[a][b] += 1;

	  degs_new = aspect[a][b];
	  theta = rads[a][b];
          ++count;

	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      || (dir == 2 && degs_new > 90 && degs_new < 270)
	      || (dir == 3 && degs_new > 180 && degs_new < 360)
	      || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

	    //DO NORMAL FLOW PATH
	    //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	    temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
	    temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
	    temp_yo2 = yi - xi*tan(theta);			// xo = 0
	    temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0

	    // can't outlet at same point as inlet
	    if (dir == 1) temp_yo2 = -1;
	    else if (dir == 2) temp_xo1 = -1;
	    else if (dir == 3) temp_yo1 = -1;
	    else if (dir == 4) temp_xo2 = -1;

	    s_local = slope[a][b];           

	    if (temp_yo1 <= 1 && temp_yo1 > 0) {              
	      xo = 1, yo = temp_yo1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 0, yi = yo,
		dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	      if (xi== 0 && yi == 0) yi = 0.00001;
	      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
	    }
	    else if (temp_xo2 <= 1 && temp_xo2 > 0) {
	      xo = temp_xo2, yo = 0;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 1,
		dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	      if (xi== 0 && yi == 1) xi = 0.00001;
	      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
	    }
	    else if (temp_yo2 <= 1 && temp_yo2 > 0) {
	      xo = 0, yo = temp_yo2;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 1, yi = yo,
		dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	      if (xi== 1 && yi == 0) yi = 0.00001;
	      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
	    }

	    else if (temp_xo1 <= 1 && temp_xo1 > 0) {
	      xo = temp_xo1, yo = 1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 0,
		dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	      if (xi == 0 && yi == 0) xi = 0.00001;
	      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
	    }           
					
          }          
	  else {

	    // ROUTE ALONG EDGES						
	    if (dir	== 1) {
	      if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north                   
		xo = 0.00001, yo = 1;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else if (degs_new > 90 && degs_new < 270) {  //secondary component is south							  
		xo = 0.00001, yo = 0;
		s_edge = abs(s_local*sin((PI/2)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }						 
	      else {
		cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 2) {
	      if 	(degs_new >= 0 && degs_new <= 180) { //secondary component is East							  
		xo = 1, yo = 1-0.00001;
		s_edge = abs(s_local*sin((2/PI)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West						
		xo = 0, yo = 1-0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 3) {
	      if 	(degs_new >= 90 && degs_new <= 270) {  //secondary component is South
		xo = 1-0.00001, yo = 0;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }
	      else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
		xo = 1-0.00001, yo = 1;
		s_edge = abs(s_local*sin((2/PI) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1- yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 4) {
	      if 	(degs_new >= 180 && degs_new <= 360) { //secondary component is West
		xo = 0, yo = 0.00001;
		s_edge = abs(s_local*sin((PI/2) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else if (degs_new >= 0 && degs_new < 180) { //secondary component is East				  
		xo = 1, yo = 0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
						
	  }

	  if (path[a][b] < 1){  // only update length on 'first slosh'
	    length += d;
          }
          else if (path[a][b] >= 3){ //update the skip trace flag so we can categorise each trace          
            skip_trace = true;          
          } 

          degs = degs_new;

	  // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
	  // The larger the counter the more convergent or divergent the trace is
	  if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)){
	    ++DivergentCountFlag;
	  }
	  else {
	    ++PlanarCountFlag;
	  }
        
	  if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
	}

        if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 ){
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
	  {
	    //if trace finished at a stream, print hillslope info.
	    if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
	      {
		path[a][b] = 1;

		++s_count;

		X = XMinimum + j*DataResolution;
		Y = YMinimum - (NRows-i)*DataResolution;
		relief = zeta[i][j] - zeta[a][b];
		mean_slope = relief/(length * DataResolution);

		// update arrays with the current metrics
		RoutedHilltops[i][j] = 1;
		HillslopeLength_Array[i][j] = (length * DataResolution);
		Slope_Array[i][j] = mean_slope;
		Relief_Array[i][j] = relief;
          
		//calculate an E* and R* Value assuming S_c of 0.8
		E_Star = (2.0 * abs(hilltops[i][j])*(length*DataResolution))/0.8;
		R_Star = relief/((length*DataResolution)*0.8);
            
		//calulate the Euclidean distance between the start and end points of the trace
		EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;

          
		if (relief > 0){
		  ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "\n";
		}
		else {
		  ++neg_count;
		}
	      }
	    else{  //unable to route using aspects
	      //this will encompass skipped traces
	      ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
	      ++ns_count;
	    }
	  }

	//This block checks the various path printing options and writes the data out accordingly
	if (print_paths_switch == true){
	  if (ht_count % thinning == 0){
	    if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.
          
	      //create stringstream object to create filename
	      ofstream pathwriter;
          
	      //create the output filename from the user supplied path
	      stringstream ss_path;
	      ss_path << trace_path << i << "_" << j << "_trace.txt";
          
	      pathwriter.open(ss_path.str().c_str());
          
	      if(pathwriter.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
		exit(EXIT_FAILURE);
	      }
          
	      for (int v = 0; v < count+1; ++v){
		if (basin_filter_switch == false){
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
		}
		else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
		}
	      }
	      pathwriter.close();
	    }
	  }
	}
	// End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;
  
  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN 
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting_RAW(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
							     LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
							     bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
							     vector<int> Target_Basin_Vector){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag = 0;
  float PI = 3.14159265;
  float degs, degs_old, degs_new, theta;
  float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;
  
  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);
  
  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData_RAW.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope,DivergentCount\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {
      
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {
        
	length = 0;
	flag = true;
	count = 1;
	path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

	++ht_count;

	degs = aspect[i][j];
	theta = rads[i][j];
	a = i;
	b = j;
	path[a][b] += 1;
	east_vec[0] = easting[b];
	north_vec[0] = northing[a];
	s_local = slope[a][b];

	//test direction, calculate outlet coordinates and update indicies
	// easterly
	if (degs >= 45 && degs < 135) {
	  //cout << "\neasterly" << endl;
	  xo = 1, yo = (1+tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 0, yi = yo;
	  dir = 1;
	  east_vec[count] = easting[b] + 0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  ++b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//southerly
	else if (degs >= 135 && degs < 225) {
	  //cout << "\nsoutherly" << endl;
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  xi = xo, yi = 1;
	  dir = 2;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] - 0.5*DataResolution;
	  ++a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	// westerly
	else if (degs >= 225 && degs < 315) {
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 1,	yi = yo;
	  dir = 3;
	  east_vec[count] = easting[b] -0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  --b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//northerly
	else if (degs >= 315 || degs < 45) {
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  xi = xo, yi = 0;
	  dir = 4;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] + 0.5*DataResolution;
	  --a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	else {
	  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

	//collect slopes and totals weighted by path length
	length += d;
	s_local = slope[a][b];

	//continue trace until a stream node is encountered
	while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the	edge of the DEM tile.
	  int a_2 = a;
          int b_2 = b;

	  path[a][b] += 1;

	  degs_old = degs;
	  degs_new = aspect[a][b];
	  theta = rads[a][b];
          ++count;

	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      || (dir == 2 && degs_new > 90 && degs_new < 270)
	      || (dir == 3 && degs_new > 180 && degs_new < 360)
	      || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

	    //DO NORMAL FLOW PATH
	    //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	    temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
	    temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
	    temp_yo2 = yi - xi*tan(theta);			// xo = 0
	    temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0

	    // can't outlet at same point as inlet
	    if (dir == 1) temp_yo2 = -1;
	    else if (dir == 2) temp_xo1 = -1;
	    else if (dir == 3) temp_yo1 = -1;
	    else if (dir == 4) temp_xo2 = -1;

	    s_local = slope[a][b];           

	    if (temp_yo1 <= 1 && temp_yo1 > 0) {              
	      xo = 1, yo = temp_yo1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 0, yi = yo,
		dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	      if (xi== 0 && yi == 0) yi = 0.00001;
	      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
	    }
	    else if (temp_xo2 <= 1 && temp_xo2 > 0) {
	      xo = temp_xo2, yo = 0;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 1,
		dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	      if (xi== 0 && yi == 1) xi = 0.00001;
	      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
	    }
	    else if (temp_yo2 <= 1 && temp_yo2 > 0) {
	      xo = 0, yo = temp_yo2;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 1, yi = yo,
		dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	      if (xi== 1 && yi == 0) yi = 0.00001;
	      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
	    }

	    else if (temp_xo1 <= 1 && temp_xo1 > 0) {
	      xo = temp_xo1, yo = 1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 0,
		dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	      if (xi == 0 && yi == 0) xi = 0.00001;
	      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
	    }           
					
          }          
	  else {

	    // ROUTE ALONG EDGES						
	    if (dir == 1) {
	      if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north                   
		xo = 0.00001, yo = 1;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else if (degs_new > 90 && degs_new < 270) {  //secondary component is south							  
		xo = 0.00001, yo = 0;
		s_edge = abs(s_local*sin((PI/2)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }						 
	      else {
		cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 2) {
	      if 	(degs_new >= 0 && degs_new <= 180) { //secondary component is East							  
		xo = 1, yo = 1-0.00001;
		s_edge = abs(s_local*sin((2/PI)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West						
		xo = 0, yo = 1-0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 3) {
	      if 	(degs_new >= 90 && degs_new <= 270) {  //secondary component is South
		xo = 1-0.00001, yo = 0;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }
	      else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
		xo = 1-0.00001, yo = 1;
		s_edge = abs(s_local*sin((2/PI) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1- yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 4) {
	      if 	(degs_new >= 180 && degs_new <= 360) { //secondary component is West
		xo = 0, yo = 0.00001;
		s_edge = abs(s_local*sin((PI/2) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else if (degs_new >= 0 && degs_new < 180) { //secondary component is East							  
		xo = 1, yo = 0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
						
	  }

	  if (path[a][b] < 1){  // only update length on 'first slosh'
	    length += d;
          }

          degs = degs_new;

          if(zeta[a][b] - zeta[a_2][b_2] > 0){
            
            length -= d;    //remove uphill length from trace

            a = a_2;
            b = b_2;

            //restart trace
            degs = aspect[a][b];
	    theta = rads[a][b];
	    path[a][b] += 1;
	    s_local = slope[a][b];

            length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

	    //test direction, calculate outlet coordinates and update indicies
	    // easterly
	    if (degs >= 45 && degs < 135) {
	      xo = 1, yo = (1+tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 0, yi = yo;
	      dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	    }
	    //southerly
	    else if (degs >= 135 && degs < 225) {
	      xo = (1-(1/tan(theta)))/2, yo = 0;
	      d = abs(1/(2*cos((PI/2)-theta)));
	      xi = xo, yi = 1;
	      dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	    }
	    // westerly
	    else if (degs >= 225 && degs < 315) {
	      xo = 0, yo = (1-tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 1,	yi = yo;
	      dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	    }
	    //northerly
	    else if (degs >= 315 || degs < 45) {
	      xo = (1+(1/tan(theta)))/2, yo = 1;
	      d = abs(1/(2*cos((PI/2) - theta)));
	      xi = xo, yi = 0;
	      dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	    }
	    else {
	      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	      exit(EXIT_FAILURE);
	    }

	    //collect slopes and totals weighted by path length

	    length += d;
	    s_local = slope[a][b];

          }

	  if (path[a][b] >= 1){  //self intersect/'slosh'
                 
            degs = aspect[a][b];
	    theta = rads[a][b];
	    path[a][b] += 1;
	    s_local = slope[a][b];

            a_2 = a;
            b_2 = b;

	    length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

	    //test direction, calculate outlet coordinates and update indicies
	    // easterly
	    if (degs >= 45 && degs < 135) {
	      xo = 1, yo = (1+tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 0, yi = yo;
	      dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	    }
	    //southerly
	    else if (degs >= 135 && degs < 225) {
	      xo = (1-(1/tan(theta)))/2, yo = 0;
	      d = abs(1/(2*cos((PI/2)-theta)));
	      xi = xo, yi = 1;
	      dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	    }
	    // westerly
	    else if (degs >= 225 && degs < 315) {
	      xo = 0, yo = (1-tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 1,	yi = yo;
	      dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	    }
	    //northerly
	    else if (degs >= 315 || degs < 45) {
	      xo = (1+(1/tan(theta)))/2, yo = 1;
	      d = abs(1/(2*cos((PI/2) - theta)));
	      xi = xo, yi = 0;
	      dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	    }
	    else {
	      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	      exit(EXIT_FAILURE);
	    }

	    //collect slopes and totals weighted by path length
	    length += d;
	    s_local = slope[a][b];                                                                                                                                                                                                                

	  }

	  // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
	  // The larger the counter the more convergent or divergent the trace is
	  if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)){
	    ++DivergentCountFlag;
	  }
	  else {
	    ++PlanarCountFlag;
	  }

	  if (path[a][b] >=3){ //update flag if a trace cannot complete, so that we can track errors.
	    skip_trace = true;
	  }

	  if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
	}

        if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 ){
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
	  {
	    //if trace finished at a stream, print hillslope info.
	    if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
	      {
		path[a][b] = 1;

		++s_count;

		X = XMinimum + j*DataResolution;
		Y = YMinimum - (NRows-i)*DataResolution;
		relief = zeta[i][j] - zeta[a][b];
		mean_slope = relief/(length * DataResolution);

		// update arrays with the current metrics
		RoutedHilltops[i][j] = 1;
		HillslopeLength_Array[i][j] = (length * DataResolution);
		Slope_Array[i][j] = mean_slope;
		Relief_Array[i][j] = relief;
          
		//calculate an E* and R* Value assuming S_c of 0.8
		E_Star = (2.0 * abs(hilltops[i][j])*(length*DataResolution))/0.8;
		R_Star = relief/((length*DataResolution)*0.8);
            
		//calulate the Euclidean distance between the start and end points of the trace
		EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;
          
		if (relief > 0){
		  ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "\n";
		}
		else {
		  ++neg_count;
		}
	      }
	    else{  //unable to route using aspects
	      //this will encompass the skipped traces
	      ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
	      ++ns_count;
	    }
	  }

	//This block checks the various path printing options and writes the data out accordingly
	if (print_paths_switch == true){
	  if (ht_count % thinning == 0){
	    if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop and has a valid trace, ie there is actually a trace to write to file.
          
	      //create stringstream object to create filename
	      ofstream pathwriter;
          
	      //create the output filename from the user supplied path
	      stringstream ss_path;
	      ss_path << trace_path << i << "_" << j << "_trace.txt";
          
	      pathwriter.open(ss_path.str().c_str());
          
	      if(pathwriter.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
		exit(EXIT_FAILURE);
	      }
          
	      for (int v = 0; v < count+1; ++v){
		if (basin_filter_switch == false){
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
		}
		else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << endl;
		}
	      }
	      pathwriter.close();
	    }
	  }
	}
	// End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;
  
  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN 
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// This version is used to generate elevation profiles of the traces. The elevation data is encoded within
// the trace files.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
//
// SWDG 25/3/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRouting_Profile(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
								 LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir, string Prefix, LSDIndexRaster Basins,
								 bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
								 vector<int> Target_Basin_Vector){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  int flag;
  int count = 0;
	
  float PI = 3.14159265;
  float degs, degs_old, degs_new, theta;
  float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);
  
  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);
  vector<float> ZetaList;
  vector<float> LengthList;

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData_RAW.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {
      
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {
        
	length = 0;
	flag = true;
	count = 1;
	path = blank.copy();
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

	++ht_count;

	degs = aspect[i][j];
	theta = rads[i][j];
	a = i;
	b = j;
	path[a][b] += 1;
	east_vec[0] = easting[b];
	north_vec[0] = northing[a];
	s_local = slope[a][b];
        ZetaList.clear();
        LengthList.clear();

	//test direction, calculate outlet coordinates and update indicies
	// easterly
	if (degs >= 45 && degs < 135) {
	  //cout << "\neasterly" << endl;
	  xo = 1, yo = (1+tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 0, yi = yo;
	  dir = 1;
	  east_vec[count] = easting[b] + 0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  ++b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//southerly
	else if (degs >= 135 && degs < 225) {
	  //cout << "\nsoutherly" << endl;
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  xi = xo, yi = 1;
	  dir = 2;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] - 0.5*DataResolution;
	  ++a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	// westerly
	else if (degs >= 225 && degs < 315) {
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 1,	yi = yo;
	  dir = 3;
	  east_vec[count] = easting[b] -0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  --b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//northerly
	else if (degs >= 315 || degs < 45) {
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  xi = xo, yi = 0;
	  dir = 4;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] + 0.5*DataResolution;
	  --a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	else {
	  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

	//collect slopes and totals weighted by path length
	length += d;
	s_local = slope[a][b];

	//continue trace until a stream node is encountered
	while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the	edge of the DEM tile.
	  int a_2 = a;
          int b_2 = b;

	  path[a][b] += 1;

	  degs_old = degs;
	  degs_new = aspect[a][b];
	  theta = rads[a][b];
          ++count;

	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      || (dir == 2 && degs_new > 90 && degs_new < 270)
	      || (dir == 3 && degs_new > 180 && degs_new < 360)
	      || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {

	    //DO NORMAL FLOW PATH
	    //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	    temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
	    temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
	    temp_yo2 = yi - xi*tan(theta);			// xo = 0
	    temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0

	    // can't outlet at same point as inlet
	    if (dir == 1) temp_yo2 = -1;
	    else if (dir == 2) temp_xo1 = -1;
	    else if (dir == 3) temp_yo1 = -1;
	    else if (dir == 4) temp_xo2 = -1;

	    s_local = slope[a][b];           

	    if (temp_yo1 <= 1 && temp_yo1 > 0) {              
	      xo = 1, yo = temp_yo1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 0, yi = yo,
		dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	      if (xi== 0 && yi == 0) yi = 0.00001;
	      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
	    }
	    else if (temp_xo2 <= 1 && temp_xo2 > 0) {
	      xo = temp_xo2, yo = 0;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 1,
		dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	      if (xi== 0 && yi == 1) xi = 0.00001;
	      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
	    }
	    else if (temp_yo2 <= 1 && temp_yo2 > 0) {
	      xo = 0, yo = temp_yo2;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 1, yi = yo,
		dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	      if (xi== 1 && yi == 0) yi = 0.00001;
	      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
	    }

	    else if (temp_xo1 <= 1 && temp_xo1 > 0) {
	      xo = temp_xo1, yo = 1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = xo, yi = 0,
		dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	      if (xi == 0 && yi == 0) xi = 0.00001;
	      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
	    }           
					
          }          
	  else {

	    // ROUTE ALONG EDGES						
	    if (dir	== 1) {
	      if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north                   
		xo = 0.00001, yo = 1;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else if (degs_new > 90 && degs_new < 270) {  //secondary component is south							  
		xo = 0.00001, yo = 0;
		s_edge = abs(s_local*sin((PI/2)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }						 
	      else {
		cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 2) {
	      if 	(degs_new >= 0 && degs_new <= 180) { //secondary component is East							  
		xo = 1, yo = 1-0.00001;
		s_edge = abs(s_local*sin((2/PI)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West						
		xo = 0, yo = 1-0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 3) {
	      if 	(degs_new >= 90 && degs_new <= 270) {  //secondary component is South
		xo = 1-0.00001, yo = 0;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }
	      else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
		xo = 1-0.00001, yo = 1;
		s_edge = abs(s_local*sin((2/PI) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1- yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 4) {
	      if 	(degs_new >= 180 && degs_new <= 360) { //secondary component is West
		xo = 0, yo = 0.00001;
		s_edge = abs(s_local*sin((PI/2) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else if (degs_new >= 0 && degs_new < 180) { //secondary component is East							  
		xo = 1, yo = 0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
                cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		skip_trace = true;
                //exit(EXIT_FAILURE);
	      }
	    }
						
	  }

	  if (path[a][b] < 1){  // only update length on 'first slosh'
	    length += d;
          }

          degs = degs_new;

          if(zeta[a][b] - zeta[a_2][b_2] > 0){
            
            length -= d;    //remove uphill length from trace

            a = a_2;
            b = b_2;

            //restart trace
            degs = aspect[a][b];
	    theta = rads[a][b];
	    path[a][b] += 1;
	    s_local = slope[a][b];

            length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

	    //test direction, calculate outlet coordinates and update indicies
	    // easterly
	    if (degs >= 45 && degs < 135) {
	      xo = 1, yo = (1+tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 0, yi = yo;
	      dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	    }
	    //southerly
	    else if (degs >= 135 && degs < 225) {
	      xo = (1-(1/tan(theta)))/2, yo = 0;
	      d = abs(1/(2*cos((PI/2)-theta)));
	      xi = xo, yi = 1;
	      dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	    }
	    // westerly
	    else if (degs >= 225 && degs < 315) {
	      xo = 0, yo = (1-tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 1,	yi = yo;
	      dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	    }
	    //northerly
	    else if (degs >= 315 || degs < 45) {
	      xo = (1+(1/tan(theta)))/2, yo = 1;
	      d = abs(1/(2*cos((PI/2) - theta)));
	      xi = xo, yi = 0;
	      dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	    }
	    else {
	      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	      exit(EXIT_FAILURE);
	    }

	    //collect slopes and totals weighted by path length

	    length += d;
	    s_local = slope[a][b];

          }

	  if (path[a][b] >= 1){  //self intersect/'slosh'
                 
            degs = aspect[a][b];
	    theta = rads[a][b];
	    path[a][b] += 1;
	    s_local = slope[a][b];

            a_2 = a;
            b_2 = b;

	    length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

	    //test direction, calculate outlet coordinates and update indicies
	    // easterly
	    if (degs >= 45 && degs < 135) {
	      xo = 1, yo = (1+tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 0, yi = yo;
	      dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	    }
	    //southerly
	    else if (degs >= 135 && degs < 225) {
	      xo = (1-(1/tan(theta)))/2, yo = 0;
	      d = abs(1/(2*cos((PI/2)-theta)));
	      xi = xo, yi = 1;
	      dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	    }
	    // westerly
	    else if (degs >= 225 && degs < 315) {
	      xo = 0, yo = (1-tan(theta))/2;
	      d = abs(1/(2*cos(theta)));
	      xi = 1,	yi = yo;
	      dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	    }
	    //northerly
	    else if (degs >= 315 || degs < 45) {
	      xo = (1+(1/tan(theta)))/2, yo = 1;
	      d = abs(1/(2*cos((PI/2) - theta)));
	      xi = xo, yi = 0;
	      dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	    }
	    else {
	      cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	      exit(EXIT_FAILURE);
	    }

	    //collect slopes and totals weighted by path length
	    length += d;
	    s_local = slope[a][b];                                                                                                                                                                                                                

	  }

	  if (path[a][b] >=3){ //update flag if a trace cannot complete, so that we can track errors.
	    skip_trace = true;
	  }

	  ZetaList.push_back(zeta[a][b]);
	  LengthList.push_back(length*DataResolution);

	  if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
	}

        if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 ){
          // avoid going out of bounds.

          // this is caused by having a hilltop on the first row or col away from the border
          // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
          // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*

          ++edge_count;

        }
        else
	  {
	    //if trace finished at a stream, print hillslope info.
	    if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
	      {
		path[a][b] = 1;

		++s_count;

		X = XMinimum + j*DataResolution;
		Y = YMinimum - (NRows-i)*DataResolution;
		relief = zeta[i][j] - zeta[a][b];
		mean_slope = relief/(length * DataResolution);

		// update arrays with the current metrics
		RoutedHilltops[i][j] = 1;
		HillslopeLength_Array[i][j] = (length * DataResolution);
		Slope_Array[i][j] = mean_slope;
		Relief_Array[i][j] = relief;
          
		if (relief > 0){
		  ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "\n";
		}
		else {
		  ++neg_count;
		}
	      }
	    else{  //unable to route using aspects
	      //this will encompass the skipped traces
	      ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
	      ++ns_count;
	    }
	  }

	//This block checks the various path printing options and writes the data out accordingly
	if (print_paths_switch == true){
	  if (ht_count % thinning == 0){
	    if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop and has a valid trace, ie there is actually a trace to write to file.
          
	      //create stringstream object to create filename
	      ofstream pathwriter;
          
	      //create the output filename from the user supplied path
	      stringstream ss_path;
	      ss_path << trace_path << i << "_" << j << "_trace.txt";
          
	      pathwriter.open(ss_path.str().c_str());
          
	      if(pathwriter.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
		exit(EXIT_FAILURE);
	      }
          
	      for (int v = 0; v < count+1; ++v){
		if (basin_filter_switch == false){
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << ZetaList[v] << " " << LengthList[v] << endl;
		}
		else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << ZetaList[v] << " " << LengthList[v] << endl;
		}
	      }
	      pathwriter.close();
	    }
	  }
	}
	// End of path printing logic
      }
    }   //for loop i,j
  }

  ofs.close();

  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);

  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;
  
  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function makes a mask of all the pixels that recieve flow (d8) 
// from a pixel that is either nodata or is on the boundary of the DEM
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDFlowInfo::find_cells_influenced_by_nodata(LSDIndexRaster& Bordered_mask,
							    LSDRaster& Topography)
{

  // set up the array
  Array2D<int> influenced_mask(NRows,NCols,int(NoDataValue));
  for(int row = 0; row <NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
	{
	  if(Topography.get_data_element(row,col) != NoDataValue)
	    { 
	      influenced_mask[row][col] = 0;
	    }
	}
    }
  
  
  int curr_node;
  int next_node;
  int next_row,next_col;
  
  // now loop through every node in the array
  for(int row = 0; row <NRows; row++)
    {
      for(int col = 0; col<NCols; col++)
	{
	  if(Topography.get_data_element(row,col) != NoDataValue)
	    {
	      // this node has data. 
	      // first see if it has already been tagged 
	      if(influenced_mask[row][col] != 1)
		{
        
		  //See if it is borderd by a NDV
		  if(Bordered_mask.get_data_element(row,col) == 1)
		    {
		      // it is bordered by nodata. Work your way down the node list
		      curr_node = retrieve_node_from_row_and_column(row, col);
		      next_node = ReceiverVector[curr_node];
          
		      influenced_mask[row][col] = 1;
		      retrieve_current_row_and_col(next_node, next_row, next_col);
            
		      //cout << "I am bordered by NDV, entering search loop" << endl;
		      //cout << "Row: " << row <<  " col: " << col << " node: " << curr_node 
		      //     << " receiver: " << next_node << " next infl mask: "
		      //     << influenced_mask[next_row][next_col] << endl;
            
		      // loop until you hit another influenced node or a baselevel node
		      while(next_node != curr_node && influenced_mask[next_row][next_col] != 1 )
			{
			  curr_node = next_node;
			  next_node = ReceiverVector[curr_node];
              
			  // the index here say next row and column but actually this is
			  // preserved from the previous loop so is the current node. 
			  influenced_mask[next_row][next_col] = 1;
              
			  // get the row and column of the receiver
			  retrieve_current_row_and_col(next_node, next_row, next_col); 
              
			  //cout << "Looping thought influence, next influenced is: " 
			  //     << influenced_mask[next_row][next_col] << endl;
                           
			}                     
		    }       
		}
	    }
	}
    }

  // now write the mask as an LSDIndexRaster
  LSDIndexRaster Influence_by_NDV(NRows,NCols,XMinimum,YMinimum,
				  DataResolution,int(NoDataValue),influenced_mask,GeoReferencingStrings);
  return Influence_by_NDV;
  
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//----------------------------------------------------------------------------------------
// get_raster_values_for_nodes
//---------------------------------------------------------------------------------------- 
// This function gets the values from a raster corresponding to the given nodes.
vector<float> LSDFlowInfo::get_raster_values_for_nodes(LSDRaster& Raster, vector<int>& node_indices)
{
  int N_nodes = node_indices.size();
  vector<float> return_values(N_nodes,float(NoDataValue));
  int row=0;
  int col=0;
  for(int i = 0; i < N_nodes; ++i)
    {
      if(node_indices[i] == NoDataValue)
	{
	  return_values[i] = NoDataValue;
	}
      else
	{
	  retrieve_current_row_and_col(node_indices[i],row,col);
	  return_values[i] = Raster.get_data_element(row,col);
	}
    }
  return return_values;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This version performs a single trace from a specified node, and routes down until it
// reaches a channel pixel
//
// THIS VERSION OF THE CODE RETAINS THE FLOODING METHOD TO ALLOW TRACES TO BE USED
// ON RAW TOPOGRPAHY TO GET EVENT SCALE HILLSLOPE LENGTHS WITH NO SMOOTHING. IN 
// MOST CASES USE THE MAIN METHOD, TO ANALYSE SEDIMENT TRANSPORT OVER GEOMORPHIC TIME.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// There are 4 outputs:
// output_trace_coordinates - the output coordinates tracing the flow path
// output_trace_metrics - the metrics derived from the flow routing
//                        (i) X
//                        (i) Y
//                        (i) mean slope
//                        (i) hillslope relief
//                        (i) hillslope length
//                        (i) channel ID @ lower boundary
// output_channel_node -the nodeindex at the bounding channel
// skip_trace - a bool object that specifies whether this trace has routed
// successfully to the channel.
//
// SWDG (adapted by DTM) 23/3/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDFlowInfo::D_Inf_single_trace_to_channel(LSDRaster Elevation, int start_node, LSDIndexRaster StreamNetwork, LSDRaster D_inf_Flowdir,
						vector< vector<float> >& output_trace_coordinates, vector<float>& output_trace_metrics,
						int& output_channel_node, bool& skip_trace)
{

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  //   int flag;
  int count = 0;
  //   int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  float PI = 3.14159265;
  float degs, degs_old, degs_new, theta;
  //   float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  //   bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  //   int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;

  float ymax = YMinimum + NRows*DataResolution;

  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = D_inf_Flowdir.get_RasterData(); //aspect
  //   Array2D<float> slope = Slope.get_RasterData(); //slope
  
  Array2D<float> rads(NRows,NCols,NoDataValue);
  Array2D<float> path(NRows, NCols,0.0);
  Array2D<float> blank(NRows,NCols,0.0);
  
  int channel_node = int(NoDataValue);
  vector<float> trace_metrics;
  vector< vector<float> > trace_coordinates;
  vector<float> empty;
  trace_coordinates.push_back(empty);
  trace_coordinates.push_back(empty);

  int vec_size = 1000000;

  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  //calculate northing and easting
  for (i=0;i<NRows;++i)
    {
      northing[i] = ymax - i*DataResolution - 0.5;
    }
  for (j=0;j<NCols;++j)
    {
      easting[j] = XMinimum + j*DataResolution + 0.5;
    }
  // find node and trace downstream
      
  // ignore edge cells and non-hilltop cells
  // route initial node by aspect and get outlet coordinates
  int start_row, start_col;
  retrieve_current_row_and_col(start_node,start_row,start_col);
  bool flag = false;
  if (zeta[start_row][start_col] != NoDataValue)
    {
      length = 0;
      flag = true;
      count = 1;
      path = blank.copy();
      //     DivergentCountFlag = 0; //initialise count of divergent cells in trace
      skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.

      //     ++ht_count;

      degs = aspect[start_row][start_col];
      theta = BearingToRad(aspect[start_row][start_col]);
      a = start_row;
      b = start_col;
      path[a][b] += 1;
      east_vec[0] = easting[b];
      north_vec[0] = northing[a];
      //     s_local = slope[a][b];

      //test direction, calculate outlet coordinates and update indicies
      // easterly
      if (degs >= 45 && degs < 135)
	{
	  //cout << "\neasterly" << endl;
	  xo = 1, yo = (1+tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 0, yi = yo;
	  dir = 1;
	  east_vec[count] = easting[b] + 0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  ++b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
      //southerly
      else if (degs >= 135 && degs < 225)
	{
	  //cout << "\nsoutherly" << endl;
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  xi = xo, yi = 1;
	  dir = 2;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] - 0.5*DataResolution;
	  ++a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
      // westerly
      else if (degs >= 225 && degs < 315)
	{
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 1,	yi = yo;
	  dir = 3;
	  east_vec[count] = easting[b] -0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  --b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
      //northerly
      else if (degs >= 315 || degs < 45)
	{
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  xi = xo, yi = 0;
	  dir = 4;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] + 0.5*DataResolution;
	  --a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
      else
	{
	  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

      //collect slopes and totals weighted by path length
      length += d;
      //     s_local = slope[a][b];
      // place coordinates into output vector
      trace_coordinates[0].push_back(east_vec[count]);
      trace_coordinates[1].push_back(north_vec[count]);
      //continue trace until a stream node is encountered 
      while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1)   //added boudary checking to catch cells which flow off the edge of the DEM tile.
	{
	  int a_2 = a;
	  int b_2 = b;

	  path[a][b] += 1;

	  degs_old = degs;
	  degs_new = aspect[a][b];
	  theta = BearingToRad(aspect[a][b]);
	  ++count;
      
	  //       cout << "TEST1" << endl;
	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      || (dir == 2 && degs_new > 90 && degs_new < 270)
	      || (dir == 3 && degs_new > 180 && degs_new < 360)
	      || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90)))
	    {
        
	      //       cout << "TEST1a" << endl;
	      //DO NORMAL FLOW PATH
	      //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	      temp_yo1 = yi + (1-xi)*tan(theta);      // xo = 1
	      temp_xo1 = xi + (1-yi)*(1/tan(theta));  // yo = 1
	      temp_yo2 = yi - xi*tan(theta);          // xo = 0
	      temp_xo2 = xi - yi*(1/tan(theta));      // yo = 0

	      // can't outlet at same point as inlet
	      if (dir == 1) temp_yo2 = -1;
	      else if (dir == 2) temp_xo1 = -1;
	      else if (dir == 3) temp_yo1 = -1;
	      else if (dir == 4) temp_xo2 = -1;

	      //         s_local = slope[a][b];
	      if (temp_yo1 <= 1 && temp_yo1 > 0)
		{              
		  xo = 1, yo = temp_yo1;
		  d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		  xi = 0;
		  yi = yo;
		  dir = 1;
		  east_vec[count] = easting[b] + 0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  ++b;
		  if (xi== 0 && yi == 0) yi = 0.00001;
		  else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
		}
	      else if (temp_xo2 <= 1 && temp_xo2 > 0)
		{
		  xo = temp_xo2, yo = 0;
		  d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		  xi = xo;
		  yi = 1;
		  dir = 2;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] - 0.5*DataResolution;
		  ++a;
		  if (xi== 0 && yi == 1) xi = 0.00001;
		  else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
		}
	      else if (temp_yo2 <= 1 && temp_yo2 > 0)
		{
		  xo = 0, yo = temp_yo2;
		  d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		  xi = 1;
		  yi = yo;
		  dir = 3;
		  east_vec[count] = easting[b] -0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  --b;
		  if (xi== 1 && yi == 0) yi = 0.00001;
		  else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
		}
	      else if (temp_xo1 <= 1 && temp_xo1 > 0)
		{
		  xo = temp_xo1, yo = 1;
		  d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		  xi = xo;
		  yi = 0;
		  dir = 4;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] + 0.5*DataResolution;
		  --a;
		  if (xi == 0 && yi == 0) xi = 0.00001;
		  else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
		}
        
	      //       cout << "TEST1_end" << endl;
	    }
	  else
	    { 
	      // ROUTE ALONG EDGES
        
	      //       cout << "TEST-" << endl;
	      if (dir == 1)
		{ 
		  //       cout << "TEST2" << endl;
		  if (degs_new <= 90 || degs_new >= 270) //secondary compenent of flow is north                   
		    {
		      xo = 0.00001;
		      yo = 1;
		      //             s_edge = abs(s_local*sin(theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = xo;
		      yi = 1-yo;
		      dir = 4;
		      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		      north_vec[count] = northing[a] + 0.5*DataResolution;
		      --a;
		    }
		  else if (degs_new > 90 && degs_new < 270)  //secondary component is south
		    {
		      xo = 0.00001;
		      yo = 0;
		      //             s_edge = abs(s_local*sin((PI/2)-theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = xo;
		      yi = 1-yo;
		      dir = 2;
		      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		      north_vec[count] = northing[a] - 0.5*DataResolution;
		      ++a;
		    } 
		  else
		    {
		      cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
		      cout << "Trace skipped.\n" << endl; 
		      skip_trace = true;
		      //exit(EXIT_FAILURE);
		    }
		}
	      else if (dir == 2)
		{
		  //       cout << "TEST3" << endl;
		  if (degs_new >= 0 && degs_new <= 180) //secondary component is East
		    {
		      xo = 1, yo = 1-0.00001;
		      //             s_edge = abs(s_local*sin((2/PI)-theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = 1-xo, yi = yo;
		      dir = 1;
		      east_vec[count] = easting[b] + 0.5*DataResolution;
		      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		      ++b;
		    }
		  else if (degs_new > 180 && degs_new <= 360)  //secondary component is West
		    {
		      xo = 0, yo = 1-0.00001;
		      //             s_edge = abs(s_local*sin(theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = 1-xo, yi = yo;
		      dir = 3;
		      east_vec[count] = easting[b] -0.5*DataResolution;
		      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		      --b;
		    }
		  else
		    {
		      cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		      cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		      skip_trace = true;
		      //exit(EXIT_FAILURE);
		    }
		}
	      else if (dir == 3)
		{
          
		  //           cout << "TEST4" << endl;
		  if(degs_new >= 90 && degs_new <= 270)   //secondary component is South
		    {
		      xo = 1-0.00001;
		      yo = 0;
		      //             s_edge = abs(s_local*sin(theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = xo;
		      yi = 1-yo;
		      dir = 2;
		      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		      north_vec[count] = northing[a] - 0.5*DataResolution;
		      ++a;
		    }
		  else if (degs_new > 270 || degs_new < 90)   //secondary component is North
		    {
		      xo = 1-0.00001, yo = 1;
		      //             s_edge = abs(s_local*sin((2/PI) - theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = xo;
		      yi = 1- yo;
		      dir = 4;
		      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		      north_vec[count] = northing[a] + 0.5*DataResolution;
		      --a;
		    }
		  else
		    {
		      cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
		      cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		      skip_trace = true;
		      //exit(EXIT_FAILURE);
		    }
		}
	      else if (dir == 4)
		{ 
		  //       cout << "TEST5" << endl;
		  if(degs_new >= 180 && degs_new <= 360) //secondary component is West
		    {
		      xo = 0, yo = 0.00001;
		      //             s_edge = abs(s_local*sin((PI/2) - theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = 1-xo;
		      yi = yo;
		      dir = 3;
		      east_vec[count] = easting[b] -0.5*DataResolution;
		      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		      --b;
		    }
		  else if (degs_new >= 0 && degs_new < 180) //secondary component is East
		    { 
		      xo = 1, yo = 0.00001;
		      //             s_edge = abs(s_local*sin(theta));
		      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		      xi = 1-xo, yi = yo;
		      dir = 1;
		      east_vec[count] = easting[b] + 0.5*DataResolution;
		      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		      ++b;
		    }
		  else
		    {
		      cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		      cout << "Trace skipped.\n" << endl; //something has gone very wrong...
		      skip_trace = true;
		      //exit(EXIT_FAILURE);
		    }
		}
        
	      //       cout << "TEST6" << endl;
	    }
      
	  if (path[a][b] < 1) length += d; // only update length on 'first slosh'

	  degs = degs_new;
      
	  if(zeta[a][b] - zeta[a_2][b_2] > 0)
	    {
	      length -= d;    //remove uphill length from trace

	      a = a_2;
	      b = b_2;
        
	      //         cout << "TEST7" << endl;
	      //restart trace
	      degs = aspect[a][b];
	      theta = BearingToRad(aspect[a][b]);
	      path[a][b] += 1;
	      //         s_local = slope[a][b];
	      length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace

	      //test direction, calculate outlet coordinates and update indices easterly
	      if (degs >= 45 && degs < 135)
		{
		  xo = 1, yo = (1+tan(theta))/2;
		  d = abs(1/(2*cos(theta)));
		  xi = 0, yi = yo;
		  dir = 1;
		  east_vec[count] = easting[b] + 0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  ++b;
		}
	      //southerly
	      else if (degs >= 135 && degs < 225) 
		{
		  xo = (1-(1/tan(theta)))/2, yo = 0;
		  d = abs(1/(2*cos((PI/2)-theta)));
		  xi = xo, yi = 1;
		  dir = 2;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] - 0.5*DataResolution;
		  ++a;
		}
	      // westerly
	      else if (degs >= 225 && degs < 315)
		{
		  xo = 0, yo = (1-tan(theta))/2;
		  d = abs(1/(2*cos(theta)));
		  xi = 1,	yi = yo;
		  dir = 3;
		  east_vec[count] = easting[b] -0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  --b;
		}
	      //northerly
	      else if (degs >= 315 || degs < 45)
		{
		  xo = (1+(1/tan(theta)))/2, yo = 1;
		  d = abs(1/(2*cos((PI/2) - theta)));
		  xi = xo, yi = 0;
		  dir = 4;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] + 0.5*DataResolution;
		  --a;
		}
	      else
		{
		  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
		  exit(EXIT_FAILURE);
		}
	      //collect slopes and totals weighted by path length
	      length += d;
	      //         s_local = slope[a][b];
	    }

	  if (path[a][b] >= 1)  //self intersect/'slosh'
	    {
	      degs = aspect[a][b];
	      theta = rads[a][b];
	      path[a][b] += 1;
	      //         s_local = slope[a][b];

	      a_2 = a;
	      b_2 = b;
        
	      length += sqrt((pow((xo-0.5),2) + pow((yo-0.5),2)));  //update length to cope with the 'jump' to the centre of the cell to restart the trace
        
	      //test direction, calculate outlet coordinates and update indices
	      // easterly
	      if (degs >= 45 && degs < 135)
		{
		  xo = 1, yo = (1+tan(theta))/2;
		  d = abs(1/(2*cos(theta)));
		  xi = 0, yi = yo;
		  dir = 1;
		  east_vec[count] = easting[b] + 0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  ++b;
		}
	      //southerly
	      else if (degs >= 135 && degs < 225)
		{
		  xo = (1-(1/tan(theta)))/2, yo = 0;
		  d = abs(1/(2*cos((PI/2)-theta)));
		  xi = xo, yi = 1;
		  dir = 2;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] - 0.5*DataResolution;
		  ++a;
		}
	      // westerly
	      else if (degs >= 225 && degs < 315)
		{
		  xo = 0, yo = (1-tan(theta))/2;
		  d = abs(1/(2*cos(theta)));
		  xi = 1,	yi = yo;
		  dir = 3;
		  east_vec[count] = easting[b] -0.5*DataResolution;
		  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		  --b;
		}
	      //northerly
	      else if (degs >= 315 || degs < 45)
		{
		  xo = (1+(1/tan(theta)))/2, yo = 1;
		  d = abs(1/(2*cos((PI/2) - theta)));
		  xi = xo, yi = 0;
		  dir = 4;
		  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		  north_vec[count] = northing[a] + 0.5*DataResolution;
		  --a;
		}
	      else
		{
		  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
		  exit(EXIT_FAILURE);
		}
	      //collect slopes and totals weighted by path length
	      length += d;
	      //         s_local = slope[a][b];
	    }
      
	  // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
	  // The larger the counter the more convergent or divergent the trace is
	  //       if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)) ++DivergentCountFlag;
	  if (path[a][b] >=3) skip_trace = true;//update flag if a trace cannot complete, so that we can track errors.
	  if (a == 0 || b == 0 || a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
      
	  // save trace coordinates for this iteration.
	  trace_coordinates[0].push_back(east_vec[count]);
	  trace_coordinates[1].push_back(north_vec[count]);
	}
    
      if (a == 0 || b == 0 || a == NRows-1 || b == NCols-1 )
	{
	  // avoid going out of bounds.
	  // this is caused by having a hilltop on the first row or col away from the border
	  // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
	  // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*
	  ++edge_count;
	}
      else
	{
	  //if trace finished at a stream, print hillslope info.
	  if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
	    {
	      path[a][b] = 1;
	      ++s_count;
	      X = XMinimum + j*DataResolution;
	      Y = YMinimum - (NRows-i)*DataResolution;
	      relief = zeta[start_row][start_col] - zeta[a][b];
	      mean_slope = relief/(length * DataResolution);

	      trace_metrics.push_back(X);
	      trace_metrics.push_back(Y);
	      trace_metrics.push_back(float(start_node));
	      trace_metrics.push_back(mean_slope);
	      trace_metrics.push_back(relief);
	      trace_metrics.push_back(length*DataResolution);
        
	      if (stnet[a][b] != NoDataValue)
		{
		  channel_node = retrieve_node_from_row_and_column(a,b);
		  trace_metrics.push_back(float(channel_node));
		}
	      // find nearest channel pixel within 1m buffer - if more than one, choose furthest downstream
	      else
		{
		  float min_elev=NoDataValue;
		  if (stnet[a-1][b-1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a-1][b-1] < min_elev)
			{
			  min_elev = zeta[a-1][b-1]; 
			  channel_node = retrieve_node_from_row_and_column(a-1,b-1);
			}
		    }
		  else if (stnet[a-1][b] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a-1][b] < min_elev)
			{
			  min_elev = zeta[a-1][b]; 
			  channel_node = retrieve_node_from_row_and_column(a-1,b);
			}
		    }
		  else if (stnet[a-1][b+1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a-1][b+1] < min_elev)
			{
			  min_elev = zeta[a-1][b+1]; 
			  channel_node = retrieve_node_from_row_and_column(a-1,b+1);
			}
		    }
		  else if (stnet[a][b-1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a][b-1] < min_elev)
			{
			  min_elev = zeta[a][b-1]; 
			  channel_node = retrieve_node_from_row_and_column(a,b-1);
			}
		    }
		  else if (stnet[a][b+1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a][b+1] < min_elev)
			{
			  min_elev = zeta[a][b+1]; 
			  channel_node = retrieve_node_from_row_and_column(a,b+1);
			}
		    }
		  else if (stnet[a+1][b-1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a+1][b-1] < min_elev)
			{
			  min_elev = zeta[a+1][b-1]; 
			  channel_node = retrieve_node_from_row_and_column(a+1,b-1);
			}
		    }
		  else if (stnet[a+1][b] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a+1][b] < min_elev)
			{
			  min_elev = zeta[a+1][b]; 
			  channel_node = retrieve_node_from_row_and_column(a+1,b);
			}
		    }
		  else if (stnet[a+1][b+1] != NoDataValue)
		    {
		      if (min_elev == NoDataValue || zeta[a+1][b+1] < min_elev)
			{
			  min_elev = zeta[a+1][b+1]; 
			  channel_node = retrieve_node_from_row_and_column(a+1,b+1);
			}
		    }
		  trace_metrics.push_back(float(channel_node));
		}
	      //         if (relief > 0) ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "\n";
	      //         else ++neg_count;
	      if (relief <= 0) ++neg_count;
	    }
	  else
	    { //unable to route using aspects
	      //this will encompass the skipped traces
	      //         ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
	      ++ns_count;
	    }
	}
    }
  output_trace_coordinates = trace_coordinates;
  output_trace_metrics = trace_metrics;
  output_channel_node = channel_node;  
}
//----------------------------------------------------------------------------------------------------------------------



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Hilltop flow routing code built around original code from Martin Hurst. Based on
// Lea (1992), with improvements discussed by Tarboton (1997) and a solution to the
// problem of looping flow paths implemented.
//
// This code is SLOW but robust, a refactored version may appear, but there may not be
// enough whisky in Scotland to support that endeavour.
//
// The algorithm now checks for local uphill flows and in the case of identifying one,
// D8 flow path is used to push the flow into the centre of the steepest downslope
// cell, at which point the trace is restarted. The same technique is used to cope
// with self intersections of the flow path. These problems are not solved in the
// original paper and I think they are caused at least in part by the high resolution
// topogrpahy we are using.
//
// The code is also now built to take a d infinity flow direction raster instead of an
// aspect raster. See Tarboton (1997) for discussions on why this is the best solution.
//
// The Basins input raster is used to code each hilltop into a basin to allow basin
// averaging to take place.
//
// The final 5 parameters are used to set up printing flow paths to files for visualisation,
// if this is not needed simply pass in false to the two boolean switches and empty variables for the
// others, and the code will run as normal.
//
// The structure of the returned vector< Array2D<float> > is as follows:
// [0] Hilltop Network coded with stream ID
// [1] Hillslope Lengths
// [2] Slope
// [3] Relief
// (4) Fraction Rock Exposure
//
// SWDG 12/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector< Array2D<float> > LSDFlowInfo::HilltopFlowRoutingBedrock(LSDRaster Elevation, LSDRaster Hilltops, LSDRaster Slope,
								LSDIndexRaster StreamNetwork, LSDRaster Aspect, string Prefix, LSDIndexRaster Basins, LSDRaster PlanCurvature,
								bool print_paths_switch, int thinning, string trace_path, bool basin_filter_switch,
								vector<int> Target_Basin_Vector, LSDRaster RockExposure){

  //Declare parameters
  int i,j;
  int a = 0;
  int b = 0;
  float X,Y;
  float mean_slope, relief;
  float length, d;
  float rock_exposure;
  int flag;
  int count = 0;
  int DivergentCountFlag = 0; //Flag used to count the number of divergent cells encountered in a trace
  int PlanarCountFlag;
  float PI = 3.14159265;
  float degs, degs_new, theta;
  float s_local, s_edge;
  float xo, yo, xi, yi, temp_yo1, temp_yo2, temp_xo1, temp_xo2;
  bool skip_trace; //flag used to skip traces where no path to a stream can be found. Will only occur on noisy, raw topography
  float E_Star = 0;
  float R_Star = 0;
  float EucDist = 0;

  //debugging counters
  int ns_count = 0;
  int s_count = 0;
  int neg_count = 0;
  int edge_count = 0;
  int ht_count = 0;

  // a direction flag numbered 1,2,3,4 for E,S,W,N respectively
  int dir;
  
  float ymax = YMinimum + NRows*DataResolution;
  
  //Get data arrays from LSDRasters
  Array2D<float> zeta = Elevation.get_RasterData(); //elevation
  Array2D<int> stnet = StreamNetwork.get_RasterData(); // stream network
  Array2D<float> aspect = Aspect.get_RasterData(); //aspect
  Array2D<float> hilltops = Hilltops.get_RasterData(); //hilltops
  Array2D<float> slope = Slope.get_RasterData(); //hilltops
  Array2D<int> basin = Basins.get_RasterData(); //basins
  Array2D<float> rock = RockExposure.get_RasterData(); // Rock Exposure

  //empty arrays for data to be stored in
  Array2D<float> rads(NRows,NCols);
  Array2D<float> path(NRows, NCols, 0.0);
  Array2D<float> blank(NRows, NCols, 0.0);
  Array2D<float> RoutedHilltops(NRows,NCols,NoDataValue);
  Array2D<float> HillslopeLength_Array(NRows,NCols,NoDataValue);
  Array2D<float> Slope_Array(NRows,NCols,NoDataValue);
  Array2D<float> Relief_Array(NRows,NCols,NoDataValue);
  Array2D<float> Rock_Array(NRows,NCols,NoDataValue);
  
  //vector to store the output data arrays in one vector that can be returned
  vector< Array2D<float> > OutputArrays;

  int vec_size = 1000000;
  
  Array1D<double> easting(NCols);
  Array1D<double> northing(NRows);
  Array1D<double> east_vec(vec_size);
  Array1D<double> north_vec(vec_size);

  ofstream ofs;

  //create the output filename from the user supplied filename prefix
  stringstream ss_filename;
  ss_filename << Prefix << "_HilltopData.csv";

  ofs.open(ss_filename.str().c_str());

  if( ofs.fail() ){
    cout << "\nFATAL ERROR: unable to write to " << ss_filename.str() << endl;
    exit(EXIT_FAILURE);
  }
  ofs << "X,Y,hilltop_id,S,R,Lh,BasinID,StreamID,HilltopSlope,DivergentCount\n";

  //calculate northing and easting
  for (i=0;i<NRows;++i){
    northing[i] = ymax - i*DataResolution - 0.5;
  }
  for (j=0;j<NCols;++j){
    easting[j] = XMinimum + j*DataResolution + 0.5;
  }

  //convert aspects to radians with east as theta = 0/2*pi
  for (i=0; i<NRows; ++i) {
    for (j=0; j<NCols; ++j) {
      //convert aspects to radians with east as theta = 0/2*pi
      if (rads[i][j] != NoDataValue) rads[i][j] = BearingToRad(aspect[i][j]);
    }
  }

  // cycle through study area, find hilltops and trace downstream
  for (i=1; i<NRows-1; ++i) {
    cout << flush <<  "\tRow: " << i << " of = " << NRows-1 << "              \r";
    for (j=1; j<NCols-1; ++j) {
      
      // ignore edge cells and non-hilltop cells
      // route initial node by aspect and get outlet coordinates
      if (hilltops[i][j] != NoDataValue) {
        
	length = 0;
	rock_exposure = 0;
	flag = true;
	count = 1;
	path = blank.copy();
        DivergentCountFlag = 0; //initialise count of divergent cells in trace
        PlanarCountFlag = 0;
        skip_trace = false; //initialise skip trace flag as false, will only be switched if no path to stream can be found. Very rare.
        
        E_Star = 0;
        R_Star = 0;
        EucDist = 0;
        
	++ht_count;

	degs = aspect[i][j];
	theta = rads[i][j];
	a = i;
	b = j;
	path[a][b] += 1;
	east_vec[0] = easting[b];
	north_vec[0] = northing[a];
	s_local = slope[a][b];

	//test direction, calculate outlet coordinates and update indicies
	// easterly
	if (degs >= 45 && degs < 135) {
	  //cout << "\neasterly" << endl;
	  xo = 1, yo = (1+tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 0, yi = yo;
	  dir = 1;
	  east_vec[count] = easting[b] + 0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  ++b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//southerly
	else if (degs >= 135 && degs < 225) {
	  //cout << "\nsoutherly" << endl;
	  xo = (1-(1/tan(theta)))/2, yo = 0;
	  d = abs(1/(2*cos((PI/2)-theta)));
	  xi = xo, yi = 1;
	  dir = 2;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] - 0.5*DataResolution;
	  ++a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	// westerly
	else if (degs >= 225 && degs < 315) {
	  xo = 0, yo = (1-tan(theta))/2;
	  d = abs(1/(2*cos(theta)));
	  xi = 1,	yi = yo;
	  dir = 3;
	  east_vec[count] = easting[b] -0.5*DataResolution;
	  north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	  --b;
	  if (yi == 0) yi = 0.00001;
	  else if (yi == 1) yi = 1 - 0.00001;
	}
	//northerly
	else if (degs >= 315 || degs < 45) {
	  xo = (1+(1/tan(theta)))/2, yo = 1;
	  d = abs(1/(2*cos((PI/2) - theta)));
	  xi = xo, yi = 0;
	  dir = 4;
	  east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	  north_vec[count] = northing[a] + 0.5*DataResolution;
	  --a;
	  if (xi == 0) xi = 0.00001;
	  else if (xi == 1) xi = 1 - 0.00001;
	}
	else {
	  cout << "FATAL ERROR, Kinematic routing algorithm enountered null aspect value" << endl;
	  exit(EXIT_FAILURE);
	}

	//collect slopes and totals weighted by path length
	length += d;
	s_local = slope[a][b];
	rock_exposure += rock[a][b]*d;

	//continue trace until a stream node is encountered
	while (flag == true && a > 0 && a < NRows-1 && b > 0 && b < NCols-1) {   //added boudary checking to catch cells which flow off the	edge of the DEM tile.
				  
	  path[a][b] += 1;
				  
	  degs_new = aspect[a][b];
	  theta = rads[a][b];
	  ++count;
				  
	  //Test for perimeter flow paths
	  if ((dir == 1 && degs_new > 0 && degs_new < 180)
	      || (dir == 2 && degs_new > 90 && degs_new < 270)
	      || (dir == 3 && degs_new > 180 && degs_new < 360)
	      || ((dir == 4 && degs_new > 270) || (dir == 4 && degs_new < 90))) {
				    
	    //DO NORMAL FLOW PATH
	    //set xo, yo to 0 and 1 in turn and test for true outlet (xi || yi == 0 || 1)
	    temp_yo1 = yi + (1-xi)*tan(theta); 		// xo = 1
	    temp_xo1 = xi + (1-yi)*(1/tan(theta)); 	// yo = 1
	    temp_yo2 = yi - xi*tan(theta);			// xo = 0
	    temp_xo2 = xi - yi*(1/tan(theta));		// yo = 0
				    
	    // can't outlet at same point as inlet
	    if (dir == 1) temp_yo2 = -1;
	    else if (dir == 2) temp_xo1 = -1;
	    else if (dir == 3) temp_yo1 = -1;
	    else if (dir == 4) temp_xo2 = -1;
				    
	    s_local = slope[a][b];           
				    
	    if (temp_yo1 <= 1 && temp_yo1 > 0) {              
	      xo = 1, yo = temp_yo1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
	      xi = 0, yi = yo,
		dir = 1;
	      east_vec[count] = easting[b] + 0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      ++b;
	      if (xi== 0 && yi == 0) yi = 0.00001;
	      else if (xi== 0 && yi == 1) yi = 1 - 0.00001;
	    }
	    else if (temp_xo2 <= 1 && temp_xo2 > 0) {
	      xo = temp_xo2, yo = 0;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));   
	      xi = xo, yi = 1,
		dir = 2;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] - 0.5*DataResolution;
	      ++a;
	      if (xi== 0 && yi == 1) xi = 0.00001;
	      else if (xi== 1 && yi == 1) xi = 1 - 0.00001;
	    }
	    else if (temp_yo2 <= 1 && temp_yo2 > 0) {
	      xo = 0, yo = temp_yo2;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));      
	      xi = 1, yi = yo,
		dir = 3;
	      east_vec[count] = easting[b] -0.5*DataResolution;
	      north_vec[count] = northing[a] + yo - 0.5*DataResolution;
	      --b;
	      if (xi== 1 && yi == 0) yi = 0.00001;
	      else if (xi== 1 && yi == 1) yi = 1 - 0.00001;
	    }
				    
	    else if (temp_xo1 <= 1 && temp_xo1 > 0) {
	      xo = temp_xo1, yo = 1;
	      d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));   
	      xi = xo, yi = 0,
		dir = 4;
	      east_vec[count] = easting[b] + xo - 0.5*DataResolution;
	      north_vec[count] = northing[a] + 0.5*DataResolution;
	      --a;
	      if (xi == 0 && yi == 0) xi = 0.00001;
	      else if (xi== 1 && yi == 0) xi = 1 - 0.00001;
	    }           
				    
	  }          
	  else {
				    
	    // ROUTE ALONG EDGES						
	    if (dir	== 1) {
	      if (degs_new <= 90 || degs_new >= 270) { //secondary compenent of flow is north                   
		xo = 0.00001, yo = 1;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else if (degs_new > 90 && degs_new < 270) {  //secondary component is south							  
		xo = 0.00001, yo = 0;
		s_edge = abs(s_local*sin((PI/2)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }						 
	      else {
		cout << "Flow unable to route N or S " << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
		//exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 2) {
	      if 	(degs_new >= 0 && degs_new <= 180) { //secondary component is East							  
		xo = 1, yo = 1-0.00001;
		s_edge = abs(s_local*sin((2/PI)-theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else if (degs_new > 180 && degs_new <= 360) {  //secondary component is West						
		xo = 0, yo = 1-0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
		//exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 3) {
	      if 	(degs_new >= 90 && degs_new <= 270) {  //secondary component is South
		xo = 1-0.00001, yo = 0;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1-yo;
		dir = 2;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] - 0.5*DataResolution;
		++a;
	      }
	      else if (degs_new > 270 || degs_new < 90) {   //secondary component is North
		xo = 1-0.00001, yo = 1;
		s_edge = abs(s_local*sin((2/PI) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = xo, yi = 1- yo;
		dir = 4;
		east_vec[count] = easting[b] + xo - 0.5*DataResolution;
		north_vec[count] = northing[a] + 0.5*DataResolution;
		--a;
	      }
	      else {
		cout << "Flow unable to route N or S" << endl;  //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
		//exit(EXIT_FAILURE);
	      }
	    }
	    else if (dir == 4) {
	      if 	(degs_new >= 180 && degs_new <= 360) { //secondary component is West
		xo = 0, yo = 0.00001;
		s_edge = abs(s_local*sin((PI/2) - theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 3;
		east_vec[count] = easting[b] -0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		--b;
	      }
	      else if (degs_new >= 0 && degs_new < 180) { //secondary component is East				  
		xo = 1, yo = 0.00001;
		s_edge = abs(s_local*sin(theta));
		d = sqrt((pow((xo-xi),2) + pow((yo-yi),2)));
		xi = 1-xo, yi = yo;
		dir = 1;
		east_vec[count] = easting[b] + 0.5*DataResolution;
		north_vec[count] = northing[a] + yo - 0.5*DataResolution;
		++b;
	      }
	      else {
		cout << "Flow unable to route E or W" << endl; //something has gone very wrong...
		cout << "Trace skipped.\n" << endl; 
		skip_trace = true;
		//exit(EXIT_FAILURE);
	      }
	    }
				    
	  }
				  
	  if (path[a][b] < 1){  // only update length on 'first slosh'
	    length += d;   
	    rock_exposure += rock[a][b]*d;
	  }
	  else if (path[a][b] >= 3){ //update the skip trace flag so we can categorise each trace          
	    skip_trace = true;          
	  } 
				  
	  degs = degs_new;
				  
	  // test for plan curvature here and set a flag if flow is divergent or convergent but continue trace regardless
	  // The larger the counter the more convergent or divergent the trace is
	  if (abs(PlanCurvature.get_data_element(a,b)) > (0.001)){
	    ++DivergentCountFlag;
	  }
	  else {
	    ++PlanarCountFlag;
	  }
				  
	  if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 || stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue || path[a][b] >= 3 || skip_trace == true) flag = false;
	}
				
	if (a == 0 || b == 0 ||	a == NRows-1 || b == NCols-1 ){
	  // avoid going out of bounds.
				  
	  // this is caused by having a hilltop on the first row or col away from the border
	  // eg i or j == 1 or nrows/ncols - 2 and flowing towards the edge.
	  // can fix with a test here for if streamnet[a][b] != NDV otherwise trace will fail *correctly*
				  
	  ++edge_count;
				  
	}
	else
	  {
	    //if trace finished at a stream, print hillslope info.
	    if (stnet[a][b] != NoDataValue || stnet[a-1][b-1] != NoDataValue || stnet[a][b-1] != NoDataValue || stnet[a+1][b-1] != NoDataValue || stnet[a+1][b] != NoDataValue || stnet[a+1][b+1] != NoDataValue || stnet[a][b+1] != NoDataValue || stnet[a-1][b+1] != NoDataValue || stnet[a-1][b] != NoDataValue)
	      {
		path[a][b] = 1;
					
		++s_count;
					
		X = XMinimum + j*DataResolution;
		Y = YMinimum - (NRows-i)*DataResolution;
		relief = zeta[i][j] - zeta[a][b];
		mean_slope = relief/(length * DataResolution);
					
		// update arrays with the current metrics
		RoutedHilltops[i][j] = 1;
		HillslopeLength_Array[i][j] = (length * DataResolution);
		Slope_Array[i][j] = mean_slope;
		Relief_Array[i][j] = relief;
		Rock_Array[i][j] = rock_exposure/length;
					
		//calculate an E* and R* Value assuming S_c of 0.8
		E_Star = (2.0 * abs(hilltops[i][j])*(length*DataResolution))/0.8;
		R_Star = relief/((length*DataResolution)*0.8);
					
		//calulate the Euclidean distance between the start and end points of the trace
		EucDist = sqrt((pow(((i+0.5)-(a+yo)),2) + pow(((j+0.5)-(b+xo)),2))) * DataResolution;
					
					
		if (relief > 0){
		  ofs << X << "," << Y << "," << hilltops[i][j] << "," << mean_slope << "," << relief << "," << length*DataResolution << "," << basin[i][j] << "," << stnet[a][b] << "," << slope[i][j] << "," << DivergentCountFlag << "," << PlanarCountFlag << "," << E_Star << "," << R_Star << "," << EucDist << "," << rock_exposure/length << "\n";
		}
		else {
		  ++neg_count;
		}
	      }
	    else{  //unable to route using aspects
	      //this will encompass skipped traces
	      ofs << "fail: " << a << " " << b << " " << i << " " << j << endl;
	      ++ns_count;
	    }
	  }
				
	//This block checks the various path printing options and writes the data out accordingly
	if (print_paths_switch == true){
	  if (ht_count % thinning == 0){
	    if (hilltops[i][j] != NoDataValue && skip_trace == false){ //check that the current i,j tuple corresponds to a hilltop, ie there is actually a trace to write to file, and check that the trace was valid.
				      
	      //create stringstream object to create filename
	      ofstream pathwriter;
				      
	      //create the output filename from the user supplied path
	      stringstream ss_path;
	      ss_path << trace_path << i << "_" << j << "_trace.txt";
				      
	      pathwriter.open(ss_path.str().c_str());
				      
	      if(pathwriter.fail() ){
		cout << "\nFATAL ERROR: unable to write to " << ss_path.str() << endl;
		exit(EXIT_FAILURE);
	      }
				      
	      for (int v = 0; v < count+1; ++v){
		if (basin_filter_switch == false){
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << "," << rock_exposure/length << endl;
		}
		else if (basin_filter_switch == true && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[a][b]) != Target_Basin_Vector.end() && find(Target_Basin_Vector.begin(), Target_Basin_Vector.end(), basin[i][j]) != Target_Basin_Vector.end()){  //is this correct? evaulating to not equal one past the end of the vector should equal that the value is found
		  pathwriter << setiosflags(ios::fixed) << setprecision(7) << east_vec[v] << " " << north_vec[v] << " " << DivergentCountFlag << " " << length << " " << PlanarCountFlag << " " << E_Star << " " << R_Star << " " << EucDist << "," << rock_exposure/length << endl;
		}
	      }
	      pathwriter.close();
	    }
	  }
	}
	// End of path printing logic
      }
    }   //for loop i,j
  }
	
  ofs.close();
	
  //add the data arrays to the output vector
  OutputArrays.push_back(RoutedHilltops);
  OutputArrays.push_back(HillslopeLength_Array);
  OutputArrays.push_back(Slope_Array);
  OutputArrays.push_back(Relief_Array);
  OutputArrays.push_back(Rock_Array);
	
  //Print debugging info to screen
  cout << endl; //push output onto new line
  cout << "Hilltop count: " << ht_count << endl;
  cout << "Stream count: " << s_count << endl;
  cout << "Fail count: " << ns_count << endl;
  cout << "Uphill count: " << neg_count << endl;
  cout << "Edge count: " << edge_count << endl;
  
  return OutputArrays;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This method removes end nodes which are not the uppermost extent of the channel network.
// SWDG 23/7/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::ProcessEndPointsToChannelHeads(LSDIndexRaster Ends){

  Array2D<int> EndArray = Ends.get_RasterData();
  vector<int> Sources;
  
  //make a map containing each nodeindex where true means it is a valid channel head, eg the top of the network
  map<int,bool> EndStatus;
  vector<int> EndNodes;
  
  for(int i=1; i<NRows-1; ++i){
    for(int j=1; j<NCols-1; ++j){
      if (EndArray[i][j] != NoDataValue){
        int nodeindex = retrieve_node_from_row_and_column (i,j);
        EndStatus[nodeindex] = true;
        EndNodes.push_back(nodeindex);
      }
    }
  }
      
  for (int q = 0; q < int(EndNodes.size());++q){
    cout << flush << q << " of " << EndNodes.size() << "\r"; 
    int CurrentNode = EndNodes[q];
    if (EndStatus[CurrentNode] == true){
    
      bool stop = false;
      
      while (stop == false){
        int DownslopeNode;
        int Downslopei;
        int Downslopej;	 
              
        //get steepest descent neighbour
        retrieve_receiver_information(CurrentNode,DownslopeNode,Downslopei,Downslopej);
      
        if (find(EndNodes.begin(), EndNodes.end(), DownslopeNode) != EndNodes.end()){
          EndStatus[DownslopeNode] = false;
          stop = true;          
        }		      
        
        //check for out of bounds
        if (Downslopei == 0 || Downslopei == NRows - 1 || Downslopej == 0 || Downslopej == NCols - 1){
          stop = true;
        }
        
        //check for a node with no downslope neughbours
        if (CurrentNode == DownslopeNode){
          stop = true;
        }
               
        CurrentNode = DownslopeNode; 
        
      }
    
    }
  }
  cout << endl; 

  for (int q = 0; q < int(EndNodes.size());++q){
    if (EndStatus[EndNodes[q]] == true){  
      Sources.push_back(EndNodes[q]);
    }
  }

  return Sources;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This method removes single pixel channels from a channel network.
// SWDG 23/7/15
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> LSDFlowInfo::RemoveSinglePxChannels(LSDIndexRaster StreamNetwork, vector<int> Sources){

  for (int q = 0; q < int(Sources.size());++q){
    
    int CurrentNode = Sources[q];
    int Currenti;
    int Currentj;
    retrieve_current_row_and_col(CurrentNode,Currenti,Currentj);
    int CurrentOrder = StreamNetwork.get_data_element(Currenti,Currentj);
    
    //get steepest descent neighbour
    int DownslopeNode;
    int Downslopei;
    int Downslopej;
    retrieve_receiver_information(CurrentNode,DownslopeNode,Downslopei,Downslopej);
    int DownslopeOrder = StreamNetwork.get_data_element(Downslopei,Downslopej);
    
    if (CurrentOrder != DownslopeOrder){
      //remove the value from the list of nodes -> Sources is passed by val, so this will not change values in sources outide this method
      Sources.erase(remove(Sources.begin(), Sources.end(), Sources[q]), Sources.end());
    }
  
  }

  return Sources;
}


#endif
