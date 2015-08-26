//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexChannel
// Land Surface Dynamics IndexChannel
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for retaining information of channels
//
// Data is mostly pointers to indeces withing the LSDFlowInfo object
//  For data about the actuall channel details such as elevation,
//  drainage area, etc. once needs to use the derivative class
//  LSDChannel
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



// Source code for the LSDChannelIndex object
// this object contains the node indices as well as the
// row and col indices for individual channel segments
// these can be arranged arbitrailiy accoridng to channel
// junctions or simply nodes downstream of a given node and upstream
// of another arbitrary node EndNode




//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------


#include <vector>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDIndexChannel.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexChannel_CPP
#define LSDIndexChannel_CPP

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Empty create routine that throws error (you need to give it parameters)
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create()
{
	//cout << "LSDIndexChannel You need to initialize with some parameters" << endl;
	//exit(EXIT_FAILURE);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// first create routine
// creates an index channel with just the node index of the starting and ending nodes
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJN, int EJN, LSDFlowInfo& FlowInfo)
{

	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();

	StartJunction = -1;
	EndJunction = -1;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;


		}
		else
		{
			curr_node = receive_node;
		}
	}

	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// second create routine
// creates an index channel with just the node index of the starting and ending nodes
// also includes junction information
// IMPORTANT
// The starting node is upstream
// the ending node is downstream
// In this create function the junction indices are left blank (this can
// describe a channel between two arbitraty points
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::create(int SJ, int SJN, int EJ, int EJN, LSDFlowInfo& FlowInfo)
{
	NRows = FlowInfo.get_NRows();
	NCols = FlowInfo.get_NCols();
	XMinimum = FlowInfo.get_XMinimum();
	YMinimum = FlowInfo.get_YMinimum();
	DataResolution = FlowInfo.get_DataResolution();
	NoDataValue = FlowInfo.get_NoDataValue();
	GeoReferencingStrings =  FlowInfo.get_GeoReferencingStrings();

	StartJunction = SJ;
	EndJunction = EJ;
	StartNode = SJN;
	EndNode = EJN;

	vector<int> RowI;
	vector<int> ColI;
	vector<int> NdI;

	int curr_node = StartNode;

	// push back the data vecotors with the starting node
	int curr_row, curr_col;
	FlowInfo.retrieve_current_row_and_col(curr_node,curr_row,curr_col);
	NdI.push_back(StartNode);
	RowI.push_back(curr_row);
	ColI.push_back(curr_col);

	int receive_node = -99;
	int receive_row, receive_col;

	// loop through receivers until you get to EndNode
	while(curr_node != EndNode)
	{
		FlowInfo.retrieve_receiver_information(curr_node, receive_node, receive_row,
                                              receive_col);

		//cout << "receive_node: " << receive_node << " and Endnode: " << EndNode << endl;

		NdI.push_back(receive_node);
		RowI.push_back(receive_row);
		ColI.push_back(receive_col);

		if (receive_node == curr_node)
		{
			EndNode = curr_node;
			cout << "Warning, the channel has come to a baselevel node before it has"
			     << endl << "reached the end node" << endl;

		}
		else
		{
			curr_node = receive_node;
		}
	}
	RowSequence = RowI;
	ColSequence = ColI;
	NodeSequence = NdI;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This gets the row and column of the end node
// row and col are replaced
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_row_column_of_end_node(LSDFlowInfo& FlowInfo, int& row, int& col)
{
  FlowInfo.retrieve_current_row_and_col(EndNode,row, col);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at one before the final pixel.
// useful for when you want to get the basin area just above the tributary
// junction, which is usually at EndNode
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_penultimate_node(LSDFlowInfo& FlowInfo)
{
	if(StartNode==EndNode)
	{
		return FlowInfo.retrieve_contributing_pixels_of_node(EndNode);
	}
	else
	{
		int n_nodes = RowSequence.size();
		return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_nodes-2]);
	}
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function removes the final node in the index channel. This is used for
// channels that extend to a downstream junction and the user
// wants to truncate the channel before it encounters a junction
//
// SMM 05/11/2013
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::truncate_final_node()
{
  EndJunction = NoDataValue;
  int n_nodes = NodeSequence.size();

  EndNode = NodeSequence[n_nodes-2];
  NodeSequence.pop_back();
  RowSequence.pop_back();
  ColSequence.pop_back();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the n contributing pixels at the node in the channel (not the node index)
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_contributing_pixels_at_node(int n_node, LSDFlowInfo& FlowInfo)
{
  int n_nodes_in_channel = int(NodeSequence.size());
  if (n_node < 0)
  {
    cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
    exit(EXIT_FAILURE);
  }
  if (n_node >= n_nodes_in_channel)
  {
    cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
    exit(EXIT_FAILURE);
  }
  return FlowInfo.retrieve_contributing_pixels_of_node( NodeSequence[n_node] );
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object
// of a node in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int LSDIndexChannel::get_node_in_channel(int n_node)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	return NodeSequence[n_node];
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the node index (from flowInfo object) and row and column indices
// into an LSDRaster or IndexRaster
// of a node in the channel
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::get_node_row_col_in_channel(int n_node, int& node, int& row, int& col)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	if (n_node < 0)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	if (n_node >= n_nodes_in_channel)
	{
		cout << "LINE 138 LSDIndexChannel Not in channel!" << endl;
		exit(EXIT_FAILURE);
	}
	node = NodeSequence[n_node];
	row = RowSequence[n_node];
	col = ColSequence[n_node];
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the channel onto an index raster. Used to test where the channel is.
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexChannel::print_index_channel_to_index_raster()
{
	int n_nodes_in_channel = int(NodeSequence.size());
	cout << "NRows: " << NRows << " NCols: " << NCols << endl;

	Array2D<int> Channel_array(NRows,NCols,NoDataValue);
	for(int i = 0; i<n_nodes_in_channel; i++)
	{
		//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
		Channel_array[RowSequence[i]][ColSequence[i]]= 1;
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array,GeoReferencingStrings);
	return Channel_loc;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this appends the channel onto an index raster. Used to test where the channel is.
//
// SMM 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexChannel::append_index_channel_to_index_raster(LSDIndexRaster& old_raster)
{
	int n_nodes_in_channel = int(NodeSequence.size());
	cout << "NRows: " << NRows << " NCols: " << NCols << endl;

	Array2D<int> Channel_array= old_raster.get_RasterData();
	for(int i = 0; i<n_nodes_in_channel; i++)
	{
		//cout << "row: " << RowSequence[i] << " col: " << ColSequence[i] << endl;
		Channel_array[RowSequence[i]][ColSequence[i]]= 1;
	}

	LSDIndexRaster Channel_loc(NRows,NCols, XMinimum, YMinimum, DataResolution, NoDataValue, Channel_array,GeoReferencingStrings);
	old_raster = Channel_loc;
}



#endif
