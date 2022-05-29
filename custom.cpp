/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

void create_cell_types(void)
{
	// set the random seed
	SeedRandom(parameters.ints("random_seed"));

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, std::vector<double> center)
{
	std::vector<std::vector<double>> cells;
	int xc = 0, yc = 0, zc = 0;
	double x_spacing = cell_radius * sqrt(3);
	double y_spacing = cell_radius * 2;
	double z_spacing = cell_radius * sqrt(3);

	std::vector<double> tempPoint(3, 0.0);
	// std::vector<double> cylinder_center(3,0.0);

	double x0 = center[0];
	double y0 = center[1];
	double z0 = center[2];

	for (double z = z0 - sphere_radius; z < z0 + sphere_radius; z += z_spacing, zc++)
	{
		for (double x = x0 - sphere_radius; x < x0 + sphere_radius; x += x_spacing, xc++)
		{
			for (double y = y0 - sphere_radius; y < y0 + sphere_radius; y += y_spacing, yc++)
			{
				tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
				tempPoint[1] = y + (xc % 2) * cell_radius;
				tempPoint[2] = z;
				double dr = sqrt((tempPoint[0] - x0) * (tempPoint[0] - x0) + (tempPoint[1] - y0) * (tempPoint[1] - y0) + (tempPoint[2] - z0) * (tempPoint[2] - z0));
				if (dr < sphere_radius)
				{
					cells.push_back(tempPoint);
				}
			}
		}
	}
	return cells;
}

void setup_tissue(void)
{
	// place a cluster of tumor cells at the center

	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;

	double tumor_radius1 =
		parameters.doubles("first_tumor_radius"); // 250.0;
	double tumor_radius2 =
		parameters.doubles("second_tumor_radius"); // 250.0;

	Cell *pCell = NULL;

	double distance_between_centers = parameters.doubles("distance_between_centers");
	std::cout << "distance between two spheroids is " << distance_between_centers << std::endl;

	if (parameters.ints("number_of_spheroids") == 1)
		distance_between_centers = 0.0;

	// first spheroid
	std::vector<double> center1 = {-distance_between_centers * 0.5, 0, 0};
	std::cout << "position of centers1 : " << center1[0] << ' ' << center1[1] << ' ' << center1[2] << std::endl;
	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius, tumor_radius1, center1);
	std::cout << "creating " << positions.size() << " closely-packed tumor cells in tissue 1... " << std::endl;

	Cell_Definition *pCD = cell_definitions_by_index[0];

	for (int i = 0; i < positions.size(); i++)
	{
		pCell = create_cell(*pCD); // tumor cell
		pCell->assign_position(positions[i]);
	}

std::cout<<"OK0"<<std::endl;
	int Ncell1=(*all_cells).size();//total number of cells in the first spheroid
	std::srand(unsigned(std::time(0)));
	std::vector<int> index1;
	for (int i=0;i<Ncell1;i++)
	{
		index1.push_back(i);
	}
	std::random_shuffle(index1.begin(),index1.end());

	double percentage1=parameters.doubles("percentage_1");
	int N1=(int)Ncell1*percentage1;////number of cells to be modified
	for (int i=0;i<N1;i++)
	{
		Cell *pCell = (*all_cells)[index1[i]];//select the cell to be modified
		pCell->phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("adhesion_strength_1");
		pCell->phenotype.motility.migration_speed= parameters.doubles("speed_1");

		pCell->type=2;

	}

	std::cout<<"OK1"<<std::endl;
	// second spheroid
	if (parameters.ints("number_of_spheroids") == 2)
	{

		std::cout << "# of cell definitions " << cell_definitions_by_index.size() << std::endl;
		Cell_Definition *pCD2;
		if (cell_definitions_by_index.size() > 1)
		{
			pCD2 = cell_definitions_by_index[1];
			std::cout << "use second type of cell for the second spheroid" << std::endl;
		}
		else
		{
			pCD2 = cell_definitions_by_index[0];
			std::cout << "only one cell type exists, so we use the same type for the second spheroid" << std::endl;
		}
		std::vector<double> center2 = {distance_between_centers * 0.5, 0, 0};
		std::cout << "position of centers2 : " << center2[0] << ' ' << center2[1] << ' ' << center2[2] << std::endl;
		positions = create_cell_sphere_positions(cell_radius, tumor_radius2, center2);
		std::cout << "creating " << positions.size() << " closely-packed tumor cells in tissue 2... " << std::endl;

		for (int i = 0; i < positions.size(); i++)
		{
			pCell = create_cell(*pCD2); // tumor cell
			pCell->assign_position(positions[i]);
		}

std::cout<<"OK2"<<std::endl;
///new

int Ncell2=(*all_cells).size()-Ncell1;//total number of cells in the second spheroid
	//std::srand(unsigned(std::time(0)));
	std::vector<int> index2;
	for (int i=0;i<Ncell2;i++)
	{
		index2.push_back(i);
	}
	std::random_shuffle(index2.begin(),index2.end());


	double percentage2=parameters.doubles("percentage_2");
	int N2=(int)Ncell2*percentage2;//number of cells to be modified
	for (int i=0;i<N2;i++)
	{
		Cell *pCell = (*all_cells)[index2[i]+Ncell1];
		pCell->phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("adhesion_strength_2");
						pCell->phenotype.motility.migration_speed= parameters.doubles("speed_2");

		pCell->type=3;
	}
std::cout<<"OK3"<<std::endl;
	////end new
std::cout<<"N1="<<N1<<" N2="<<N2<<std::endl;
std::cout<<"Nc1="<<Ncell1<<" Nc2="<<Ncell2<<std::endl;

	}

	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml();
	return;
}

std::vector<std::string> my_coloring_function(Cell *pCell)
{
	return paint_by_number_cell_coloring(pCell);
}

void phenotype_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}

void custom_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}

void contact_function(Cell *pMe, Phenotype &phenoMe, Cell *pOther, Phenotype &phenoOther, double dt)
{
	return;
}

void move_tissues()
{

	// find center and radius of the tissues
	std::vector<double> center1 = {0, 0, 0};
	std::vector<double> center2 = {0, 0, 0};
	int Np1 = 0;
	int Np2 = 0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pCell = (*all_cells)[i];
//		std::cout << "type:" << pCell->type << std::endl;

		if (pCell->type == 0||pCell->type == 2)
		{
			center1[0] += pCell->position[0];
			center1[1] += pCell->position[1];
			center1[2] += pCell->position[2];
			Np1 += 1;
		}
		else if (pCell->type == 1||pCell->type == 3)
		{
			center2[0] += pCell->position[0];
			center2[1] += pCell->position[1];
			center2[2] += pCell->position[2];
			Np2 += 1;
		}
	}

	center1[0] /= (double)Np1;
	center1[1] /= (double)Np1;
	center1[2] /= (double)Np1;

	center2[0] /= (double)Np2;
	center2[1] /= (double)Np2;
	center2[2] /= (double)Np2;

	std::cout << "center of tissue 1 : " << center1[0] << ' ' << center1[1] << ' ' << center1[2] << std::endl;
	std::cout << "center of tissue 2 : " << center2[0] << ' ' << center2[1] << ' ' << center2[2] << std::endl;

	double distance = sqrt((center1[0] - center2[0]) * (center1[0] - center2[0]) + (center1[1] - center2[1]) * (center1[1] - center2[1]) + (center1[2] - center2[2]) * (center1[2] - center2[2]));
	std::cout << "distance of the two tissues is :" << distance << std::endl;

	// we use original tumor radius which should be modified if the radius of the spheroids changed
	double tumor_radius1 =
		parameters.doubles("first_tumor_radius"); // 250.0;
	double tumor_radius2 =
		parameters.doubles("second_tumor_radius"); // 250.0;

	std::vector<double> R1, R2;
	double R1max, R2max;
	R1max = 0;
	R2max = 0;
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pCell = (*all_cells)[i];
		// std::cout << "type:" << pCell->type << std::endl;

		if (pCell->type == 0||pCell->type == 2)
		{
			double tmpr = sqrt((pCell->position[0] - center1[0]) * (pCell->position[0] - center1[0]) + (pCell->position[1] - center1[1]) * (pCell->position[1] - center1[1]) + (pCell->position[2] - center1[2]) * (pCell->position[2] - center1[2]));
			R1.push_back(tmpr);
			if (tmpr > R1max)
				R1max = tmpr;
		}
		else if (pCell->type == 1||pCell->type == 3)
		{
			double tmpr = sqrt((pCell->position[0] - center2[0]) * (pCell->position[0] - center2[0]) + (pCell->position[1] - center2[1]) * (pCell->position[1] - center2[1]) + (pCell->position[2] - center2[2]) * (pCell->position[2] - center2[2]));
			R2.push_back(tmpr);
			if (tmpr > R2max)
				R2max = tmpr;
		}
	}
	tumor_radius1 = R1max;
	tumor_radius2 = R2max;

	std::cout << "LUO: radius 1:" << tumor_radius1 << std::endl;
	std::cout << "LUO: radius 2:" << tumor_radius2 << std::endl;

	double dx = 0.5 * (distance - tumor_radius1 - tumor_radius2);
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pCell = (*all_cells)[i];
		if (pCell->type == 0||pCell->type == 2)
		{
			double newx = pCell->position[0] + dx;
			pCell->assign_position(newx, pCell->position[1], pCell->position[2]);
		}
		else if (pCell->type == 1||pCell->type == 3)
		{
			double newx = pCell->position[0] - dx;
			pCell->assign_position(newx, pCell->position[1], pCell->position[2]);
		}
	}

	return;
}