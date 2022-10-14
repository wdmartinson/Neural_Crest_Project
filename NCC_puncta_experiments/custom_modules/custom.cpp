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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

// declare cell definitions here

Cell_Definition follower_cell;
Cell_Definition fibronectin_puncta;
Cell_Definition edge_agent;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	// SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
	SeedRandom(); // Completely random initial condition

	// housekeeping

	initialize_default_cell_definition();
	// Name the default cell type

	cell_defaults.type = 0;
	cell_defaults.name = "leader_cell";
	// set default cell cycle model

 	cell_defaults.functions.cycle_model = live;

	// set default_cell_functions
	cell_defaults.functions.update_migration_bias = NULL;
	// cell_defaults.functions.update_migration_bias = leader_cell_migration_bias_and_speed_update_rule;
	cell_defaults.functions.contact_function = NULL;
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;

	cell_defaults.functions.update_velocity = simpler_function_leader_follower_update_cell_velocity;

	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.volume_update_function = NULL;
	// For follower cells and leader cells, this is the boundary condition
 	cell_defaults.functions.custom_cell_rule = record_total_velocity;
	std::vector<double> zeros(3,0.0);
	cell_defaults.custom_data.add_vector_variable("total_velocity", "um/min", zeros);
	cell_defaults.custom_data.add_vector_variable("haptotaxis_direction", "dimensionless", zeros);
	cell_defaults.custom_data.add_vector_variable("contact_guidance_direction", "dimensionless", zeros);
	cell_defaults.custom_data.add_vector_variable("volume_exclusion_direction", "dimensionless", zeros);
	cell_defaults.custom_data.add_variable("did_you_sense_fibronectin", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("are_you_covered_by_a_cell", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("time_until_next_filopodia_drop", "min", 0);
	cell_defaults.custom_data.add_variable("time_since_last_filopodia_drop", "min", 0);
	cell_defaults.custom_data.add_variable("b_mot", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("b_vol_exclusion", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("fibronectin_ID", "dimensionless", -1);

	/* grab code from heterogeneity */

	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;

	cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance = 1.0;

	// make sure the defaults are self-consistent.

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	static int fibronectin_index = microenvironment.find_density_index("fibronectin");
	// static int VEGF_index = microenvironment.find_density_index("VEGF");
	// static int DAN_index = microenvironment.find_density_index("DAN");

	cell_defaults.phenotype.secretion.secretion_rates[fibronectin_index] = parameters.doubles("leader_cell_fibronectin_secretion_rate");
	// cell_defaults.phenotype.secretion.secretion_rates[VEGF_index] = 0.0;
	// cell_defaults.phenotype.secretion.secretion_rates[DAN_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[fibronectin_index] = 99.9;
	// cell_defaults.phenotype.secretion.saturation_densities[VEGF_index] = 0.0;
	// cell_defaults.phenotype.secretion.saturation_densities[DAN_index] = 0.0;
	cell_defaults.phenotype.secretion.uptake_rates[fibronectin_index] = 0.0;
	// cell_defaults.phenotype.secretion.uptake_rates[VEGF_index] = parameters.doubles("leader_cell_VEGF_uptake_rate");
	// cell_defaults.phenotype.secretion.uptake_rates[DAN_index] = parameters.doubles("leader_cell_DAN_uptake_rate");
	// cell_defaults.phenotype.secretion.set_all_secretion_to_zero();
	// cell_defaults.phenotype.secretion.set_all_uptake_to_zero();

	cell_defaults.phenotype.cycle.data.transition_rate( 0 , 0 ) = 0.0;
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0;
	cell_defaults.phenotype.motility.is_motile = true;

	cell_defaults.phenotype.motility.persistence_time = parameters.doubles("cell_persistence_time"); // how long on average it takes for cells to change their current velocity vector, min
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("default_cell_speed"); // um/min;
	cell_defaults.phenotype.motility.migration_bias_direction = {1.0, 0.0, 0.0};
	cell_defaults.phenotype.motility.migration_bias = 1.0;

	cell_defaults.phenotype.mechanics.cell_BM_adhesion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("leaders_cell_cell_repulsion_strength");
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_BM_repulsion_strength = 0.0;

	cell_defaults.phenotype.geometry.radius = parameters.doubles("filopodia_sensing_radius");
	cell_defaults.phenotype.volume.total = 4.0/3.0*PhysiCell_constants::pi*pow(cell_defaults.phenotype.geometry.radius, 3);
	cell_defaults.phenotype.geometry.nuclear_radius = parameters.doubles("default_cell_radius");
	cell_defaults.phenotype.volume.nuclear = 4.0/3.0*PhysiCell_constants::pi*pow(cell_defaults.phenotype.geometry.nuclear_radius, 3);
	cell_defaults.parameters.pReference_live_phenotype = &(cell_defaults.phenotype);

	// Now, let's define another cell type.
	// It's best to just copy the default and modify it.

	follower_cell = cell_defaults;
	follower_cell.type = 1;
	follower_cell.name = "follower_cell";

	follower_cell.functions.update_migration_bias = NULL;

	// follower_cell.phenotype.motility.migration_bias = parameters.doubles("follower_cell_migration_bias");
	follower_cell.phenotype.motility.persistence_time = parameters.doubles("cell_persistence_time"); // how long it takes for cells to change their current velocity vector
	follower_cell.phenotype.secretion.secretion_rates[fibronectin_index] = parameters.doubles("follower_cell_fibronectin_secretion_rate");
	// follower_cell.phenotype.secretion.secretion_rates[VEGF_index] = 0.0;
	// follower_cell.phenotype.secretion.secretion_rates[DAN_index] = 0.0;
	follower_cell.phenotype.secretion.saturation_densities[fibronectin_index] = 99.9;
	// follower_cell.phenotype.secretion.saturation_densities[VEGF_index] = 0.0;
	// follower_cell.phenotype.secretion.saturation_densities[DAN_index] = 0.0;
	follower_cell.phenotype.secretion.set_all_uptake_to_zero();
	// follower_cell.phenotype.secretion.uptake_rates[DAN_index] = parameters.doubles("follower_cell_DAN_uptake_rate");
	follower_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("followers_cell_cell_repulsion_strength");
	follower_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
	// make sure the new cell type has its own reference live phenotype
	follower_cell.parameters.pReference_live_phenotype = &( follower_cell.phenotype );

	// Create a new "cell", which is really just a fibronectin puncta that can't interact with cells, but which cells may sense:
	// In this way, possible for leaders to "secrete" FN by dividing and putting a new cell at their current location.
	fibronectin_puncta = cell_defaults;
	fibronectin_puncta.type = 2;
	fibronectin_puncta.name = "fibronectin_puncta";

	fibronectin_puncta.functions.update_migration_bias = NULL;
	fibronectin_puncta.functions.update_velocity = NULL;
	fibronectin_puncta.functions.custom_cell_rule = NULL;

	fibronectin_puncta.phenotype.motility.migration_speed = 0.0;
	fibronectin_puncta.phenotype.motility.migration_bias_direction = zeros;
	fibronectin_puncta.phenotype.geometry.radius = 0.5*parameters.doubles("default_cell_radius");
	fibronectin_puncta.phenotype.volume.total = 4.0/3.0*PhysiCell_constants::pi*pow(fibronectin_puncta.phenotype.geometry.radius, 3);
	fibronectin_puncta.phenotype.geometry.nuclear_radius = fibronectin_puncta.phenotype.geometry.radius;
	fibronectin_puncta.phenotype.volume.nuclear = fibronectin_puncta.phenotype.volume.total;
	fibronectin_puncta.phenotype.secretion.set_all_uptake_to_zero();
	fibronectin_puncta.phenotype.secretion.set_all_secretion_to_zero();
	fibronectin_puncta.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
	fibronectin_puncta.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
	fibronectin_puncta.parameters.pReference_live_phenotype = &(fibronectin_puncta.phenotype);

	edge_agent = fibronectin_puncta;
	edge_agent.type = 3;
	edge_agent.name = "edge_agent";
	edge_agent.phenotype.geometry.radius = parameters.doubles("default_cell_radius");
	edge_agent.phenotype.volume.total = 4.0/3.0*PhysiCell_constants::pi*pow(edge_agent.phenotype.geometry.radius, 3);
	edge_agent.phenotype.geometry.nuclear_radius = edge_agent.phenotype.geometry.radius;
	edge_agent.phenotype.volume.nuclear = edge_agent.phenotype.volume.total;
	edge_agent.parameters.pReference_live_phenotype = &(edge_agent.phenotype);

	// Construct all cell definitions:
	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

/* now this is in XML
	default_microenvironment_options.X_range = {-1000, 1000};
	default_microenvironment_options.Y_range = {-1000, 1000};
	default_microenvironment_options.simulate_2D = true;
*/

	// make sure to override and go back to 2D
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << "\n";
		default_microenvironment_options.simulate_2D = true;
	}

/* now this is in XML
	// no gradients need for this example

	default_microenvironment_options.calculate_gradients = false;

	// set Dirichlet conditions

	default_microenvironment_options.outer_Dirichlet_conditions = true;

	// if there are more substrates, resize accordingly
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;

	// set initial conditions
	default_microenvironment_options.initial_condition_vector = { 38.0 };
*/

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	// ***************************************************************************
	// Custom initial condition (non homogeneous)
	// if(parameters.bools("VEGF_strip_initial_condition"))
	// {
	// 	int VEGF_substrate_index = microenvironment.find_density_index( "VEGF" );
	// 	double Y_strip = parameters.doubles("y_length_of_VEGF_strip");
	// 	for(int i=0; i<microenvironment.number_of_voxels(); i++)
	// 	{
	// 		if(fabs(microenvironment.voxels(i).center[1])>Y_strip)
	// 		{
	// 			microenvironment.density_vector(i)[VEGF_substrate_index] = 0.0;
	// 		}
	// 	}
	// }
	// if(parameters.bools("DAN_is_in_domain"))
	// {
	// 	double X_left = microenvironment.mesh.bounding_box[0];
	// 	int DAN_substrate_index = microenvironment.find_density_index("DAN");
	// 	double width_of_DAN_rectangle = parameters.doubles("width_of_DAN_rectangle");
	// 	for(int i = 0; i<microenvironment.number_of_voxels(); i++)
	// 	{
	// 		if(microenvironment.voxels(i).center[0]>=(X_left + width_of_DAN_rectangle))
	// 		{
	// 			microenvironment.density_vector(i)[DAN_substrate_index] = 0.0;
	// 		}
	// 	}
	// }

	double X_left = microenvironment.mesh.bounding_box[0];
	double Y_left = microenvironment.mesh.bounding_box[1];
	double X_length = (microenvironment.mesh.bounding_box[3]-X_left);
	double Y_length = (microenvironment.mesh.bounding_box[4]-Y_left);
	int fibronectin_index = microenvironment.find_density_index("fibronectin");
	// int number_of_fibronectin_columns = int(X_length/parameters.doubles("fibronectin_puncta_wavelength"));
	// int number_of_fibronectin_columns = int((microenvironment.mesh.bounding_box[3] - (X_left + 2*parameters.doubles("default_cell_radius") + 0.25*parameters.doubles("filopodia_sensing_radius")) )/parameters.doubles("fibronectin_puncta_wavelength"));
	int number_of_fibronectin_rows = int(Y_length/parameters.doubles("fibronectin_puncta_wavelength"));
	std::vector<double> position(3,0.0);
	// find nearest voxel to the current position
	for(int n = 0; n<number_of_fibronectin_rows; n++)
	{
		position[1] = Y_left + (n+0.5)*parameters.doubles("fibronectin_puncta_wavelength");
		double left_most_puncta = X_left + parameters.doubles("default_cell_radius") + 0.5*parameters.doubles("filopodia_sensing_radius");
		position[0] = left_most_puncta;
		int k = 0;
		while(position[0] <= microenvironment.mesh.bounding_box[3])
		{
			// position[0] = X_left + (k + 0.5)*parameters.doubles("fibronectin_puncta_wavelength");
			position[0] = left_most_puncta + k*parameters.doubles("fibronectin_puncta_wavelength");
			int i = microenvironment.nearest_voxel_index(position);
			microenvironment.density_vector(i)[fibronectin_index] = 30.0;
			if(parameters.bools("fibronectin_strip_initial_condition") && fabs(microenvironment.voxels(i).center[1]) > parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
			}
			else if (parameters.ints("scenario") == 78 && fabs(microenvironment.voxels(i).center[0]-100.0) <= 5.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
			}
			else if (parameters.ints("scenario") == 79 && fabs(microenvironment.voxels(i).center[0]-100.0) <= 25.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
			}
			else if (parameters.ints("scenario") == 80 && fabs(microenvironment.voxels(i).center[1]) <= 5.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
			}
			else if (parameters.ints("scenario") == 81 && fabs(microenvironment.voxels(i).center[1]) <= 25.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
			}
			k++;
		}
	}
	// ***************************************************************************

	microenvironment.display_information( std::cout );
	return;
}

void setup_tissue( void )
{
	// get domain size
	double X_left  = microenvironment.mesh.bounding_box[0];
	double X_length = microenvironment.mesh.bounding_box[3] - X_left;

	double Y_left  = microenvironment.mesh.bounding_box[1];
	double Y_length =  microenvironment.mesh.bounding_box[4] - Y_left;

	// int number_of_columns = int(X_length/2.0/parameters.doubles("default_cell_radius"));
	int number_of_rows = int(Y_length/2.0/parameters.doubles("default_cell_radius"));
	Cell* pC;
	for(int n = 0; n<number_of_rows; n++)
	{
		std::vector<double> position(3,0.0);
		std::vector<double> orientation(3,0.0);
		position[1] = Y_left + (2*n+1)*parameters.doubles("default_cell_radius");

		if (fabs(position[1] - (Y_left + Y_length/2.0)) <= parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
		{
			// Put a follower cell in the first column, with random orientation in [0, 2*pi]:
			position[0] = X_left + parameters.doubles("default_cell_radius");
			double angle = 2.0*PhysiCell_constants::pi*UniformRandom();
			orientation[0] = cos(angle);
			orientation[1] = sin(angle);

			#pragma omp critical(initialize_cells)
			{
				// pC = create_cell( follower_cell );
				// pC->assign_position( position );
				// pC->phenotype.motility.migration_bias_direction = orientation;

				// // Now create a leader cell next to the followers
				// position[0] = X_left + 2*parameters.doubles("default_cell_radius");
				// pC = create_cell();
				// pC->assign_position( position );

				// Now create a leader cell at the entrance of the domain:
				pC = create_cell();
				pC->assign_position( position );
				int fibronectin_ID_index = pC->custom_data.find_variable_index("fibronectin_ID");
				pC->custom_data.variables[fibronectin_ID_index].value = 1.0*pC->ID;
			}

			//NB: Get rid of the following 2 lines of code when you have leaders adopt the FN migration rule!
			int did_you_sense_fibronectin_index = pC->custom_data.find_variable_index("did_you_sense_fibronectin");
			pC->custom_data.variables[did_you_sense_fibronectin_index].value = 1;

			int time_until_next_filopodia_drop_index = pC->custom_data.find_variable_index("time_until_next_filopodia_drop");
			pC->custom_data.variables[time_until_next_filopodia_drop_index].value = -parameters.doubles("average_time_until_next_filopodia_drop")*log(UniformRandom());
			if(parameters.bools("leader_cell_initially_has_random_orientation"))
			{
				angle = 2.0*PhysiCell_constants::pi*UniformRandom();
				orientation[0] = cos(angle);
				orientation[1] = sin(angle);
				pC->phenotype.motility.migration_bias_direction = orientation;
			}// end if

			// Add in an edge agent:
			position[0] = X_left;
			#pragma omp critical (initialize_cells)
			{
				pC = create_cell(edge_agent);
				pC->assign_position(position);
			}
			// These boundary agents are already covered by cells:
			int cover_index = pC->custom_data.find_variable_index("are_you_covered_by_a_cell");
			pC->custom_data.variables[cover_index].value = 1;

		}// end if
	} // end for n

	// Add fibronectin_puncta wherever the initial condition of FN is greater than 0:
	// int number_of_fibronectin_columns = int(X_length/parameters.doubles("fibronectin_puncta_wavelength"));
	// int number_of_fibronectin_columns = int((microenvironment.mesh.bounding_box[3] - (X_left + 2*parameters.doubles("default_cell_radius") + 0.25*parameters.doubles("filopodia_sensing_radius")) )/parameters.doubles("fibronectin_puncta_wavelength"));
	int number_of_fibronectin_rows = int(Y_length/parameters.doubles("fibronectin_puncta_wavelength"));
	std::vector<double> position(3, 0.0);
	for(int n = 0; n < number_of_fibronectin_rows; n++)
	{
		position[1] = Y_left + (n+0.5)*parameters.doubles("fibronectin_puncta_wavelength");
		double left_most_puncta = X_left + parameters.doubles("default_cell_radius") + 0.5*parameters.doubles("filopodia_sensing_radius");
		position[0] = left_most_puncta;
		int k = 0;
		while(position[0] <= microenvironment.mesh.bounding_box[3])
		{
			position[0] = left_most_puncta + k*parameters.doubles("fibronectin_puncta_wavelength");
			k++;
			if(parameters.bools("fibronectin_strip_initial_condition") && fabs(microenvironment.voxels(microenvironment.nearest_voxel_index(position)).center[1]) > parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
			{
				continue;
			}
			else if (parameters.ints("scenario") == 78 && fabs(microenvironment.voxels(microenvironment.nearest_voxel_index(position)).center[0]-100.0) <= 5.0)
			{
				continue;
			}
			else if (parameters.ints("scenario") == 79 && fabs(microenvironment.voxels(microenvironment.nearest_voxel_index(position)).center[0]-100.0) <= 25.0)
			{
				continue;
			}
			else if (parameters.ints("scenario") == 80 && fabs(microenvironment.voxels(microenvironment.nearest_voxel_index(position)).center[1]) <= 5.0)
			{
				continue;
			}
			else if (parameters.ints("scenario") == 81 && fabs(microenvironment.voxels(microenvironment.nearest_voxel_index(position)).center[1]) <= 25.0)
			{
				continue;
			}
			#pragma omp critical (initialize_cells)
			{
				pC = create_cell(fibronectin_puncta);
				pC->assign_position(position);
				int passed_over_index = pC->custom_data.find_variable_index("time_since_last_filopodia_drop");
				int puncta_angle_index = pC->custom_data.find_variable_index("time_until_next_filopodia_drop");
				if(parameters.ints("scenario") == 8)
				{
					double theta = 0;
					pC->custom_data.variables[passed_over_index].value = 1; // For posterity's sake
					pC->custom_data.variables[puncta_angle_index].value = theta;
				}
				else if(parameters.ints("scenario") == 9)
				{
					double theta = PhysiCell_constants::pi/2;
					pC->custom_data.variables[passed_over_index].value = 1; // For posterity's sake
					pC->custom_data.variables[puncta_angle_index].value = theta;
				}
			}
		}
	}
	return;
}

void apply_reflective_boundary_condition(Cell* pCell, Phenotype& phenotype, double dt)
{
	// NOTE: In the main function, you already screen for leaders and followers.
	// If you're still in the domain and movable, do nothing:
	if(pCell->is_movable && !(pCell->is_out_of_domain))
	{
		return;
	}

	// Else, figure out which boundary the cell has crossed:
	static double Xmin  = microenvironment.mesh.bounding_box[0];
	static double Xmax = microenvironment.mesh.bounding_box[3];
	static double Ymin  = microenvironment.mesh.bounding_box[1];
	static double Ymax = microenvironment.mesh.bounding_box[4];

	std::vector<double> new_position = pCell->position;
	// Check to see where it has left:
	// Left-right:
	if (pCell->position[0] < Xmin)
	{
		new_position[0] = Xmin + (Xmin - pCell->position[0]);
		phenotype.motility.migration_bias_direction[0] *= -1;
		phenotype.motility.motility_vector[0] *= -1;
	}
	else if (pCell->position[0] > Xmax)
	{
		new_position[0] = Xmax - (pCell->position[0] - Xmax);
		phenotype.motility.migration_bias_direction[0] *= -1;
		phenotype.motility.motility_vector[0] *= -1;
	}
	// Top-bottom
	if (pCell->position[1]  < Ymin)
	{
		new_position[1] = Ymin + (Ymin-pCell->position[1]);
		phenotype.motility.migration_bias_direction[1] *= -1;
		phenotype.motility.motility_vector[1] *= -1;
	}
	else if (pCell->position[1] > Ymax)
	{
		new_position[1] = Ymax - (pCell->position[1] - Ymax);
		phenotype.motility.migration_bias_direction[1] *= -1;
		phenotype.motility.motility_vector[1] *= -1;
	}
	#pragma omp critical(remove_agent_from_vectors)
	{
		// Adding the pragma omp critical to these two lines fixes the segfault.
		// I guess what was going on was that you were fiddling with the vector
		// sizes of the agent grids in the parallel loop, and this sometimes
		// caused a crash when you were trying to add/delete cells simulateously.
		pCell->get_container()->remove_agent(pCell); // This ensures you don't have duplicate cells in the agent grid vector.
		pCell->assign_position(new_position); // This is where bug is occuring! (Specific location: register_agent() in assign_position.
	}
	pCell->is_active = true; // controls secretion/uptake
	pCell->is_movable = true;
	pCell->is_out_of_domain = false;

	return;
}

void record_total_velocity(Cell* pCell, Phenotype& phenotype, double dt)
{
	// NOTE: In cell definitions, you have already screened out the edge agents and FN puncta.

	// If you're out of the domain/not movable, do nothing:
	if(pCell->is_out_of_domain || !(pCell->is_movable) )
	{
		return;
	}
	// If you're in the domain, record your total velocity vector
	static int total_velocity_index = pCell->custom_data.find_vector_variable_index("total_velocity");
	pCell->custom_data.vector_variables[total_velocity_index].value = pCell->velocity;

	// Set the migration bias direction to be along current trajectory, if the velocity is large enough:
	if(norm(pCell->velocity) > 1e-16)
	{
		phenotype.motility.migration_bias_direction = pCell->velocity;
		normalize( &(phenotype.motility.migration_bias_direction) );
	}
	return;
}

void add_in_new_cells(double dt)
{
	double probability_of_cell_entrance = dt*parameters.doubles("rate_of_cell_entrance");
	if(UniformRandom()<probability_of_cell_entrance)
	{
		static double X_left  = microenvironment.mesh.bounding_box[0];
		static double X_length = microenvironment.mesh.bounding_box[3] - X_left;
		static double Y_left  = microenvironment.mesh.bounding_box[1];
		static double Y_length = microenvironment.mesh.bounding_box[4] - Y_left;
		static double Y_center = Y_left + Y_length/2.0;

		std::vector<double> position(3,0.0);
		std::vector<double> orientation(3,0.0);
		position[0] = X_left + parameters.doubles("default_cell_radius");
		position[1] = Y_center + parameters.doubles("y_length_of_cell_entrance_strip")*( -0.5 + UniformRandom() );
		double angle = 2.0*PhysiCell_constants::pi*UniformRandom();
		orientation[0] = cos(angle);
		orientation[1] = sin(angle);
		Cell* pC = create_cell(follower_cell);
		pC->assign_position(position);
		pC->phenotype.motility.migration_bias_direction = orientation;
		static int time_until_next_filopodia_drop_index = pC->custom_data.find_variable_index("time_until_next_filopodia_drop");
		pC->custom_data.variables[time_until_next_filopodia_drop_index].value = -parameters.doubles("average_time_until_next_filopodia_drop")*log(UniformRandom());
		static int fibronectin_ID_index = pC->custom_data.find_variable_index("fibronectin_ID");
		pC->custom_data.variables[fibronectin_ID_index].value = 1.0*pC->ID;
	}
	return;
}

void add_in_new_cells_if_free_space( void )
{
	static int cover_index = (*all_cells)[0]->custom_data.find_variable_index("are_you_covered_by_a_cell");
	// Check boundary agents for cell covering:
	#pragma omp parallel for
	for(int i = 0; i < (*all_cells).size(); i++)
	{
		if((*all_cells)[i]->type != 3)
		{
			continue;
		}
		if((*all_cells)[i]->custom_data.variables[cover_index].value < 1e-14)
		{
			std::vector<double> position = (*all_cells)[i]->position;
			position[0] = microenvironment.mesh.bounding_box[0];
			double angle = 2*PhysiCell_constants::pi*UniformRandom();
			std::vector<double> orientation = {cos(angle), sin(angle), 0.0};
			#pragma omp critical (add_in_a_new_cell)
			{
				Cell* pC = create_cell(follower_cell);
				// Initial Velocity is drawn uniformly at random
				pC->assign_position(position);
				pC->phenotype.motility.migration_bias_direction = orientation;
				int time_until_next_filopodia_drop_index = pC->custom_data.find_variable_index("time_until_next_filopodia_drop");
				pC->custom_data.variables[time_until_next_filopodia_drop_index].value = -parameters.doubles("average_time_until_next_filopodia_drop")*log(UniformRandom());
				static int fibronectin_ID_index = pC->custom_data.find_variable_index("fibronectin_ID");
				pC->custom_data.variables[fibronectin_ID_index].value = 1.0*pC->ID;
			}
		}
	}

	return;
}

std::vector<double> sample_from_Von_Mises_distribution(double mu, double kappa)
{
		if(std::isinf(kappa))
		{
			std::vector<double> direction = {cos(mu),sin(mu),0};
			return direction;
		}
    // Sample from non-normalised density with rejection sampling:
		// Find the maximum value of the non-normalized VM distribution
		// 		non_normalised_von_mises_distribution = exp(kappa * cos(theta - mu));
    // double max_of_von_mises_distribution = non_normalised_von_mises_distribution(mu,mu,kappa);
		double max_of_von_mises_distribution = exp(kappa);

		// Rejection sampling:
    double theta = PhysiCell_constants::pi*( 2*UniformRandom() - 1 ); // theta in [-pi,pi]
    double pdf_value = exp(kappa*cos(theta-mu)); // non_normalised_von_mises_distribution(theta,mu,kappa);
    double u = UniformRandom()*max_of_von_mises_distribution;

    while (u > pdf_value)
		{
        theta = PhysiCell_constants::pi*( 2*UniformRandom() - 1 );
        pdf_value = exp(kappa*cos(theta-mu)); // non_normalised_von_mises_distribution(theta,mu,kappa);
        u = UniformRandom()*max_of_von_mises_distribution;
    }
		std::vector<double> sample = {cos(theta), sin(theta), 0.0};
    return sample;
}

void add_new_fibronectin_puncta(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	static int time_until_next_filopodia_drop_index = pCell->custom_data.find_variable_index("time_until_next_filopodia_drop");
	static int time_since_last_filopodia_drop_index = pCell->custom_data.find_variable_index("time_since_last_filopodia_drop");
	static int fibronectin_index = microenvironment.find_density_index("fibronectin");
	static int fibronectin_ID_index = pCell->custom_data.find_variable_index("fibronectin_ID");

	if( (pCell->custom_data.variables[time_since_last_filopodia_drop_index].value > pCell->custom_data.variables[time_until_next_filopodia_drop_index].value) && (pCell->custom_data.variables[did_you_sense_fibronectin_index].value > 1e-14) )
	{
		// Drop a new filopodia puncta:
		#pragma omp critical(remove_agent_from_vectors)
		{
			Cell* pC = create_cell(fibronectin_puncta);
			pC->assign_position(pCell->position);
			pC->custom_data.variables[fibronectin_ID_index].value = 1.0*pCell->ID;
		}
		int voxel_index = microenvironment.nearest_voxel_index(pCell->position);
		microenvironment.density_vector(voxel_index)[fibronectin_index] = 30.0; // This is the initial FN concentration at puncta sites.
		// Draw a new time from an exponential distribution:
		pCell->custom_data.variables[time_until_next_filopodia_drop_index].value = -parameters.doubles("average_time_until_next_filopodia_drop")*log(UniformRandom());
		pCell->custom_data.variables[time_since_last_filopodia_drop_index].value = 0;
		return;
	}

	pCell->custom_data.variables[time_since_last_filopodia_drop_index].value += dt;
	return;
}

void add_new_fibronectin_puncta_followers(Cell* pCell, Phenotype& phenotype, double dt, double average_interval)
{
	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	static int time_until_next_filopodia_drop_index = pCell->custom_data.find_variable_index("time_until_next_filopodia_drop");
	static int time_since_last_filopodia_drop_index = pCell->custom_data.find_variable_index("time_since_last_filopodia_drop");
	static int fibronectin_index = microenvironment.find_density_index("fibronectin");
	static int fibronectin_ID_index = pCell->custom_data.find_variable_index("fibronectin_ID");

	if( (pCell->custom_data.variables[time_since_last_filopodia_drop_index].value > pCell->custom_data.variables[time_until_next_filopodia_drop_index].value) && (pCell->custom_data.variables[did_you_sense_fibronectin_index].value > 1e-14) )
	{
		// Drop a new filopodia puncta:
		#pragma omp critical(remove_agent_from_vectors)
		{
			Cell* pC = create_cell(fibronectin_puncta);
			pC->assign_position(pCell->position);
			pC->custom_data.variables[fibronectin_ID_index].value = 1.0*pCell->ID; // Should you record IDs for follower secretion? Probably not?
		}
		int voxel_index = microenvironment.nearest_voxel_index(pCell->position);
		microenvironment.density_vector(voxel_index)[fibronectin_index] = 30.0; // This is the initial FN concentration at puncta sites.
		// Draw a new time from an exponential distribution:
		pCell->custom_data.variables[time_until_next_filopodia_drop_index].value = -average_interval*log(UniformRandom());
		pCell->custom_data.variables[time_since_last_filopodia_drop_index].value = 0;
		return;
	}

	pCell->custom_data.variables[time_since_last_filopodia_drop_index].value += dt;
	return;
}

void check_puncta_for_cell_cover(Cell* pCell, Phenotype& phenotype, double dt)
{
	// This function checks if a puncta or edge_agent is covered by a cell so that you can use this information in the next time step.
	// It also resets the FN concentration at the puncta lattice site if the puncta is sensed and is not covered.
	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	static int cover_index = pCell->custom_data.find_variable_index("are_you_covered_by_a_cell");
	static int fibronectin_index = microenvironment.find_density_index("fibronectin");
	static int passed_over_index = pCell->custom_data.find_variable_index("time_since_last_filopodia_drop");
	static int puncta_angle_index = pCell->custom_data.find_variable_index("time_until_next_filopodia_drop");
	static int total_velocity_index = pCell->custom_data.find_vector_variable_index("total_velocity");
	static bool followers_can_update_FN_orientation = ((parameters.ints("scenario") != 6) && (parameters.ints("scenario") != 7) && (parameters.ints("scenario") != 20) && (parameters.ints("scenario") != 22) && (parameters.ints("scenario") != 23) && (parameters.ints("scenario") != 24) && (parameters.ints("scenario") != 33) && (parameters.ints("scenario") != 34) && (parameters.ints("scenario") != 35) && (parameters.ints("scenario") != 36) && (parameters.ints("scenario") != 37) && (parameters.ints("scenario") != 38) && (parameters.ints("scenario") != 39) && (parameters.ints("scenario") != 49) && (parameters.ints("scenario") != 50) && (parameters.ints("scenario") != 51) && (parameters.ints("scenario") != 52) && (parameters.ints("scenario") != 53) && (parameters.ints("scenario") != 54) && (parameters.ints("scenario") != 55) && (parameters.ints("scenario") != 56) && (parameters.ints("scenario") != 57) && (parameters.ints("scenario") != 58) && (parameters.ints("scenario") != 59) && (parameters.ints("scenario") != 60));
	static bool leaders_can_update_FN_orientation = ((parameters.ints("scenario") != 15) && (parameters.ints("scenario") != 16) && (parameters.ints("scenario") != 20) && (parameters.ints("scenario") != 22) && (parameters.ints("scenario") != 23) && (parameters.ints("scenario") != 25) && (parameters.ints("scenario") != 33) && (parameters.ints("scenario") != 34) && (parameters.ints("scenario") != 35) && (parameters.ints("scenario") != 36) && (parameters.ints("scenario") != 37) && (parameters.ints("scenario") != 38) && (parameters.ints("scenario") != 39) && (parameters.ints("scenario") != 49) && (parameters.ints("scenario") != 50) && (parameters.ints("scenario") != 51) && (parameters.ints("scenario") != 52) && (parameters.ints("scenario") != 53) && (parameters.ints("scenario") != 54) && (parameters.ints("scenario") != 55) && (parameters.ints("scenario") != 56) && (parameters.ints("scenario") != 57) && (parameters.ints("scenario") != 58) && (parameters.ints("scenario") != 59) && (parameters.ints("scenario") != 60));

	// First determine if you are covered by a cell:
	pCell->custom_data.variables[cover_index].value = 0;

	// Check all cells within current and neighboring mechanics voxels:
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// Other puncta / edge_agents don't count
		if((*neighbor)->type == 2 || (*neighbor)->type == 3)
		{
			continue;
		}
		std::vector<double> displacement = (*neighbor)->position - pCell->position;
		double distance_to_cell = norm(displacement);
		if((pCell->type == 2 && distance_to_cell <= ((*neighbor)->phenotype.geometry.nuclear_radius)) || (pCell->type == 3 && distance_to_cell < (phenotype.geometry.nuclear_radius + (*neighbor)->phenotype.geometry.nuclear_radius)))
		{
			// You can't be sensed by a filopodia if you are covered
			pCell->custom_data.variables[did_you_sense_fibronectin_index].value = 0;
			pCell->custom_data.variables[cover_index].value = 1;

			if(((*neighbor)->type == 1) && !followers_can_update_FN_orientation)
			{
				return; // Don't update the FN orientations if you aren't in the right scenario
			}
			if (((*neighbor)->type == 0) && !leaders_can_update_FN_orientation)
			{
				return; // Don't update the FN orientations if you aren't in the right scenario
			}
			// Store/update the angle of the fibronectin at the puncta location:
			double theta = atan2((*neighbor)->custom_data.vector_variables[total_velocity_index].value[1], (*neighbor)->custom_data.vector_variables[total_velocity_index].value[0]);
			if(pCell->custom_data.variables[passed_over_index].value < 1e-14)
			{
				// You've just been passed over for the first time --
				// 		Record the direction in which the leader/follower cell is going:
				pCell->custom_data.variables[passed_over_index].value = 1; // For posterity's sake
				pCell->custom_data.variables[puncta_angle_index].value = theta;
			}
			else
			{
				// Assume half-life of difference between angles = T_half
				pCell->custom_data.variables[puncta_angle_index].value += dt*log(2)/parameters.doubles("half_life_of_puncta_orientation_change")*sin(theta - pCell->custom_data.variables[puncta_angle_index].value);
			}
			// Stop the function immediately once you know you are covered:
			return;
		}
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end;
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			if((*neighbor)->type == 2 || (*neighbor)->type == 3)
			{
				continue;
			}
			std::vector<double> displacement = (*neighbor)->position - pCell->position;
			double distance_to_cell = norm(displacement);
			if((pCell->type == 2 && distance_to_cell <= ((*neighbor)->phenotype.geometry.nuclear_radius)) || (pCell->type == 3 && distance_to_cell < (phenotype.geometry.nuclear_radius + (*neighbor)->phenotype.geometry.nuclear_radius)))
			{
				pCell->custom_data.variables[did_you_sense_fibronectin_index].value = 0;
				pCell->custom_data.variables[cover_index].value = 1;

				if(((*neighbor)->type == 1) && !followers_can_update_FN_orientation)
				{
					return; // Don't update the FN orientations if you aren't in the right scenario
				}
				if (((*neighbor)->type == 0) && !leaders_can_update_FN_orientation)
				{
					return; // Don't update the FN orientations if you aren't in the right scenario
				}
				// Store/update the angle of the fibronectin at the puncta location:
				double theta = atan2((*neighbor)->custom_data.vector_variables[total_velocity_index].value[1], (*neighbor)->custom_data.vector_variables[total_velocity_index].value[0]);
				if(pCell->custom_data.variables[passed_over_index].value < 1e-14)
				{
					// You've just been passed over for the first time --
					// 		Record the direction in which the leader/follower cell is going:
					pCell->custom_data.variables[passed_over_index].value = 1; // For posterity's sake
					pCell->custom_data.variables[puncta_angle_index].value = theta;
				}
				else
				{
					// Assume half-life of difference between angles = 5 min:
					pCell->custom_data.variables[puncta_angle_index].value += dt*log(2)/parameters.doubles("half_life_of_puncta_orientation_change")*sin(theta - pCell->custom_data.variables[puncta_angle_index].value);
				}
				// Stop the function immediately once you know you are covered:
				return;
			}
		}
	}

	if(pCell->type == 2 && pCell->custom_data.variables[did_you_sense_fibronectin_index].value > 1e-14)
	{
		// Secrete FN into microenvironment:
		int voxel_index = microenvironment.nearest_voxel_index(pCell->position);
		microenvironment.density_vector(voxel_index)[fibronectin_index] = 30.0;
	}
	// For all puncta, reset the did_you_sense_fibronectin variable to FALSE:
	pCell->custom_data.variables[did_you_sense_fibronectin_index].value = 0;
	return;
}

void calculate_surface_gradients(Cell* pCell, Cell* neighbor, Phenotype& phenotype, std::vector<double>& FN_gradient, std::vector<double>& direction_from_puncta_angle, std::vector<double>& contact_surface_gradient, std::vector<double>& spring_force_PhysiCell, std::vector<double>& spring_force_FN, int& number_of_fibronectin_neighbors, int& number_of_contact_neighbors, int& number_of_puncta_contributors)
{
	static int fibronectin_index = microenvironment.find_density_index("fibronectin");
	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	static int cover_index = pCell->custom_data.find_variable_index("are_you_covered_by_a_cell");
	static int passed_over_index = pCell->custom_data.find_variable_index("time_since_last_filopodia_drop");
	static int puncta_angle_index = pCell->custom_data.find_variable_index("time_until_next_filopodia_drop");
	static bool have_larger_follower_FN_surface_height = (parameters.ints("scenario") == 44 || parameters.ints("scenario") == 45 || parameters.ints("scenario") == 46 || parameters.ints("scenario") == 47 || parameters.ints("scenario") == 53 || parameters.ints("scenario") == 54 || parameters.ints("scenario") == 55 || parameters.ints("scenario") == 56 || parameters.ints("scenario") == 59 || parameters.ints("scenario") == 74 || parameters.ints("scenario") == 75 || parameters.ints("scenario") == 76 || parameters.ints("scenario") == 77);

	int neighbor_type = neighbor->type;
	std::vector<double> displacement = (neighbor->position);
	double neighbor_cell_radius = neighbor->phenotype.geometry.nuclear_radius;
	bool puncta_is_covered = (neighbor->custom_data.variables[cover_index].value > 0.5);
	bool neighbor_is_pCell = (neighbor == pCell);

	// Don't count yourself or edge agents:
	if( neighbor_type == 3 || neighbor_is_pCell )
	{
		return;
	}

	int voxel_index = microenvironment.nearest_voxel_index(displacement);
	double fibronectin_concentration = microenvironment.density_vector(voxel_index)[fibronectin_index];

	static double one = 1.0;
	naxpy( &(displacement), one, pCell->position);

	double distance = norm(displacement);
	double distance_squared = norm_squared(displacement);
	bool can_sense_fibronectin = true;
	std::vector<double> temp_gradient(3,0.0);
	double FN_sigma_squared = pow(phenotype.geometry.radius-phenotype.geometry.nuclear_radius,2);
	// double FN_sigma_squared = pow(2.0*((*neighbor)->phenotype.geometry.radius),2);
	double contact_sigma_squared = pow(neighbor_cell_radius,2);
	// double contact_sigma_squared = pow(2.0*((*neighbor)->phenotype.geometry.nuclear_radius),2);

	if(neighbor_type == 2 && parameters.bools("sense_fibronectin_along_movement_direction"))
	{
		double scalar_product = std::inner_product(displacement.begin(), displacement.end(), pCell->velocity.begin(), 0.0);
		if(scalar_product < 0)
		{
			can_sense_fibronectin = false;
		}
	}
	// Recall that the "cell radius in Physicell" is actually the fibronectin sensing radius, and the nuclear radius is the actual cell!
	if((neighbor_type == 2) && (distance <= phenotype.geometry.radius) && can_sense_fibronectin && (!puncta_is_covered) && (fibronectin_concentration >= parameters.doubles("fibronectin_speed_threshold")))
	{
		// Calculate the velocity you would have gotten from the built in PhysiCell cell-cell potential:
		double temp_r = -distance; // -d
		temp_r /= parameters.doubles("filopodia_sensing_radius"); // -d/R_filo
		temp_r += 1.0; // 1- d/R_filo
		temp_r *= temp_r; // (1-d/R)^2
		// Add on to the spring force
		temp_r /= distance; // This is for the unit vector of displacement
		axpy(&(spring_force_FN), temp_r, displacement);

		number_of_fibronectin_neighbors++;
		double gamma_FN = parameters.doubles("fibronectin_kernel_max_value");
		if (pCell->type == 1 && have_larger_follower_FN_surface_height)
		{
			gamma_FN *= 10; // Set the height of the gaussian kernel to be 10x higher when you are conducting the relevant rescue expts that only affect follower cells
		}
		temp_gradient[0] = -fibronectin_concentration*gamma_FN;
		temp_gradient[0] /= FN_sigma_squared;
		temp_gradient[0] *= exp(-distance_squared/2.0/FN_sigma_squared);
		temp_gradient[1] = temp_gradient[0];
		temp_gradient[0] *= displacement[0]; // Note: this gives the negative gradient! (Positive is (x_pCell - x_i))
		temp_gradient[1] *= displacement[1]; // Note: this gives the negative gradient! (Positive is (y_pCell - y_i))
		// FN_gradient -= temp_gradient;
		naxpy( &(FN_gradient), one, temp_gradient);

		// Set state variable of the puncta to be TRUE (state variable = you have been sensed by leader cell. Use this in case you want to use secretion methods to calculate FN increase)
		neighbor->custom_data.variables[did_you_sense_fibronectin_index].value = 1;
		// Set state variable of the cell to be TRUE (state variable = are you sensing FN puncta?)
		pCell->custom_data.variables[did_you_sense_fibronectin_index].value = 1;

		bool does_the_puncta_have_an_angle = false;
		double puncta_angle = 0.0;
		does_the_puncta_have_an_angle = (neighbor->custom_data.variables[passed_over_index].value > 0.5);
		puncta_angle = neighbor->custom_data.variables[puncta_angle_index].value;

		// If the puncta has already been passed over, record the puncta angle in your vector:
		if(does_the_puncta_have_an_angle)
		{
			number_of_puncta_contributors++;
			direction_from_puncta_angle[0] += cos(puncta_angle);
			direction_from_puncta_angle[1] += sin(puncta_angle);
		}
		// Move to the next cell, since you already know that this neighbor is a puncta:
		return;
	} // end if the puncta is sensed by the cell

	if( (neighbor_type != 2) && (distance <= phenotype.geometry.radius) )
	{
		double temp_r = -distance; // -d
		double R_min = parameters.doubles("filopodia_sensing_radius");

		temp_r /= R_min; // -d/R_filo
		temp_r += 1.0; // 1- d/R_filo
		temp_r *= temp_r; // (1-d/R)^2
		// Have to multiply by -1 because force is repulsive along this displacement vector (force is from neighbor to pCell, not other way around)
		temp_r *= -sqrt(phenotype.mechanics.cell_cell_repulsion_strength * neighbor->phenotype.mechanics.cell_cell_repulsion_strength);
		// Also multiply by the singularity, if it occurs
		temp_r *= pow(parameters.doubles("scaling_of_singularity"), parameters.doubles("power_of_repulsive_singularity"))/pow(distance, parameters.doubles("power_of_repulsive_singularity")); // This is NOT for the unit vector!
		// Add on to the spring force
		temp_r /= distance; // This is for the unit vector of displacement
		// spring_force_PhysiCell += temp_r*displacement;
		axpy(&(spring_force_PhysiCell), temp_r, displacement);

		// Compute the contact_surface_gradient:
		number_of_contact_neighbors++;
		// std::vector<double> contact_gradient = {0,0,0};
		temp_gradient[0] = -parameters.doubles("contact_guidance_kernel_max_value");
		temp_gradient[0] /= contact_sigma_squared;
		temp_gradient[0] *= exp(-distance_squared/2.0/contact_sigma_squared);
		temp_gradient[1] = temp_gradient[0];
		temp_gradient[0] *= displacement[0]; // Note: this gives the negative gradient! (Positive is (x_pCell - x_i))
		temp_gradient[1] *= displacement[1]; // Note: this gives the negative gradient! (Positive is (y_pCell - y_i))
		// contact_surface_gradient += temp_gradient;
		axpy( &(contact_surface_gradient), one, temp_gradient );
	}// end if you had a cell-cell contact

	return;
}

void compute_directions(Cell* pCell, Phenotype& phenotype, double dt, std::vector<double>& motility_direction, std::vector<double>& contact_velocity, int& number_of_contact_neighbors)
{
	std::vector<double> FN_gradient(3,0.0);
	std::vector<double> direction_from_puncta_angle(3,0.0);
	std::vector<double> contact_surface_gradient(3,0.0);
	std::vector<double> spring_force_PhysiCell(3,0.0);
	std::vector<double> spring_force_FN(3,0.0);
	int number_of_fibronectin_neighbors = 0;
	int number_of_puncta_contributors = 0;

	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	static int volume_exclusion_direction_index = pCell->custom_data.find_vector_variable_index("volume_exclusion_direction");
	static int contact_guidance_direction_index = pCell->custom_data.find_vector_variable_index("contact_guidance_direction");
	static int haptotaxis_direction_index = pCell->custom_data.find_vector_variable_index("haptotaxis_direction");
	static int b_mot_index = pCell->custom_data.find_variable_index("b_mot");
	static int b_vol_exclusion_index = pCell->custom_data.find_variable_index("b_vol_exclusion");
	static bool follower_has_different_ECM_bias = (parameters.ints("scenario") == 40 || parameters.ints("scenario") == 41 || parameters.ints("scenario") == 42 || parameters.ints("scenario") == 43 || parameters.ints("scenario") == 49 || parameters.ints("scenario") == 50 || parameters.ints("scenario") == 51 || parameters.ints("scenario") == 52 || parameters.ints("scenario") == 58 || parameters.ints("scenario") == 64 || parameters.ints("scenario") == 65 || parameters.ints("scenario") == 66 || parameters.ints("scenario") == 70 || parameters.ints("scenario") == 74);

	// Set cell state variable did_you_sense_fibronectin to be FALSE (haven't sensed puncta yet)
	pCell->custom_data.variables[did_you_sense_fibronectin_index].value = 0;

	// Use the filopodia_puncta/"cells" to efficiently sample FN concentrations.
	//First check the neighbors in the current mechanics voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		calculate_surface_gradients(pCell, (*neighbor), phenotype, FN_gradient, direction_from_puncta_angle, contact_surface_gradient, spring_force_PhysiCell, spring_force_FN, number_of_fibronectin_neighbors, number_of_contact_neighbors, number_of_puncta_contributors);
	}

	// // Now check cells within Moore neighborhood of the mechanics grid:
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor_voxel_index =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end;
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
		{
			continue;
		}
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			calculate_surface_gradients(pCell, (*neighbor), phenotype, FN_gradient, direction_from_puncta_angle, contact_surface_gradient, spring_force_PhysiCell, spring_force_FN, number_of_fibronectin_neighbors, number_of_contact_neighbors, number_of_puncta_contributors);
		}
	}
	std::vector<double> FN_direction(3,0.0);
	if(number_of_fibronectin_neighbors == 0)
	{
		// If no neighbors, move towards current velocity:
		FN_direction = phenotype.motility.migration_bias_direction;
	}
	else
	{
		double mu_FN = atan2(FN_gradient[1], FN_gradient[0]);
		double kappa_FN = norm(FN_gradient);
		FN_direction = sample_from_Von_Mises_distribution(mu_FN, kappa_FN);
	}
	pCell->custom_data.vector_variables[haptotaxis_direction_index].value = FN_direction;
	if(number_of_puncta_contributors != 0)
	{
		// Find the average direction (so far you have sum of cos(phi_i) and sin(phi_i))
		direction_from_puncta_angle *= 1.0/number_of_puncta_contributors;
		normalize( &(direction_from_puncta_angle) );
		pCell->custom_data.vector_variables[contact_guidance_direction_index].value = direction_from_puncta_angle;
		double p1 = parameters.doubles("bias_towards_VM_direction");
		if (pCell->type == 1 && follower_has_different_ECM_bias)
		{
			p1 = 0.33; // Set rho to be equal to 0.33 if you are a follower cell and you are conducting rescue expts that only affect followers
		}
		double one_minus_p1 = 1-p1;
		// motility_direction = p1*FN_direction + (1-p1)*direction_from_puncta_angle;
		axpy( &(motility_direction), p1, FN_direction);
		axpy( &(motility_direction), one_minus_p1, direction_from_puncta_angle);
		normalize( &(motility_direction) );
	}
	else
	{
		pCell->custom_data.vector_variables[contact_guidance_direction_index].value = direction_from_puncta_angle; // This vector is just a zero vector
		static double one = 1.0;
		axpy(&(motility_direction), one, FN_direction);
	}
	double speed_spring_force_FN = norm(spring_force_FN); // NOTE: if no FN puncta in neighborhood, this is equal to 0
	pCell->custom_data.variables[b_mot_index].value = speed_spring_force_FN;
	motility_direction *= speed_spring_force_FN;
	// std::cout<<"Done!\n\n";
	// std::cout<<"Computing contact direction...\n\n";
	if(number_of_contact_neighbors != 0)
	{
		double mu_contact = atan2(contact_surface_gradient[1], contact_surface_gradient[0]);
		double kappa_contact = norm(contact_surface_gradient);
		// With height of kernel = 100, kappa is about 2-3

		// Random sampling biased toward direction of steepest descent:
		contact_velocity = sample_from_Von_Mises_distribution(mu_contact, kappa_contact);

		// Cell-cell strength = 5.0: contact speed about 0.5 um/min
		// Cell-cell strength = 1.0: contact speed about 0.05 um/min
		// double contact_speed = norm(spring_force_PhysiCell); // ORIGINAL FORCE!
		double dot_prod = std::inner_product(contact_velocity.begin(), contact_velocity.end(), spring_force_PhysiCell.begin(), 0.0);
		double contact_speed = norm(spring_force_PhysiCell)*( (dot_prod >0) - (dot_prod < 0)); // Take the sign of the dot product to figure out where the force points
		pCell->custom_data.variables[b_vol_exclusion_index].value = contact_speed;
		pCell->custom_data.vector_variables[volume_exclusion_direction_index].value = contact_velocity;
		// contact speed is magnitude of spring_force_PhysiCell dotted with VM unit vector (|Physicell_speed*cos(angle_between_VM_direction_and_Phys)|)
		// double scalar_product = std::inner_product(contact_velocity.begin(), contact_velocity.end(), spring_force_PhysiCell.begin(), 0.0);
		// double contact_speed = fabs(scalar_product);

		contact_velocity *= contact_speed;
	}
	else
	{
		pCell->custom_data.variables[b_vol_exclusion_index].value = 0;
		pCell->custom_data.vector_variables[volume_exclusion_direction_index].value = contact_velocity;
	}
	return;
}

void simpler_function_leader_follower_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt)
{
	/* DETERMINE NEW VELOCITIES */
	// 1. Compute bias direction from fibronectin, velocity from cell-cell contact:
	std::vector<double> motility_direction(3,0.0); // NOT A UNIT VECTOR (will use to compute resultant velocities)
	std::vector<double> contact_velocity(3,0.0); // NOT A UNIT VECTOR
	int number_of_contact_neighbors = 0;
	compute_directions( pCell, phenotype, dt, motility_direction, contact_velocity, number_of_contact_neighbors);
	std::vector<double> FN_direction = motility_direction;
	// If you don't have any cell-cell contact and you would get a zero motility direction after normalizing, maintain your current course.
	if(norm(motility_direction) <= 1e-16 && number_of_contact_neighbors == 0)
	{
		// normalize(&(phenotype.motility.migration_bias_direction));
		motility_direction = phenotype.motility.migration_bias_direction;
		FN_direction = phenotype.motility.migration_bias_direction;
	}

	// Leader cells don't change their path if you don't have them moving randomly:
	if(pCell->type == 0 && !(parameters.bools("leader_cell_initially_has_random_orientation")))
	{
		pCell->velocity = phenotype.motility.migration_speed*phenotype.motility.migration_bias_direction;
		return;
	}

	// 2. Update the migration bias and migration_speed:
			// Always update the speed, even if you have nonzero persistence_time:
			// If you didn't sense FN, then set your velocity to be the cell speed off fibronectin:
	static int did_you_sense_fibronectin_index = pCell->custom_data.find_variable_index("did_you_sense_fibronectin");
	if (pCell->custom_data.variables[did_you_sense_fibronectin_index].value < 1e-14)
	{
		phenotype.motility.migration_speed = parameters.doubles("cell_speed_off_fibronectin");
	}
	else
	{
		phenotype.motility.migration_speed = parameters.doubles("default_cell_speed");
	}
	// Check to see if you actually are changing your motility vector/Fibronectin velocity:
	// 		NOTE: dt = mechanics_dt (default = 0.1)
	if( UniformRandom() < dt / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt )
	{
		phenotype.motility.migration_bias_direction = FN_direction;
		// std::cout<<"New migration_bias_direction = "<<phenotype.motility.migration_bias_direction<<"\n"<<"\n";

		phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motility = FN_bias

		phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
	}
	else
	{
		double bias = norm(motility_direction);
		motility_direction = phenotype.motility.migration_bias_direction;
		motility_direction *= bias;
		phenotype.motility.motility_vector = phenotype.motility.migration_speed*phenotype.motility.migration_bias_direction;
	}


	// 3. Compute the resultant cell velocity:
	static double one = 1.0;
	std::vector<double> sum_of_vectors(3,0.0);
	axpy( &(sum_of_vectors), one, motility_direction);
	axpy( &(sum_of_vectors), one, contact_velocity);
	pCell->velocity = sum_of_vectors; // If you don't have contact neighbors and didn't sense FN, then the norm of this vector is 1 > default_cell_speed!
	double speed = norm(pCell->velocity);

	// If you have a global signal directing leaders/followers to the right site, apply it now
	bool implement_global_signal_leaders = (parameters.bools("leaders_have_global_signal") && pCell->type == 0);
	bool implement_global_signal_followers = (parameters.bools("followers_have_global_signal") && pCell->type == 1);
	if( implement_global_signal_leaders || implement_global_signal_followers )
	{
		normalize( &(pCell->velocity));
		static double z = parameters.doubles("weight_towards_chemotaxis_cues");
		static std::vector<double> right = {1.0, 0.0, 0.0};
		pCell->velocity *= (1.0-z);
		axpy( &(pCell->velocity), z, right); // v_final = (1-z)*v + z*right
		normalize( &(pCell->velocity));
		pCell->velocity *= speed;
	}

	if( (pCell->custom_data.variables[did_you_sense_fibronectin_index].value < 1e-14) && number_of_contact_neighbors == 0)
	{
		normalize( &(pCell->velocity));
		pCell->velocity *= parameters.doubles("cell_speed_off_fibronectin");
	}
	// std::cout<<"New velocity = "<<sum_of_vectors<<"\n"<<"\n";


	// Speed is always less than or equal to the maximum cell speed possible,
	// 	which is the cell speed on fibronectin:
	if(speed > parameters.doubles("default_cell_speed") + 1e-14)
	{
		speed = parameters.doubles("default_cell_speed");
		normalize( &(pCell->velocity) );
		pCell->velocity *= speed;
	}

	return;
}

void inject_new_fibronectin( double left_most_position, double right_most_position )
{
	// double X_left = microenvironment.mesh.bounding_box[0];
	double Y_left = microenvironment.mesh.bounding_box[1];
	double Y_length = (microenvironment.mesh.bounding_box[4]-Y_left);
	double injected_fibronectin_spacing = 5; // microns
	int fibronectin_index = microenvironment.find_density_index("fibronectin");
	int number_of_fibronectin_rows = int(Y_length/injected_fibronectin_spacing);
	std::vector<double> position(3,0.0);
	// find nearest voxel to the current position
	for(int n = 0; n<number_of_fibronectin_rows; n++)
	{
		position[1] = Y_left + (n+0.5)*injected_fibronectin_spacing;
		double left_most_puncta = left_most_position + 0.5*injected_fibronectin_spacing;
		position[0] = left_most_puncta;
		int k = 0;
		while(position[0] <= right_most_position)
		{
			// position[0] = X_left + (k + 0.5)*parameters.doubles("fibronectin_puncta_wavelength");
			position[0] = left_most_puncta + k*injected_fibronectin_spacing;
			int i = microenvironment.nearest_voxel_index(position);
			microenvironment.density_vector(i)[fibronectin_index] = 30.0;
			if(parameters.bools("fibronectin_strip_initial_condition") && fabs(microenvironment.voxels(i).center[1]) > parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
				k++;
				continue;
			}
			Cell* pC = create_cell(fibronectin_puncta);
			pC->assign_position(position);
			k++;
		}
	}
	return;
}

void inject_new_fibronectin_outside_corridor( double right_most_position )
{
	double X_left = microenvironment.mesh.bounding_box[0];
	double Y_left = microenvironment.mesh.bounding_box[1];
	double Y_length = (microenvironment.mesh.bounding_box[4]-Y_left);
	double injected_fibronectin_spacing = 5; // microns
	int fibronectin_index = microenvironment.find_density_index("fibronectin");
	int number_of_fibronectin_rows = int(Y_length/injected_fibronectin_spacing);
	std::vector<double> position(3,0.0);
	// find nearest voxel to the current position
	for(int n = 0; n<number_of_fibronectin_rows; n++)
	{
		position[1] = Y_left + (n+0.5)*injected_fibronectin_spacing;
		if(fabs(position[1]) <= parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
		{
			continue;
		}
		double left_most_puncta = X_left + 0.5*injected_fibronectin_spacing;
		position[0] = left_most_puncta;
		int k = 0;
		while(position[0] <= right_most_position)
		{
			// position[0] = X_left + (k + 0.5)*parameters.doubles("fibronectin_puncta_wavelength");
			position[0] = left_most_puncta + k*injected_fibronectin_spacing;
			int i = microenvironment.nearest_voxel_index(position);
			microenvironment.density_vector(i)[fibronectin_index] = 30.0;
			if(parameters.bools("fibronectin_strip_initial_condition") && fabs(microenvironment.voxels(i).center[1]) > parameters.doubles("y_length_of_cell_entrance_strip")/2.0)
			{
				microenvironment.density_vector(i)[fibronectin_index] = 0.0;
				continue;
			}
			Cell* pC = create_cell(fibronectin_puncta);
			pC->assign_position(position);
			k++;
		}
	}
	return;
}

void delete_puncta( int scenario, std::vector<double> right_most_cell_position)
{
	int fibronectin_index = microenvironment.find_density_index("fibronectin");
	std::vector<int> indices_of_dead_puncta;

	// Loop over all cells/puncta, and record indices of all puncta to be deleted:
	for (int j = 0; j < (*all_cells).size(); j++)
	{
		// Scenario 82: Delete all puncta immediately to the right of the right most cell location (just delete up to 20 microns/fibronectin spacing):
		if ((*all_cells)[j]->type == 2 && scenario == 82 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) >= 0.0 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) < parameters.doubles("fibronectin_puncta_wavelength"))
		{
			int i = microenvironment.nearest_voxel_index((*all_cells)[j]->position);
			microenvironment.density_vector(i)[fibronectin_index] = 0.0; // Reset microenvironment to be zero as well
			indices_of_dead_puncta.push_back(j); // Record index of the punctum to be deleted
		}
		// Scenario 82: Delete all puncta immediately to the right of the right most cell location (just delete up until 50 microns):
		else if ((*all_cells)[j]->type == 2 && scenario == 83 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) >= 0.0 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) < 50.0)
		{
			int i = microenvironment.nearest_voxel_index((*all_cells)[j]->position);
			microenvironment.density_vector(i)[fibronectin_index] = 0.0; // Reset microenvironment to be zero as well
			indices_of_dead_puncta.push_back(j); // Record index of the punctum to be deleted
		}
		// Scenario 84: Delete all puncta along the y = 0 axis (plus or minus 5 microns) located to the right of the furthest cell:
		else if ((*all_cells)[j]->type == 2 && scenario == 84 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) >= 0.0 && fabs((*all_cells)[j]->position[1]) < 5.0)
		{
			int i = microenvironment.nearest_voxel_index((*all_cells)[j]->position);
			microenvironment.density_vector(i)[fibronectin_index] = 0.0; // Reset microenvironment to be zero as well
			indices_of_dead_puncta.push_back(j); // Record index of the punctum to be deleted
		}
		// Scenario 85: Delete all puncta along the y = 0 axis (plus or minus 25 microns) located to the right of the furthest cell:
		else if ((*all_cells)[j]->type == 2 && scenario == 85 && ((*all_cells)[j]->position[0] - right_most_cell_position[0]) >= 0.0 && fabs((*all_cells)[j]->position[1]) < 25.0)
		{
			int i = microenvironment.nearest_voxel_index((*all_cells)[j]->position);
			microenvironment.density_vector(i)[fibronectin_index] = 0.0; // Reset microenvironment to be zero as well
			indices_of_dead_puncta.push_back(j); // Record index of the punctum to be deleted
		}
	}

	// Now delete the puncta in order of descending index, so you don't accidentally remove a leader/follower:
	std::sort(indices_of_dead_puncta.rbegin(), indices_of_dead_puncta.rend());
	#pragma omp critical (delete_fibronectin)
	{
		for (int k = 0; k < indices_of_dead_puncta.size(); k++)
		{
			delete_cell(indices_of_dead_puncta[k]);
		}
	}
	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring

	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);

	if( pCell->phenotype.death.dead == false && pCell->type == 0 )
	{
		 output[0] = "black";
		 output[2] = "black";
	}
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "red";
		 output[2] = "red";
	}
	if( pCell->phenotype.death.dead == false && pCell->type == 2 )
	{
		 output[0] = "blue";
		 output[2] = "blue";
	}

	return output;
}
