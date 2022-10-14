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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"

// put custom code modules here!

#include "./custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)

	bool XML_status = false;
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// time setup
	std::string time_units = "min";

	/* Microenvironment setup */
	setup_microenvironment(); // modify this in the custom code

	/* PhysiCell setup */

	// set mechanics voxel size, and match the data structure to BioFVM
	// double mechanics_voxel_size = 30;
	// Maximum possible filopodia length = 30 um
	double mechanics_voxel_size = 50; // CAUTION: IF THE FILOPODIA SENSING RADIUS IS GREATER THAN mechanics_voxel_size, YOU MAY NOT COMPUTE ALL OF THE NEIGHBORS A CELL CAN SENSE
	// double mechanics_voxel_size = parameters.doubles("filopodia_sensing_radius");
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	/* Users typically start modifying here. START USERMODS */

	create_cell_types();

	setup_tissue();

	if (parameters.ints("scenario") == 11)
	{
		double extra_fibronectin_end_distance = 100; // microns
		inject_new_fibronectin(microenvironment.mesh.bounding_box[0], extra_fibronectin_end_distance);
	}

	if (parameters.ints("scenario") == 14)
	{
		double extra_fibronectin_end_distance = 150; // microns, in x-direction
		inject_new_fibronectin_outside_corridor(extra_fibronectin_end_distance);
	}


	/* Users typically stop modifying here. END USERMODS */

	// set MultiCellDS save options

	set_save_biofvm_mesh_as_matlab( true );
	set_save_biofvm_data_as_matlab( true );
	set_save_biofvm_cell_data( true );
	set_save_biofvm_cell_data_as_custom_matlab( true );

	// save a simulation snapshot

	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

	// save a quick SVG cross section through z = 0, after setting its
	// length bar to 200 microns

	PhysiCell_SVG_options.length_bar = 200;

	// for simplicity, set a pathology coloring function

	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;

	// sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
	// SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	display_citations();

	// set the performance timers

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() );

		report_file.open(filename); 	// create the data log file
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}

	// main loop
	double last_mechanics_time = 0.0;
	bool initialized = false;
	static double mechanics_tolerance = 0.001 * mechanics_dt;
	bool experiment_allows_FN_secretion = (parameters.ints("scenario") != 4 && parameters.ints("scenario") != 13 && parameters.ints("scenario") != 15 && parameters.ints("scenario") != 16 && parameters.ints("scenario") != 20 && parameters.ints("scenario") != 22 && parameters.ints("scenario") != 23 && parameters.ints("scenario") != 25 && parameters.ints("scenario") != 26 && parameters.ints("scenario") != 34 && parameters.ints("scenario") != 35 && parameters.ints("scenario") != 36 && parameters.ints("scenario") != 38 && parameters.ints("scenario") != 49 && parameters.ints("scenario") != 50 && parameters.ints("scenario") != 51 && parameters.ints("scenario") != 52 && parameters.ints("scenario") != 53 && parameters.ints("scenario") != 54 && parameters.ints("scenario") != 55 && parameters.ints("scenario") != 56 && parameters.ints("scenario") != 57);
	bool followers_can_secrete_FN = (parameters.ints("scenario") == 16 || parameters.ints("scenario") == 17 || parameters.ints("scenario") == 18 || parameters.ints("scenario") == 19 || parameters.ints("scenario") == 21 || parameters.ints("scenario") == 25 || parameters.ints("scenario") == 26 || parameters.ints("scenario") == 27 || parameters.ints("scenario") == 29 || parameters.ints("scenario") == 34 || parameters.ints("scenario") == 37 || parameters.ints("scenario") == 41 || parameters.ints("scenario") == 45 || parameters.ints("scenario") == 48 || parameters.ints("scenario") == 50 || parameters.ints("scenario") == 54 || parameters.ints("scenario") == 57 || parameters.ints("scenario") == 67 || parameters.ints("scenario") == 68 || parameters.ints("scenario") == 69);
	bool can_inject_new_FN_in_domain = (parameters.ints("scenario") == 5 || parameters.ints("scenario") == 7 || parameters.ints("scenario") == 13 || parameters.ints("scenario") == 22);
	bool inject_FN_ahead_of_cells = (parameters.ints("scenario") == 13) || (parameters.ints("scenario") == 22);
	bool have_you_removed_cells = false;
	bool have_you_deleted_puncta = false;

	try
	{
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time.
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout );
				if( PhysiCell_settings.enable_legacy_saves == true )
				{
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}

				if( PhysiCell_settings.enable_full_saves == true )
				{
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );
				}

				PhysiCell_globals.full_output_index++;
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}

			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

					PhysiCell_globals.SVG_output_index++;
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
		 	// microenvironment.simulate_diffusion_decay( diffusion_dt );

			// run PhysiCell
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			/*
			  Custom add-ons could potentially go here.
			*/
			// Once you update the mechanics in PhysiCell, make sure you apply the reflective boundary conditions:
			if(fabs((PhysiCell_globals.current_time-last_mechanics_time) - mechanics_dt)< mechanics_tolerance || !initialized)
			{
				double time_since_last_mechanics = PhysiCell_globals.current_time-last_mechanics_time;
				if(!initialized)
				{
					time_since_last_mechanics = mechanics_dt;
					initialized = true;
				}
				// std::cout<<"Applying boundary conditions... \n\n";
				#pragma omp parallel for
				for( int i=0; i < (*all_cells).size(); i++ )
				{
					if((*all_cells)[i]->type == 2 || (*all_cells)[i]->type == 3)
					{
						continue;
					}
					// Only leader cells can secrete FN puncta:
					if((*all_cells)[i]->type == 0 && experiment_allows_FN_secretion)
					{
						// Add in new fibronectin puncta, so long as the experiment allows it:
						add_new_fibronectin_puncta((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
					}
					if((*all_cells)[i]->type == 1 && followers_can_secrete_FN)
					{
						double average_interval = parameters.doubles("average_time_until_next_filopodia_drop"); // min
						if (parameters.ints("scenario") == 18 || parameters.ints("scenario") == 26)
						{
							average_interval = 10.0;
						}
						else if (parameters.ints("scenario") == 19)
						{
							average_interval = 5.0;
						}
						else if (parameters.ints("scenario") == 21)
						{
							average_interval = 30.0;
						}
						else if (parameters.ints("scenario") == 29 || parameters.ints("scenario") == 34 || parameters.ints("scenario") == 37 || parameters.ints("scenario") == 41 || parameters.ints("scenario") == 45 || parameters.ints("scenario") == 48 || parameters.ints("scenario") == 50 || parameters.ints("scenario") == 54 || parameters.ints("scenario") == 57 || parameters.ints("scenario") == 67 || parameters.ints("scenario") == 68 || parameters.ints("scenario") == 69)
						{
							average_interval = 90.0;
						}
						add_new_fibronectin_puncta_followers((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics, average_interval);
					}
				 	// Apply any and all boundary conditions:
				 	apply_reflective_boundary_condition((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
				}
				if (!have_you_removed_cells && parameters.ints("scenario")==12 && PhysiCell_globals.current_time>720)
				{
					have_you_removed_cells = true;
					#pragma omp parallel for
					for( int i=0; i < (*all_cells).size(); i++ )
					{
						if((*all_cells)[i]->type == 2 || (*all_cells)[i]->type == 3)
						{
							int cover_index = (*all_cells)[i]->custom_data.find_variable_index("are_you_covered_by_a_cell");
							(*all_cells)[i]->custom_data.variables[cover_index].value = 0;
							continue;
						}
						#pragma omp critical(delete_cells)
						{
							delete_cell((*all_cells)[i]);
						}
					}
				}
				// std::cout<<"Done!\n\n";
				// std::cout<<"Adding in new cells...\n\n";
				if(parameters.bools("cell_enters_regardless_of_space"))
				{
					// Add in new cells at the left-boundary at a constant rate:
					add_in_new_cells(time_since_last_mechanics);
				}
				else
				{
					add_in_new_cells_if_free_space();
				}
				// std::cout<<"Done!\n\n";
				// std::cout<<"Checking puncta and edge agents to see if they're covered...\n\n";
				// Make sure you cover all puncta:
				#pragma omp parallel for
				for(int j = 0; j < (*all_cells).size(); j++)
				{
					if((*all_cells)[j]->type == 0 || (*all_cells)[j]->type == 1)
					{
						continue;
					}
					check_puncta_for_cell_cover((*all_cells)[j], (*all_cells)[j]->phenotype, time_since_last_mechanics);
				}
				// If you can inject new FN and the simulation time is over 4 hours, then
				if(can_inject_new_FN_in_domain && fabs(PhysiCell_globals.current_time-240) < time_since_last_mechanics-1e-5)
				{
					double left_most_position = microenvironment.mesh.bounding_box[0];
					double right_most_position = 50; // microns
					if(inject_FN_ahead_of_cells)
					{
						// left_most_position = 50; // This is the original value for Scenario 13
						left_most_position = 70; // There is less chance that you place FN puncta where the cells are
						right_most_position = microenvironment.mesh.bounding_box[3];
					}
					inject_new_fibronectin(left_most_position, right_most_position);
				}

				if(!have_you_deleted_puncta && PhysiCell_globals.current_time > 240 && (parameters.ints("scenario") == 82 || parameters.ints("scenario") == 83 || parameters.ints("scenario") == 84 || parameters.ints("scenario") == 85))
				{
					have_you_deleted_puncta = true;
					std::vector<double> right_most_position(3, 0.0);
					for (int j = 0; j < (*all_cells).size(); j++)
					{
						if (((*all_cells)[j]->type == 0 || (*all_cells)[j]->type == 1) && (*all_cells)[j]->position[0] > right_most_position[0])
						{
							right_most_position = (*all_cells)[j]->position;
						}
					}
					delete_puncta(parameters.ints("scenario"), right_most_position);
				}

				// std::cout<<"Done!\n\n";
				last_mechanics_time = PhysiCell_globals.current_time;
			}
			// Advance the time forward
			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}

	// save a final simulation snapshot

	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

	// sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() );
	// SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );

	return 0;
}
