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

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <algorithm>

using namespace BioFVM;
using namespace PhysiCell;

// any additional cell types (beyond cell_defaults)
extern Cell_Definition follower_cell;
extern Cell_Definition fibronectin_puncta;
extern Cell_Definition edge_agent;

// custom cell phenotype functions could go here
void apply_reflective_boundary_condition(Cell* pCell, Phenotype& phenotype, double dt);
void record_total_velocity(Cell* pCell, Phenotype& phenotype, double dt);
std::vector<double> sample_from_Von_Mises_distribution(double mu, double kappa);
void add_new_fibronectin_puncta(Cell* pCell, Phenotype& phenotype, double dt);
void add_new_fibronectin_puncta_followers(Cell* pCell, Phenotype& phenotype, double dt, double average_interval);
void check_puncta_for_cell_cover(Cell* pCell, Phenotype& phenotype, double dt);

void calculate_surface_gradients(Cell* pCell, Cell* neighbor, Phenotype& phenotype, std::vector<double>& FN_gradient, std::vector<double>& direction_from_puncta_angle, std::vector<double>& contact_surface_gradient, std::vector<double>& spring_force_PhysiCell, std::vector<double>& spring_force_FN, int& number_of_fibronectin_neighbors, int& number_of_contact_neighbors, int& number_of_puncta_contributors);
void compute_directions(Cell* pCell, Phenotype& phenotype, double dt, std::vector<double>& motility_direction, std::vector<double>& contact_velocity, int& number_of_contact_neighbors);
void simpler_function_leader_follower_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt);
void inject_new_fibronectin( double left_most_position, double right_most_position );
void inject_new_fibronectin_outside_corridor( double right_most_position );
void delete_puncta( int scenario, std::vector<double> right_most_cell_position);

// setup functions to help us along

void create_cell_types( void );
void setup_tissue( void );

// set up the BioFVM microenvironment
void setup_microenvironment( void );
void add_in_new_cells(double dt);
void add_in_new_cells_if_free_space( void );

// custom pathology coloring function
std::vector<std::string> my_coloring_function( Cell* );
