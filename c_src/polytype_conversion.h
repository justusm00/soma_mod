/**
 * @Author: Your name
 * @Date:   2022-12-09 16:16:32
 * @Last Modified by:   Your name
 * @Last Modified time: 2023-01-19 16:32:41
 */
/* Copyright (C) 2016-2021 Ludwig Schneider

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOMA_POLYTYPE_CONVERSION_H
#define SOMA_POLYTYPE_CONVERSION_H

#include "soma_config.h"
#include "soma_util.h"

//! Top level struct for polymer type conversions.
//! controls the execution frequency. All other fields are only valid, if deltaMC != 0
typedef struct PolyConversion {
    unsigned int deltaMC;       //!< control execution frequency of the conversion
    uint8_t *array;             //!< Array that contains the reaction start index of the conversion list.

    unsigned int *input_type;   //!< Array that contains the input poly type for each reaction (educt)
    unsigned int *output_type;  //!< Array that contains the output poly type for each reaction (product)
    unsigned int *reaction_end; //!< Array indicating if this is the last reaction in the list. (boolean)
    soma_scalar_t *rate;        //!< control execution probability of the conversion
    unsigned int *dependency_ntype;     //!<Array that contains the number of  dependency indices
    unsigned int *dependency_type_offset;       //!<Array that contains the start/offset of dependency indices
    unsigned int *dependency_type;      //!<Array that contains the dependency types
    unsigned int len_reactions; //!< length of the reaction related arrays input_type, output_type and reaction_end
    unsigned int len_dependencies;      //!< length of the density dependency array dependency_type (=sum over dependency_ntype)
    unsigned int *num_conversions; //!< n_poly_type * n_poly_type Matrix containing the corresponding conversions counters


} PolyConversion;

//! Helper function to copy the pc data to the device
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyin_poly_conversion(struct Phase *p);

//! Helper function delete the pc data from the device and copy it to the CPU memory
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyout_poly_conversion(struct Phase *p);

//! Helper function to update the host with the pc data
//! \private
//! \param p Fully initialized Phase struct
//! \return Errorcode
int update_self_poly_conversion(const struct Phase *const p);

/*! Helper function to read the conversion array from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_poly_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the polyconversion array to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_poly_conversion_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to free the CPU memory resources of the pc struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \return Errorcode
*/
int free_poly_conversion(struct Phase *p);

/*! Convert polymer types according to the reaction description of the PolyConversion struct.
  This updates the center of mass of the polymers and chooses between full or partial (with rates) conversions.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int convert_polytypes(struct Phase *p);

/*! Fully convert polymer types according to the reaction description of the PolyConversion struct.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int fully_convert_polytypes(struct Phase *p);

/*! Partially Convert polymer types according to the reaction description of the PolyConversion struct.
 This converts the polymer only with a probability given by the rate which may depend (linearly) on the normalized density of some type (for reactions involving multiple types).
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int partially_convert_polytypes(struct Phase *p);


/*! Convert polymer types using a mixture of simulated annealing and random flipping.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int simulated_annealing(struct Phase *p);


/*! Convert polymer types without simulated annealing.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int convert_target(struct Phase *p);


/*! Compute cost function based on composition. Close to an area51 it makes more sense to use composition instead of density.
  \param p Phase struct describing the simulation
  \param delta_fields_unified Array containing changes in density fields.
  \return Errorcode
*/
soma_scalar_t get_composition_cost(struct Phase *p, int64_t * delta_fields_unified);


/*! Helper function to compare elements for quicksort.
  \param elem1 First element
  \param elem2 Second element
*/
int comp (const void * elem1, const void * elem2);

/*! Helper function to flip polymer type.
  \param initial_type Initial type.
*/
int flip(int initial_type);

/*! Get all information about flip candidates and where they have how many monomers.
  \param p Phase struct describing the simulation
  \param poly_isflippable Arrays are 0 if polymer is not flippable, 1 if it is
  \param poly_cell_indices For flippable polymers, contains indices in which cells the polymer has monomers
  \param poly_cell_num Number of monomers corresponding to the cells in poly_cell_indices
*/
void get_flip_candidates(struct Phase * p, int64_t * poly_isflippable, int64_t * poly_cell_indices, int64_t * poly_cell_num);


/*! Calculate cost difference based on composition if a polymer is flipped.
  \param p Phase struct describing the simulation
  \param poly Polymer index
  \param initial_type Type before flip
  \param final_type Type after flip
  \param poly_cell_indices For flippable polymers, contains indices in which cells the polymer has monomers
  \param poly_cell_num Number of monomers corresponding to the cells in poly_cell_indices
  \param delta_fields_unified Changes in density fields caused by flips
*/
soma_scalar_t get_composition_flip_cost(struct Phase * p, uint64_t poly, unsigned int initial_type, unsigned int final_type, int64_t * poly_cell_indices, int64_t * poly_cell_num,int64_t * delta_fields_unified);

/*! Update differences in density fields cause by flips
  \param p Phase struct describing the simulation
  \param poly Polymer index
  \param initial_type Type before flip
  \param final_type Type after flip
  \param poly_cell_indices For flippable polymers, contains indices in which cells the polymer has monomers
  \param poly_cell_num Number of monomers corresponding to the cells in poly_cell_indices
  \param delta_fields_unified Changes in density fields caused by flips
*/
void update_delta_fields(struct Phase * p, uint64_t poly, unsigned int initial_type, unsigned int final_type, int64_t * poly_cell_indices, int64_t * poly_cell_num,int64_t * delta_fields_unified);


/*! Run simulated annealing
  \param p Phase struct describing the simulation
  \param total_cost Cost value before function call
  \param num_poly_flippable Polymers available for flip
  \param Tmin Minimum temperature after which annealing will be stopped
  \param Tmax Starting temperature for SA run
  \param alpha Temperature decreasing factor
  \param sa_buffer_size Maximum number of polymers in single SA run
  \param poly_cell_indices For flippable polymers, contains indices in which cells the polymer has monomers
  \param poly_cell_num Number of monomers corresponding to the cells in poly_cell_indices
  \param poly_flippable_indices Sequentially stores indices of flippable polymers for faster iteration
  \param delta_fields_unified Changes in density fields caused by flips
  \param delta_fields_unified_best Changes in density fields caused by flips corresponding to current best value of cost function
  \param poly_types Polymer types of flippable polymers
  \param poly_types_best Polymer types of flippable polymers corresponding to current best value of cost function
*/


soma_scalar_t anneal_polytypes(struct Phase * p,soma_scalar_t total_cost, uint64_t num_poly_flippable, uint64_t * total_flip_attempts, uint64_t * total_flips_accepted,  int64_t * poly_cell_indices, int64_t * poly_cell_num, int64_t * poly_flippable_indices,  int64_t * delta_fields_unified, int64_t * delta_fields_unified_best,unsigned int * poly_types,unsigned int * poly_types_best,uint64_t num_target_cells);


/*! Simulated annealing at T=0
  \param p Phase struct describing the simulation
  \param total_cost Cost value before function call
  \param num_poly_flippable Polymers available for flip
  \param acc_rate_target Target acceptance rate after which iteration will be stopped
  \param poly_cell_indices For flippable polymers, contains indices in which cells the polymer has monomers
  \param poly_cell_num Number of monomers corresponding to the cells in poly_cell_indices
  \param poly_flippable_indices Sequentially stores indices of flippable polymers for faster iteration
  \param delta_fields_unified Changes in density fields caused by flips
  \param delta_fields_unified_best Changes in density fields caused by flips corresponding to current best value of cost function
  \param poly_types Polymer types of flippable polymers
  \param poly_types_best Polymer types of flippable polymers corresponding to current best value of cost function
*/

soma_scalar_t flip_polytypes(struct Phase * p,soma_scalar_t total_cost, uint64_t num_poly_flippable, uint64_t * total_flip_attempts,uint64_t * total_flips_accepted, int64_t * poly_cell_indices, int64_t * poly_cell_num, int64_t * poly_flippable_indices, int64_t * delta_fields_unified, int64_t * delta_fields_unified_best,unsigned int * poly_types ,unsigned int * poly_types_best);

#endif                          //SOMA_POLYTYPE_CONVERSION_H
