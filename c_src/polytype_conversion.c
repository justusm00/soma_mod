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

#include "polytype_conversion.h"

#include <stdlib.h>
#include <hdf5.h>
#include "phase.h"
#include "io.h"
#include "mesh.h"

//! \file polytype_conversion.c
//! \brief Implementation of polytype_conversion.h

int read_poly_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{

    p->pc.deltaMC = 0;
    p->pc.array = NULL;
    p->pc.len_reactions = 0;
    p->pc.len_dependencies = 0;
    p->pc.num_conversions = NULL;
    p->pc.rate = NULL;
    p->pc.input_type = NULL;
    p->pc.output_type = NULL;
    p->pc.reaction_end = NULL;
    p->pc.dependency_ntype = NULL;
    p->pc.dependency_type = NULL;
    p->pc.dependency_type_offset = NULL;
    //Quick exit if no poly conversion is present in the file
    if (!(H5Lexists(file_id, "/polyconversion", H5P_DEFAULT) > 0))
        return 0;

    herr_t status;
    unsigned int tmp_deltaMC;
    status = read_hdf5(file_id, "polyconversion/deltaMC", H5T_STD_U32LE, plist_id, &tmp_deltaMC);
    HDF5_ERROR_CHECK(status);

    //Quick exit if no conversion updates are scheduled
    if (tmp_deltaMC == 0)
        return 0;

#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    fprintf(stderr,
            "ERROR: %s: %d, Monotype Conversions are activated so polytype conversions do not work in the current implementation. Switch off the this feature and rerun.",
            __FILE__, __LINE__);
    return -1;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    //Read the conversion list before reading the spatial array.
    const hid_t dset_input = H5Dopen(file_id, "/polyconversion/input_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_input);
    const hid_t dspace_input = H5Dget_space(dset_input);
    HDF5_ERROR_CHECK(dspace_input);
    const unsigned int ndims_input = H5Sget_simple_extent_ndims(dspace_input);
    if (ndims_input != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/input_type");
            return -1;
        }
    hsize_t dim_input;
    status = H5Sget_simple_extent_dims(dspace_input, &dim_input, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_output = H5Dopen(file_id, "/polyconversion/output_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_output);
    const hid_t dspace_output = H5Dget_space(dset_output);
    HDF5_ERROR_CHECK(dspace_output);
    const unsigned int ndims_output = H5Sget_simple_extent_ndims(dspace_output);
    if (ndims_output != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/output_type");
            return -1;
        }
    hsize_t dim_output;
    status = H5Sget_simple_extent_dims(dspace_output, &dim_output, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_end = H5Dopen(file_id, "/polyconversion/end", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_end);
    const hid_t dspace_end = H5Dget_space(dset_end);
    HDF5_ERROR_CHECK(dspace_end);
    const unsigned int ndims_end = H5Sget_simple_extent_ndims(dspace_end);
    if (ndims_end != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/end_type");
            return -1;
        }
    hsize_t dim_end;
    status = H5Sget_simple_extent_dims(dspace_end, &dim_end, NULL);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_output)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input and output type for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_output);
            return -3;
        }
    if (dim_input != dim_end)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and end for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_end);
            return -3;
        }

    p->pc.input_type = (unsigned int *)malloc(dim_input * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.input_type, dim_input * sizeof(unsigned int));
    p->pc.output_type = (unsigned int *)malloc(dim_output * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.output_type, dim_output * sizeof(unsigned int));
    p->pc.reaction_end = (unsigned int *)malloc(dim_end * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.reaction_end, dim_end * sizeof(unsigned int));
    p->pc.num_conversions = (unsigned int *)malloc(p->n_poly_type * p->n_poly_type * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.num_conversions,p->n_poly_type * p->n_poly_type * sizeof(unsigned int));

    for(int i = 0; i < p->n_poly_type * p->n_poly_type;i++)
        {
            p->pc.num_conversions[i]=0;
        }


    p->pc.len_reactions = dim_input;

    //Actuablly read the polyconversion list data
    status = H5Dread(dset_input, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.input_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_output, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.output_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_end, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.reaction_end);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_end);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_end);
    HDF5_ERROR_CHECK(status);

    //Read the array information
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/polyconversion/array", H5P_DEFAULT);
    p->pc.array =
        (uint8_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(uint8_t));
    if (p->pc.array == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->pc.array + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

#if ( ENABLE_MPI == 1 )
    const int left_neigh_rank =
        (((my_domain - 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;
    const int right_neigh_rank =
        (((my_domain + 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;

    MPI_Request req[4];
    MPI_Status stat[4];

    uint8_t *ptr = p->pc.array + ghost_buffer_size;
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->pc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->pc.array;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->pc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 3);

    MPI_Waitall(4, req, stat);
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    //If rate is defined in the hdf5, partial conversions are activated and "rate", "n_density_dependencies", "density_dependencies" are read.
    if (!(H5Lexists(file_id, "/polyconversion/rate", H5P_DEFAULT) > 0))
        {
            //Enable the updat only if everything worked fine so far
            p->pc.deltaMC = tmp_deltaMC;
            return 0;
        }

    const hid_t dset_rate = H5Dopen(file_id, "/polyconversion/rate", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_rate);
    const hid_t dspace_rate = H5Dget_space(dset_rate);
    HDF5_ERROR_CHECK(dspace_rate);
    const unsigned int ndims_rate = H5Sget_simple_extent_ndims(dspace_rate);
    if (ndims_rate != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/rate");
            return -1;
        }
    hsize_t dim_rate;
    status = H5Sget_simple_extent_dims(dspace_rate, &dim_rate, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_ndependency = H5Dopen(file_id, "/polyconversion/n_density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_ndependency);
    const hid_t dspace_ndependency = H5Dget_space(dset_ndependency);
    HDF5_ERROR_CHECK(dspace_ndependency);
    const unsigned int ndims_ndependency = H5Sget_simple_extent_ndims(dspace_ndependency);
    if (ndims_ndependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/n_density_dependencies");
            return -1;
        }
    hsize_t dim_ndependency;
    status = H5Sget_simple_extent_dims(dspace_ndependency, &dim_ndependency, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_dependency = H5Dopen(file_id, "/polyconversion/density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_dependency);
    const hid_t dspace_dependency = H5Dget_space(dset_dependency);
    HDF5_ERROR_CHECK(dspace_dependency);
    const unsigned int ndims_dependency = H5Sget_simple_extent_ndims(dspace_dependency);
    if (ndims_dependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/density_dependencies");
            return -1;
        }
    hsize_t dim_dependency;
    status = H5Sget_simple_extent_dims(dspace_dependency, &dim_dependency, NULL);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_rate)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }
    if (dim_input != dim_ndependency)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }

    p->pc.rate = (soma_scalar_t *) malloc(dim_rate * sizeof(soma_scalar_t));
    MALLOC_ERROR_CHECK(p->pc.rate, dim_rate * sizeof(soma_scalar_t));
    p->pc.dependency_ntype = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_ntype, dim_ndependency * sizeof(unsigned int));
    p->pc.dependency_type_offset = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_type_offset, dim_ndependency * sizeof(unsigned int));
    p->pc.dependency_type = (unsigned int *)malloc(dim_dependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_type, dim_dependency * sizeof(unsigned int));


    p->pc.len_dependencies = dim_dependency;

    //Read rate and dependencies:
    status = H5Dread(dset_rate, H5T_SOMA_NATIVE_SCALAR, H5S_ALL, H5S_ALL, plist_id, p->pc.rate);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_ndependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.dependency_ntype);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_ndependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_ndependency);
    HDF5_ERROR_CHECK(status);
    p->pc.dependency_type_offset[0] = 0;
    for (unsigned int i = 1; i < dim_ndependency; i++)
        {
            p->pc.dependency_type_offset[i] = p->pc.dependency_type_offset[i - 1] + p->pc.dependency_ntype[i - 1];
        }
    if (p->pc.len_dependencies > 0)
        {
            status = H5Dread(dset_dependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.dependency_type);
            HDF5_ERROR_CHECK(status);
        }
    status = H5Sclose(dspace_dependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_dependency);
    HDF5_ERROR_CHECK(status);

    //Enable the updat only if everything worked fine so far
    p->pc.deltaMC = tmp_deltaMC;

    return 0;
}

int write_poly_conversion_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Quick exit for no poly conversions
    if (p->pc.deltaMC == 0)
        return 0;
    hid_t group = H5Gcreate2(file_id, "/polyconversion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status;
    const hsize_t one = 1;
    status =
        write_hdf5(1, &one, file_id, "/polyconversion/deltaMC", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   &(p->pc.deltaMC));
    HDF5_ERROR_CHECK(status);

    const hsize_t list_len = p->pc.len_reactions;
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/input_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.input_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/output_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.output_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/end", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.reaction_end);
    HDF5_ERROR_CHECK(status);

    if (p->pc.rate != NULL)
        {
            status =
                write_hdf5(1, &list_len, file_id, "/polyconversion/rate", H5T_SOMA_NATIVE_SCALAR,
                           H5T_SOMA_NATIVE_SCALAR, plist_id, p->pc.rate);
            HDF5_ERROR_CHECK(status);
            status =
                write_hdf5(1, &list_len, file_id, "/polyconversion/n_density_dependencies", H5T_STD_U32LE,
                           H5T_NATIVE_UINT, plist_id, p->pc.dependency_ntype);
            HDF5_ERROR_CHECK(status);
            const hsize_t dep_list_len = p->pc.len_dependencies;
            if (dep_list_len > 0)
                {
                    status =
                        write_hdf5(1, &dep_list_len, file_id, "/polyconversion/density_dependencies", H5T_STD_U32LE,
                                   H5T_NATIVE_UINT, plist_id, p->pc.dependency_type);
                    HDF5_ERROR_CHECK(status);
                }
            else
                {               //no dependencies --> empty dataset, so only create, don't write.
                    herr_t status;
                    const hid_t dataspace = H5Screate_simple(1, &dep_list_len, NULL);
                    HDF5_ERROR_CHECK(dataspace);
                    const hid_t dataset =
                        H5Dcreate2(file_id, "/polyconversion/density_dependencies", H5T_STD_U32LE, dataspace,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    HDF5_ERROR_CHECK(dataset);
                    status = H5Sclose(dataspace);
                    HDF5_ERROR_CHECK(status);
                    status = H5Dclose(dataset);
                    HDF5_ERROR_CHECK(status);
                }
        }

    //Write the array to disk
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset =
        H5Dcreate2(file_id, "/polyconversion/array", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI

    if ((status =
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->pc.array + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Gclose(group)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    return 0;
}

int copyin_poly_conversion(struct Phase *p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
            //The pc struct itself is part of the phase struct and is already present of the device
#pragma acc enter data copyin(p->pc.array[0:p->n_cells_local])
#pragma acc enter data copyin(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.reaction_end[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.num_conversions[0:p->n_poly_type*p->n_poly_type])
            if (p->pc.rate != NULL)
                {
#pragma acc enter data copyin(p->pc.rate[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int copyout_poly_conversion(struct Phase *p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc exit data copyout(p->pc.array[0:p->n_cells_local])
#pragma acc exit data copyout(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.reaction_end[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.num_conversions[0:p->n_poly_type*p->n_poly_type])
            if (p->pc.rate != NULL)
                {
#pragma acc exit data copyout(p->pc.rate[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int update_self_poly_conversion(const struct Phase *const p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc update self(p->pc.array[0:p->n_cells_local])
#pragma acc update self(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc update self(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc update self(p->pc.reaction_end[0:p->pc.len_reactions])
#pragma acc update self(p->pc.num_conversions[0:p->n_poly_type*p->n_poly_type])
            if (p->pc.rate != NULL)
                {
#pragma acc update self(p->pc.rate[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int free_poly_conversion(struct Phase *p)
{
    free(p->pc.array);
    free(p->pc.input_type);
    free(p->pc.output_type);
    free(p->pc.reaction_end);
    free(p->pc.rate);
    free(p->pc.dependency_ntype);
    free(p->pc.dependency_type_offset);
    free(p->pc.dependency_type);
    free(p->pc.num_conversions);

    return 0;
}


int convert_polytypes(struct Phase *p)
{
    //Quick exit for
    static unsigned int last_time = 0;
    if (last_time >= p->time)
        return 0;
    last_time = p->time;
    update_polymer_rcm(p);

    if (p->umbrella_field != NULL)
        {
            //do polymer conversion only with root rank
            if (p->info_MPI.sim_rank == 0) 
                {
                    return simulated_annealing(p);
                }
            
        }
    else
        {
            printf("ERROR: No umbrella field available \n");
            return -1;
        }

/*     if (p->pc.rate == NULL)
        {
            return fully_convert_polytypes(p);
        }
    else
        {
            return partially_convert_polytypes(p);
        } */
}

int fully_convert_polytypes(struct Phase *p)
{
//Iterate all polymers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            const Monomer rcm = p->polymers[poly].rcm;
            const uint64_t cell = coord_to_index(p, rcm.x, rcm.y, rcm.z);


            if (p->pc.array[cell] != 0)
                {
                    //Minus 1 because the index in array are shifted by 1
                    int i = p->pc.array[cell] - 1;
                    if (p->polymers[poly].type == p->pc.input_type[i])
                        {
#pragma acc atomic update
                            p->pc.num_conversions[p->polymers[poly].type * p->n_poly_type + p->pc.output_type[i]]++;                            
                            p->polymers[poly].type = p->pc.output_type[i];
                        }
                    for (i++; !p->pc.reaction_end[i - 1]; i++)
                        {
                            if (p->polymers[poly].type == p->pc.input_type[i])
                            {
#pragma acc atomic update
                                p->pc.num_conversions[p->polymers[poly].type * p->n_poly_type + p->pc.output_type[i]]++;                                     
                                p->polymers[poly].type = p->pc.output_type[i];

                            }
                        }
                }
        }
    return 0;
}

int partially_convert_polytypes(struct Phase *p)
{
    //Iterate all polymers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            Polymer *mypoly = p->polymers + poly;
            const Monomer rcm = mypoly->rcm;
            const uint64_t cell = coord_to_index(p, rcm.x, rcm.y, rcm.z);

            if (p->pc.array[cell] != 0)
                {
                    soma_scalar_t probability = 0.;
                    int i = p->pc.array[cell] - 2;
                    do
                        {
                            i++;
                            if (mypoly->type == p->pc.input_type[i])
                                {
                                    soma_scalar_t norm = 1 - probability;
                                    probability = p->pc.rate[i];
                                    for (unsigned int j = 0; j < p->pc.dependency_ntype[i]; j++)
                                        {
                                            unsigned int type_offset = p->pc.dependency_type_offset[i];
                                            unsigned int dependency_type = p->pc.dependency_type[type_offset + j];
                                            probability *=
                                                p->fields_unified[dependency_type * p->n_cells_local +
                                                                  cell] * p->field_scaling_type[dependency_type];
                                        }
                                    probability /= norm;
                                    soma_scalar_t random_number = soma_rng_soma_scalar(&(mypoly->poly_state), p);
                                    if (random_number < probability)
                                        {
#pragma acc atomic update
                                            p->pc.num_conversions[p->polymers[poly].type * p->n_poly_type + p->pc.output_type[i]]++;                                             
                                            p->polymers[poly].type = p->pc.output_type[i];      //got modifiable lvalue compile error when using mypoly->type = ... and was not able to fix this otherwise.
                                            break;      //to continue with next polymer if conversion has taken place.
                                        }
                                    else
                                        {
                                            probability += (1 - norm);
                                        }
                                }
                    } while (!p->pc.reaction_end[i]);
                }
        }

    return 0;
}


int simulated_annealing(struct Phase *p)
{
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }
 

    // adjustable parameters
    uint64_t flip_buffer_size = p->n_polymers; //maximum number of flippable polymers, need to find optimal value


    // initialize other parameters
    uint64_t num_target_cells =0; // Total number of cells for which target density is available * polytypes
    uint64_t num_poly_flippable = 0; // number of flippable polymers
    uint64_t total_flip_attempts = 0; //total number of flip attempts
    uint64_t total_flips_accepted = 0; //total number of flip attempts
    soma_scalar_t total_cost = 0.0; 

 

    // arrays and constants
    int64_t * poly_isflippable = (int64_t *)malloc( p->n_polymers* sizeof(int64_t)); //Boolean array that stores whether or not polymer has monomers in target density area.
    int64_t * poly_cell_indices = (int64_t *)malloc(p->n_polymers * p->reference_Nbeads * sizeof(int64_t)); //Array of length p->n_polymers * p->reference_Nbeads that stores cell indices in which a polymer has monomers.
    int64_t * poly_cell_num = (int64_t *)malloc(p->n_polymers  * p->n_poly_type* p->n_types * p->reference_Nbeads * sizeof(int64_t)); //Array of length p->n_polymers * p->reference_Nbeads * p->n_types * p->n_poly_type that stores number of monomers in cells corresponding to poly_cell_indices for each possible polymer type. 
    int64_t * poly_flippable_indices = (int64_t *)malloc( flip_buffer_size * sizeof(int64_t)); //array that contains indices of flippable polymers
    int64_t * delta_fields_unified = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); //array that stores changes in density
    int64_t * delta_fields_unified_best = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); 
    unsigned int * poly_types=(unsigned int *)malloc(flip_buffer_size * sizeof(unsigned int)); //array that stores polymer types
    unsigned int * poly_types_best=(unsigned int *)malloc(flip_buffer_size * sizeof(unsigned int)); //array that stores best polymer types


    //initialize poly_isflippable
    for (uint64_t poly = 0; poly < p->n_polymers; poly++) poly_isflippable[poly]=0;

    //initialize poly_cell_indices with -1's
    for (uint64_t i = 0; i < p->n_polymers * p->reference_Nbeads; i++) poly_cell_indices[i] = -1;

    //initialize poly_cell_num
    for (uint64_t i = 0; i < p->n_polymers * p->reference_Nbeads * p->n_poly_type* p->n_types; i++) poly_cell_num[i] = 0;

    //get flippable polymers
    get_flip_candidates(p, poly_isflippable, poly_cell_indices, poly_cell_num);





    //get number of target cells * polymer types (not actually needed for optimization, just for normalization of the mean squared error)
    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            for(uint64_t type = 0; type < p->n_types; type++)
                {
                    if(p->umbrella_field[type*p->n_cells_local + cell] > 0) num_target_cells++;
                }
        }


    //save flippable polymer indices sequentially to new array for quicker access
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            if(poly_isflippable[poly]==1)
                {
                    poly_flippable_indices[num_poly_flippable]=poly;
                    num_poly_flippable++;

                }
        }


/*     //print cell information of some polymer
    uint64_t some_poly = poly_flippable_indices[0];
    unsigned int N = p->reference_Nbeads;
    unsigned int polytype = p->polymers[some_poly].type;
    unsigned int monotype = polytype;
    for(unsigned int mono= 0; mono < p->reference_Nbeads; mono++)
        {
            unsigned int bla = 0;
            if(poly_cell_indices[some_poly * N + mono] < 0) break;
            printf("%llu\n",poly_cell_indices[some_poly * N + mono]);
            bla = poly_cell_num[some_poly * p->n_poly_type * p->n_types * N + polytype * p->n_types * N +monotype * N + mono];
            printf("%u\n",bla);
        } */
    
    //check if there are more flippable polymers than the buffer allows
    if(num_poly_flippable>flip_buffer_size)
        {
            printf("ERROR: FLIP BUFFER SIZE IS TOO SMALL\n");
            return -1;
        }


    // intialize delta_fields_unified

    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            for(uint64_t type = 0; type < p->n_types; type++)
                {
                    delta_fields_unified[type*p->n_cells_local + cell] = 0;
                    delta_fields_unified_best[type*p->n_cells_local + cell] = 0;
                }
        }


    // inititalize poly_types
    
    for (uint64_t poly = 0; poly < num_poly_flippable; poly++)
        {
            poly_types[poly]=p->polymers[poly_flippable_indices[poly]].type;
            poly_types_best[poly]=poly_types[poly];
        }
    
    //initialize cost
    total_cost=get_composition_cost(p, delta_fields_unified);

    printf("Start configuration optimization at t=%d on testing branch\n",p->time);
    printf("MSE before annealing: %f \n",total_cost/(soma_scalar_t)num_target_cells);   

    //get new cost value from simulated annealing
    total_cost = anneal_polytypes(p,total_cost, num_poly_flippable, &total_flip_attempts, &total_flips_accepted,  poly_cell_indices, poly_cell_num, poly_flippable_indices,  delta_fields_unified, delta_fields_unified_best, poly_types, poly_types_best,num_target_cells);

    printf("MSE after annealing: %f \n",total_cost/(soma_scalar_t)num_target_cells);  

    //do some more flips at T=0
    total_cost = flip_polytypes(p,total_cost, num_poly_flippable, &total_flip_attempts, &total_flips_accepted, poly_cell_indices, poly_cell_num, poly_flippable_indices,  delta_fields_unified,delta_fields_unified_best, poly_types, poly_types_best);

    printf("MSE after flips at T=0: %f \n",total_cost/(soma_scalar_t)num_target_cells);

    printf("Polymers flipped: %llu\n",total_flip_attempts);
    printf("Accepted flips: %llu\n",total_flips_accepted);


    //update polymer types
    for(uint64_t poly = 0; poly < num_poly_flippable; poly++) p->polymers[poly_flippable_indices[poly]].type=poly_types_best[poly];
    
    free(poly_isflippable);
    free(poly_cell_indices);
    free(poly_cell_num);
    free(poly_types);
    free(poly_types_best);
    free(delta_fields_unified);
    free(delta_fields_unified_best);
    free(poly_flippable_indices);
    return 0;
}


void get_flip_candidates(struct Phase * p, int64_t * poly_isflippable, int64_t * poly_cell_indices, int64_t * poly_cell_num)
{
    const unsigned int N = p->reference_Nbeads; //polymer length (only if all polymers have the same length)
    int64_t * mono_cells=(int64_t *)malloc( N *  sizeof(int64_t)); //Array of length p->reference_Nbeads that contains monomer cell indices of a polymer. Values are -1 if no target density available in that cell.
    //loop over polymers to identify the ones that may be flipped
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            unsigned int initial_poly_type = p->polymers[poly].type;
            unsigned int target_count = 0; //counts number of available target densities for polymer accounting for all possible types after potential flips
            //loop over monomers to get cell information while disregarding the monomer type for now
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //Read Monomer cell
                    mono_cells[mono]=mono_cell; //write to mono_cells array
                    //check if umbrella field is available for some (or more) types
                    for (unsigned int monotype = 0; monotype < p->n_types; monotype ++)
                        {
                            if(p->umbrella_field[monotype*p->n_cells_local + mono_cell] > 0)
                                {
                                    target_count++;
                                    break;
                                }

                        }
                }

            
            //check for target density availability
            if(target_count != 0)
                {
                    poly_isflippable[poly]=1; //1 means that polymer has monomers in target density area
                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp); //sort cells to get unique ones
                    unsigned int k=0;
                    //loop over monomers to get unique cells
                    for(unsigned int mono = 0; mono < N - 1 ; mono++)
                        {
                            //update array only for new cells
                            if(mono_cells[mono]!= mono_cells[mono+1])
                                {
                                    poly_cell_indices[poly * N + k]=mono_cells[mono]; //unique monomer cell indices
                                    k++;
                                }
                        }

                    poly_cell_indices[poly * N + k]=mono_cells[N-1];


                    //loop over polymer types to fill poly_cell_num array
                    for(unsigned int polytype = 0; polytype < p->n_poly_type ; polytype++)
                        {
                            //temporarily change polytype
                            p->polymers[poly].type = polytype;
                            //loop over monomers
                            for(unsigned int mono = 0; mono < N ; mono++)
                                {
                                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //Read Monomer cell
                                    unsigned int monotype = get_particle_type(p, poly, mono); //Read Monomer type
                                    unsigned int mono_cell_offset=0; //offset of monomer cell in poly_cell_indices
                                    //Get cell offset for poly_cell_indices
                                    for(unsigned int cell_offset = 0; cell_offset < N; cell_offset++)
                                        {
                                            if(poly_cell_indices[poly * N + cell_offset] == mono_cell)
                                                {
                                                    mono_cell_offset = cell_offset;
                                                    break;
                                                }
                                        }
                                    poly_cell_num[poly * p->n_poly_type * p->n_types * N + polytype * p->n_types * N +monotype * N + mono_cell_offset]++; //increment according entry in poly_cells_num (tedious but I don't say another way)
                                }

                        }

                    //reset polymer type to original one
                    p->polymers[poly].type = initial_poly_type;

                }
            
        }
    free(mono_cells);
    return;
}



void update_delta_fields(struct Phase * p, uint64_t poly, unsigned int initial_type, unsigned int final_type, int64_t * poly_cell_indices, int64_t * poly_cell_num,int64_t * delta_fields_unified)
{
    
    const Polymer *polymer = p->polymers + poly;
    const unsigned int N = p->reference_Nbeads;

    // loop over cell indices
    for(unsigned int mono_cell_offset = 0; mono_cell_offset < N; mono_cell_offset++)
        {
            if(poly_cell_indices[poly  * N + mono_cell_offset] < 0) break; //check if end has been reached
            int64_t cell = poly_cell_indices[poly * N + mono_cell_offset];
            uint64_t beads_in_cell = 0; //total number of beads in the cell
            //loop over monotypes to update delta_fields
            for(unsigned int monotype = 0; monotype < p->n_types; monotype++)
                {
                    if(p->umbrella_field[monotype*p->n_cells_local + cell] > 0)
                        {
                            //number of monomers of current type in current cell before and after flip
                            int64_t num_mono_initial = poly_cell_num[poly * p->n_poly_type * p->n_types * N + initial_type * p->n_types * N + monotype * N + mono_cell_offset ];
                            int64_t num_mono_final = poly_cell_num[poly * p->n_poly_type * p->n_types * N + final_type * p->n_types * N + monotype * N + mono_cell_offset ];
                            //difference of the two
                            int64_t delta_num_mono = num_mono_final - num_mono_initial;
                            delta_fields_unified[monotype*p->n_cells_local + cell] += delta_num_mono;
                        }
                }
        }
        
    return;
}



soma_scalar_t anneal_polytypes(struct Phase * p,soma_scalar_t total_cost, uint64_t num_poly_flippable, uint64_t * total_flip_attempts, uint64_t * total_flips_accepted, int64_t * poly_cell_indices, int64_t * poly_cell_num, int64_t * poly_flippable_indices,  int64_t * delta_fields_unified, int64_t * delta_fields_unified_best,unsigned int * poly_types,unsigned int * poly_types_best,uint64_t num_target_cells)
{
    soma_scalar_t total_cost_old = total_cost;
    soma_scalar_t total_cost_best = total_cost;
    soma_scalar_t T = p->Tmax;
    soma_scalar_t acc_rate = 0;
    uint64_t flip_counter = 0; //counts flips per temperature for acceptance rate calculation
    uint64_t flips_acc = 0; //counts accepted flips per temperature for acceptance rate calculation
    uint64_t num_poly_flipped=0; //number of accepted flips before improvement in cost function
    int64_t * poly_flipped_indices = (int64_t *)malloc( num_poly_flippable * sizeof(int64_t)); //array that contains indices of flipped polymers. the same polymer can be flipped multiple times


    uint64_t temperature_count=0; //used for temperature calculation
    uint64_t n_cycles = 50; //number of temperatures considered

    const unsigned int N = p->reference_Nbeads; //maximum polymer length
    
    while(T > p->Tmin)
        {
            temperature_count++;
            total_cost=total_cost_best;
            total_cost_old=total_cost;
            //do num_poly_flippable random flips at current temperature
            for(uint64_t i = 0; i< num_poly_flippable; i++)
                {
                    flip_counter++;
                    //choose random polymer to flip
                    uint64_t random_index = (uint64_t)(soma_rng_soma_scalar(&((p->polymers)->poly_state), p) * (soma_scalar_t)(num_poly_flippable - 1));
                    uint64_t poly = poly_flippable_indices[random_index];
                    Polymer *mypoly = p->polymers + poly;
                    unsigned int initial_type = poly_types[random_index];
                    unsigned int final_type = flip(p,poly);
                    *total_flip_attempts=*total_flip_attempts+1;
                    
                    //calculate cost 

                    total_cost += get_composition_flip_cost(p, poly, initial_type, final_type, poly_cell_indices, poly_cell_num, delta_fields_unified);

                    soma_scalar_t random_number = soma_rng_soma_scalar(&(mypoly->poly_state), p);
                    //check acceptance criterion
                    if((total_cost < total_cost_old) || (random_number < exp(-(total_cost - total_cost_old)/T)) )
                        {
                            //accept flip
                            flips_acc++;
                            poly_flipped_indices[num_poly_flipped]=random_index;
                            num_poly_flipped++;
                            *total_flips_accepted=*total_flips_accepted+1;
                            poly_types[random_index]=final_type;
                            total_cost_old=total_cost;
                            //update delta fields unified
                            update_delta_fields(p, poly, initial_type, final_type, poly_cell_indices, poly_cell_num, delta_fields_unified);

                        }
                    else
                        {
                            //reject
                            total_cost = total_cost_old;
                        }



                    //update best solution so far
                    if (total_cost < total_cost_best)
                        {
                            total_cost_best=total_cost;
                            //loop over flipped polymers
                            //note that the same polymer could have been flipped multiple times
                            for(uint64_t j = 0; j < num_poly_flipped; j++)
                                {
                                    poly_types_best[poly_flipped_indices[j]]=poly_types[poly_flipped_indices[j]];
                                    poly=poly_flippable_indices[poly_flipped_indices[j]];
                                    for(unsigned int i = 0; i < N; i++)
                                        {
                                            if(poly_cell_indices[poly * N + i] < 0) break;
                                            unsigned int cell = poly_cell_indices[poly * N + i];
                                            for(unsigned int type = 0; type < p->n_types; type++) delta_fields_unified_best[type*p->n_cells_local + cell] = delta_fields_unified[type*p->n_cells_local + cell];
                                        }
                                }
                            //reset flip counter   
                            num_poly_flipped=0;
                        }

                }
            //Revert to best solution if bad flips remain
            for(uint64_t j = 0; j < num_poly_flipped; j++)
                {
                    poly_types[poly_flipped_indices[j]]=poly_types_best[poly_flipped_indices[j]];
                    uint64_t poly=poly_flippable_indices[poly_flipped_indices[j]];
                    for(unsigned int i = 0; i < N; i++)
                        {
                            if(poly_cell_indices[poly * N + i] < 0) break;
                            unsigned int cell = poly_cell_indices[poly * N + i];
                            for(unsigned int type = 0; type < p->n_types; type++) delta_fields_unified[type*p->n_cells_local + cell] = delta_fields_unified_best[type*p->n_cells_local + cell];
                        }
                }
                
            num_poly_flipped=0;
            printf("T: %f\n",T);
            printf("acceptance_rate: %f\n",(soma_scalar_t)(flips_acc)/(soma_scalar_t)(flip_counter));
            printf("MSE: %f\n",(soma_scalar_t)total_cost_best/num_target_cells);
            flip_counter=0;
            flips_acc=0;

            //update temperature
            T *= p->alpha; //exponential cooling
            //T=p->Tmax/(1+p->alpha*log(1.0+(soma_scalar_t)temperature_count)); //logarithmic cooling
            //T=p->Tmax/(1.0+p->alpha*(soma_scalar_t)temperature_count); //linear multiplicative cooling
            //T=p->Tmax/(1.0+p->alpha*(soma_scalar_t)temperature_count*(soma_scalar_t)temperature_count); //quadratic multiplicative cooling
            //T=p->Tmin+(p->Tmax-p->Tmin)*pow((soma_scalar_t)(n_cycles-temperature_count)/(soma_scalar_t)n_cycles,2); //quadratic additive cooling
            //T-=0.0001; //linear cooling
        }

    free(poly_flipped_indices);
    return total_cost_best;
}


soma_scalar_t flip_polytypes(struct Phase * p,soma_scalar_t total_cost, uint64_t num_poly_flippable, uint64_t * total_flip_attempts,uint64_t * total_flips_accepted, int64_t * poly_cell_indices, int64_t * poly_cell_num, int64_t * poly_flippable_indices, int64_t * delta_fields_unified, int64_t * delta_fields_unified_best,unsigned int * poly_types ,unsigned int * poly_types_best)
{
    const unsigned int N = p->reference_Nbeads; //monomers per polymer
    soma_scalar_t total_cost_old = total_cost;
    soma_scalar_t acc_rate = 1.0;
    uint64_t flip_counter =0; //counts number of flips 
    uint64_t flip_counter_acc =0; //counts number of accepted flips 

    //do num_poly_flippable random flips at T=0
    //for(uint64_t i = 0; i < num_poly_flippable; i++)
    while(acc_rate > 0.1)
        {
            total_cost=total_cost_old;
            //choose random polymer to flip
            flip_counter++;
            uint64_t random_index = rand() % (num_poly_flippable - 1);
            uint64_t poly = poly_flippable_indices[random_index];
            Polymer *mypoly = p->polymers + poly;
            unsigned int initial_type = poly_types_best[random_index];
            unsigned int final_type = flip(p,poly);
            *total_flip_attempts=*total_flip_attempts+1;
            //calculate cost 
            total_cost += get_composition_flip_cost(p, poly, initial_type, final_type, poly_cell_indices, poly_cell_num, delta_fields_unified);
            //check acceptance criterion
            if((total_cost < total_cost_old))
                {
                    flip_counter_acc++;
                    *total_flips_accepted=*total_flips_accepted+1;
                    //accept flip
                    poly_types[random_index]=final_type;
                    poly_types_best[random_index]=final_type;
                    total_cost_old=total_cost;
                    //update delta fields unified
                    update_delta_fields(p, poly, initial_type, final_type, poly_cell_indices, poly_cell_num, delta_fields_unified);
                    update_delta_fields(p, poly, initial_type, final_type, poly_cell_indices, poly_cell_num, delta_fields_unified_best);

                }
            else
                {
                    //reject
                    total_cost = total_cost_old;
                }
            if(flip_counter % 10 == 0) acc_rate = (soma_scalar_t)(flip_counter_acc)/(soma_scalar_t)(flip_counter);


        }

 //   printf("%d\n",flip_counter);
    return total_cost;
}


soma_scalar_t get_composition_cost(struct Phase *p, int64_t * delta_fields_unified)
{   
    soma_scalar_t total_cost=0.0;
    //loop over cells
    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            //update cost
            for(unsigned int type = 0; type < p->n_types; type++)
                {
                    if(p->umbrella_field[type*p->n_cells_local + cell] > 0)
                        {
                            //get number of beads in cell
                            uint16_t beads_in_cell = 0; 
                            for(uint64_t type = 0; type < p->n_types; type++) beads_in_cell += p->fields_unified[type*p->n_cells_local + cell];
                            total_cost+=powl((soma_scalar_t)p->umbrella_field[type*p->n_cells_local + cell]-(soma_scalar_t)( p->fields_unified[type*p->n_cells_local + cell] + delta_fields_unified[type*p->n_cells_local + cell]) /beads_in_cell,2.0);
                        }
                }
        }
    return total_cost;
}


soma_scalar_t get_composition_flip_cost(struct Phase * p, uint64_t poly, unsigned int initial_type, unsigned int final_type, int64_t * poly_cell_indices, int64_t * poly_cell_num,int64_t * delta_fields_unified)
{
    
    soma_scalar_t delta_cost = 0.0; //change in cost function
    const unsigned int N = p->reference_Nbeads; //maximum polymer length

    // loop over cell indices
    for(unsigned int mono_cell_offset = 0; mono_cell_offset < N; mono_cell_offset++)
        {
            if(poly_cell_indices[poly  * N + mono_cell_offset] < 0) break; //check if end has been reached
            int64_t cell = poly_cell_indices[poly * N + mono_cell_offset];
            uint64_t beads_in_cell = 0; //total number of beads in the cell
            //get beads in cell
            for(unsigned int monotype = 0; monotype < p->n_types; monotype++) beads_in_cell += p->fields_unified[monotype*p->n_cells_local + cell];
            //loop over monotypes to get delta_cost
            for(unsigned int monotype = 0; monotype < p->n_types; monotype++)
                {
                    if(p->umbrella_field[monotype*p->n_cells_local + cell] > 0)
                        {
                            //number of monomers of current type in current cell before and after flip
                            int64_t num_mono_initial = poly_cell_num[poly * p->n_poly_type * p->n_types * N + initial_type * p->n_types * N + monotype * N + mono_cell_offset];
                            int64_t num_mono_final = poly_cell_num[poly * p->n_poly_type * p->n_types * N + final_type * p->n_types * N + monotype * N + mono_cell_offset];
                            //difference of the two
                            int64_t delta_num_mono = num_mono_final - num_mono_initial;
                            //subtract old value
                            delta_cost-=powl((soma_scalar_t)p->umbrella_field[monotype*p->n_cells_local + cell]-(soma_scalar_t)(p->fields_unified[monotype*p->n_cells_local + cell] + delta_fields_unified[monotype*p->n_cells_local + cell]) / beads_in_cell,2.0);
                            //add new value
                            delta_cost+=powl((soma_scalar_t)p->umbrella_field[monotype*p->n_cells_local + cell]-(soma_scalar_t)(p->fields_unified[monotype*p->n_cells_local + cell] + delta_fields_unified[monotype*p->n_cells_local + cell] + delta_num_mono)/beads_in_cell,2.0);
                        }
                   
                }
        }
        
        
    return delta_cost;
}

int comp (const void * elem1, const void * elem2) 
{
    int f = *((uint64_t*)elem1);
    int s = *((uint64_t*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

unsigned int flip(struct Phase * p, uint64_t poly)
{
    unsigned int initial_poly_type = p->polymers[poly].type;
    unsigned int * target_types=(unsigned int *)malloc( (p->n_poly_type - 1) *  sizeof(unsigned int));
    unsigned int k=0;
    Polymer *mypoly = p->polymers + poly; //needed only for rng
    //get array with every type except current one
    for(unsigned int polytype=0; polytype < p->n_poly_type;polytype++)
        {
            if(polytype != initial_poly_type)
                {
                    target_types[k]=polytype;
                    k++;
                }
        }
    unsigned int random_idx = (unsigned int)soma_rng_soma_scalar(&(mypoly->poly_state), p)* (soma_scalar_t)(k - 1);
    unsigned int final_type=target_types[random_idx];
    free(target_types);
/*     printf("Initial type : %u\n",initial_poly_type);
    printf("Final type : %u\n",final_type); */
    return final_type;
}

