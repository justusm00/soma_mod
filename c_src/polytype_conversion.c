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
                    return convert_target_4(p);
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

int convert_target_old(struct Phase *p)
{
    
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }

    soma_scalar_t diff=0; //cumulative deviation after swap
    soma_scalar_t diff_old=0; //cumulative deviation before swap
    const unsigned int N = p->reference_Nbeads;
    //compute cumulative deviation from target density
    for (uint64_t type = 0; type < p->n_types; type++)
        {
            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                {
                if (p->umbrella_field[type*p->n_cells_local + cell] > 0)
                    {
                        diff_old=diff_old+fabs(p->umbrella_field[type*p->n_cells_local + cell]-p->fields_unified[type*p->n_cells_local + cell] * p->field_scaling_type[type]);
                    }
                }
        }

    //loop over polymers
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono += 1)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if target density is available
                    if (p->umbrella_field[mono_cell] > 0)
                        {
                            //flip polymer type
                            if(p->polymers[poly].type==0) p->polymers[poly].type=1;
                            else p->polymers[poly].type=0;
                            //check if this was smart move
                            //calculate new density
                            update_density_fields(p);
                            for (uint64_t type = 0; type < p->n_types; type++)
                                {
                                    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                                        {
                                        if (p->umbrella_field[type*p->n_cells_local + cell] > 0)
                                            {
                                                diff=diff+fabs(p->umbrella_field[type*p->n_cells_local + cell]-p->fields_unified[type*p->n_cells_local + cell] * p->field_scaling_type[type]);
                                            }
                                        }
                                }
                            if (diff>diff_old)
                                {
                                    //revert swap
                                    if(p->polymers[poly].type==0) p->polymers[poly].type=1;
                                    else p->polymers[poly].type=0;
                                    update_density_fields(p);
                                    diff=diff_old;
                                }
                            diff_old=diff;
                            diff=0;
                            break;
                        }


                        
                }

        }  
    return 0;
}


int convert_target_1(struct Phase *p)
{
    
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }

    soma_scalar_t diff=0; //cumulative deviation 
    soma_scalar_t diff_flip=0; //cumulative deviation with flip
    unsigned int mono_counter=0; //count monomers of polymer in specific cell
    const unsigned int N = p->reference_Nbeads;
    unsigned int type = 0;
    unsigned int mono_offset = 0; //number of monomers in cells for which target density is not avaiable
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //array that contains monomer cell indices. Values are -1 if no target density available in that cell



    //loop over polymers
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            mono_offset=0;
            diff=0;
            diff_flip=0;
            // get polymer type
            type = p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }


            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {

                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);

                    //compute cumulative density deviation in those cells and check how it changes if monomers are flipped
                    for (unsigned int mono = mono_offset; mono < N - 1; mono ++)
                        {
                            
                            mono_counter++;
                            //update difference only for new cells
                            if(mono_cells[poly*N+mono]!= mono_cells[mono+1])
                                {

                                    for(int types =0; types < p->n_types; types++)
                                        {
                                            diff=diff+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-p->fields_unified[types*p->n_cells_local + mono_cells[mono]] * p->field_scaling_type[types]);
                                            if(types == type)
                                                {
                                                    diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] - mono_counter )* p->field_scaling_type[types]);
                                                }
                                            else
                                                {
                                                    diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] + mono_counter )* p->field_scaling_type[types]);
                                                }
                                            mono_counter=0;
                                        }
                                }
                        }
                    //final monomer
                    if (mono_counter != 0)
                        {
                            mono_counter++;
                        }
                    else
                        {
                            mono_counter=1;
                        }
                    for(int types =0; types < p->n_types; types++)
                        {
                            diff=diff+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] * p->field_scaling_type[types]);
                            if(types == type)
                                {
                                    diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[type*p->n_cells_local + mono_cells[N-1]] - mono_counter )* p->field_scaling_type[types]);
                                }
                            else
                                {
                                    diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] + mono_counter )* p->field_scaling_type[types]);
                                }
                        }


                    mono_counter=0;


                    //finally, check if flip leads to improvement and update density fields
                    if(diff_flip < diff)
                        {
                            if(type==0) p->polymers[poly].type=1;
                            else p->polymers[poly].type=0;

                            //update density fields
                            for (unsigned int mono = mono_offset; mono < N - 1; mono ++)
                                {
                                    mono_counter++;
                                    //update difference only for new cells
                                    if(mono_cells[poly*N+mono]!= mono_cells[mono+1])
                                        {

                                            for(int types =0; types < p->n_types; types++)
                                                {
                                                    if(types == type)
                                                        {
                                                            p->fields_unified[types*p->n_cells_local + mono_cells[mono]] -= mono_counter;
                                                        }
                                                    else
                                                        {
                                                            p->fields_unified[types*p->n_cells_local + mono_cells[mono]] += mono_counter;
                                                        }
                                                    mono_counter=0;
                                                }
                                        }
                                }
                            //final monomer
                            if (mono_counter != 0)
                                {
                                    mono_counter++;
                                }
                            else
                                {
                                    mono_counter=1;
                                }
                            for(int types =0; types < p->n_types; types++)
                                {
                                    if(types == type)
                                        {
                                            p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] -= mono_counter;
                                        }
                                    else
                                        {
                                            p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] += mono_counter;
                                        }
                                }




                        }

                }
                    
        
        }



    free(mono_cells);
    return 0;
}


int convert_target_2(struct Phase *p)
{
    
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }

    soma_scalar_t diff=0; //cumulative deviation 
    soma_scalar_t diff_flip=0; //cumulative deviation with flip
    unsigned int mono_counter=0; //count monomers of polymer in specific cell
    unsigned int k=0; //counter variable
    const unsigned int N = p->reference_Nbeads;
    unsigned int type = 0;
    unsigned int mono_offset = 0; //number of monomers in cells for which target density is not avaiable
    unsigned int mono_cell =0;
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //monomer cell indices. Values are -1 if no target density available in that cell
    int64_t * mono_cells_unique=(int64_t *)malloc( N* sizeof(int64_t)); //unique monomer cell indices
    int64_t * num_mono=(int64_t *)malloc( N* sizeof(int64_t)); //number of monomers corresponding to unique cells



    //loop over polymers
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            mono_offset=0;
            mono_cell=0;
            diff=0;
            diff_flip=0;
            mono_counter=0;
            //get polymer type
            type=p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }

            

            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {
             
                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);
                    //get unique cells and number of monomers in them
                    k=0;
                    for(unsigned int mono = mono_offset; mono < N-1 ; mono++)
                        {
                            mono_counter++;
                            //update array only for new cells
                            if(mono_cells[mono]!= mono_cells[mono+1])
                            
                                {
                                    mono_cells_unique[k]=mono_cells[mono];
                                    num_mono[k]=mono_counter;
                                    mono_counter=0;
                                    k++;
                                }
                        }
   
                    //final monomer
                    if (mono_counter == 0) num_mono[k]=1;
               
                    else 
                        {
                            mono_counter++;
                            num_mono[k]=mono_counter;
                        }
                    mono_cells_unique[k]=mono_cells[N-1];
                    k++;
                    //set end of arrays
                    if(k<N-1) 
                        {
                            mono_cells_unique[k]=-1;
                            num_mono[k]=-1;
                        } 


                    //compute cumulative density deviation in unique cells and check how it changes if monomers are flipped
                    for (unsigned int mono = 0; mono < N; mono ++)
                        {
                            if(mono_cells_unique[mono] < 0) break;
                            mono_cell=mono_cells_unique[mono];
                            mono_counter=num_mono[mono];
                            //update diff and diff_flip
                            for(int types =0; types < p->n_types; types++)
                                {
                                    diff=diff+fabs(p->umbrella_field[types*p->n_cells_local + mono_cell]-p->fields_unified[types*p->n_cells_local + mono_cell] * p->field_scaling_type[types]);
                                    if(types == type)
                                        {
                                            diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cell]-( p->fields_unified[types*p->n_cells_local + mono_cell] - mono_counter )* p->field_scaling_type[types]);
                                        }
                                    else
                                        {
                                            diff_flip=diff_flip+fabs(p->umbrella_field[types*p->n_cells_local + mono_cell]-( p->fields_unified[types*p->n_cells_local + mono_cell] + mono_counter )* p->field_scaling_type[types]);
                                        }
                                }
                        }


                    //finally, check if flip leads to improvement and update density fields
                    if(diff_flip < diff)
                        {
                            if(type==0) p->polymers[poly].type=1;
                            else p->polymers[poly].type=0;
                            

                            //update density fields
                            for (unsigned int mono = 0; mono < N; mono++)
                                {
                                    if(mono_cells_unique[mono] < 0) break;
                                    mono_cell=mono_cells_unique[mono];
                                    mono_counter=num_mono[mono];
                                    for(int types =0; types < p->n_types; types++)
                                        {
                                            if(types == type)
                                                {
                                                    p->fields_unified[types*p->n_cells_local + mono_cell] -= mono_counter;
                                                }
                                            else
                                                {
                                                    p->fields_unified[types*p->n_cells_local + mono_cell] += mono_counter;
                                                }
                                        }
                                        
                                }




                        }

                }
        }


    
    free(mono_cells);
    free(mono_cells_unique);
    free(num_mono);
    return 0;
}


int print_info(struct Phase *p)
{
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }

    const unsigned int N = p->reference_Nbeads; //monomers per polymer 
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //monomer cell indices. Values are -1 if no target density available in that cell
    int64_t * poly_flippable = (int64_t *)malloc( p->n_polymers* sizeof(int64_t)); //boolean array that stores whether or not polymer has monomers in target density area
    int64_t * poly_cell_indices = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores in which cells a given polymer has monomers
    int64_t * poly_cell_num = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores number of monomers in cells. Values correspond to cells specified in poly_cell_indices



    //loop over polymers to identify the ones that may be flipped
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            unsigned int mono_offset=0;
            unsigned int mono_cell=0;
            unsigned int mono_counter=0;
            //get polymer type
            unsigned int type = p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }

            

            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {
                    poly_flippable[poly]=1; //1 means that polymer has monomers in target density area
                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);
                    //get unique cells and number of monomers in them
                    unsigned int k=0;
                    for(unsigned int mono = mono_offset; mono < N-1 ; mono++)
                        {
                            mono_counter++;
                            //update array only for new cells
                            if(mono_cells[mono]!= mono_cells[mono+1])
                            
                                {
                                    poly_cell_indices[poly * N + k]=mono_cells[mono]; //unique monomer cell indices
                                    poly_cell_num[poly * N + k]=mono_counter; //corresponding number of monomers
                                    mono_counter=0;
                                    k++;
                                }
                        }
   
                    //final monomer
                    if (mono_counter == 0) poly_cell_num[poly * N + k]=1;
               
                    else 
                        {
                            mono_counter++;
                            poly_cell_num[poly * N + k]=mono_counter; 
                        }
                    poly_cell_indices[poly * N + k]=mono_cells[N-1];
                    k++;

                    //set end of arrays 
                    if(k<N-1) 
                        {
                            poly_cell_indices[poly * N + k]=-1;
                            poly_cell_num[poly * N + k]=-1;
                        } 
                }
            else poly_flippable[poly]=0; //0 means that polymer has no monomers in target density area
        }


    //print density field for use outside soma
    for(uint64_t type = 0; type<p->n_poly_type; type++)
        {
            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                {
                    printf("%f ", p->fields_unified[type*p->n_cells_local + cell] * p->field_scaling_type[type]);
                }
            printf("\n");
        }

    //print flippable polymer info for use outside soma
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            
            if(poly_flippable[poly]==1)
                {
                    printf("%d\n",p->polymers[poly].type);
                    for(int i = 0; i < N ; i++)
                        {
                            if(poly_cell_indices[poly * N + i] < 0) break;
                            printf("%d ",poly_cell_indices[poly * N + i]);
                        }   
                    printf("\n");
                    for(int i = 0; i < N ; i++)
                        {
                            if(poly_cell_indices[poly * N + i] < 0) break;
                            printf("%d ",poly_cell_num[poly * N + i]);
                        }   
                    printf("\n");
                }
        }





    free(mono_cells);
    free(poly_flippable);
    free(poly_cell_indices);
    free(poly_cell_num);
    return 0;
}



int convert_target_3(struct Phase *p)
{
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }

    


    // ARRAY AND VARIABLE DECLARATION
    
    const unsigned int N = p->reference_Nbeads; //monomers per polymer
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //monomer cell indices. Values are -1 if no target density available in that cell
    int64_t * poly_flippable = (int64_t *)malloc( p->n_polymers* sizeof(int64_t)); //boolean array that stores whether or not polymer has monomers in target density area
    int64_t * poly_flippable_indices = (int64_t *)malloc( p->n_polymers* sizeof(int64_t)); //array that contains indices of flippable polymers
    int64_t * poly_cell_indices = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores in which cells a given polymer has monomers
    int64_t * poly_cell_num = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores number of monomers in cells. Values correspond to cells specified in poly_cell_indices
    int64_t * delta_fields_unified = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); //array that stores changes in density
    int64_t * delta_fields_unified_best = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); 
    unsigned int * poly_types=(uint64_t *)malloc(p->n_polymers * sizeof(uint64_t)); //array that stores polymer types
    unsigned int * poly_types_best=(uint64_t *)malloc(p->n_polymers * sizeof(uint64_t));
    uint64_t num_poly_flippable = 0; // number of flippable polymers
    soma_scalar_t acc_rate= 1.0; //initial flip acceptance rate
    uint64_t num_iter=0; 
    uint64_t num_acc = 0;
    soma_scalar_t total_cost = 0.0;
    soma_scalar_t total_cost_old = 0.0;
    soma_scalar_t total_cost_best = 0.0;


    // SIMULATED ANNEALING PARAMETERS 

    soma_scalar_t Tmax = 0.001;
    soma_scalar_t Tmin = 0.0001;
    soma_scalar_t alpha = 0.85;
    uint64_t max_iter = 10000; //maximum number of flips 
    soma_scalar_t acc_rate_target=0.02; //flip acceptance rate after which converison will be stopped

    srand(time(0)); //rng for selecting polymers to flip



    // INITIALIZE delta_fields_unified

    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            for(uint64_t type = 0; type < p->n_types; type++)
                {
                    delta_fields_unified[type*p->n_cells_local + cell] = 0;
                    delta_fields_unified_best[type*p->n_cells_local + cell] = 0;
                }
        }

    // INTIALIZE poly_types
    
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            poly_types[poly]=p->polymers[poly].type;
            poly_types_best[poly]=p->polymers[poly].type;
        }


    // LOOP OVER POLYMERS TO IDENTIFY THE ONES THAT MAY BE FLIPPED
    
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            
            unsigned int mono_offset=0;
            unsigned int mono_cell=0;
            unsigned int mono_counter=0;
            //get polymer type
            unsigned int type = p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }

            

            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {
                    poly_flippable_indices[num_poly_flippable]=poly;
                    poly_flippable[poly]=1; //1 means that polymer has monomers in target density area
                    num_poly_flippable ++;
                    
                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);
                    //get unique cells and number of monomers in them
                    unsigned int k=0;
                    for(unsigned int mono = mono_offset; mono < N-1 ; mono++)
                        {
                            mono_counter++;
                            //update array only for new cells
                            if(mono_cells[mono]!= mono_cells[mono+1])
                            
                                {
                                    poly_cell_indices[poly * N + k]=mono_cells[mono]; //unique monomer cell indices
                                    poly_cell_num[poly * N + k]=mono_counter; //corresponding number of monomers
                                    mono_counter=0;
                                    k++;
                                }
                        }
   
                    //final monomer
                    if (mono_counter == 0) poly_cell_num[poly * N + k]=1;
               
                    else 
                        {
                            mono_counter++;
                            poly_cell_num[poly * N + k]=mono_counter; 
                        }
                    poly_cell_indices[poly * N + k]=mono_cells[N-1];
                    k++;

                    //set end of arrays 
                    if(k<N-1) 
                        {
                            poly_cell_indices[poly * N + k]=-1;
                            poly_cell_num[poly * N + k]=-1;
                        } 
                }
            else poly_flippable[poly]=0; //0 means that polymer has no monomers in target density area
        }



    //initialize cost
    total_cost=get_cost(p, delta_fields_unified);
    total_cost_old = total_cost;
    total_cost_best= total_cost;

    // SIMULATED ANNEALING

    //printf("Total cost before : %f \n",total_cost);
    while((acc_rate > acc_rate_target) && (num_iter < max_iter))
        {
            soma_scalar_t T = Tmax;
            while(T > Tmin)
                {
                    num_iter++;
                    //choose random polymer to flip
                    uint64_t random_index = rand() % (num_poly_flippable - 1);
                    uint64_t poly = poly_flippable_indices[random_index];
                    Polymer *mypoly = p->polymers + poly;
                    unsigned int initial_type = poly_types[poly];
                    unsigned int final_type = flip(initial_type);
                    total_cost=total_cost_old;
                    //calculate cost (only need to update it for the cells in which the polymer has monomers)
                    for(unsigned int i = 0; i < N; i++)
                        {
                            if(poly_cell_indices[poly * N + i] < 0) break;
                            unsigned int cell = poly_cell_indices[poly * N + i];
                            unsigned int num_mono = poly_cell_num[poly * N + i];
                            for(unsigned int type = 0; type < p->n_types; type++)
                                {
                                    total_cost-=pow(p->umbrella_field[type*p->n_cells_local + cell]-( ( p->fields_unified[type*p->n_cells_local + cell] + delta_fields_unified[type*p->n_cells_local + cell])* p->field_scaling_type[type]),2);
                                }
                            total_cost+=pow(p->umbrella_field[initial_type*p->n_cells_local + cell]-( ( p->fields_unified[initial_type*p->n_cells_local + cell] + delta_fields_unified[initial_type*p->n_cells_local + cell] - num_mono)* p->field_scaling_type[initial_type]),2);
                            total_cost+=pow(p->umbrella_field[final_type*p->n_cells_local + cell]-( ( p->fields_unified[final_type*p->n_cells_local + cell] + delta_fields_unified[final_type*p->n_cells_local + cell] + num_mono)* p->field_scaling_type[final_type]),2);
                        }
                    if(total_cost < total_cost_old)
                        {
                            //accept flip
                            num_acc++;
                            poly_types[poly]=final_type;
                            total_cost_old=total_cost;
                            //update delta fields unified
                            for(unsigned int i = 0; i < N; i++)
                                {
                                    if(poly_cell_indices[poly * N + i] < 0) break;
                                    unsigned int cell = poly_cell_indices[poly * N + i];
                                    unsigned int num_mono = poly_cell_num[poly * N + i];
                                    delta_fields_unified[initial_type*p->n_cells_local + cell]-=num_mono;
                                    delta_fields_unified[final_type*p->n_cells_local + cell]+=num_mono;
                                }
                        }
                    else
                        {
                            //accept flip with probability 
                            soma_scalar_t random_number = soma_rng_soma_scalar(&(mypoly->poly_state), p);
                            if(random_number < exp(-(total_cost - total_cost_old)/T))
                                {
                                    //accept
                                    num_acc++;
                                    poly_types[poly]=final_type;
                                    total_cost_old=total_cost;
                                    //update delta fields unified
                                    for(unsigned int i = 0; i < N; i++)
                                        {
                                            if(poly_cell_indices[poly * N + i] < 0) break;
                                            unsigned int cell = poly_cell_indices[poly * N + i];
                                            unsigned int num_mono = poly_cell_num[poly * N + i];
                                            delta_fields_unified[initial_type*p->n_cells_local + cell]-=num_mono;
                                            delta_fields_unified[final_type*p->n_cells_local + cell]+=num_mono;
                                        }
                                    }
                            //else reject
                            else total_cost = total_cost_old;
                        }

                    //update best solution so far
                    if (total_cost < total_cost_best)
                        {
                            total_cost_best=total_cost;
                            poly_types_best[poly]=poly_types[poly];
                            for(unsigned int i = 0; i < N; i++)
                                {
                                    if(poly_cell_indices[poly * N + i] < 0) break;
                                    unsigned int cell = poly_cell_indices[poly * N + i];
                                    for(unsigned int type = 0; type < p->n_types; type++) delta_fields_unified_best[type*p->n_cells_local + cell] = delta_fields_unified[type*p->n_cells_local + cell];
                                }
                        }
                    //update temperature
                    T *= alpha;
                }
            acc_rate=(float)(num_acc)/(float)(num_iter);
            //set everything to its best values
            for (uint64_t poly = 0; poly < p->n_polymers; poly++) poly_types[poly]=poly_types_best[poly];
            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                {
                    for(uint64_t type = 0; type < p->n_types; type++)
                        {
                            delta_fields_unified[type*p->n_cells_local + cell] = delta_fields_unified_best[type*p->n_cells_local + cell];
                        }
                }
            total_cost=total_cost_best;
            

        }


        //printf("Total cost after: %f \n",total_cost_best);
        //update density fields
        for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
            {
                for(uint64_t type = 0; type < p->n_types; type++)
                    {
                        
                        p->fields_unified[type*p->n_cells_local + cell] += delta_fields_unified_best[type*p->n_cells_local + cell];
                    }
            }

        //update polymer types
        for (uint64_t poly = 0; poly < p->n_polymers; poly++) p->polymers[poly].type=poly_types_best[poly];


    free(mono_cells);
    free(poly_flippable);
    free(poly_cell_indices);
    free(poly_cell_num);
    free(poly_types);
    free(poly_types_best);
    free(delta_fields_unified);
    free(delta_fields_unified_best);
    free(poly_flippable_indices);
    return 0;
}

int convert_target_4(struct Phase *p)
{
    if (p->n_types != 2)
        {
            printf("ERROR: convert_target requires exactly two polymer types.\n");
            return -1;
        }
 

    // ADJUSTABLE PARAMETERS

    soma_scalar_t Tmax = 0.001; // maximum SA temperature
    soma_scalar_t Tmin = 0.0001; // minimum SA temperature
    soma_scalar_t alpha = 0.85; // SA temperature decrease factor
    uint64_t max_sa_runs = 10000; //maximum simulated annealing runs
    soma_scalar_t sa_acc_rate_target=0.5; //sa acceptance rate after which converison will be stopped
    int64_t sa_buffer_size = p->n_polymers; //maximum number of polymers flipped in one SA run. needed to update best values in the end.
    int64_t flip_buffer_size = p->n_polymers; //maximum number of flippable polymers, need to find optimal value

    // INITIALIZE PARAMETERS

    uint64_t num_poly_flippable = 0; // number of flippable polymers
    soma_scalar_t sa_acc_rate= 1.0; // SA acceptance rate
    uint64_t num_sa_runs=0; // counts number of simulated annealing runs
    uint64_t num_sa_runs_acc = 0; // number of simulated annealing runs that lead to an improvement
    uint64_t flip_counter =0; //counts number of accepted flips within one SA iteration

    // cost values used inside SA loop
    soma_scalar_t total_cost = 0.0; // cost 
    soma_scalar_t total_cost_old = 0.0;

    // cost values used outside SA loop
    soma_scalar_t total_cost_best = 0.0; //best value after SA run
    soma_scalar_t total_cost_best_old = 0.0; //best value before SA run


    // ARRAYS
    
    const unsigned int N = p->reference_Nbeads; //monomers per polymer
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //monomer cell indices. Values are -1 if no target density available in that cell
    int64_t * poly_flippable = (int64_t *)malloc( p->n_polymers* sizeof(int64_t)); //boolean array that stores whether or not polymer has monomers in target density area
    int64_t * poly_cell_indices = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores in which cells a given polymer has monomers
    int64_t * poly_cell_num = (int64_t *)malloc(p->n_polymers * N * sizeof(int64_t)); //array that stores number of monomers in cells. Values correspond to cells specified in poly_cell_indices
    int64_t * delta_fields_unified = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); //array that stores changes in density
    int64_t * delta_fields_unified_best = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); 
    int64_t * poly_flippable_indices = (int64_t *)malloc( flip_buffer_size * sizeof(int64_t)); //array that contains indices of flippable polymers
    int64_t * poly_flipped_indices = (int64_t *)malloc(sa_buffer_size * sizeof(int64_t)); //contains indices of flipped polymers after sa run is done
    unsigned int * poly_types=(uint64_t *)malloc(flip_buffer_size * sizeof(uint64_t)); //array that stores polymer types
    unsigned int * poly_types_best=(uint64_t *)malloc(flip_buffer_size * sizeof(uint64_t));


    // RNG FOR POLYMER FLIP SELECTION
    srand(time(0));

    // LOOP OVER POLYMERS TO IDENTIFY THE ONES THAT MAY BE FLIPPED
    // probably the most expensive part, needs to run in parallel
    

    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            
            unsigned int mono_offset=0;
            unsigned int mono_cell=0;
            unsigned int mono_counter=0;
            //get polymer type
            unsigned int type = p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }

            

            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {
                    poly_flippable[poly]=1; //1 means that polymer has monomers in target density area
                    
                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);
                    //get unique cells and number of monomers in them
                    unsigned int k=0;
                    for(unsigned int mono = mono_offset; mono < N-1 ; mono++)
                        {
                            mono_counter++;
                            //update array only for new cells
                            if(mono_cells[mono]!= mono_cells[mono+1])
                            
                                {
                                    poly_cell_indices[poly * N + k]=mono_cells[mono]; //unique monomer cell indices
                                    poly_cell_num[poly * N + k]=mono_counter; //corresponding number of monomers
                                    mono_counter=0;
                                    k++;
                                }
                        }
   
                    //final monomer
                    if (mono_counter == 0) poly_cell_num[poly * N + k]=1;
               
                    else 
                        {
                            mono_counter++;
                            poly_cell_num[poly * N + k]=mono_counter; 
                        }
                    poly_cell_indices[poly * N + k]=mono_cells[N-1];
                    k++;

                    //set end of arrays 
                    if(k<N-1) 
                        {
                            poly_cell_indices[poly * N + k]=-1;
                            poly_cell_num[poly * N + k]=-1;
                        } 
                }
            else poly_flippable[poly]=0; //0 means that polymer has no monomers in target density area
        }

    //save flippable polymer indices to new array, this must be done sequentially
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            if(poly_flippable[poly]==1)
                {
                    poly_flippable_indices[num_poly_flippable]=poly;
                    num_poly_flippable++;

                }
        }

    //check if there are more flippable polymers than the buffer allows
    if(num_poly_flippable>flip_buffer_size)
        {
            printf("ERROR: FLIP BUFFER SIZE IS TOO SMALL\n");
            return -1;
        }


    // INITIALIZE delta_fields_unified

    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            for(uint64_t type = 0; type < p->n_types; type++)
                {
                    delta_fields_unified[type*p->n_cells_local + cell] = 0;
                    delta_fields_unified_best[type*p->n_cells_local + cell] = 0;
                }
        }



    // INTIALIZE poly_types
    
    for (uint64_t poly = 0; poly < num_poly_flippable; poly++)
        {
            poly_types[poly]=p->polymers[poly_flippable_indices[poly]].type;
            poly_types_best[poly]=poly_types[poly];
        }
    
    //initialize cost
    total_cost=get_cost(p, delta_fields_unified);
    total_cost_old = total_cost;
    total_cost_best= total_cost;

    // SIMULATED ANNEALING
    printf("Total cost before: %f \n",total_cost_best);
    while((sa_acc_rate > sa_acc_rate_target) && (num_sa_runs < max_sa_runs))
        {
            soma_scalar_t T = Tmax;
            num_sa_runs++; //increment SA counter
            total_cost_best_old=total_cost_best;




            // ############# START SA LOOP ##############

            while(T > Tmin)
                {
                    if(flip_counter >= sa_buffer_size)
                        {
                            printf("ERROR: SA BUFFER TOO SMALL\n");
                            return -1;
                        }
                    //choose random polymer to flip
                    uint64_t random_index = rand() % (num_poly_flippable - 1);
                    uint64_t poly = poly_flippable_indices[random_index];
                    Polymer *mypoly = p->polymers + poly;
                    unsigned int initial_type = poly_types[random_index];
                    unsigned int final_type = flip(initial_type);
                    total_cost=total_cost_old;
                    //calculate cost (only need to update it for the cells in which the polymer has monomers)
                    for(unsigned int i = 0; i < N; i++)
                        {
                            if(poly_cell_indices[poly * N + i] < 0) break;
                            unsigned int cell = poly_cell_indices[poly * N + i];
                            unsigned int num_mono = poly_cell_num[poly * N + i];
                            for(unsigned int type = 0; type < p->n_types; type++)
                                {
                                    total_cost-=pow(p->umbrella_field[type*p->n_cells_local + cell]-( ( p->fields_unified[type*p->n_cells_local + cell] + delta_fields_unified[type*p->n_cells_local + cell])* p->field_scaling_type[type]),2);
                                }
                            total_cost+=pow(p->umbrella_field[initial_type*p->n_cells_local + cell]-( ( p->fields_unified[initial_type*p->n_cells_local + cell] + delta_fields_unified[initial_type*p->n_cells_local + cell] - num_mono)* p->field_scaling_type[initial_type]),2);
                            total_cost+=pow(p->umbrella_field[final_type*p->n_cells_local + cell]-( ( p->fields_unified[final_type*p->n_cells_local + cell] + delta_fields_unified[final_type*p->n_cells_local + cell] + num_mono)* p->field_scaling_type[final_type]),2);
                        }
                    if(total_cost < total_cost_old)
                        {
                            //accept flip
                            poly_flipped_indices[flip_counter]=random_index;
                            flip_counter++;
                            poly_types[random_index]=final_type;
                            total_cost_old=total_cost;
                            //update delta fields unified
                            for(unsigned int i = 0; i < N; i++)
                                {
                                    if(poly_cell_indices[poly * N + i] < 0) break;
                                    unsigned int cell = poly_cell_indices[poly * N + i];
                                    unsigned int num_mono = poly_cell_num[poly * N + i];
                                    delta_fields_unified[initial_type*p->n_cells_local + cell]-=num_mono;
                                    delta_fields_unified[final_type*p->n_cells_local + cell]+=num_mono;
                                }
                        }
                    else
                        {
                            //accept flip with probability 
                            soma_scalar_t random_number = soma_rng_soma_scalar(&(mypoly->poly_state), p);
                            if(random_number < exp(-(total_cost - total_cost_old)/T))
                                {
                                    //accept
                                    poly_flipped_indices[flip_counter]=random_index;
                                    flip_counter++;
                                    poly_types[random_index]=final_type;
                                    total_cost_old=total_cost;
                                    //update delta fields unified
                                    for(unsigned int i = 0; i < N; i++)
                                        {
                                            if(poly_cell_indices[poly * N + i] < 0) break;
                                            unsigned int cell = poly_cell_indices[poly * N + i];
                                            unsigned int num_mono = poly_cell_num[poly * N + i];
                                            delta_fields_unified[initial_type*p->n_cells_local + cell]-=num_mono;
                                            delta_fields_unified[final_type*p->n_cells_local + cell]+=num_mono;
                                        }
                                    }
                            //else reject
                            else total_cost = total_cost_old;
                        }


                    //update best solution so far
                    if (total_cost < total_cost_best)
                        {
                            total_cost_best=total_cost;
                            for(int64_t polyy = 0; polyy < flip_counter; polyy++)
                                {
                                    poly=poly_flipped_indices[polyy]; //index in flippable arrays
                                    poly_types_best[poly]=poly_types[poly];
                                    poly = poly_flippable_indices[poly]; //index in polymers array
                                    for(unsigned int i = 0; i < N; i++)
                                        {
                                            if(poly_cell_indices[poly * N + i] < 0) break;
                                            unsigned int cell = poly_cell_indices[poly * N + i];
                                            for(unsigned int type = 0; type < p->n_types; type++) delta_fields_unified_best[type*p->n_cells_local + cell] = delta_fields_unified[type*p->n_cells_local + cell];
                                        }
                                }
                            //reset flip counter
                            flip_counter=0;

                        }
                    //update temperature
                    T *= alpha;
                }

                // ############# END  SA LOOP ##############
            



            //set everything to its best values
            if(total_cost_best<total_cost_best_old)
                {
                    num_sa_runs_acc++;
                    for (uint64_t poly = 0; poly < num_poly_flippable; poly++) poly_types[poly]=poly_types_best[poly];
                    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                        {
                            for(uint64_t type = 0; type < p->n_types; type++)
                                {
                                    delta_fields_unified[type*p->n_cells_local + cell] = delta_fields_unified_best[type*p->n_cells_local + cell];
                                }
                        }
                    total_cost=total_cost_best;
                }

            //update sa acceptance rate every 10 steps 
            if((num_sa_runs % 10) == 0) sa_acc_rate=(float)(num_sa_runs_acc)/(float)(num_sa_runs);

        }

    printf("Simulated annealing iterations: %d\n",num_sa_runs);
    printf("Total cost after: %f \n",total_cost_best);


    //update polymer types
    for(int64_t polyy = 0; polyy < num_poly_flippable; polyy++) p->polymers[poly_flippable_indices[polyy]].type=poly_types_best[polyy];


    free(mono_cells);
    free(poly_flippable);
    free(poly_cell_indices);
    free(poly_cell_num);
    free(poly_types);
    free(poly_types_best);
    free(delta_fields_unified);
    free(delta_fields_unified_best);
    free(poly_flippable_indices);
    free(poly_flipped_indices);
    return 0;
}



int comp (const void * elem1, const void * elem2) 
{
    int f = *((uint64_t*)elem1);
    int s = *((uint64_t*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int flip(int initial_type)
{
    if(initial_type == 0) return 1;
    else return 0;
}


void get_flip_candidates(struct Phase *p, int64_t * poly_flippable_indices)

soma_scalar_t get_cost(struct Phase *p, int64_t * delta_fields_unified)
{   
    soma_scalar_t total_cost=0;
    //loop over cells
    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        {
            //loop over types
            for(uint64_t type = 0; type < p->n_types; type++)
                {
                    if(p->umbrella_field[type*p->n_cells_local + cell] > 0)
                        {
                            total_cost+=pow(p->umbrella_field[type*p->n_cells_local + cell]-( ( p->fields_unified[type*p->n_cells_local + cell] + delta_fields_unified[type*p->n_cells_local + cell])* p->field_scaling_type[type]),2);
                        }
                }
        }
    return total_cost;
}

void get_flip_prob(struct Phase *p, soma_scalar_t * poly_flip)
{
    soma_scalar_t diff=0; //cumulative deviation 
    soma_scalar_t diff_flip=0; //cumulative deviation with flip
    unsigned int mono_counter=0; //count monomers of polymer in specific cell
    const unsigned int N = p->reference_Nbeads;
    unsigned int type = 0;
    unsigned int mono_offset = 0; //number of monomers in cells for which target density is not avaiable
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //array that contains monomer cell indices. Values are -1 if no target density available in that cell


    //loop over polymers
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            mono_offset=0;
            diff=0;
            diff_flip=0;
            // get polymer type
            type = p->polymers[poly].type;
            //loop over monomers
            for (unsigned int mono = 0; mono < N; mono ++)
                {
                    const Monomer pos = ((Monomer *) p->ph.beads.ptr)[p->polymers[poly].bead_offset + mono];       //Read Monomer position
                    const uint64_t mono_cell = coord_to_index(p, pos.x, pos.y, pos.z);    //cell of current monomer
                    //check if umbrella field exists in monomer cell
                    if(p->umbrella_field[type*p->n_cells_local + mono_cell] > 0)
                        {
                            mono_cells[mono]=mono_cell;
                        }
                    else 
                        {
                            mono_cells[mono]=-1;
                            mono_offset++;
                        }
                }


            //check if there are monomers in cells with target density
            if(mono_offset != N)
                {

                    //sort mono_cells array
                    qsort(mono_cells,N,sizeof(int64_t),comp);

                    //compute cumulative density deviation in those cells and check how it changes if monomers are flipped
                    for (unsigned int mono = mono_offset; mono < N - 1; mono ++)
                        {
                            
                            mono_counter++;
                            //update difference only for new cells
                            if(mono_cells[poly*N+mono]!= mono_cells[mono+1])
                                {

                                    for(int types =0; types < p->n_types; types++)
                                        {
                                            diff=diff+(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-p->fields_unified[types*p->n_cells_local + mono_cells[mono]] * p->field_scaling_type[types])*(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-p->fields_unified[types*p->n_cells_local + mono_cells[mono]] * p->field_scaling_type[types]);
                                            if(types == type)
                                                {
                                                    diff_flip=diff_flip+(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] - mono_counter )* p->field_scaling_type[types])*(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] - mono_counter )* p->field_scaling_type[types]);
                                                }
                                            else
                                                {
                                                    diff_flip=diff_flip+(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] + mono_counter )* p->field_scaling_type[types])*(p->umbrella_field[types*p->n_cells_local + mono_cells[mono]]-( p->fields_unified[types*p->n_cells_local + mono_cells[mono]] + mono_counter )* p->field_scaling_type[types]);
                                                }
                                            mono_counter=0;
                                        }
                                }
                        }
                    //final monomer
                    if (mono_counter != 0)
                        {
                            mono_counter++;
                        }
                    else
                        {
                            mono_counter=1;
                        }
                    for(int types =0; types < p->n_types; types++)
                        {
                            diff=diff+(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] * p->field_scaling_type[types])*(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] * p->field_scaling_type[types]);
                            if(types == type)
                                {
                                    diff_flip=diff_flip+(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[type*p->n_cells_local + mono_cells[N-1]] - mono_counter )* p->field_scaling_type[types]) * (p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[type*p->n_cells_local + mono_cells[N-1]] - mono_counter )* p->field_scaling_type[types]);
                                }
                            else
                                {
                                    diff_flip=diff_flip+(p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] + mono_counter )* p->field_scaling_type[types]) * (p->umbrella_field[types*p->n_cells_local + mono_cells[N-1]]-( p->fields_unified[types*p->n_cells_local + mono_cells[N-1]] + mono_counter )* p->field_scaling_type[types]);
                                }
                        }


                    mono_counter=0;


                    //finally, calculate difference
                    if(diff_flip < diff)
                        {
                            poly_flip[poly]=(diff-diff_flip)/diff;
                        }
                    else
                        {
                            poly_flip[poly]=-1;
                        }

                }
                    
        
        }
    free(mono_cells);
    return;
}
    