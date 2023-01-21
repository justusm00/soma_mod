    // ARRAY AND VARIABLE DECLARATION
    
    const unsigned int N = p->reference_Nbeads; //monomers per polymer
    int64_t * mono_cells=(int64_t *)malloc( N* sizeof(int64_t)); //monomer cell indices. Values are -1 if no target density available in that cell
    int64_t * delta_fields_unified = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); //array that stores changes in density
    int64_t * delta_fields_unified_best = (int64_t *)malloc(p->n_types * p->n_cells_local * sizeof(int64_t)); 
    uint64_t num_poly_flippable = 0; // number of flippable polymers
    soma_scalar_t acc_rate= 1.0; //initial flip acceptance rate
    uint64_t num_iter=0; 
    uint64_t num_acc = 0;
    soma_scalar_t total_cost = 0.0;
    soma_scalar_t total_cost_old = 0.0;
    soma_scalar_t total_cost_best = 0.0;

    // ARRAYS THAT CONTAIN FLIP CANDIDATE INFO
    // BUFFER SIZE MAY NEED TO BE ADJUSTED
    uint64_t flip_buffer_size = (uint64_t)(p->n_polymers)/(uint64_t)(10);
    unsigned int * poly_types=(uint64_t *)malloc(flip_buffer_size * sizeof(uint64_t)); //array that stores polymer types
    unsigned int * poly_types_best=(uint64_t *)malloc(flip_buffer_size * sizeof(uint64_t));
    int64_t * poly_flippable_indices = (int64_t *)malloc( flip_buffer_size * sizeof(int64_t)); //array that contains indices of flippable polymers
    int64_t * poly_flippable_types = (int64_t *)malloc( flip_buffer_size * sizeof(int64_t)); //array that contains types of flippable polymers
    int64_t * poly_flippable_cell_indices = (int64_t *)malloc(flip_buffer_size * N * sizeof(int64_t)); //array that stores in which cells a given polymer has monomers
    int64_t * poly_flippable_cell_num = (int64_t *)malloc(flip_buffer_size * N * sizeof(int64_t)); //array that stores number of monomers in cells. Values correspond to cells specified in poly_flippable_cell_indices

    // SIMULATED ANNEALING PARAMETERS 

    soma_scalar_t Tmax = 0.001;
    soma_scalar_t Tmin = 0.0001;
    soma_scalar_t alpha = 0.85;
    uint64_t max_iter = 10000; //maximum number of flips 
    soma_scalar_t acc_rate_target=0.02; //flip acceptance rate after which converison will be stopped

    srand(time(0)); //rng for selecting polymers to flip



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
                    poly_flippable_types[num_poly_flippable]=p->polymers[poly].type;
                    num_poly_flippable ++;
                    if(num_poly_flippable > flip_buffer_size)
                        {
                            printf("ERROR: POLYMER FLIP BUFFER TOO SMALL");
                            return -1;
                        }
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
                                    poly_flippable_cell_indices[num_poly_flippable * N + k]=mono_cells[mono]; //unique monomer cell indices
                                    poly_flippable_cell_num[num_poly_flippable * N + k]=mono_counter; //corresponding number of monomers
                                    mono_counter=0;
                                    k++;
                                }
                        }
   
                    //final monomer
                    if (mono_counter == 0) poly_flippable_cell_num[num_poly_flippable * N + k]=1;
               
                    else 
                        {
                            mono_counter++;
                            poly_flippable_cell_num[num_poly_flippable * N + k]=mono_counter; 
                        }
                    poly_flippable_cell_indices[num_poly_flippable * N + k]=mono_cells[N-1];
                    k++;

                    //set end of arrays 
                    if(k<N-1) 
                        {
                            poly_flippable_cell_indices[num_poly_flippable * N + k]=-1;
                            poly_flippable_cell_num[num_poly_flippable * N + k]=-1;
                        } 
                }
        }



    //initialize cost
    total_cost=get_cost(p, delta_fields_unified);
    total_cost_old = total_cost;
    total_cost_best= total_cost;


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
            poly_types[poly]=poly_flippable_types[poly];
            poly_types_best[poly]=poly_flippable_types[poly];
        }

    // SIMULATED ANNEALING

    //printf("Total cost before : %f \n",total_cost);
    while((acc_rate > acc_rate_target) && (num_iter < max_iter))
        {
            soma_scalar_t T = Tmax;
            while(T > Tmin)
                {
                    num_iter++;
                    //choose random polymer to flip
                    uint64_t poly = rand() % (num_poly_flippable - 1);
                    Polymer *mypoly = p->polymers + 1; //only needed for rng
                    unsigned int initial_type = poly_types[poly];
                    unsigned int final_type = flip(initial_type);
                    total_cost=total_cost_old;
                    //calculate cost (only need to update it for the cells in which the polymer has monomers)
                    for(unsigned int i = 0; i < N; i++)
                        {
                            if(poly_flippable_cell_indices[poly * N + i] < 0) break;
                            unsigned int cell = poly_flippable_cell_indices[poly * N + i];
                            unsigned int num_mono = poly_flippable_cell_num[poly * N + i];
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
                                    if(poly_flippable_cell_indices[poly * N + i] < 0) break;
                                    unsigned int cell = poly_flippable_cell_indices[poly * N + i];
                                    unsigned int num_mono = poly_flippable_cell_num[poly * N + i];
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
                                            if(poly_flippable_cell_indices[poly * N + i] < 0) break;
                                            unsigned int cell = poly_flippable_cell_indices[poly * N + i];
                                            unsigned int num_mono = poly_flippable_cell_num[poly * N + i];
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
                                    if(poly_flippable_cell_indices[poly * N + i] < 0) break;
                                    unsigned int cell = poly_flippable_cell_indices[poly * N + i];
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
    free(poly_flippable_indices);
    free(poly_flippable_cell_indices);
    free(poly_flippable_cell_num);
    free(poly_types);
    free(poly_types_best);
    free(delta_fields_unified);
    free(delta_fields_unified_best);
    return 0;