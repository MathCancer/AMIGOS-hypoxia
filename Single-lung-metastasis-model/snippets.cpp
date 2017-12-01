//
//  snippets.cpp
//  
//
//  Created by John Metzcar on 10/5/17.
//

#include "snippets.hpp"

//    std::cout<<
//    breast_cancer.phenotype_ID << "\n" << breast_cancer.name << "\n" <<
//    breast_cancer.name << "\n" <<
//    breast_cancer.time_units << "\n" <<
//    breast_cancer.spatial_units << "\n" <<
//    breast_cancer.birth_rate << "\n" <<
//    breast_cancer.apoptotic_death_rate << "\n" <<
//    breast_cancer.apoptotic_clearance_rate << "\n" <<
//    breast_cancer.necrotic_death_rate << "\n" <<
//    breast_cancer.necrotic_clearance_rate << "\n" <<
//    breast_cancer.motility << "\n" <<
//    breast_cancer.spatial_mechanical_factor << "\n" <<
//    breast_cancer.spatial_proliferation_factor << "\n" <<
//    breast_cancer.hypoxic_o2_threshold << "\n" <<
//    breast_cancer.critical_o2_threshold << "\n" <<
//    breast_cancer.cell_volume << "\n" <<
//    breast_cancer.max_cells<<std::endl;
//
//    breast_cancer.name = "breast_cancer";
//
//    std::cout<<breast_cancer.name << std::endl;
//
//    breast_cancer.max_cells =  test_voxel.volume / breast_cancer.cell_volume;
//
//    std::cout<<breast_cancer.max_cells << std::endl;
//
//    //breast_cancer.birth_rate = 1/18;
//
//    //Testing voxel_populaton_vector members/things
//
//    Voxel_Population_Vector CB1;
//
//    CB1.add_population("breast_cancer_default", breast_cancer);
//
//    CB1.live_cell_counts[0] = 6;

//    for (int i=0; i<8; i++)  Links work for multiple links among voxels
//    {
//        //int j = i+1;
//        //std::cout<<"Voxel1 is "<<tissue.tissue_mesh.voxels[i].mesh_index<<std::endl;
//        //std::cout<<"Voxel2 is "<<tissue.tissue_mesh.voxels[i+1].mesh_index<<std::endl;
//        tissue.tissue_mesh.link_voxels(tissue.tissue_mesh.voxels[i].mesh_index,tissue.tissue_mesh.voxels[i+2].mesh_index,10);
//        //std::cout<<"Voxel1 is "<<tissue.tissue_mesh.voxel_links[i].pVoxel1->mesh_index<<std::endl;
//        //std::cout<<"Voxel2 is "<<tissue.tissue_mesh.voxel_links[i].pVoxel2->mesh_index<<std::endl;
//
//        //std::cout<<tissue.tissue_mesh.voxel_links.size()<<std::endl;
//    }

//    std::ofstream fileoutput1 ( "voxel_link_indices1.txt", std::ofstream::out);
//
//    for( int i=0; i<pseudo_tissue.mesh.voxel_links.size(); i++ )
//    {
//
//        fileoutput1 << i << "\t" << pseudo_tissue.mesh.voxel_links[i].ID << "\n" << std::endl;
//
//    }
//
//    fileoutput1.close();
//
//    double dx = 60;
//    pseudo_tissue.resize_space( 0 , 1000.0 , 0, 1000.0 , 0.0 , 1000.0 , dx, dx, dx );
//
//    pseudo_tissue.display_information( std::cout );
//
//    //int j = pseudo_tissue.mesh.voxels.size();
//
//    std::ofstream fileoutput2 ("voxel_indices.txt", std::ofstream::out);
//
//    for( int i=0; i<pseudo_tissue.mesh.voxels.size(); i++ )
//     {
//         fileoutput2 << i << "\t" << pseudo_tissue.mesh.voxels[i].mesh_index << "\n" << std::endl;
//
//     }
//
//    fileoutput2.close();
//
//    std::ofstream fileoutput3 ("voxel_link_indices2.txt", std::ofstream::out);
//
//    for( int i=0; i<pseudo_tissue.mesh.voxel_links.size(); i++ )
//    {
//
//        fileoutput3 << i << "\t" << pseudo_tissue.mesh.voxel_links[i].ID << "\t" << pseudo_tissue.mesh.voxel_links[i].surface_area << pseudo_tissue.mesh.voxel_links[i].pVoxel1->mesh_index << "\t" << pseudo_tissue.mesh.voxel_links[i].pVoxel2->mesh_index << std::endl;
//
//    }
//
//    fileoutput3.close();



//             fileoutput4.open( filename );

//             for # of voxels 0 to i -1 - loop over all the voxels, which will eventually be the voxel_population_vectors
//                 for # of phenotypes per voxel - loop over all the phenotypes

// might need turned into a more standard/standardizable function at some point, but gets the job done now ...
//             for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
//              {
//                  fileoutput4 << t << "\t" << tissue.voxel_population_vectors[i].live_cell_counts << "\t" << tissue.voxel_population_vectors[i].apoptotic_cell_counts << "\t" << tissue.voxel_population_vectors[i].necrotic_cell_counts << std::endl;
//              }
//             fileoutput4.close();
//std::cout << "Done!\n";

//delete [] filename;  Why did you do this?

//    std::ofstream fileoutput("cell_shells_output.txt", std::ofstream::out);  // needs header info ... somewhere ...
//    while (t < t_max)
//     {
//
//         // if it's time, save the simulation
//         if( fabs( t - next_output_time ) < dt/2.0 )
//         {
//             std::ofstream fileoutput4;
//             char* filename;
//             filename = new char [1024];
//
//             sprintf( filename, "output_%06.0f.txt" , next_output_time );
//
//             fileoutput4.open( filename );
//
////             for # of voxels 0 to i -1 - loop over all the voxels, which will eventually be the voxel_population_vectors
////                 for # of phenotypes per voxel - loop over all the phenotypes
//
//             // might need turned into a more standard/standardizable function at some point, but gets the job done now ...
//             for(int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
//              {
//                  fileoutput4 << t << "\t" << tissue.voxel_population_vectors[i].live_cell_counts << "\t" << tissue.voxel_population_vectors[i].apoptotic_cell_counts << "\t" << tissue.voxel_population_vectors[i].necrotic_cell_counts << std::endl;
//              }
//             fileoutput4.close();
//             //std::cout << "Done!\n";
//
//             //delete [] filename;  Why did you do this?
//             next_output_time += output_interval;
//         }
//
//      for (int i=0; i<tissue.tissue_mesh.voxels.size(); i++)
//       {
//           tissue.voxel_population_vectors[i].update_populations( dt );
//       }
//
//         t += dt;
//     }


/*    // create a microenvironment, and set units
 
 Microenvironment M;
 M.name = "microenvironment";
 M.time_units = "min";
 M.spatial_units = "micron";
 M.mesh.units = M.spatial_units;
 
 // set up and add all the densities you plan
 
 M.set_density( 0 , "oxygen" , "mmHg" );
 
 // here's how you add a new substrate
 M.add_density( "crayons" , "Megacrayola" );
 
 
 // set the properties of the diffusing substrates
 
 M.diffusion_coefficients[0] = 1e5;
 M.decay_rates[0] = 10; // 100 micron length scale
 M.diffusion_coefficients[1] = 1e4;
 M.decay_rates[1] = 1.0/9.0; // 300 micron length scale
 
 // set the mesh size
 
 double dx = 20; //
 M.resize_space( 0 , 1000.0 , 0, 1000.0 , 0.0 , 1000.0 , dx, dx, dx );
 
 // display summary information
 
 M.display_information( std::cout );
 
 // set up metadata
 
 BioFVM_metadata.program.user.surname = "Kirk";
 BioFVM_metadata.program.user.given_names = "James T.";
 BioFVM_metadata.program.user.email = "Jimmy.Kirk@starfleet.mil";
 BioFVM_metadata.program.user.organization = "Starfleet";
 BioFVM_metadata.program.user.department = "U.S.S. Enterprise (NCC 1701)";
 
 BioFVM_metadata.program.creator.surname = "Roykirk";
 BioFVM_metadata.program.creator.given_names = "Jackson";
 BioFVM_metadata.program.creator.organization = "Yoyodyne Corporation";
 
 BioFVM_metadata.program.program_name = "Nomad";
 BioFVM_metadata.program.program_version = "MK-15c";
 BioFVM_metadata.program.program_URL = "";
 
 // set initial conditions
 
 // use this syntax to create a zero vector of length 3
 // std::vector<double> zero(3,0.0);
 
 // use this syntax for a parallelized loop over all the
 // voxels in your mesh:
 #pragma omp parallel for
 for( int i=0 ; i < M.number_of_voxels() ; i++ )
 {
 // use this syntax to access the coordinates (as a vector) of
 // the ith voxel;
 // M.mesh.voxels[i].center
 
 // use this access the jth substrate at the ith voxel
 // M.density_vector(i)[j]
 
 }
 
 // save the initial profile
 
 // M.write_to_matlab( "initial.mat" ); // barebones
 save_BioFVM_to_MultiCellDS_xml_pugi( "initial" , M , 0.0 ); // MultiCellDS digital snapshot
 
 // set up the diffusion solver, sources and sinks
 
 M.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D;
 
 M.bulk_supply_rate_function = supply_function;
 M.bulk_supply_target_densities_function = supply_target_function;
 M.bulk_uptake_rate_function = uptake_function;
 
 double t     = 0.0;
 double t_max = 100.0;
 double dt    = 0.1;
 
 double output_interval  = 10.0;  // how often you save data
 double next_output_time = t;     // next time you save data
 
 while( t < t_max )
 {
 // if it's time, save the simulation
 if( fabs( t - next_output_time ) < dt/2.0 )
 {
 std::cout << "simulation time: " << t << " " << M.time_units << " (" << t_max << " " << M.time_units << " max)" << std::endl;
 
 char* filename;
 filename = new char [1024];
 
 // sprintf( filename, "output_%6f.mat" , next_output_time );
 // M.write_to_matlab( filename );
 
 sprintf( filename, "output_%6f" , next_output_time );
 save_BioFVM_to_MultiCellDS_xml_pugi( filename , M , 0.0 ); // MultiCellDS digital snapshot
 
 delete [] filename;
 next_output_time += output_interval;
 }
 
 M.simulate_bulk_sources_and_sinks( dt );
 M.simulate_diffusion_decay( dt );
 M.simulate_cell_sources_and_sinks( dt );
 
 t += dt;
 }
 
 // M.write_to_matlab( "final.mat"); // barebones
 save_BioFVM_to_MultiCellDS_xml_pugi( "final" , M , 0.0 ); // MultiCellDS digital snapshot  */

//    tissue.voxel_population_vectors[251].live_cell_counts[0] = 20.0;
//    tissue.voxel_population_vectors[252].live_cell_counts[0] = 20.0;
//    tissue.voxel_population_vectors[253].live_cell_counts[0] = 20.0;
//    tissue.voxel_population_vectors[254].live_cell_counts[0] = 20.0;
//    tissue.voxel_population_vectors[255].live_cell_counts[0] = 20.0;
//    tissue.voxel_population_vectors[256].live_cell_counts[0] = 20.0;
