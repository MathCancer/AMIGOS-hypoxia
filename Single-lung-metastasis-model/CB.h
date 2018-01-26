//
//  CB.h
//  
//
//  Created by John Metzcar on 8/25/17.
//
//

#ifndef __CB_h__
#define __CB_h__

#include "./BioFVM.h"

namespace BioFVM{


class Phenotype
 {
 private:
     // Screen output stuff
     
 public:
     
     int phenotype_ID;
     std::string name; // Do we need this here?  Does it really get stored as the actual name of the instantiation?  Like "Phenotype MCF-7;"???
     std::string time_units;
     std::string spatial_units;
     double birth_rate;
     double apoptotic_death_rate;
     double apoptotic_clearance_rate;
     double necrotic_death_rate;
     double necrotic_clearance_rate;
     double motility;
     double spatial_mechanical_factor;
     double spatial_proliferation_factor;
     double base_secretion_rate;
     double hypoxic_o2_threshold;
     double critical_o2_threshold;
     double cell_volume;
     double max_cells; // Voxel carrying capacity expressed as cells - a bit of a voluem measurement
     std::vector<std::string> other_properties_names;
     std::vector<double> other_properties;  // Will include extra normal phenotype stuff - like a vascular cut off
     
     // These vectors store the base rates.  Environmental modifications made on the fly, just like the birth and death rates.  There is one three member vector (live, apop, necr) per substrate
     std::vector <std::vector<double>> secretion_rate; // 1/time  Vectors allow for multiple substrate uptake and secrection
     std::vector <std::vector<double>> substrate_target; // 1/time (at times perhaps belongs in substrate/properties??
     std::vector <std::vector<double>> uptake_rate; // 1/time

     
    
     Phenotype();  // Done
    
 };

class Voxel_Population_Vector  // Population in a tissue voxel
 {
 private:
    // Screen output stuff
 public:
    int tissue_voxel_index;  // Index of tissue voxel associated with a particular voxel population vector
    std::vector<int> microenvironment_voxel_indices;
    Microenvironment* pMicroenvironment;
    std::vector<std::string> names_vector;  // Will be moved out eventually
    std::vector<Phenotype> phenotypes_vector;
    std::vector<double> live_cell_counts;  // Vector of live cell counts for each population (on per voxel level)
    std::vector<double> apoptotic_cell_counts;  // Vector of apoptotic cell counts for each population (on per voxel basis) - Could see this going in as just a population instead of subpopulation as needed
    std::vector<double> necrotic_cell_counts;  // Vector of necrotic_cell_counts for each population (on per voxel basis) - Could see this going in as just a population instead of subpopulation as needed
    
    // Changes the live cell count and name of the population given by the population_number to the given population size and name
    void set_population(int population_number, double live_cell_count, std::string name);  // Done
    
    // Changes the live cell count etc and name of the population given by the population_number to the passed counts and name
    void set_population(int population_number, double live_cell_count, double apoptotic_cell_count, double necrotic_cell_count, std::string name); // Done
    
    // Adds a population of the given population name, initilizizes all the populations of the added phenotype to 0 counts, and adds the passed phenotype to that tissue voxel's phenotype vector
    void add_population(std::string name, Phenotype phenotype);  // done
     
    void pair_to_microenvironment( Microenvironment* pME );  // Done
     
    void setup_microenvironment_voxel_indices( void );
     
    std::vector<double> sample_microenvironment( void ); // Done
     double update_AF_secretion_rate ( double oxygen );  // Done This is hardcoded.  Is there more general way to do this?
     double update_vascular_creep_rate ( double AF );
     void update_vascular_population ( double dt, double AF );
     
    // Updates cell populations within each voxel, based on cell phenotype
     
    void update_populations( double dt, double oxygen ); 
    
    // Updates cell phenotypes within each voxel, based on uE - is this really the same as update populations.  Yes ...
     
    void update_phenotypes (Microenvironment& chemical_microenvironment, double dt);  // Not using currently
    
    Voxel_Population_Vector();  //done
    
    /* Possible additions
    void cell_density_vector( std::vector<double> max_cell_counts, std::vector<double> total_cell_counts);
    */
 };

class Substrate_Vector
{
    public:
    
    int tissue_voxel_index; // Index of tissue voxel associated with a particular voxel population vector
    std::vector<std::string> names_vector;  // Will be moved out eventually
    std::vector<std::string> units_vector;  // Will be moved out eventually
    std::vector<double> substrate_quantity;  // in units of "stuff"
    std::vector<double> substrate_decay_constant;  // 1/time
    
    // Changes the live cell count and name of the population given by the population_number to the given population size and name
    void set_substrate(int substrate_number, double quantity, std::string name, std::string unit);  // Done
    
    // Changes the live cell count etc and name of the population given by the population_number to the passed counts and name
    void set_population(int population_number, double live_cell_count, double apoptotic_cell_count, double necrotic_cell_count, std::string name); // Done
    
    // Adds a population of the given population name, initilizizes all the populations of the added phenotype to 0 counts, and adds the passed phenotype to that tissue voxel's phenotype vector
    void add_substrate(std::string name, std::string unit);  // done
    void set_decay_constant (int substrate_number, double decay_constant);
    
    
    Substrate_Vector();
    
};
    
class Substrate_Properties
{
    public:
//    int voxel_link_index; // Index of voxel link associated with a particular
    std::vector<double> diffusion_coefficients;
    void add_diffusion_coefficient(double diffusion_coefficient);
    
    Substrate_Properties();
};
    
class Vasculature  // Make it its on thing!  See the notes/photo.  This iwll be a different thing - "bulk" source - after all it is blood.
 {
     double vascular_density;
     std::string time_units;
     std::string spatial_units;
     double vascular_birth_rate;
     double vascular_death_rate;
     double apoptotic_clearance_rate;
     double necrotic_death_rate;
     double necrotic_clearance_rate;
     double extension_rate;
     double spatial_mechanical_factor;
     double spatial_proliferation_factor;
     double hypoxic_o2_threshold;
     double critical_o2_threshold;
     double vascular_cutoff_threshold;
     double cell_volume;
     double max_cells;
     
     
    /* Will be moved into here from the voxel population vectors
     */
 };
    
class Tissue
 {
 public:
     
     // phenotype names and units eventually goes in tissue?
     
     Microenvironment chemical_microenvironment;  // This will hold the diffusive microenvironment fields  (unless it doesn't ... see below)
     General_Mesh tissue_mesh;  // In BIOFVM, BUT W/O resize functions.  // This holds the cell bucket voxels - it is the mesh for the semi-discrete cell movement method AND for the vascular FVM.  
     //Spherical_Mesh spherical_mesh; // NOT YET IN bioTissueBox
     
     std::vector< Voxel_Population_Vector > voxel_population_vectors;
     std::vector< Substrate_Vector > substrate_vectors;
     std::vector< Substrate_Properties> substrate_properties_vectors;
     std::vector<double> vascular_densities;
     std::vector<double> max_vascular_densities;
//     std::vector<double> center_to_center_voxel_distances; // Holds center to center distances between all linked voxels.
     
     void make_tumor_spheroid_2D( double tumor_radius , double voxel_radius );  // Hexagonal lattice.  Close to working for 3D spheroid but needs work.
     void link_tumor_spheroid( double test_distance, double shared_surface_area ); // Linking the hexoganol lattice.  doe I neeesd to pass it a reference to the Voxels vector?
     void make_voxel_links( void ); // Done
     void spherical_geometry_linker (double dr, double tumor_radius ); // Makes and links a domain of concentric spheres.
     void create_population_vectors (void);  // Done
     void create_substrate_vectors (void);  // Done
     void create_substrate_properties_vectors (void);  // Done
     void resize_space ( double dx, double dy, double dz);
     void update_substrates( double dt );  // Done
     void run_diffusion ( double dt, int substrate_index); // Done
     void flux_vascular_density ( double dt );  // Done
     void flux_cell_populations ( double dt );  // Done.  What about when have a bunch of tissues?  Then what?  Will I need to pass it a specific voxel vector?  Probably ...
     void flux_cell_populations (std::vector< Voxel >* voxels, std::vector< Voxel_Link >* voxel_links);
     
     void write_voxel_populations_to_matlab( std::string filename );  // In development ...
     void write_substrates_to_matlab( std::string filename );  // done
     void write_substrate_properties_to_matlab( std::string filename );  // Deon
     void write_all_to_matlab( std::string filename, double time );  // Done
     void write_population_svg (int population_number, double time, std::string filename, double height, double width);  // done
     void write_vasculature_svg(std::string filename, double height, double width);
     void write_tissue_apoptotic_svg(std::string filename, double height, double width);
     void write_tissue_necrotic_svg(std::string filename, double height, double width);
     void write_subsrate_svg (int substrate_number, double time, std::string filename, double height, double width);
     void write_oxygen_svg(std::string filename, double height, double width);
     void write_AF_svg(std::string filename, double height, double width);
     
// Need this block to write text out to the file.  The whole system keeps giving me a type definition error, but I don't know why.  It seems perferctly well defined ...
     
//     struct Tissue_SVG_options_struct {
//         //bool plot_nuclei = true;
//
//         std::string simulation_time_units = "min";
//         std::string mu = "&#956;";
//         std::string simulation_space_units = "&#956;m";
//
//         std::string label_time_units = "days";
//
//         double font_size = 200;
//         std::string font_color = "black";
//         std::string font = "Arial";
//
//         double length_bar = 100;
//     };
//
//     extern Tissue_SVG_options_struct Tissue_SVG_options; // Why make this an external variable?
//
     std::vector<std::string> population_blue_to_yellow_coloring( int population_number, BioFVM::Voxel_Population_Vector* pVPV);  // Done
     std::vector<std::string> substrate_blue_to_yellow_coloring( int substrate_number, BioFVM::Substrate_Vector* pSV);  // Done
     std::string formatted_minutes_to_DDHHMM( double minutes );  // Done, took straight from PC
     
     void SVG_plot( std::string filename, Microenvironment& TM, double z_slice, double time, std::vector<std::string> (*voxel_coloring_function)(Voxel_Population_Vector*));
//     void read_from_matlab( std::string filename );
    
    Tissue();
 };
     
class Cartesian_Tissue : public Tissue
 {
     
     
     Cartesian_Mesh mesh;  // Could I put in another mesh constructor defined only in the tissue/CB scope?  So I could get the voxel_links?
     
     void create_voxel_population_vectors ( void );
     void resize_space ( double dx, double dy, double dz);
     
     Cartesian_Tissue();
     
 };
    
};

#endif 

/* CB_h */
