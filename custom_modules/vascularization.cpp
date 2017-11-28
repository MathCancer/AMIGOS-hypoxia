#include "./vascularization.h" 

namespace PhysiCell{

Vascularization_Options vascular_options; 
	
Vascularization_Options::Vascularization_Options()
{
	vascular_mesh_multiplier = 4; 
	vascular_extension_speed = 1.0; 
	vascular_birth_rate = 1.0/18.0; 
	
	vascular_substrate_densities.resize( get_default_microenvironment()->number_of_densities() , 1.0 ); 
	return; 
}
	
	
void Vascularization_Options::sync_to_BioFVM( void )
{
	// make sure it has the right number of densities. 
	// set them equal to the boundary conditions for BioFVM 
	
	vascular_substrate_densities = default_microenvironment_options.Dirichlet_condition_vector; 
	
	return; 
}

Vascular_Densities::Vascular_Densities()
{
	functional = 0.0; 
	total = 0.0; 
	return; 
}

Coarse_Vasculature::Coarse_Vasculature()
{
	
	
	
	return; 
}

/*
extern 

class Vascular_Densities
{
 private:
 public: 
	double functional; 
	double total; 
	
	Vascular_Densities(); 
}; 

class Coarse_Vasculature
{
 private:
 public: 
	General_Mesh mesh; 
	std::vector<Vascular_Densities> vascular_densities; 
	
	Vascular_Densities& operator()( std::vector<double> position );
	Vascular_Densities& operator()( int i );
	Vascular_Densities& operator()( int i );
	
	Coarse_Vasculature(); 
	
	std::vector<double> vascular_substrate_densities; 
	
	// set to size of BioFVM, add VEGF to TME, check for O2 
	// 
	void sync_to_BioFVM( void ); 
}; 

extern Coarse_Vasculature coarse_vasculature; 

void add_VEGF_to_BioFVM( void ); 
void add_VEGF_to_cells( void ); 





*/

};
