#include "./vascularization.h" 

namespace PhysiCell{

Vascularization_Options default_vascular_options; 
Coarse_Vasculature coarse_vasculature; 
	
Vascularization_Options::Vascularization_Options()
{
	vascular_mesh_multiplier = 4; 
	vascular_extension_speed = 1.0; 
	vascular_birth_rate = 1.0/18.0; 
	
	std::cout << __FUNCTION__ << std::endl; 
	
	vascular_substrate_densities.resize( 1 , 1.0 ); 

	std::cout << __FUNCTION__ << std::endl; 

	return; 
}
	
	
void Vascularization_Options::sync_to_BioFVM( void )
{
	// make sure it has the right number of densities. 
	// set them equal to the boundary conditions for BioFVM 
	
	std::cout << __FUNCTION__ << std::endl; 

	vascular_substrate_densities = default_microenvironment_options.Dirichlet_condition_vector; 
	
	std::cout << __FUNCTION__ << std::endl; 

	return; 
}

Vascular_Densities::Vascular_Densities()
{
	functional = 0.0; 
	total = 0.0; 
	return; 
}

Coarse_Vasculature::Coarse_Vasculature()
{ return; }

void Coarse_Vasculature::sync_to_BioFVM( void )
{
	// first, resize the mesh 
	
	Microenvironment* pME = get_default_microenvironment(); 
	
	
	int Xnodes = pME->mesh.x_coordinates.size(); 
	int Ynodes = pME->mesh.y_coordinates.size();  
	int Znodes = pME->mesh.z_coordinates.size(); 
	
	if( Xnodes > 1 )
	{ Xnodes /= default_vascular_options.vascular_mesh_multiplier; }
	if( Ynodes > 1 )
	{ Ynodes /= default_vascular_options.vascular_mesh_multiplier; }
	if( Znodes > 1 )
	{ Znodes /= default_vascular_options.vascular_mesh_multiplier; }

	// next, make sure the microenvironment has oxygen 
	
	int oxygen_i = pME->find_density_index( "oxygen" ); 
	if( oxygen_i < 0 )
	{
		std::cout << "Adding oxygen to the microenvironment ... " << std::endl; 
		pME->add_density( "oxygen", "mmHg" , 1e5 , 0.1 ); 
		oxygen_i = pME->find_density_index( "oxygen" ); 
		
		default_microenvironment_options.Dirichlet_condition_vector[oxygen_i] = 38.0;  
		default_microenvironment_options.Dirichlet_activation_vector[oxygen_i] = true; 		
	}
	
	// next, make sure the microenvironment has VEGF 
	
	int VEGF_i = pME->find_density_index( "VEGF" ); 
	{
	if( VEGF_i < 0 )
		// 5.8 × 10−11 m2 s−1. // https://www.nature.com/articles/nprot.2012.051
		
		// decay : http://walter.deback.net/old/media/KohnLuque_etal_PhysBiol_2013.pdf 
		
		
		std::cout << "Adding VEGF to the microenvironment ... " << std::endl; 
		
		
		
		pME->add_density( "VEGF", "dimensionless" , 3.5e3 , 0.035 ); 
		oxygen_i = pME->find_density_index( "VEGF" ); 
		
		default_microenvironment_options.Dirichlet_condition_vector[oxygen_i] = 38.0;  
		default_microenvironment_options.Dirichlet_activation_vector[oxygen_i] = true; 		
	}

	// next, make sure the 


	
	std::cout << __FUNCTION__ << std::endl; 
	
	pME->display_information( std::cout ); 

	mesh.resize( default_microenvironment_options.X_range[0] , default_microenvironment_options.X_range[1] ,
		default_microenvironment_options.Y_range[0] , default_microenvironment_options.Y_range[1] ,
		default_microenvironment_options.Z_range[0] , default_microenvironment_options.Z_range[1] ,
		Xnodes, Ynodes, Znodes ); 
		
	mesh.units = default_microenvironment_options.spatial_units; 
		
	std::cout << __FUNCTION__ << std::endl; 
		
	get_default_microenvironment()->mesh.display_information( std::cout ); 
	mesh.display_information( std::cout ); 

	std::cout << __FUNCTION__ << std::endl; 
	
	return; 
}

/*
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
