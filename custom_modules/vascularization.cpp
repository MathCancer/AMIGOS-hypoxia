#include "./vascularization.h" 

namespace PhysiCell{

Vascular_Options default_vascular_options; 
Coarse_Vasculature coarse_vasculature; 
	
Vascular_Options::Vascular_Options()
{
	vascular_mesh_multiplier = 4; 
	vascular_extension_speed = 1.0; 
	vascular_birth_rate = 1.0/18.0; 
	
	angiogenesis_dt = 60.0; 
	
	std::cout << __FUNCTION__ << std::endl; 
	
	vascular_substrate_densities.resize( 1 , 1.0 ); 

	std::cout << __FUNCTION__ << std::endl; 

	return; 
}
	
	
void Vascular_Options::sync_to_BioFVM( void )
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
	mesh.resize(1,1,1) ; 
	vascular_densities.resize(1); 
	vascular_substrate_densities = default_vascular_options.vascular_substrate_densities; 
	
	return; 
}

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
		
		// decay 72h half-life : http://walter.deback.net/old/media/KohnLuque_etal_PhysBiol_2013.pdf 
		// ECM binding: 1.5e-3 s-1 ~ 0.09 min-1 
		
		std::cout << "Adding VEGF to the microenvironment ... " << std::endl; 
		
		pME->add_density( "VEGF", "dimensionless" , 3.5e3 , 0.09 ); 
		VEGF_i = pME->find_density_index( "VEGF" ); 
		
		default_microenvironment_options.Dirichlet_condition_vector[VEGF_i] = 0.0;  
		default_microenvironment_options.Dirichlet_activation_vector[VEGF_i] = false; 		
	}
	
	// next, resize the vascular mesh 
	
	mesh.resize( default_microenvironment_options.X_range[0] , default_microenvironment_options.X_range[1] ,
		default_microenvironment_options.Y_range[0] , default_microenvironment_options.Y_range[1] ,
		default_microenvironment_options.Z_range[0] , default_microenvironment_options.Z_range[1] ,
		Xnodes, Ynodes, Znodes ); 
	mesh.units = default_microenvironment_options.spatial_units; 
		
	// set the substrate densities to the correct values 
	
	vascular_substrate_densities = default_vascular_options.vascular_substrate_densities;  

	// next, make sure the vascular densities are of the right size 

	vascular_densities.resize( mesh.voxels.size() ); 

	return; 
}

Vascular_Densities& Coarse_Vasculature::operator()( std::vector<double> position )
{ return vascular_densities[ mesh.nearest_voxel_index( position ) ]; }

Vascular_Densities& Coarse_Vasculature::operator()( int n )
{ return vascular_densities[n]; } 

Vascular_Densities& Coarse_Vasculature::operator()( int i, int j, int k )
{ return vascular_densities[ mesh.voxel_index(i,j,k) ]; }

Vascular_Densities& Coarse_Vasculature::operator()( Cell* pCell )
{ return this->operator()( pCell->position ); }

void coarse_vasculature_setup( void )
{
	// sync the options to BioFVM 
	default_vascular_options.sync_to_BioFVM();
	
	// USER EDITS TO default_vascular_options GO HERE!!!
	
	
	
	// END OF USER EDITS TO default_vascular_options
	
	// sync the environment to BioFVM
	
	coarse_vasculature.sync_to_BioFVM(); 
	
	// now, set bulk source and sink functions 
	
	std::cout << "need to set up source functions!!!" << std::endl; 
	
	return; 
}

void update_coarse_vasculature( double dt )
{
	static double t_angio = 0.0; 
	static double t_last_angio_update = 0.0; 
	static double t_next_angio_update = 0.0 + default_vascular_options.angiogenesis_dt; 
	
	// simulate bulk sources and sinks 
	// microenvironment.simulate_bulk_sources_and_sinks( dt ); 
	
	if( t_angio > t_next_angio_update )
	{
		std::cout << "angio update!" << std::endl; 
		double dt_temp = t_angio - t_last_angio_update; 
		
		
		// birth and death in the vasculature 
		
		// fluxes 		
		
		// custom function 
		
		
		t_last_angio_update = t_angio; 
		t_next_angio_update = t_angio + default_vascular_options.angiogenesis_dt; 
	}
	

	
	
	t_angio += dt; 
	
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
