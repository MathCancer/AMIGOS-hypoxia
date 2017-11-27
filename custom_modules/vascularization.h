#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#ifndef __PhysiCell_vascularization__
#define __PhysiCell_vascularization__

class Vascularization_Options
{
 private:
 public:
	int vascular_mesh_multiplier; // a multiple of the BioFVM mesh size 
	double vascular_extension_speed; 
	double vascular_birth_rate; 
	
	std::vector<double> vascular_substrate_densities; 
	
	Vascularization_Options(); 
	
	void sync_to_BioFVM( void ); 
};

extern Vascularization_Options vascular_options; 

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
	
	void sync_to_BioFVM( void ); 
}; 

extern Coarse_Vasculature coarse_vasculature; 

void setup( void ); 



#endif
