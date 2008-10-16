// $Id$

#include "PhononDispersion.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "mesh.h"
#include "elem.h"
#include "PhononModel.h"
 


#include "Material.h"
#include "Boundary.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"



//To Cut
//extern "C" void zheev_(char *jobz, char *uplo, int *n, complex<double> *a,int *lda,
// double *w, complex<double> *work, int *lwork, double *rwork, int *info);
//
extern "C" void dsyev_(char* jobz,  char* uplo,int* n, double* a,int* lda, double* w,double* work,int* lwork, int* info );


using namespace std;
Device* PhononDispersion::_device;

//-----------------------------------------------------------------//


void PhononDispersion::parse_options( )
{ 
 
  //const ModelOptions& sim_opt = SimulationInterface::get_options();

  const ModelOptions& sim_opt = get_options();


}

void  PhononDispersion::do_init( ) 
{

   SimulationEnvironment& si = get_environment();   

   _device = &( si.get_device() );

  mesh = & (_device->get_mesh());
 

}


//--------------------------------------------------------------------------------//
 PhononDispersion::~PhononDispersion()
{


}
//---------------------------------------------------------------------------------//
 PhononDispersion::PhononDispersion()
{
  

}
//----------------------------------------------------------------------------------//
PhysicalModel*
PhononDispersion::create_physical_model(const ModelOptions &options,
    const Material* mat) const throw (ModelErrorException)
{
  
    
    PhononModel* model = dynamic_cast<PhononModel*> ( PhysicalModelInterface::create("phonon",options) );
    
   if (model == NULL) 
     throw ModelErrorException("PhononModel: PhononDispersion physical model is not created" );
    return model;      
   

}
//----------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------//
PhononDispersion*   PhononDispersion::create (void)
{
  return new  PhononDispersion;
}



BoundaryProperties* PhononDispersion::create_boundary_model (const ModelOptions &options) const 
                    throw (ModelErrorException)

{
 

}









void
PhononDispersion::build_elemental_results(const std::set<std::string>& variables,
					  std::vector<double>& results, std::vector<std::string>& legend)
{
  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;

  const set<string>::const_iterator varend(variables.end());

  
  vector<ID> ids;
  unsigned int nm; 
   
  unsigned int n_vars = 0;  
  
  std::vector<unsigned int> W;

  const unsigned int nn  = mesh->n_active_elem();
  const unsigned int dim = mesh->mesh_dimension();
  legend.resize(variables.size());
  const Device& device = *(_device);

  
  int PD = -1;
  if (variables.count("RamanShift") ||
      variables.count("PhononVariables")  )
  {
    PD = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_1";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_2";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="dw_3";
    n_vars++;
  }

  
  //To cut
  int PFS = -1;
  if (variables.count("PhononEnergy") ||
      variables.count("PhononVariables") )
  {
    PFS = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_1";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_2";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="w_3";
    n_vars++;
  }
  ///

  int EV = -1;
  int EV1 = -1;
  int EV2 = -1;
  int EV3 = -1;
  if (variables.count("PhononPolarization")  ||
      variables.count("PhononVariables")  )
  {
    EV = 1;
    legend.resize(legend.size() + dim);
    EV1 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV1 + 2] = "E1_z";
        n_vars++;
      case 2:
        legend[EV1 + 1] = "E1_y";
        n_vars++;
        legend[EV1 + dim] = "modE1";
        n_vars++;
      default:
        legend[EV1] = "E1_x";
        n_vars++;
    }

    legend.resize(legend.size() + dim);
    EV2 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV2 + 2] = "E2_z";
        n_vars++;
      case 2:
        legend[EV2 + 1] = "E2_y";
        n_vars++;
        legend[EV2 + dim] = "modE2";
        n_vars++;
      default:
        legend[EV2] = "E2_x";
        n_vars++;
    }

    legend.resize(legend.size() + dim+1);
    EV3 = n_vars;
    switch (dim)
    {
      case 3:
        legend[EV3 + 2] = "E3_z";
        n_vars++;
      case 2:
        legend[EV3 + 1] = "E3_y";
        n_vars++;
        legend[EV3 + dim] = "modE3";
        n_vars++;
      default:
        legend[EV3] = "E3_x";
        n_vars++;
    }

  }

 //To cut
  int OA = -1;
  if (variables.count("OverAll") ||
      variables.count("PhononVariables") )
  {
    OA = n_vars;
    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelOverAll";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossOverAll";
    n_vars++;
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolOverAll";
    n_vars++;
    
  }
  ///

  //To cut
  int RI = -1;
  if (variables.count("RamanIntensity") ||
      variables.count("PhononVariables")  )
  {
    RI = n_vars;

    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="ParallelRamanIntensity_3";
    n_vars++;
  
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="CrossRamanIntensity_3";
    n_vars++;
  
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_1";
    n_vars++;

    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_2";
    n_vars++;
    
    legend.resize(legend.size() + 1);
    legend[n_vars]="NoPolRamanIntensity_3";
    n_vars++;

  }
  ///


    legend.resize(n_vars);
    
    results.resize(nn * n_vars,0.0);
 
  
  MeshBase::const_element_iterator it =    mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     mesh->active_local_elements_end();
  
  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
  { 
   unsigned int id = n_vars * elem_number;
    const Elem* elem = *it;

    if (PFS != -1)
    {
      std::vector<double> sol = PD_full_sol[elem];
      results[id + PFS]   = sol[0];
      results[id + PFS+1] = sol[1]; 
      results[id + PFS+2] = sol[2]; 
    }

    if (PD != -1)
    {
      std::vector<double> sol = PD_sol[elem];
      results[id + PD]   = sol[0];
      results[id + PD+1] = sol[1]; 
      results[id + PD+2] = sol[2]; 
    }

    if (EV != -1)
    {
      std::vector< std::vector <double> >  sol = sol_eignvectors[elem];

      //std::cout<< sol[0][0]<<std::endl;
 
      RealGradient eign(0);




      eign(0) = sol[0][0];
      eign(1) = sol[1][0];
      eign(2) = sol[2][0];

      switch (dim)
      {
        case 3:
          results[id + EV1 + 2] = eign(2);
        case 2:
          results[id + EV1 + 1] = eign(0);
          results[id + EV1 + dim] = eign.size();
        default:
          results[id + EV1] = eign(1); 
      
      }

 
      eign(0) = sol[0][1];
      eign(1) = sol[1][1];
      eign(2) = sol[2][1];

      switch (dim)
      {
        case 3:
          results[id + EV2 + 2] = eign(2);
        case 2:
          results[id + EV2 + 1] = eign(0);
          results[id + EV2 + dim] = eign.size();
        default:
          results[id + EV2] = eign(1); 
      
      }

 
      eign(0) = sol[0][2];
      eign(1) = sol[1][2];
      eign(2) = sol[2][2];

      switch (dim)
      {
        case 3:
          results[id + EV3 + 2] = eign(2);
        case 2:
          results[id + EV3 + 1] = eign(0);
          results[id + EV3 + dim] = eign.size();
        default:
          results[id + EV3] = eign(1); 
      
      }




    }

    if (OA != -1)
    {
      std::vector<double> sol = OverAll[elem];
      results[id + OA]    = sol[0];
      results[id + OA+1]  = sol[1]; 
      results[id + OA+2]  = sol[2];
          }
  

    if (RI != -1)
    {
      std::vector<std::vector <double> > sol = Intensity[elem];
      results[id + RI]   = sol[0][0];
      results[id + RI+1] = sol[0][1];
      results[id + RI+2] = sol[0][2];

      results[id + RI +3] = sol[1][0];
      results[id + RI +4] = sol[1][1];
      results[id + RI +5] = sol[1][2];
 
      results[id + RI +6] = sol[2][0];
      results[id + RI +7] = sol[2][1];
      results[id + RI +8] = sol[2][2];
     
 
    }
  

    elem_number++;
  } //over element

  results.resize(elem_number * n_vars);
}


//-------------------------------------------------------------------------------//
void  PhononDispersion::do_solve()
{
   
    solve_phonon_dispersion();		 

}

// void 
// PhononDispersion::diagonalize_complex(void)
//  {


//    std::cout<<"Diagonalize...start"<<std::endl;


//    int  n = 3;
   
//    char jobz = 'V';
//    char uplo = 'U';
//    int info;
//    int lwork = 3*n-1;
   
//    double *b;
//    complex<double> *work = new complex<double>[lwork]; //The work array to be used by zheev an its size
//    double *rwork = new double[3*n-2];
//    complex<double> *a = new complex<double>[n*n];

//    for (int j=0; j < n; j++){
//      for (int i=0; i < n; i++){
//         double value;
//        if (i < j) 
// 	 value = dynamical_matrix(j+1, i+1);
//        else
// 	 value = dynamical_matrix(i+1, j+1);

   
//        complex<double> s(value,0.0);
//        a[n*j+i] =s;//passed as rods
//      }
//    }
    
//    zheev_(&jobz, &uplo, &n, a, &n, b, work, &lwork, rwork, &info);//perform the diagonalization
     
//    delete [] work;
//    delete [] rwork;
   
//    //sorting
//    {  
//      complex<double> ctemp;
     
//      for(int j = 0; j < n; j++){
//        int it = j;
//        double temp = b[j];
      
//        for(int i = j; i < n; i++){
// 	 if(b[i] > temp){
// 	   temp = b[i];
// 	   it = i;
// 	 }
//        }

//        b[it] = b[j];
//        b[j] = temp;
      
//        for(int k = 0; k < n; k++){
// 	 ctemp = a[n*j+k];
// 	 a[n*j+k] = a[n*it+k];
// 	 a[n*it+k] = ctemp;
//        }
//      }
     
//    }
//       //empty the autovector matrix
//    _eignvectors.resize(n);
//    E.resize(n);
//    for (unsigned int k = 0; k<n; k++)
//       _eignvectors[k].resize(n,0.0);
   
//    //Put the value in the new style matrix
//    for(int i=0; i < n; i++)
//    {
//      for(int j = 0; j < n; j++)
//      {
//         _eignvectors[i][j] = a[n*i+j].real();//passed as lines
//      }
//      E[i] = b[i];
//    }
   
//    delete [] a;
//    delete [] b;

// }          
//

void 
PhononDispersion::diagonalize_double(void)
 {


   // std::cout<<"Diagonalize...start"<<std::endl;

   //int n = 2;   
    //  D.resize(n);
//   for (unsigned int i = 0; i<n; i++)
//   {
//    D[i].resize(n,0.0);
//   }
 
  
//   D[0][0] = 1.0;
//    D[0][1] = 0.0;
//    D[1][0] = 0.0;
//    D[1][1] = 2.0; 

   //int  n = D[0].size();


   int n = 3;
   
   char jobz = 'V';
   char uplo = 'U';
   int info;
   int lwork = 3*n-1;

   double *b = new double[n];;
   double *work = new double[lwork]; //The work array to be used by zheev an its size
   double *a = new double[n*n];

   // std::cout<<"DM"<<std::endl;
   for (int j=0; j < n; j++){
    for (int i=0; i < n; i++){
       
     double value;
      if (i < j) 
        value = dynamical_matrix(j+1, i+1);
      else
        value = dynamical_matrix(i+1, j+1);

       a[n*j+i] = (double) value;
       //std::cout<< a[n*j+i] <<std::endl;
       }
    }

    dsyev_(&jobz, &uplo, &n, a, &n, b,work,&lwork,&info);

    //  std::cout<<"eignvalue:"<<std::endl;
    //std::cout<<std::sqrt(b[0])*8065.6<<std::endl;
    //std::cout<<std::sqrt(b[1])*8065.6<<std::endl;
    //std::cout<<std::sqrt(b[2])*8065.6<<std::endl;

   delete [] work;
    
   _eignvectors.resize(n);
   E.resize(n);
   for (unsigned int k = 0; k<n; k++)
      { _eignvectors[k].resize(n,0.0);}
   
   //Put the value in the new style matrix
   for(int i=0; i < n; i++)
    {
    for(int j = 0; j < n; j++)
     {
          _eignvectors[i][j] = a[n*i+j];//passed as lines
    }
     E[i] = std::sqrt(b[i]) *8065.6;
   }
   
   delete [] a;
   delete [] b;


}          
//

void PhononDispersion::solve_phonon_dispersion(void)
{

     std::vector<Tensor2Sym> raman_tensor;
     std::vector<Tensor1> light_polarization;
  

     mesh = & (_device->get_mesh());

     MeshBase::const_element_iterator       el     = mesh->active_elements_begin();
     const MeshBase::const_element_iterator end_el = mesh->active_elements_end();

      for ( ; el != end_el ; ++el)
     {   //loop over elements
       const Elem* elem = *el;

       ID subdomain = elem->subdomain_id();

       const Material* mat = _device->get_material(subdomain);

       PhononModel* phonon_model =  (  dynamic_cast<PhononModel*> (  mat -> get_model(get_id()) )  );

       phonon_model->set_element(elem);     

       phonon_model->re_init(); 

       //Full
       phonon_model->get_full_dynamical_matrix(dynamical_matrix);
       diagonalize_tensor();
       std::vector<double> full_sol  =  E;
       sol_eignvectors[elem] = _eignvectors;
          
        

       //Free
       phonon_model->get_free_dynamical_matrix(dynamical_matrix);
       diagonalize_tensor();
       std::vector<double> free_sol  =  E;
       
       std::vector<double> diff(3);
       for (unsigned int k = 0; k<3; k++)       
       {diff[k] = full_sol[k]-free_sol[k];}

       PD_full_sol[elem] = full_sol;
       PD_sol[elem] = diff; 

       phonon_model->get_raman_tensor(raman_tensor);

     phonon_model->get_light_polarization(light_polarization);
     Tensor1 e0       = light_polarization[0];
     Tensor1 es_paral = light_polarization[1];
     Tensor1 es_cross = light_polarization[2];
     Tensor1 es_nopol = light_polarization[3];

       std::vector<double> I_paral(3);
       std::vector<double> I_cross(3);
       std::vector<double> I_nopol(3);

       for(unsigned int nm =0; nm<3;nm ++)
       {
         
          Tensor2Gen D =  raman_tensor[0]*_eignvectors[0][nm] +
                          raman_tensor[1]*_eignvectors[1][nm] +
                          raman_tensor[2]*_eignvectors[2][nm];
       
         

          Tensor1 tens = D*e0;


          //Parallel
          double int_temp = es_paral*tens;
          I_paral[nm] = std::abs(int_temp) * std::abs(int_temp);

          //Cross
          int_temp = es_cross*tens;
          I_cross[nm] = std::abs(int_temp) * std::abs(int_temp);

          //NoPol
          int_temp = es_nopol*tens;
          I_nopol[nm] = std::abs(int_temp) * std::abs(int_temp);

       }
  
       //Paral
       double temp_paral = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_paral += diff[k] * I_paral[k];}
        
        double sum_paral =  I_paral[0] +  I_paral[1] + I_paral[2];
        if (sum_paral < 1e-13)
       {
            temp_paral = 0.0;
         }
       else
       {
        temp_paral /=sum_paral; 
        }


       //Cross
       double temp_cross = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_cross += diff[k] * I_cross[k];}

            
       double sum_cross =  I_cross[0] +  I_cross[1] + I_cross[2];
       if (sum_cross < 1e-13)
       {
            temp_cross = 0.0;
         }
       else
       {
        temp_cross /=sum_cross; 
        }

       //NoPol
       double temp_nopol = 0;
       for (unsigned int k = 0; k<3; k++)
         {temp_nopol += diff[k] * I_nopol[k];}

            
       double sum_nopol =  I_nopol[0] +  I_nopol[1] + I_nopol[2];
       if (sum_nopol < 1e-13)
       {
            temp_nopol = 0.0;
         }
       else
       {
        temp_nopol /=sum_nopol; 
       }


        //Write Output
        OverAll[elem].resize(3,0.0);
        OverAll[elem][0] = temp_paral;
        OverAll[elem][1] = temp_cross;
        OverAll[elem][2] = temp_nopol;

        //Intensity
        Intensity[elem].resize(3);
        Intensity[elem][0].resize(3,0.0);
        Intensity[elem][0] = I_paral;
        Intensity[elem][1].resize(3,0.0);
        Intensity[elem][1] = I_cross;
        Intensity[elem][2].resize(3,0.0);
        Intensity[elem][2] = I_nopol;

        
     }
     
}

void 
PhononDispersion::diagonalize_tensor(void)
 {

   int n = 3;
   
   Tensor2Gen  V;
   double landa1;
   double landa2;
   double landa3;

   dynamical_matrix.eigen(&landa1,&landa2,&landa3,&V);

      _eignvectors.resize(n);
      E.resize(n);
      for (unsigned int k = 0; k<n; k++)
      { _eignvectors[k].resize(n,0.0);}

    

    for(unsigned int i=0; i < n; i++)
    {
    for(unsigned int j = 0; j < n; j++)
     {
          
       _eignvectors[i][j] = V(i+1,j+1);

    }
    
  }

      double conv = 8065.6;

      E[0]=sqrt(landa1)*conv;
      E[1]=sqrt(landa2)*conv;
      E[2]=sqrt(landa3)*conv;
   

}          
