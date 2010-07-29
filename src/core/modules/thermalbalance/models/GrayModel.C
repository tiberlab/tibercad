// $Id$

#include "GrayModel.h"
#include "SimulationOptions.h"


TIBER_MODULE(GrayModel,heat_transport, gray)


using namespace std;


GrayModel::GrayModel(const ModelOptions& options):HeatTransportModel(options)
{
 
 set_type(HeatTransportModel::Gray);

}

void
GrayModel::do_init(void)
{

  //const ModelOptions& options = get_options();
   //myopts.max_error    =  options.get_option("max_error",1e-3);
   //myopts.max_iter     =  options.get_option("max_iter",1);
  //myopts.phi_slices   = options.get_option("phi_slices",0);
  //myopts.theta_slices = options.get_option("theta_slices",0);



//   ID phi_slices   = options.get_option("phi_slices",0);
//   ID theta_slices = options.get_option("theta_slices",0);

//   //Compute_directions
//   double min_theta, max_theta, min_phi, max_phi; 
//   double d_theta;
//   double d_phi;
//   double total_angle;
//   double weight;
//   ID n_slices;

//   switch (SimulationOptions::dim)
//     {
      
//     case 1 :
      
//       theta_slices = 1;
//       phi_slices = 2;
      
//       min_theta = 0.0;
//       max_theta = M_PI;

//       min_phi = M_PI * 0.5;
//       max_phi = M_PI * 2.0 + M_PI * 0.5;
      
//       weight = 2.0 * M_PI;
      
//       n_slices = theta_slices * phi_slices;
//       myopts.spec.resize(n_slices);
      
//       break;
      
//     case 2 :
      
//       //std::cout<<"2D"<<std::endl;
     
//       if (phi_slices == 0)
//       	phi_slices = 4;
      
      
//       min_theta = 0.0;
//       max_theta = M_PI;
      
//       min_phi = M_PI * 0.5;
//       max_phi = M_PI * 2.0 + M_PI * 0.5;

      
//       weight =  1.0;
      
//       n_slices = theta_slices * phi_slices;
//       myopts.spec.resize(n_slices);
    
      
//       break;
      
//     case 3 :
      
//       //std::cout<<"3D"<<std::endl;
      
//       if (theta_slices == 0)
// 	theta_slices = 2;
      
//       if (phi_slices == 0)
// 	phi_slices = 4;
      
      
//       min_theta = 0.0;
//       max_theta = M_PI;
      
//       min_phi = 0.0;
//       max_phi = M_PI * 2.0;
      
//       weight = 1.0;
      
      
//       n_slices = theta_slices * phi_slices;

//       myopts.spec.resize(n_slices);

    
//       break;
//     }
    

//     total_angle = 4.0 * M_PI;
//     myopts.directions.resize(n_slices);
//     myopts.d_omega.resize(n_slices);
//     myopts.dir.resize(n_slices);
//     myopts.theta_vec.resize(n_slices);
//     myopts.phi_vec.resize(n_slices);
    
//     d_theta =  (max_theta - min_theta) / theta_slices;
//     d_phi =  (max_phi - min_phi) / phi_slices;
//     double theta, phi;


//     ID k = 0;
//     for (ID n_phi = 0; n_phi < phi_slices; n_phi++)
//     {
//       // phi = min_phi + d_phi * 0.5 + d_phi * n_phi;
//       phi = min_phi + d_phi * n_phi;
      
//       for (ID n_theta = 0; n_theta < theta_slices; n_theta++)
//       {
// 	theta = min_theta + d_theta * 0.5 + d_theta * n_theta;
	
// 	myopts.d_omega[k] =  2.0 * sin(theta) * sin (0.5 * d_theta) * d_phi;
	
// 	myopts.dir[k](0) = sin(theta) * sin(phi);
// 	myopts.dir[k](1) = sin(theta) * cos(phi);
// 	myopts.dir[k](2) = cos(theta);
	
// 	myopts.directions[k](0) =  weight * sin(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
// 	myopts.directions[k](1) =  weight * cos(phi) * sin(0.5 * d_phi) * (d_theta - cos(2.0 * theta) * sin(d_theta) );
// 	myopts.directions[k](2) =  weight * 0.5 * d_phi * sin(2.0 * theta) * sin(d_theta);
	
// 	myopts.theta_vec[k] = theta;
// 	myopts.phi_vec[k] = phi;
	
// 	k++;
//       }
//     }

//     //Spec vectors
//     for (ID k1 = 0;k1<n_slices; k1++)
//       for (ID k2 = 0;k2<n_slices; k2++)
//       {
// 	Point sum = myopts.directions[k1] + myopts.directions[k2];
//         if (sum.size() < 1e-4)
// 	  myopts.spec[k1]=k2;
	
//       }
      


//    //cout<<_test<<endl;
}




 
