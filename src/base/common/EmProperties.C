
#include "EmProperties.h"
//#include "DataBaseCall.h"

// EmProperties::EmProperties()
// {
// }


// EmProperties::~EmProperties()
// {
// }


void EmProperties::read_database_bowing_parameters(const Dummy& db)
{}

//void EmProperties::read_database(DataBaseCall& db1)
void EmProperties::read_database(Dummy& db1)

{

//   //DataBaseCall  db1;
//   string test;
//  //template <typename T> 
// //       T  db_class::get_data(string&  label) 
// //  OR overload

// //const Material*  PhysicalProperties::get_material(void) const;
// // get material  ,  structure
//   const string  sim_class = "1";
//   const string  model = "1";
//  // vector_query.push_back(crystal_structure);
// //    vector_query.push_back(sim_class);
//  //   vector_query.push_back(model);
//     vector<string> vector_query;

//     vector_query.push_back( get_material()->get_name() );
//     vector_query.push_back(get_material()->get_structure() );
//     vector_query.push_back(sim_class);
//     vector_query.push_back(model);
//     db1.set_query(vector_query);
//     test=db1.get_data("lattice_const_a");
//     cout << "TEST+-+-+-+-+-+-DB:   "  << test <<endl ;
    
//       cout << " MATERIAL = "<< (get_material()->get_name())<< endl;
//   //  (this->latt_const_a) = atof(test.c_str());
//     latt_const_a = atof(test.c_str());
    
//     if (test == "")
//     {latt_const_a = 0.0;}
    
// //dbcall.set_query(vector_query);
// // latt_const_a = dbcall.get_data<double>("lattice_const_a")


//   // get_material -> get_structure
// //  if ( ( get_material()->get_structure() == "zb"  )  &&
// //       ( get_material()->get_name() == "GaAs") )
// //  {
// //    cout << " MATERIAL = GaAs  " << endl;
// //    (this -> latt_const_a) = 0.1543;
// //  }
// //  
// //  
// //  else if ( ( get_material()->get_structure() == "zb"  )  &&
// //            ( get_material()->get_name() == "AlAs") )
// //  {
// //    cout << " MATERIAL = AlAs " << endl;
// //    (this->latt_const_a) = 0.743;
// //  }
// //  else {latt_const_a = 0.0;}


//   //   const Material* _material;
//   //   latt_const_a =
//   //dbcall.set_query(vector_query);
//   //latt_const_a = dbcall.get_data<double>("lattice_const_a")


}

void EmProperties::set_properties_alloy(const  PhysicalProperties* prop_comp1,
                                        const  PhysicalProperties* prop_comp2,
                                        double molar_fraction)
{
  
  cout << "  molar_fraction,  ((EmProperties*) prop_comp2) ->latt_const_a " <<
    molar_fraction << "   "  <<  ((EmProperties*) prop_comp2) ->latt_const_a << endl ;
  
  latt_const_a = molar_fraction *( ((EmProperties*) prop_comp1) ->latt_const_a) + (1-molar_fraction)*
    (  ((EmProperties*) prop_comp2) ->latt_const_a);
  
  
}



double  EmProperties::get_lattice_constant_a()
{

  return latt_const_a;

}









