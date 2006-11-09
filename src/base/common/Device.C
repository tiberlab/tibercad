
#include <iostream>  
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include "Device.h"






Device::Device()
{
  material_regions.clear();	
}

Device::~Device()
{

}


void
Device::init_device(const  vector<RegionDefinition>& dev_reg)
{
  
  set_material_regions( dev_reg);
  set_materials ();
  set_map_ID_material_region();
  makes_set_of_materials();
}

//***************************************

void Device::set_map_alloy_model(const map <unsigned int, AlloyModel*>&  map_alloy_model)
{
  
  reg_alloy_model_map = map_alloy_model;
  
}
//***************************************




  //  public  method to get  material regions description (from input  file)
void 
Device::set_material_regions( const  vector<RegionDefinition>& dev_reg)
{
  // device_regions = dev_reg;
  material_regions = dev_reg;    //   vector  of  RegionDefinition
  
  
  //  TEST
  //   set_materials ();
  //   set_map_ID_material_region();
  //   
}

//DeviceRegion*  dev_region_point;

//vector<DeviceRegion> dev_region_vec;





// //   crea  DeviceRegion ,  p. es. una  Material Reg.   con  material_regions[i] 
// //  +  Material*
// void Device::init_device()
// {
//  for  (int i = 0; i < number_material_regions; i++)
//     {




//       dev_region_point = &dev_region_vec[i];

//       ((Materialregion *)dev_region_point)->init(  material_regions[i])

//     }

//  for  (int i = 0; i < number_doping_regions; i++)
//     {
//       // init DopingRegion

//     }

// for  (int i = 0; i < number_model_regions; i++)
//     {
//       // init ModelRegion

//     }


// }



// NEW 24/2/06
//
// void
// Device::set_materials()
// {
//   string mat_name;
//   string  crystal_name;
//   unsigned int reg_id;

//   int number_material_regions = material_regions.size ();
//   map < string, Material* >::iterator p;
//   map < unsigned int, Alloy* > ::iterator p_alloy;  //  reg_id_alloy_map;  
  
//   map <unsigned int, AlloyModel*>::iterator p_alloy_model;

//   for (int i = 0; i < number_material_regions; i++)
//   {
//     mat_name = material_regions[i].get_material_name ();  // get   material name and  crystal structure from region definition data
//     // inline string  RegionDefinition::get_crystal_name() const
//     crystal_name = material_regions[i].get_crystal_name();
//     reg_id = material_regions[i].get_region_number();
//     //if call_db(database,mat_name,crystal_struct, alloy ) !=  fail
//  //   if ( ( reg_id == 1) || ( reg_id == 2) )
//     if (false)
//     //**********************************************
//       //if  mat_name is  alloy  !!!
//       // ****************************************
//       // makes an  Alloy  object for  each material region
//       // and pass also pointer to AlloyModel
//       // then put ALL Alloy objects in  map <id,Alloy>
//     { 
//       // find  alloy model associated to present material region
//       p_alloy_model = reg_alloy_model_map.find(reg_id);
      
//       if (p_alloy_model != reg_alloy_model_map.end())
//       {alloy_model_point = p_alloy_model->second ;}
//       //
//       // test
//       const double  x = alloy_model_point -> get_x_min();
//       cout << "molar frac********* = " << x << endl;
//       //
            
//       //else error
//       //
//       matpoint_alloy = new Alloy (mat_name,crystal_name,alloy_model_point );
//       // insert in map <region id, pointer to alloy>
//       // Alloys  are  inserted in a special map  which associate them to region ID,
//       // because there can be  2 or  more Alloy objects with the  same name 
//       // but  different AlloyModel properties,
//       // while 2 physical regions with the  same SIMPLE material point to the same 
//       // Material object (so there is only  one  Material object  with a given name)
//       // Instead 2 Alloy object can have the  same  name but different AlloyModel
//       // 
//       reg_id_alloy_map.insert(make_pair (reg_id, matpoint_alloy));
//       // call_db -> vector<string> = vector(alloy components)
//       // read from database:
//       alloy_components.push_back("AlN"); // AlAs");
//       alloy_components.push_back("GaN"  );  //As");
                      
//       for (int i = 0; i < alloy_components.size (); i++)
//       {
//         p = name_mat_map.find (alloy_components[i]);
//         //  cout << " alloy_components[i]  " << alloy_components[i] << endl;

//         if (p == name_mat_map.end ())
//           //  not yet present
//         {
//           matpoint = new Material (alloy_components[i],crystal_name  );
//           cout << "***alloy_components[i]**** " << alloy_components[i] << endl;
//           name_mat_map.insert (make_pair (alloy_components[i], matpoint));

//           // put mat.  components in  alloy  object
//           matpoint_alloy->set_components (matpoint);

//         }

//         else
//           // comp. mat.. already  present   in  map
//         {
//           // put Material*  in alloy  as  a  component
//           matpoint_alloy->set_components (p->second);

//         }
//       }
//     }
              
//     else 
//       // not  alloy !!!
//       //  for simple  materials we  make an  object for EACH distinct
//       // material name  and  put Material in  map <name, Material>
//     { 
//       p = name_mat_map.find (mat_name);
//       if (p == name_mat_map.end ())
//         // not found in  material map   
//       {
//         matpoint = new Material (mat_name,crystal_name );
//         name_mat_map.insert (make_pair (mat_name, matpoint)); 
//       }
            
//     }
           
           
//   }
  
    
// }
              
              
            
            
            
void
Device::makes_set_of_materials() 
//
//makes set with all the  Material* (both simple material and alloy  material,
// included material components of  each alloy material)
//  to  be  visited from outside with an  iterator . So we can add properties and 
// init ALL the materials in a  single  loop (both Material and  Alloy objects)
//

{
  
  ////in "main":
  //  
  //  
  //  set<Material *>::iterator it = set_all_materials.begin():
  //  const set<Material *>::const_iterator end = set_all_materials.end();
  //  for ( ; it != end; ++it)  
  //    {
  //      // it points to the element of  set (Material* or  Alloy*)
  //      ( it->add_properties   ) ;
  //      it->init(db);
  //      
  //    }
  //  
  //  
  //map < string, Material* >name_mat_map;
  //map < unsigned int, Alloy* >reg_id_alloy_map;
   
   
  //set<Material *>  set_all_materials;
   
  map < string, Material* >::iterator it = name_mat_map.begin();
  const map < string, Material* >::const_iterator end = name_mat_map.end();
   
  for ( ; it != end; ++it)  
  {
    set_all_materials.insert( it->second) ;
  }
    
  map < unsigned int, Alloy* >::iterator it_alloy = reg_id_alloy_map.begin();
  const map < unsigned int, Alloy* >::const_iterator end_alloy = reg_id_alloy_map.end();
    
  for ( ; it_alloy != end_alloy; ++it_alloy)  
  {
    set_all_materials.insert( it_alloy->second) ;
  }

  
  
}       
         
         
set<Material *> &  
Device::get_set_all_materials()
{
  
  return  set_all_materials;
  
}         
          
    
    
    
//
//
//// old 14/2/06
//// ********************************
//
////in  .h
////map < string ,Material* > name_mat_map;
////Material* matpoint, matpoint_alloy 
//
//
//void
//Device::set_materials ()
//{
//  string mat_name;
//  string  crystal_name;
//
//  int number_material_regions = material_regions.size ();
//  map < string, Material * >::iterator p;
//
//
//  for (int i = 0; i < number_material_regions; i++)
//  {
//    mat_name = material_regions[i].get_material_name ();  // get   material name and  crystal structure from region definition data
//   // inline string  RegionDefinition::get_crystal_name() const
//    crystal_name = material_regions[i].get_crystal_name();
//    
//    
//    // first : if is  alloy 
//    //  {
//    //    put in map_alloy = map <reg_id, matpoint_alloy>
//   //  id = get_region_number()  ,  find id in map <reg_id , AlloyModel*> -> &AlloyModel
//   //  new Alloy (mat_name,crystal_name, AlloyModel*);
//   //   }
//   // else  //not alloy !!
//   //........
//    
//    p = name_mat_map.find (mat_name);
//    if (p == name_mat_map.end ())
//	      // not found in  map
//    {
//	         
//	         //if call_db(database,mat_name,crystal_struct, alloy ) !=  fail
//      if (false)
//	         //mat_name is  alloy  !!!
//      {
//
//        matpoint_alloy = new Alloy (mat_name,crystal_name);
//        name_mat_map.insert (make_pair (mat_name, matpoint_alloy));
//
//	             // call_db -> vector<string> = vector(alloy components)
//        for (int i = 0; i < alloy_components.size (); i++)
//        {
//          p = name_mat_map.find (alloy_components[i]);
//
//          if (p == name_mat_map.end ())
//		               //  not yet present
//          {
//            matpoint = new Material (alloy_components[i],crystal_name  );
//            name_mat_map.insert (make_pair (mat_name, matpoint));
//
//		                   // put mat.  components in  alloy  object
//            matpoint_alloy->set_components (matpoint);
//
//          }
//
//          else
//		                 // comp. mat.. already  present   in  map
//          {
//		                   // put Material*  in alloy  as  a  component
//            matpoint_alloy->set_components (p->second);
//
//          }
//
//        } 
//
//      }
//
//      else
//
//      {
//	             // mat_name is  not alloy 
//        matpoint = new Material (mat_name,crystal_name );
//        name_mat_map.insert (make_pair (mat_name, matpoint)); 
//                             
//      }
//
//    }  
//
//
//  }
//
//
//
//}
//
//
//

// *******************************



// NEW  14/2/06
// ********************************
// instatiate an object MaterialRegion for  each  RegionDefinition present in vector 
//  material_regions, then makes a map <MaterialRegion * , 
// fai  una  mappa <reg ID,  MaterialRegion>  
void
Device::set_map_ID_material_region()
{
  Material*  mat_point;	
  string mat_name;
  int number_material_regions = material_regions.size();
  map < string ,Material* >::iterator p; // map <region name, simple Material>
  
  map < unsigned int, Alloy* >::iterator p_alloy; //reg_id_alloy_map;
  unsigned int reg_id;
  
  
  // *************************************************

  //  for  (int i = 0; i < number_material_regions; i++)
  //  {
  //      // pointer to a MaterialRegion object 
  //      
  //      // pass  also  Material* !!
  //    mat_name = material_regions[i].get_material_name();
  //    	
  //    p = name_mat_map.find(mat_name);
  //    	
  //    if (p != name_mat_map.end())
  //    {mat_point = p->second ;}
  //      //  		else
  //      //  		{ //  ERROR}}
  //      //  		}
  //  		
  //      //  mat_reg_point  = new MaterialRegion (material_regions[i]);
  //      // makes object MaterialRegion which points to Material*
  //    mat_reg_point  = new MaterialRegion( material_regions[i], mat_point);
  //
  //    reg_ID = mat_reg_point -> get_region_number();
  //
  //    ID_mat_reg_map.insert(make_pair(reg_ID, mat_reg_point));
  //
  //  }
  
  // *****************************************************
  
  
  // for  alloy :
  // search in map <reg_id, matpoint_alloy> 
  // 

  // NEW  24/02/06

  for  (int i = 0; i < number_material_regions; i++)
  {

    mat_name = material_regions[i].get_material_name();
    reg_id = material_regions[i].get_region_number();
    cout << "reg_id =  " <<  reg_id <<  endl;
    
    p_alloy = reg_id_alloy_map.find(reg_id);
    // search in alloy map
    if (p_alloy != reg_id_alloy_map.end())
    {mat_point = p_alloy->second ;}
    else
    {  //  not alloy,  search  in  material  map
      p = name_mat_map.find(mat_name);
      if (p != name_mat_map.end())
      {mat_point = p->second ;}
    }
    //else  // error
    
    // makes object MaterialRegion which points to Material*
    mat_reg_point  = new MaterialRegion( material_regions[i], mat_point);
    ID_mat_reg_map.insert(make_pair(reg_id, mat_reg_point));


  }


 
 
  //  for  regions with the SAME  alloy (eg. AlGaAs) : 
  // one  MaterialRegion for each  region which  points 
  //  to a DIFFERENT Material object with the  same  material 
  //name but with DIFFERENT (in general) AlloyModel
  //
  // reg_id = material_regions[i].get_region_number();
  // p_alloy = alloy_map.find(reg_id);
  // if p_alloy != alloy_map.end())
  //  {mat_point = p_alloy->second ;}
  //  else  // not  alloy
  //       p = name_mat_map.find(mat_name);
  //      
  //        if (p != name_mat_map.end())
  //        {mat_point = p->second ;}
 
 

}








// ********************************************************************************************************



void
Device::set_map_ID_doping_region()
{

  // ..............................

  //phys_ID_dev_reg_map.insert(make_pair(reg_ID, ...));


}


void 
Device::set_map_ID_model_region()
{

  // ...........

}





MaterialRegion* 
Device::get_material_region(unsigned int  ID)
{

  map < unsigned int,MaterialRegion* >::iterator p;

  p = ID_mat_reg_map.find(ID);
  if (p != ID_mat_reg_map.end())
    return p->second ;
  else
  {
    cerr << "ID not found in Material Regions";
    exit(1);
  } 

}


// ************************************************
//  TO BE IMPLEMENTED 
// ************************************************


// DopingRegion* Device::get_doping_region(unsigned int  ID)
// {

//   map < unsigned int,DopingRegion* >::iterator p;

//   p = ID_dop_reg_map.find(ID);
//   if (p != ID_dop_reg_map.end())
//   return p->second ;
//   else
//     {
//       cerr << "ID not found in Doping Regions";
//       exit(1);
//     } 


// }


// ModelRegion* Device::get_model_region(unsigned int  ID)
// {

//   map < unsigned int,ModelRegion* >::iterator p;

//   p = ID_model_reg_map.find(ID);
//   if (p != ID_model_reg_map.end())
//   return p->second ;
//   else
//     {
//       cerr << "ID not found in Model Regions";

//       return NULL;
     
//     } 


// }


// **********************************************************




// OLD  OBSOLETE ????
// void Device::set_device_structure(vector<string>&  region_name_v, vector<unsigned int>& region_number_v,
// 				  vector<string>&  material_name_v, vector<double>& doping_concentration_v, 
// 				  vector<string>& doping_type_v   )

// {

//   // RegionDefinition current_region;

//   //  string r_name, r_num, m_name, dop_c, dop_t;
//   //   r_name = "reg_name";
//   //   r_num = "reg_numb";
//   //   m_name = "mat_name";
//   //   dop_c = "dop_conc";
//   //   dop_t = "dop_type";

//   device_regions.resize(region_number_v.size()); // allocate  correct  number of  objects RegionDefinition

//   for (int i =0; i< region_number_v.size();++i)
//     {

//       device_regions[i].set_region_name(region_name_v[i]);
//       device_regions[i].set_region_number(region_number_v[i]);
//       device_regions[i].set_material_name(material_name_v[i]);
//       device_regions[i].set_doping_concentration(doping_concentration_v[i]);
//       device_regions[i].set_doping_type(doping_type_v[i]);



//       //    cout <<  "here............."<< endl;
//       //       device_regions[i].set_field(r_name, reg_name_v[i]);
//       //  cout <<  "here............."<< endl;
//       //       device_regions[i].set_field(r_num, reg_numb_v[i]);
//       //       device_regions[i].set_field(m_name, mat_name_v[i]);
//       //       device_regions[i].set_field(dop_c, dop_conc_v[i]);
//       //       device_regions[i].set_field(dop_t, dop_type_v[i]);




//       //RegionDefinition.set_field


//       //  current_region.reg_name = reg_name_v[i];
//       //       current_region.reg_numb = reg_numb_v[i];
//       //       current_region.mat_name = mat_name_v[i];
//       //       current_region.dop_conc = dop_conc_v[i];
//       //       current_region.dop_type = dop_type_v[i];

//       //    device_regions.push_back(current_region);

//       //       current_region.reg_name = "" ;
//       //       current_region.reg_numb = 0;
//       //       current_region.mat_name = "";
//       //       current_region.dop_conc = 0.0;
//       //       current_region.dop_type = "";



//     }


//   // string file_name,section_name ;
//   //    file_name = "options.in";
//   //    section_name = "Regions";

//   //    device_regions.clear();

//   //   InputParser  input_device(file_name,section_name  );   
 
//   //   //input_device.get_device_data( vector<region_definition>& device_regions )
//   //   input_device.get_device_data( device_regions );

//   //   cout << "device_regions[0].reg_name " << device_regions[0].reg_name;


// }




// ?????????????????????????
// ***************************************************************
// in main : ref_reg_struct = device_object.get_device_data(reg_query)
//  value = ref_reg_struct.get_reg_name();


const  RegionDefinition& 
Device::get_device_data(unsigned int region_query)
{ 
  unsigned int numb;

  for (int i =0; i< material_regions.size();++i)
  {

    numb = material_regions[i].get_region_number();

    if  (numb == region_query) 
    {
      return   (material_regions[i]) ;

    }
  }
}
//????????????????????????????????????

// **********************************************************************


//  BC regions:   TO  BE  PUT  IN  "SolverEnvironment"  object (?)

void 
Device::set_device_boundary_cond(vector<string>&  BC_region_name_v,
                                 vector<unsigned int>& BC_region_number_v,
                                 vector<string>&  BC_type_v, 
                                 vector<double>& BC_value_v   )
{

  BC_device_regions.resize(BC_region_number_v.size()); // allocate  correct  number of  objects BcRegionDefinition

  for (int i =0; i< BC_region_number_v.size();++i)
  {

    BC_device_regions[i].set_BC_region_name(BC_region_name_v[i]);
    BC_device_regions[i].set_BC_region_number(BC_region_number_v[i]);
    BC_device_regions[i].set_BC_type(BC_type_v[i]);
    BC_device_regions[i].set_BC_value(BC_value_v[i]);
     
  }

}





const  BcRegionDefinition& 
Device::get_device_boundary_cond(unsigned int region_query)
{ 

  unsigned int numb;

  for (int i =0; i< BC_device_regions.size();++i)
  {

    numb = BC_device_regions[i].get_BC_region_number();

    if  (numb == region_query) 
    {	 
      return   (BC_device_regions[i]) ;

    }
  }
}


// *************************************************************






//  OLD *****************************************************
//************************************************************

  // string  Device::get_device_data(unsigned int reg_query, string& field)

  // // const  RegionDefinition&                                                           

  // {

  //   unsigned int  numb;
  //   string field_value, null;
  //   string reg_numb = "reg_numb";

  //   null = "";

  //   for (int i =0; i< device_regions.size();++i)
  //     {

  //       numb = device_regions[i].find_field2(reg_numb);

  //       // if (device_regions[i].reg_numb == reg_query)
  //       if  (numb == reg_query) 
  // 	{

  // 	  field_value = device_regions[i].find_field(field);
  // 	  return   field_value;

  // 	}
	 
   	  
  //     }

  //    cerr <<  "Error : wrong  value  of query "; 
  //    return null;

  // }



  // double  Device::get_device_data2(unsigned int reg_query, string& field)

  // {

  //   double  field_value, null;
  //   null = 0.0;

  //   for (int i =0; i< device_regions.size();++i)
  //     {

  //       if (device_regions[i].reg_numb == reg_query)
  // 	{

  // 	  field_value = device_regions[i].find_field(field);
  // 	  return   field_value;

  // 	}
	 
   	  
  //     }

  //    cerr <<  "Error : wrong  value  of query "; 
  //    return null;

  // }
