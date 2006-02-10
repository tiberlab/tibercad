#ifndef __DEVICEREGION_H__
#define __DEVICEREGION_H__

#include <iostream>  
#include <sstream>
#include <vector>
#include <string>



using namespace std;

// base  class
class DeviceRegion{



 public:

  DeviceRegion();

 ~DeviceRegion();

 //RegionDefinition  material_regions;




 // virtual void set_regions_data() = 0;

 // virtual void init() = 0;

 virtual unsigned int get_region_number() = 0;


};


#endif //  
