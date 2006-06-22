#ifndef __DEVICEREGION_H__
#define __DEVICEREGION_H__




// base  class
class DeviceRegion{



 public:

  DeviceRegion();

 ~DeviceRegion();

 //RegionDefinition  material_regions;




 // virtual void set_regions_data() = 0;

 // virtual void init() = 0;

 virtual unsigned int get_region_number() const= 0;


};


#endif //  
