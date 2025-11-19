// $Id: GrayModel.C 2457 2011-03-06 23:52:12Z gromano $

#include "GrayModel.h"

#include "tibercad/module/TiberModule.h"


using namespace std;


GrayModel::GrayModel(const ModelOptions& options):HeatTransportModel(options)
{
 
 set_type(HeatTransportModel::Gray);

}

void
GrayModel::do_init(void)
{

 
}




 
