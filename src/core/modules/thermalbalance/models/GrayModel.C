// $Id$

#include "GrayModel.h"


TIBER_MODULE(GrayModel,HeatTransport, gray)


using namespace std;


GrayModel::GrayModel(const ModelOptions& options):HeatTransportModel(options)
{
 
 set_type(HeatTransportModel::Gray);

}

void
GrayModel::do_init(void)
{

 
}




 
