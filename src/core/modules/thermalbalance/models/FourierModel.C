// $Id$

#include "FourierModel.h"

TIBER_MODULE(FourierModel, HeatTransport, fourier)


using namespace std;


FourierModel::FourierModel(const ModelOptions& options):HeatTransportModel(options)
{
  set_type(HeatTransportModel::Fourier);
}


void
FourierModel::do_init(void)
{
  //get_parameter("test", _test);
  //  cout<<_test<<endl;
}





