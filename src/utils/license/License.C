// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "License.h"
#include "Utils.h"
#include "Messages.h"
#include "TiberCad.h"

#ifdef LICENSE_CHECK
# include "lmx.h"
#endif

#include <iostream>
#include <sstream>
#include <fstream>


using namespace std;


namespace
{
#ifdef LICENSE_CHECK
  LMX_HANDLE LmxHandle;

  /* Callbacks which prints out various states when connection is lost to license server. */
  void LMX_CALLBACK UserConnectionLostRoutine(void *pVendorData, const char *szHost, int nPort, int nFailedHeartbeats)
  {
    /* Only only print out first time we get a notification of a lost server. */
    if (nFailedHeartbeats == 1)
      printf("Lost connection to server %s : %d\n", szHost, nPort);
  }

  void LMX_CALLBACK UserCheckoutFailureRoutine(void *pVendorData, const char *szFeatureName, int nUsedLicCount, LMX_STATUS LmxStat)
  {
    printf("Server could not re-checkout feature because:\n%s\n", LMX_GetErrorMessageSimple(LmxStat));
  }

  void LMX_CALLBACK UserCheckoutSuccessRoutine(void *pVendorData, const char *szFeatureName, int nUsedLicCount)
  {
    printf("Successful regain %d\n", nUsedLicCount);
  }

  void LMX_CALLBACK UserRetryRoutine(void *pVendorData, const char *szFeatureName, int nUsedLicCount)
  {
    printf("About to retry regaining feature %d\n", nUsedLicCount);
  }

  void LMX_CALLBACK UserExitRoutine(void *pVendorData)
  {
    printf("License server is down or unavailable. Your work will be saved and this application will shut down.\n");

    /* Remember to save user's work. */
    /* Shutdown in a very cruel way. */
    exit(1);
  }
#endif

}


void
License::init(void)
{
#ifdef LICENSE_CHECK
  if (LMX_Init(&LmxHandle) != LMX_SUCCESS)
  {
    printf("Unable to initialize license system!\n");
    exit(1);
  }

  /* Enable the callback functions with license server failure behavior. */
  LMX_SetOption(LmxHandle, LMX_OPT_HEARTBEAT_CHECKOUT_FAILURE_FUNCTION, (LMX_OPTION) UserCheckoutFailureRoutine);
  LMX_SetOption(LmxHandle, LMX_OPT_HEARTBEAT_CHECKOUT_SUCCESS_FUNCTION, (LMX_OPTION) UserCheckoutSuccessRoutine);
  LMX_SetOption(LmxHandle, LMX_OPT_HEARTBEAT_RETRY_FEATURE_FUNCTION, (LMX_OPTION) UserRetryRoutine);
  LMX_SetOption(LmxHandle, LMX_OPT_HEARTBEAT_CONNECTION_LOST_FUNCTION, (LMX_OPTION) UserConnectionLostRoutine);
  LMX_SetOption(LmxHandle, LMX_OPT_HEARTBEAT_EXIT_FUNCTION, (LMX_OPTION) UserExitRoutine);
  
  /* Enable automatic heartbeats. */
  LMX_SetOption(LmxHandle, LMX_OPT_AUTOMATIC_HEARTBEAT_ATTEMPTS, (LMX_OPTION) 5);

  /* Enable automatic server discovery. */
  LMX_Putenv("LMX_AUTOMATIC_SERVER_DISCOVERY=1");
#endif

}


void
License::close(void)
{
#ifdef LICENSE_CHECK
  LMX_Free(LmxHandle);
#endif
}


void
License::check_out(const std::string& feature, int major, int minor, int count)
{
#ifdef LICENSE_CHECK
  ostringstream os;
  os << "Trying to check out " << count << " license(s) for "
   << feature;
  Messages::info(os.str());

  if (LMX_Checkout(LmxHandle, feature.c_str(), major, minor, count) != LMX_SUCCESS)
  {
    printf("Unable to checkout:\n");
    printf("%s\n", LMX_GetErrorMessage(LmxHandle));
    
    LMX_Free(LmxHandle);
    exit(1);
  }
#endif
}


void
License::check_in(const std::string& feature, int count)
{
#ifdef LICENSE_CHECK
  if (count <= 0)
    count = LMX_ALL_LICENSES;

  LMX_Checkin(LmxHandle, feature.c_str(), count);
#endif
}


