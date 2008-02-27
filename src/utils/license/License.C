// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "License.h"
#include "Utils.h"
#include "TiberCad.h"

#include "hex.h"
#include "randpool.h"
#include "rsa.h"
#include "validate.h"

#include <iostream>
#include <sstream>
#include <fstream>


using namespace CryptoPP;
using namespace std;


namespace
{
  RandomPool randomPool;
}



const char*
License::_public_key = "30819D300D06092A864886F70D010101050003818B00308187028181008B2D541C099BD67293B49C95220379A4B50DAA85897C733D61C8E5892E173F34166C92D7D0365CCFFA76F7BA72E84DC97002AD39F5E1351AA68DC23F22EBFD0F2B494B7DB0088ADFE688B3F8967CB43917B5DAA4419B3C7A277671AD1CC01FC1E08EC76549C8A60369B4FCA4ECB256F342285195AA52B3870DF8D948D7F99B01020111";




bool
License::read_license_file(const std::string& licensefile, LicenseData& data)
{

  fstream file;
  file.open(licensefile.c_str(), ios::in);
  if (file.fail())
  {
#ifdef DEBUG
    cerr << "Cannot open licensefile " << licensefile << "!" << endl;
#endif
    return false;
  }

  const int bufsize = 512;
  char linebuf[bufsize];
  char strbuf[bufsize];
  
  // get rid of comments
  do {
    file.getline(linebuf, bufsize);
  } while (!file.eof() && (linebuf[0] == '#'));

  string dummy, holder, expiry, signature;

  // now we should be on the line with label HOLDER
  {
    istringstream istr(linebuf);
    istr >> dummy;
    boost::algorithm::trim(dummy);
    if (dummy != "HOLDER")
    {
      cerr << "Your license file seems to be corrupt!" << endl;
      return false;
    }

    istr.getline(strbuf, bufsize);
    holder = strbuf;
    boost::algorithm::trim(holder);
  }

  {
    file.getline(linebuf, bufsize);

    istringstream istr(linebuf);
    istr >> dummy;
    boost::algorithm::trim(dummy);
    if (dummy != "EXPIRY")
    {
      cerr << "Your license file seems to be corrupt!" << endl;
      return false;
    }

    istr.getline(strbuf, bufsize);
    expiry = strbuf;
    boost::algorithm::trim(expiry);
  }

  {
    file.getline(linebuf, bufsize);

    istringstream istr(linebuf);
    istr >> dummy;
    boost::algorithm::trim(dummy);
    if (dummy != "LICENSE")
    {
      cerr << "Your license file seems to be corrupt!" << endl;
      return false;
    }

    istr.getline(strbuf, bufsize);
    signature = strbuf;
    boost::algorithm::trim(signature);
  }

  data.holder = holder;
  data.expiry = expiry;
  data.signature = signature;

  return true;
}


bool
License::create_license(std::string& licensefile, const std::string& private_key)
{
  bool result = false;

  LicenseData data;
  if (read_license_file(licensefile, data))
  {

    const char* header =
        "#\n" \
        "# TiberCAD license file\n" \
        "#\n" \
        "# Do not touch the lines below!\n" \
        "#\n";

    string to_sign = data.holder + data.expiry;

    string signed_string;
    sign_string(private_key, to_sign, signed_string);

    ofstream out(licensefile.c_str());
    out << header
      << "HOLDER " << data.holder << endl
      << "EXPIRY " << data.expiry << endl
      << "LICENSE " << signed_string << endl
      << "END";

    result = true;
  }

  return result;
}


bool
License::check_license(const std::string& tool)
{
  bool result = false;

  string default_licfile("tibercad.lic");
  string licensefile(default_licfile);

  LicenseData data;

  string tiberroot;
  char* root = getenv("TIBERCADROOT");
  if (root != NULL)
    tiberroot = string(root);
  
  // we could have specified it in the environment
  char* licfile = getenv("TIBERLICENSEFILE");
  if (licfile != NULL)
  {
    licensefile = licfile;
    if (read_license_file(licensefile, data))
      result = verify_string(_public_key, data.holder + data.expiry, data.signature);
  }
  
  // if not, we check a few default locations
  if (!result)
  {
    licensefile = default_licfile;
    if (read_license_file(licensefile, data))
      result = verify_string(_public_key, data.holder + data.expiry, data.signature);
  }

  if (!result)
  {
    licensefile = tiberroot + "/" + default_licfile;
    if (read_license_file(licensefile, data))
      result = verify_string(_public_key, data.holder + data.expiry, data.signature);
  }
  
  if (!result)
  {
    licensefile = tiberroot + "/license/" + default_licfile;
    if (read_license_file(licensefile, data))
      result = verify_string(_public_key, data.holder + data.expiry, data.signature);
  }


  if (result)
  {
    cout << "Found license file: " << licensefile << ". "
      << "Validate ... ";
    // now we check the system date
    time_t t = time(NULL);
    tm* now = gmtime(&t);

    vector<string> tokens(3);
    size_t index1 = data.expiry.find_first_of('-');
    size_t index2 = data.expiry.find_last_of('-');
    assert(index1 > 0);
    assert(index1 < data.expiry.size());

    tokens[0] = data.expiry.substr(0, index1);
    tokens[1] = data.expiry.substr(index1 + 1, index2 - index1 - 1);
    tokens[2] = data.expiry.substr(index2 + 1, data.expiry.size() - index2 - 1);

    int exp_year = atoi(tokens[0].c_str());
    int exp_month = atoi(tokens[1].c_str());
    int exp_day = atoi(tokens[2].c_str());

    int now_year = now->tm_year + 1900;
    int now_month = now->tm_mon + 1;
    int now_day = now->tm_mday;

    if (now_year > exp_year)
      result = false;
    else if (now_year == exp_year)
      if (now_month > exp_month)
        result = false;
      else if (now_month == exp_month)
        if (now_day > exp_day)
          result = false;

    if (result)
      cout << "OK" << endl;
    else
      cout << "expired license" << endl;
  }

  return result;
}


bool
License::verify_string(const string& pubkey, const string& message,
    const string& signature)
{
  StringSource pubString(pubkey, true, new HexDecoder);
  RSASS<PKCS1v15, SHA>::Verifier pub(pubString);

  StringSource signatureString(signature, true, new HexDecoder);
  if (signatureString.MaxRetrievable() != pub.SignatureLength())
    return false;
  SecByteBlock sig(pub.SignatureLength());
  signatureString.Get(sig, sig.size());

  VerifierFilter *verifierFilter = new VerifierFilter(pub);
  verifierFilter->Put(sig, pub.SignatureLength());
  StringSource f(message, true, verifierFilter);

  return verifierFilter->GetLastResult();
}


void
License::sign_string(const string& privkey, const string& message,
    string& signature)
{
  StringSource privString(privkey, true, new HexDecoder);
  RSASS<PKCS1v15, SHA>::Signer priv(privString);
  StringSource f(message, true, new SignerFilter(randomPool, priv,
        new HexEncoder(new StringSink(signature))));
}



void
License::hex_encode(const string& in, string& out)
{
  StringSource(in, true, new HexEncoder(new StringSink(out)));
}




void
License::hex_decode(const string& in, string& out)
{
  StringSource(in, true, new HexDecoder(new StringSink(out)));
}


