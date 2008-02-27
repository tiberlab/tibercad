// $Id$

#include <string>


//! Our license handling
class License
{

  public:

    struct LicenseData
    {
      std::string holder;
      std::string expiry;
      std::string signature;
    };


    //! Returns true if license is ok
    static bool check_license(const std::string& tool = "");
    
    //! Create (sign) the license file
    static bool create_license(std::string& licensefile,
        const std::string& private_key);


  private:

    License(void);

    ~License(void);


    //! The public key
    static const char* _public_key;


    static bool read_license_file(const std::string& licensefile,
        LicenseData& data);

    static bool verify_string(const std::string& pubkey,
        const std::string& message, const std::string& signature);

    static void sign_string(const std::string& privkey,
        const std::string& message, std::string& signature);
    
    static void hex_encode(const std::string& in, std::string& out);

    static void hex_decode(const std::string& in, std::string& out);

};
