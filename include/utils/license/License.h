// $Id$

#include <string>


//! Our license handling
class License
{

  public:

    static void init(void);

    static void close(void);

    static void check_out(const std::string& feature, int major, int minor, int count = 1);

    static void check_in(const std::string& feature, int count = 0);

  private:

    License(void);

    ~License(void);

};
