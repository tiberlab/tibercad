// $Id$

#ifndef _DATABASE_H_
#define _DATABASE_H_

#include <string>

class Database
{

  public:

    //! Constructor
    Database(void) : _path("."), _material("Si") {};

    void set_search_path(const std::string& path)
      { _path = path; };

    void set_material(const std::string& material)
      { _material = material; };

    const std::string get_data_file(void) const;

    bool is_alloy(const std::string& material) const;

    void get_alloy_components(const std::string& alloy,
        std::string& comp_A, std::string& comp_B) const;

  private:

    std::string _path;
    std::string _material;

};


inline
const std::string
Database::get_data_file(void) const
{
  std::string s(_path);
  s += "/" + _material + ".dat";

  return s;
}


#endif // _DATABASE_H_
