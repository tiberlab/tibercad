// $Id$

#ifndef _HASHMAP_H_
#define _HASHMAP_H_
    
#include "tiber_config.h"

#if defined(HAVE_UNORDERED_MAP)
# include <unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_map>
#elif defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#else
# include <map>
#endif


namespace TiberCad
{

  template <typename K, typename V>
  struct HashMap
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_map<K, V> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_map<K, V> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_map<K, V> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_map<K, V> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_map<K, V> Type;
# else
    typedef std::map<K, V> Type;
# endif
#else
    typedef std::map<K, V> Type;
#endif
  };


  template <typename K, typename V>
  struct HashMultiMap
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_multimap<K, V> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_multimap<K, V> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_multimap<K, V> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_multimap<K, V> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<K, V> Type;
# else
    typedef std::multimap<K, V> Type;
# endif
#else
    typedef std::multimap<K, V> Type;
#endif
  };

}

#endif // _HASHMAP_H_
