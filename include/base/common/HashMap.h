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

#if   defined(HAVE_UNORDERED_MAP)
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
#elif defined(HAVE_TR1_UNORDERED_MAP)
  template <typename Key, typename Value, typename Hash = std::tr1::hash<Key> >
#elif defined(HAVE_HASH_MAP)
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
  template <typename Key, typename Value, typename Hash = __gnu_cxx::hash<Key> >
# else
  template <typename Key, typename Value, typename Hash = Key >
# endif
#else
  template <typename Key, typename Value, typename Hash = Key >
#endif
  struct HashMap
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_map<Key, Value, Hash> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_map<Key, Value, Hash> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_map<Key, Value, Hash> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_map<Key, Value, Hash> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_map<Key, Value, Hash> Type;
# else
    typedef std::map<Key, Value> Type;
# endif
#else
    typedef std::map<Key, Value> Type;
#endif
  };



#if   defined(HAVE_UNORDERED_MAP)
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
#elif defined(HAVE_TR1_UNORDERED_MAP)
  template <typename Key, typename Value, typename Hash = std::tr1::hash<Key> >
#elif defined(HAVE_HASH_MAP)
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
  template <typename Key, typename Value, typename Hash = std::hash<Key> >
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
  template <typename Key, typename Value, typename Hash = __gnu_cxx::hash<Key> >
# else
  template <typename Key, typename Value, typename Hash = Key >
# endif
#else
  template <typename Key, typename Value, typename Hash = Key >
#endif
  struct HashMultiMap
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_multimap<Key, Value, Hash> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_multimap<Key, Value, Hash> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_multimap<Key, Value, Hash> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_multimap<Key, Value, Hash> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<Key, Value, Hash> Type;
# else
    typedef std::multimap<Key, Value> Type;
# endif
#else
    typedef std::multimap<Key, Value> Type;
#endif
  };

}

#endif // _HASHMAP_H_
