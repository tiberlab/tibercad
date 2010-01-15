// $Id$

#ifndef _HASHSET_H_
#define _HASHSET_H_

#include "tiber_config.h"

#if defined(HAVE_UNORDERED_MAP)
# include <unordered_set>
#elif defined(HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_set>
#elif defined(HAVE_HASH_MAP)
# include <hash_set>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_set>
#else
# include <set>
#endif


namespace TiberCad
{

  template <typename K>
  struct HashSet
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_set<K> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_set<K> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_set<K> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_set<K> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_set<K> Type;
# else
    typedef std::set<K> Type;
# endif
#else
    typedef std::set<K> Type;
#endif
  };


  template <typename K>
  struct HashMultiSet
  {
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_multiset<K> Type;
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_multiset<K> Type;
#elif defined(HAVE_HASH_MAP)
    typedef std::hash_multiset<K> Type;
#elif defined(HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_multiset<K> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multiset<K> Type;
# else
    typedef std::multiset<K> Type;
# endif
#else
    typedef std::multiset<K> Type;
#endif
  };

}

#endif // _HASHSET_H_
