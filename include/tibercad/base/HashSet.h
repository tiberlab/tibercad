/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file HashSet.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_HASHSET_H
#define TC_HASHSET_H

#include "tibercad/base/tiber_config.h"

#if defined(TC_HAVE_UNORDERED_MAP)
# include <unordered_set>
#elif defined(TC_HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_set>
#elif defined(TC_HAVE_HASH_MAP)
# include <hash_set>
#elif defined(TC_HAVE_EXT_HASH_MAP)
# include <ext/hash_set>
#else
# include <set>
#endif


#if   defined(TC_HAVE_UNORDERED_MAP)
  template <typename Key, typename Hash = std::hash<Key> >
#elif defined(TC_HAVE_TR1_UNORDERED_MAP)
  template <typename Key, typename Hash = std::tr1::hash<Key> >
#elif defined(TC_HAVE_HASH_MAP)
  template <typename Key, typename Hash = std::hash<Key> >
#elif defined(TC_HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
  template <typename Key, typename Hash = std::hash<Key> >
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
  template <typename Key, typename Hash = __gnu_cxx::hash<Key> >
# else
  template <typename Key, typename Hash = Key >
# endif
#else
  template <typename Key, typename Hash = Key >
#endif
  struct HashSet
  {
#if   defined(TC_HAVE_UNORDERED_MAP)
    typedef std::unordered_set<Key, Hash> Type;
#elif defined(TC_HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_set<Key, Hash> Type;
#elif defined(TC_HAVE_HASH_MAP)
    typedef std::hash_set<Key, Hash> Type;
#elif defined(TC_HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_set<Key, Hash> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_set<Key, Hash> Type;
# else
    typedef std::set<K> Type;
# endif
#else
    typedef std::set<K> Type;
#endif
  };



#if   defined(TC_HAVE_UNORDERED_MAP)
  template <typename Key, typename Hash = std::hash<Key> >
#elif defined(TC_HAVE_TR1_UNORDERED_MAP)
  template <typename Key, typename Hash = std::tr1::hash<Key> >
#elif defined(TC_HAVE_HASH_MAP)
  template <typename Key, typename Hash = std::hash<Key> >
#elif defined(TC_HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
  template <typename Key, typename Hash = std::hash<Key> >
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
  template <typename Key, typename Hash = __gnu_cxx::hash<Key> >
# else
  template <typename Key, typename Hash = Key >
# endif
#else
  template <typename Key, typename Hash = Key >
#endif
  struct HashMultiSet
  {
#if   defined(TC_HAVE_UNORDERED_MAP)
    typedef std::unordered_multiset<Key, Hash> Type;
#elif defined(TC_HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_multiset<Key, Hash> Type;
#elif defined(TC_HAVE_HASH_MAP)
    typedef std::hash_multiset<Key, Hash> Type;
#elif defined(TC_HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_multiset<Key, Hash> Type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multiset<Key, Hash> Type;
# else
    typedef std::multiset<K> Type;
# endif
#else
    typedef std::multiset<K> Type;
#endif
  };


#endif // TC_HASHSET_H
