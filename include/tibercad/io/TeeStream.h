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
 * \file TeeStream.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _TEESTREAM_H_
#define _TEESTREAM_H_

#include <iostream>
#include <fstream>

template<typename Elem, typename Traits = std::char_traits<Elem> >
struct basic_TeeStream : std::basic_ostream<Elem,Traits>
{
   typedef std::basic_ostream<Elem,Traits> SuperType;

   basic_TeeStream(std::ostream& o1, std::ostream& o2)
      :  SuperType(o1.rdbuf()), o1_(o1), o2_(o2) { }

   basic_TeeStream& operator<<(SuperType& (*manip)(SuperType& ))
   {
      o1_ << manip;
      o2_ << manip;
      return *this;
   }

   template<typename T>
   basic_TeeStream& operator<<(const T& t)
   {
      o1_ << t;
      o2_ << t;
      return *this;
   }

private:
   std::ostream& o1_;
   std::ostream& o2_;
};


typedef basic_TeeStream<char> TeeStream;

#endif // _TEESTREAM_H_
