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
 * \file InitializerBase.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_INITIALIZERBASE_H
#define TC_INITIALIZERBASE_H

class TiberModelObject;

//! The base class for initializer functors
template <typename T>
class InitializerBase
{

  public:

    //! Destructor
    virtual ~InitializerBase(void) {};

    //! The operator to be overloaded
    virtual void operator()(T& val) = 0;


  protected:


};

#endif /* _INITIALIZERBASE_H_ */
