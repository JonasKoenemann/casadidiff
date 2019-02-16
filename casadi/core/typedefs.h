/*
 *
 *    Copyright (C) 2019 Jonas Koenemann
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with this program;
 *    if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_TYPEDEFS_H
#define CASADI_TYPEDEFS_H

namespace casadi
{
  class SXElem;
  typedef Matrix<SXElem> SX;
  typedef std::vector<SX> SXVector;
  typedef std::initializer_list<SX> SXIList;
  typedef std::vector<SXVector> SXVectorVector;
  typedef std::map<std::string, SX> SXDict;

  typedef std::vector<MX> MXVector;
  typedef std::initializer_list<MX> MXIList;
  typedef std::vector<MXVector> MXVectorVector;
  typedef std::map<std::string, MX> MXDict;

  typedef Matrix<double> DM;
  typedef std::vector<DM> DMVector;
  typedef std::vector<DMVector> DMVectorVector;
  typedef std::map<std::string, DM> DMDict;

  typedef std::map<std::string, Sparsity> SpDict;

}

#endif // CASADI_TYPEDEFS_H
