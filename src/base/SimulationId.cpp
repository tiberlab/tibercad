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
 * \file SimulationId.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/base/SimulationId.h"

#include <random>
#include <sstream>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>


std::string
SimulationId::_uuid;

std::time_t
SimulationId::_timestamp = 0;

uint64_t
SimulationId::_seed = 0;

bool
SimulationId::_seed_set = false;


void
SimulationId::initialize(void)
{
  boost::uuids::random_generator gen;
  _uuid = boost::uuids::to_string(gen());

  _timestamp = std::time(NULL);

  if (!_seed_set)
  {
    std::random_device rd;
    _seed = (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
  }
}


const std::string&
SimulationId::uuid(void)
{
  return _uuid;
}


std::time_t
SimulationId::timestamp(void)
{
  return _timestamp;
}


std::string
SimulationId::start_time(void)
{
  std::tm utc = *std::gmtime(&_timestamp);

  char datetime[32];
  std::strftime(datetime, sizeof(datetime), "%Y-%m-%d %H:%M:%SZ", &utc);

  return std::string(datetime);
}


uint64_t
SimulationId::seed(void)
{
  return _seed;
}


void
SimulationId::set_seed(uint64_t seed)
{
  _seed = seed;
  _seed_set = true;
}


std::string
SimulationId::run_id(void)
{
  std::tm utc = *std::gmtime(&_timestamp);

  char datetime[32];
  std::strftime(datetime, sizeof(datetime), "%Y%m%dT%H%M%SZ", &utc);

  std::ostringstream os;
  os << datetime << "-" << _uuid.substr(0, 8);
  return os.str();
}
