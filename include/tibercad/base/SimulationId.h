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
 * \file SimulationId.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_SIMULATIONID_H
#define TC_SIMULATIONID_H

#include "tibercad/base/tiber_dll.h"

#include <string>
#include <cstdint>
#include <ctime>

//! Unique identifying information for one simulation run
/*!
 * Collects the pieces of information needed to tell apart and reproduce a
 * particular simulation run: a UUID minted for that run, the time it
 * started, the master random seed it used, and (via TiberCad) the exact
 * code revision it ran. This is meant to be attached to output data, so
 * that data sets can later be traced back to the run, code version and
 * seed that produced them.
 */
class SimulationId
{

  public:

    //! Mint a new run identifier
    /*!
     * Generates a fresh UUID and start timestamp for this run, and a
     * master random seed unless one was already given via set_seed().
     * Called once by Control::init(), at the start of each simulation.
     */
    static void initialize(void) TC_DLLOCAL;

    //! The UUID minted for this run, as a string
    static const std::string& uuid(void);

    //! The time this run started, in seconds since the epoch (UTC)
    static std::time_t timestamp(void);

    //! The time this run started, as a human-readable UTC string
    /*!
     * ISO 8601 extended format, e.g. "2026-11-16 09:04:00Z".
     */
    static std::string start_time(void);

    //! The master random seed for this run
    static uint64_t seed(void);

    //! Explicitly set the master random seed for this run
    /*!
     * Must be called before initialize() to have an effect; otherwise
     * initialize() picks a seed from a non-deterministic source.
     */
    static void set_seed(uint64_t seed) TC_DLLOCAL;

    //! A short, human-readable identifier for this run
    /*!
     * Combines the run's start timestamp (UTC, ISO 8601 basic format)
     * with the first 8 hex characters of its UUID into a compact,
     * sortable, filename-safe string, e.g. "20261113T154230Z-3f9a2b1c"
     * -- meant for quick visual
     * identification (log banners, output file names), where the full
     * uuid() would be needlessly long.
     */
    static std::string run_id(void);

  private:

    SimulationId(void) TC_DLLOCAL {};

    //! The UUID minted for this run
    static std::string _uuid;

    //! The time this run started, in seconds since the epoch (UTC)
    static std::time_t _timestamp;

    //! The master random seed for this run
    static uint64_t _seed;

    //! Whether \c _seed was set explicitly, rather than generated
    static bool _seed_set;
};

#endif // TC_SIMULATIONID_H
