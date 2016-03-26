/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/bind.hpp>
#include <netcdfcpp.h>
#include <giss/ncutil.hpp>
#include <galahad/zd11_c.hpp>

namespace galahad {

static void netcdf_write(zd11_c const *mat, NcFile *nc, std::string const &vname)
{
    auto rowVar = nc->get_var((vname + ".row").c_str());
    rowVar->put(mat->row, mat->ne);
    auto colVar = nc->get_var((vname + ".col").c_str());
    colVar->put(mat->col, mat->ne);
    auto valVar = nc->get_var((vname + ".val").c_str());
    valVar->put(mat->val, mat->ne);
}


boost::function<void()> zd11_c::netcdf_define(NcFile &nc, std::string const &vname)
{
    auto oneDim = giss::get_or_add_dim(nc, "one", 1);
    auto infoVar = nc.add_var((vname + ".info").c_str(), ncDouble, oneDim);
        infoVar->add_att("m", this->m);
        infoVar->add_att("n", this->n);

    auto neDim = nc.add_dim((vname + ".ne").c_str(), this->ne);

    nc.add_var((vname + ".row").c_str(), ncInt, neDim);
    nc.add_var((vname + ".col").c_str(), ncInt, neDim);
    nc.add_var((vname + ".val").c_str(), ncDouble, neDim);

    return boost::bind(&netcdf_write, this, &nc, vname);
}

}
