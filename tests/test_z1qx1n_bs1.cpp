/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/fortranio.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace ibmisc;
using namespace icebin::modele;
using namespace blitz;

// The fixture for testing class Foo.
class Z1qx1n_Bs1Test : public ::testing::Test {
protected:
    // You can do set-up work for each test here.
    Z1qx1n_Bs1Test() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~Z1qx1n_Bs1Test() {}

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() {
    }

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() {
        // remove("z1qx1n_bs14_common");
    }

//    // The mock bar library shaed by all tests
//    MockBar m_bar;
};


TEST_F(Z1qx1n_Bs1Test, read_inputs)
{
    auto topo_inputs(make_topo_inputs());

    NcIO ncio("topo_inputs.nc", 'w');

    // Need to read the files!!!
    EnvSearchPath locator("MODELE_FILE_PATH");
    read_raw(topo_inputs, locator);

    topo_inputs.bundle.ncio(ncio, "", "double");
    ncio.close();
}





int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
