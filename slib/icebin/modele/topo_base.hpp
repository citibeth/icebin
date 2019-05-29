#ifndef ICEBIN_MODELE_MAKE_TOPOO_HPP
#define ICEBIN_MODELE_MAKE_TOPOO_HPP

#include <vector>
#include <string>
#include <ibmisc/filesystem.hpp>
#include <ibmisc/bundle.hpp>
#include <icebin/modele/hntr.hpp>

/** Generate the basic TOPOO file, based on Gary Russell's old TOPO generator.
This is converted from FORTRAN, and uses 1-based indexing. */

namespace icebin {
namespace modele {

struct MakeTopoO {
    HntrSpec hspec;    // Describes the grid the variables use
    ibmisc::ArrayBundle<double,2> bundle;

    MakeTopoO(
        ibmisc::FileLocator const &files,
        std::vector<std::string> const &_varinputs);
};


}}
#endif     // guard
