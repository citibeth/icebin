#pragma once

#include <memory>
#include <glint2/Grid.hpp>

namespace glint2 {

/** @param grid2 Put in an RTree */
extern std::unique_ptr<Grid> compute_exchange_grid(Grid &grid1, Grid &grid2);

};	// namespace glint2
