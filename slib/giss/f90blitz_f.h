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

#define GET_SPEC_DOUBLE_1(arr, spec) get_spec_double_1(arr, lbound(arr,1), spec)
#define GET_SPEC_DOUBLE_2(arr, spec) get_spec_double_2(arr, lbound(arr,1),lbound(arr,2), spec)
#define GET_SPEC_DOUBLE_3(arr, spec) get_spec_double_3(arr, lbound(arr,1),lbound(arr,2),lbound(arr,3), spec)
#define GET_SPEC_INT_1(arr, spec) get_spec_int_1(arr, lbound(arr,1), spec)
#define GET_SPEC_INT_2(arr, spec) get_spec_int_2(arr, lbound(arr,1),lbound(arr,2), spec)
#define GET_SPEC_INT_3(arr, spec) get_spec_int_3(arr, lbound(arr,1),lbound(arr,2),lbound(arr,3), spec)
