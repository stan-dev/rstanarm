// This file is part of RStanArm
// Copyright (C) 2017 Trustees of Columbia University
//
// RStan is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// RStan is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

/*
 * To register the functions implemented in C++, see
 * http://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
 *
 * But it seems not to work as it is supposed to be in that
 * they are still working if not registered.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>

static const R_CallMethodDef CallEntries[] = {
  {NULL, NULL, 0}
};

void attribute_visible R_init_rstanarm(DllInfo *dll) {
  // next line is necessary to avoid a NOTE from R CMD check
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE); // necessary for .onLoad() to work
}
