/* Copyright (c) 1997-2018
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

/** @file dense_linalg.h
    @brief Linear Algebra algorithms for dense vector and matrix types
 */

#ifndef MIN_CONVEX_HULL_NEW_INV
#define MIN_CONVEX_HULL_NEW_INV

#include "polymake/vector"
#include "polymake/Vector.h"
#include "polymake/GenericStruct.h"
#include "polymake/internal/linalg_exceptions.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"

namespace polymake {
namespace polytope {

using namespace pm;

/// matrix inversion
template <typename E>
std::enable_if_t<is_field<E>::value, Matrix<E>>
new_inv(Matrix<E> M)
{
   const int dim=M.rows();
   std::vector<int> row_index(dim);
   copy_range(entire(sequence(0, dim)), row_index.begin());
   Matrix<E> u=unit_matrix<E>(dim);

   for (int c=0; c<dim; ++c) {
      int r=c;
      while (is_zero(M(row_index[r],c))) {
         if (++r==dim) throw degenerate_matrix();
      }
      E *ppivot=&M(row_index[r],c);
      const E pivot=*ppivot;
      E *urow=&u(row_index[r],0);
      if (r!=c) std::swap(row_index[r],row_index[c]);
      if (!is_one(pivot)) {
         E *e=ppivot;
         for (int i=c+1; i<dim; ++i) (*++e)/=pivot;
         for (int i=0; i<=c; ++i) urow[row_index[i]]/=pivot;
      }
      for (r=0; r<dim; ++r) {
         if (r==c) continue;
         E *e2=&M(row_index[r],c);
         const E factor=*e2;
         if (!is_zero(factor)) {
            E *e=ppivot;
            for (int i=c+1; i<dim; ++i) (*++e2)-=(*++e)*factor;
            E *urow2=&u(row_index[r],0);
            for (int i=0; i<=c; ++i) urow2[row_index[i]]-=urow[row_index[i]]*factor;
         }
      }
   }
   return Matrix<E>(dim, dim, select(rows(u),row_index).begin());
}

template <typename TMatrix, typename E>
typename std::enable_if<is_field<E>::value, typename TMatrix::persistent_nonsymmetric_type>::type
new_inv(const GenericMatrix<TMatrix, E>& m)
{
   if (POLYMAKE_DEBUG || is_wary<TMatrix>()) {
      if (m.rows() != m.cols())
         throw std::runtime_error("inv - non-square matrix");
   }
   return new_inv(typename TMatrix::persistent_nonsymmetric_type(m));
}

template <typename TMatrix, typename E>
typename std::enable_if<!std::is_same<E, typename algebraic_traits<E>::field_type>::value,
                        typename GenericMatrix<TMatrix, typename algebraic_traits<E>::field_type>::persistent_nonsymmetric_type>::type
new_inv(const GenericMatrix<TMatrix, E>& m)
{
   if (POLYMAKE_DEBUG || is_wary<TMatrix>()) {
      if (m.rows() != m.cols())
         throw std::runtime_error("inv - non-square matrix");
   }
   return new_inv(typename GenericMatrix<TMatrix, typename algebraic_traits<E>::field_type>::persistent_nonsymmetric_type(m));
}


} // end namespace polymake
} // end namespace polytope


#endif // MIN_CONVEX_HULL_NEW_INV

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
