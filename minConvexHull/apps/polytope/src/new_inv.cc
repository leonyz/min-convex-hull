#include "polymake/polytope/new_inv.h"
#include "polymake/client.h"
#include "polymake/PuiseuxFraction.h"

namespace polymake {
namespace polytope {

using namespace pm;

// This just introduces an abbreviation so we can write 'pf' instead of
// 'PuiseuxFraction<Min,Rational,int>' in the following.
typedef PuiseuxFraction<Min,Rational,int> pf;

/// matrix inversion
Matrix<pf> new_inv2(Matrix<pf> M)
{
   const int dim=M.rows();
   std::vector<int> row_index(dim);
   copy_range(entire(sequence(0, dim)), row_index.begin());
   Matrix<pf> u=unit_matrix<pf>(dim);

   for (int c=0; c<dim; ++c) {
      int r=c;
      while (is_zero(M(row_index[r],c))) {
         if (++r==dim) throw degenerate_matrix();
      }
      pf *ppivot=&M(row_index[r],c);
      const pf pivot=*ppivot;
      pf *urow=&u(row_index[r],0);
      if (r!=c) std::swap(row_index[r],row_index[c]);
      if (!is_one(pivot)) {
         pf *e=ppivot;
         for (int i=c+1; i<dim; ++i) (*++e)/=pivot;
         for (int i=0; i<=c; ++i) urow[row_index[i]]/=pivot;
      }
      for (r=0; r<dim; ++r) {
         if (r==c) continue;
         pf *e2=&M(row_index[r],c);
         const pf factor=*e2;
         if (!is_zero(factor)) {
            pf *e=ppivot;
            for (int i=c+1; i<dim; ++i) (*++e2)-=(*++e)*factor;
            pf *urow2=&u(row_index[r],0);
            for (int i=0; i<=c; ++i) urow2[row_index[i]]-=urow[row_index[i]]*factor;
         }
      }
   }
   return Matrix<pf>(dim, dim, select(rows(u),row_index).begin());
}

Function4perl(&new_inv2, "new_inv2( $ )");

} // end namespace polytope
} // end namespace polymake

