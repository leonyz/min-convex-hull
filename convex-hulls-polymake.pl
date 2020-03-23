#-------------#
# Example 2.20
#-------------#
$M = new Matrix<TropicalNumber<Min,Rational>>(4,3);
$M->elem(0,0)=0;
$M->elem(0,1)=0;
$M->elem(1,0)=0;
$M->elem(2,0)=0;
$M->elem(3,0)=0;
$M->elem(1,1)=1;
$M->elem(2,1)=2;
$M->elem(3,1)=3;
$M->elem(0,2)=-1;
$M->elem(1,2)=-2;
$M->elem(2,2)=1;
$M->elem(3,2)=-1;
print $M;
# 0 0 -1
# 0 1 -2
# 0 2 1
# 0 3 -1
$P = new tropical::Polytope<Min>(POINTS=>$M);
$P->VISUAL;






#--------------#
# Example 4.8
#--------------#
$t=monomials<Rational,Rational>(1);
$M1 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1,1,1],[1,$t,$t**2],[1,$t**(-2),$t]]);
$M2 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1,1,1],[$t,$t**2,$t**3],[$t**(-2),$t,$t**5]]);
$M3 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1,1,1],[$t**2,$t**3,$t**4],[$t,$t**5,$t**8]]);
$M4 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1,1,1],[$t**3,$t**4,$t**5],[$t**5,$t**8,$t**12]]);
@list_of_matrices=($M1,$M2,$M3,$M4);
print @list_of_matrices;
# (1) (1) (1)
# (1) (t) (t^2)
# (1) (t^-2) (t)
# (1) (1) (1)
# (t) (t^2) (t^3)
# (t^-2) (t) (t^5)
# (1) (1) (1)
# (t^2) (t^3) (t^4)
# (t) (t^5) (t^8)
# (1) (1) (1)
# (t^3) (t^4) (t^5)
# (t^5) (t^8) (t^12)
$M = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1,1,1,1,1,1],[1,$t,$t**2,$t**3,$t**4,$t**5],[1,$t**(-2),$t,$t**5,$t**8,$t**12]]);
print $M;
# (1) (1) (1) (1) (1) (1)
# (1) (t) (t^2) (t^3) (t^4) (t^5)
# (1) (t^-2) (t) (t^5) (t^8) (t^12)
$P1 = tropical_row_sum_of_matrix(val_matrix(inv($M1)*$M));
$P2 = tropical_row_sum_of_matrix(val_matrix(inv($M2)*$M));
$P3 = tropical_row_sum_of_matrix(val_matrix(inv($M3)*$M));
$P4 = tropical_row_sum_of_matrix(val_matrix(inv($M4)*$M));
@list_of_tropical_points = ($P1,$P2,$P3,$P4);
$tropical_points_matrix = matrix_from_rows(@list_of_tropical_points);
$tropical_points_matrix = remove_duplicate_columns_tropical_matrix($tropical_points_matrix);
$tropical_points_matrix = transpose($tropical_points_matrix);
print $tropical_points_matrix;
# 0 -2 -3 -6
# 0 0 -4 -8
# 0 0 0 -5
# 0 0 0 0
$P = new tropical::Polytope<Min>(POINTS=>$tropical_points_matrix);
$P->VISUAL;
$MPrime = enveloping_membrane(@list_of_matrices);
$Q1 = tropical_row_sum_of_matrix(val_matrix(inv($M1)*$MPrime));
$Q2 = tropical_row_sum_of_matrix(val_matrix(inv($M2)*$MPrime));
$Q3 = tropical_row_sum_of_matrix(val_matrix(inv($M3)*$MPrime));
$Q4 = tropical_row_sum_of_matrix(val_matrix(inv($M4)*$MPrime));
@Q_list_of_tropical_points = ($Q1, $Q2, $Q3, $Q4);
$Q_tropical_points_matrix = matrix_from_rows(@Q_list_of_tropical_points);
$Q_tropical_points_matrix = remove_duplicate_columns_tropical_matrix($Q_tropical_points_matrix);
$Q_tropical_points_matrix = transpose($Q_tropical_points_matrix);
print $Q_tropical_points_matrix;
# 0 -2 -3 -6
# 0 0 -4 -8
# 0 0 0 -5
# 0 0 -1 -5
# 0 -2 -2 -7
# 0 0 0 0
# 0 0 0 -4
# 0 0 -1 -1
# 0 -2 -3 -3
$Q = new tropical::Polytope<Min>(POINTS=>$Q_tropical_points_matrix);
$Q->VISUAL;






#--------------#
# Example 4.9
#--------------#
$t=monomials<Rational,Rational>(1);
$M1 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1, 0, 0, 0],[0, 1, 0, 0],[0,0,1,0],[0,0,0,1]]);
$M2 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1, 0, 0, 0],[0, $t**3, 0, 0],[0,0,$t**(-2),0],[0,0,0,$t**(-2)]]);
$M3 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[$t**(-3)-$t**2, 1-$t**2, -$t**(-2)+1, $t**(-2)-$t],[$t**2-$t**3,-$t**(-3)+$t,1-$t,0],[0,-1+$t,$t**(-3)-$t**3,$t**(-3)-1],[-$t+$t**2,-$t**(-1)+1, 0, -$t**(-1)+$t**2]]);
$M4 = new Matrix<PuiseuxFraction<Min,Rational,Rational>>([[1-$t**3, $t**(-1)-1, 1-$t**2, 1-$t**3],[$t**(-3)-1, 1-$t, 1-$t**2, 1-$t**3],[-$t**(-3)+$t, -$t**(-2)+1, -$t**(-3)+$t**(-1), -1+$t],[$t**(-3)-$t**(-2),-$t**(-1)+1,-$t**(-1)+1,$t**(-1)-1]]);
@list_of_matrices=($M1,$M2,$M3,$M4);
print @list_of_matrices;
# (1) (0) (0) (0)
# (0) (1) (0) (0)
# (0) (0) (1) (0)
# (0) (0) (0) (1)
# (1) (0) (0) (0)
# (0) (t^3) (0) (0)
# (0) (0) (t^-2) (0)
# (0) (0) (0) (t^-2)
# (t^-3 - t^2) (1 - t^2) (- t^-2 + 1) (t^-2 - t)
# (t^2 - t^3) (- t^-3 + t) (1 - t) (0)
# (0) (- 1 + t) (t^-3 - t^3) (t^-3 - 1)
# (- t + t^2) (- t^-1 + 1) (0) (- t^-1 + t^2)
# (1 - t^3) (t^-1 - 1) (1 - t^2) (1 - t^3)
# (t^-3 - 1) (1 - t) (1 - t^2) (1 - t^3)
# (- t^-3 + t) (- t^-2 + 1) (- t^-3 + t^-1) (- 1 + t)
# (t^-3 - t^-2) (- t^-1 + 1) (- t^-1 + 1) (t^-1 - 1)
$P = convex_hull_with_visual(@list_of_matrices);
$P->VISUAL;







#--------------#
# Example 5.3
#--------------#

$spanning_points_matrix = new Matrix<TropicalNumber<Min,Rational>>([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],[7, 12, 8, 12, 10, 12, 7, 13, 7, 18, 7, 13, 12, 18, 12, 18, 13, 15, 7, 
15, 12, 19, 13, 19],[20, 20, 16, 15, -8, 20, 20, 16, 20, 15, 20, 16, 20, 15, 20, 15, 16, 
-8, 12, -8, -1, -8, 16, -8]]);
$spanning_points_matrix = transpose($spanning_points_matrix);
print $spanning_points_matrix;
# 0 7 20
# 0 12 20
# 0 8 16
# 0 12 15
# 0 10 -8
# 0 12 20
# 0 7 20
# 0 13 16
# 0 7 20
# 0 18 15
# 0 7 20
# 0 13 16
# 0 12 20
# 0 18 15
# 0 12 20
# 0 18 15
# 0 13 16
# 0 15 -8
# 0 7 12
# 0 15 -8
# 0 12 -1
# 0 19 -8
# 0 13 16
# 0 19 -8
$P = new tropical::Polytope<Min>(POINTS=>$spanning_points_matrix);
$P->VISUAL;
$smaller_points_matrix = new Matrix<TropicalNumber<Min,Rational>>([[0,0,0,0,0],[7,12,13,18,19],[20,20,16,15,-8]]);
$smaller_points_matrix = transpose($smaller_points_matrix);
print $smaller_points_matrix;
$Q = new tropical::Polytope<Min>(POINTS=>$smaller_points_matrix);
$Q->VISUAL;

