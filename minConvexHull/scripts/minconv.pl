set_var_names<UniPolynomial<Rational,Rational>>(qw(t));
$t=monomials<Rational,Rational>(1);
$I=new Matrix<PuiseuxFraction<Min>>(diag(new Vector<PuiseuxFraction<Min,Rational,Rational>>(1,1,1,1)));

$a1=new PuiseuxFraction<Min>($t**(-1));
$a2=new PuiseuxFraction<Min>($t**(-1));
$b3=new PuiseuxFraction<Min>($t*(3*$t-2)/(2*$t-1));
$b4=new PuiseuxFraction<Min>(1);

print (($a1+$a2)*($b3+$b4) - (3+$a1*$a2*$b3*$b4));

# matrices act on the left; i.e., vectors are columns
$A=new Matrix<PuiseuxFraction<Min,Rational,Rational>>(4,4);
$A->elem(0,0)=1; $A->elem(1,1)=1; $A->elem(3,2)=1; $A->elem(2,3)=-1; $A->elem(3,3)=-1; $A->elem(0,3)=$a1; $A->elem(1,3)=$a2; 
print val_matrix($A);
@v = min_valuation_of_columns($A);
print @v;
@column_pivots = (1,3);
@merged_array = (\@v, \@column_pivots);
print min_entry_of_subarray(@merged_array);
print next_pivot($A, ());
print next_pivot($A, ([0,3]));

declare $B=new Matrix<PuiseuxFraction<Min,Rational,Rational>>(4,4);
$B->elem(0,0)=-1; $B->elem(0,1)=1; $B->elem(1,0)=-1; $B->elem(2,2)=1; $B->elem(3,3)=1; $B->elem(2,0)=$b3; $B->elem(3,0)=$b4;
@sigma = [0,1,2,3];
@sigma2= [1,3,2,0];
@sigma3 = [2,1,0,3];
print sa_basis_with_prescribed_permutation($A, $B, @sigma);
print sa_basis_with_prescribed_permutation($A, $B, @sigma2);
print sa_basis($A, $B);
print sa_basis_with_prescribed_permutation($A, $B, @sigma3);