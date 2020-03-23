$A = dense(unit_matrix<PuiseuxFraction<Min,Rational,Int>>(3));
print test1($A);
print test2($A, $A);
print new_inv($A);
print new_inv2($A);

set_var_names<UniPolynomial<Rational,Int>>(qw(t));
$t=monomials<Rational,Int>(1);
$a1=new PuiseuxFraction<Min,Rational,Int>($t**(-1));
$a2=new PuiseuxFraction<Min,Rational,Int>($t**(-1));
$b3=new PuiseuxFraction<Min,Rational,Int>($t*(3*$t-2)/(2*$t-1));
$b4=new PuiseuxFraction<Min,Rational,Int>(1);
$A=new Matrix<PuiseuxFraction<Min,Rational,Int>>(4,4);
$A->elem(0,0)=1; $A->elem(1,1)=1; $A->elem(3,2)=1; $A->elem(2,3)=-1; $A->elem(3,3)=-1; $A->elem(0,3)=$a1; $A->elem(1,3)=$a2; 
print $A;
$B=new Matrix<PuiseuxFraction<Min,Rational,Int>>(4,4);
$B->elem(0,0)=-1; $B->elem(0,1)=1; $B->elem(1,0)=-1; $B->elem(2,2)=1; $B->elem(3,3)=1; $B->elem(2,0)=$b3; $B->elem(3,0)=$b4;
print $B;

print new_inv2($B);
print new_inv2($B) * $B;
