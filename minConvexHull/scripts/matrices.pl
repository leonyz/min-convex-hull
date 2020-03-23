# Michael Joswig & Stephan Tillmann, August 2010
# producing a partial Cayley graph of a matrix group

use application "graph";

set_var_names<UniPolynomial<Rational,Int>>(qw(t));
declare $t=monomials<Rational,Int>(1);
declare $I=new Matrix<PuiseuxFraction<Min,Rational,Int>>(diag(new Vector<PuiseuxFraction<Min,Rational,Int>>(1,1,1,1)));

declare $a1=new PuiseuxFraction<Min,Rational,Int>($t**(-1));
declare $a2=new PuiseuxFraction<Min,Rational,Int>($t**(-1));
declare $b3=new PuiseuxFraction<Min,Rational,Int>($t*(3*$t-2)/(2*$t-1));
declare $b4=new PuiseuxFraction<Min,Rational,Int>(1);

print (($a1+$a2)*($b3+$b4) - (3+$a1*$a2*$b3*$b4));

# matrices act on the left; i.e., vectors are columns
declare $A=new Matrix<PuiseuxFraction<Min,Rational,Int>>(4,4);
$A->elem(0,0)=1; $A->elem(1,1)=1; $A->elem(3,2)=1; $A->elem(2,3)=-1; $A->elem(3,3)=-1; $A->elem(0,3)=$a1; $A->elem(1,3)=$a2; 
declare $B=new Matrix<PuiseuxFraction<Min,Rational,Int>>(4,4);
$B->elem(0,0)=-1; $B->elem(0,1)=1; $B->elem(1,0)=-1; $B->elem(2,2)=1; $B->elem(3,3)=1; $B->elem(2,0)=$b3; $B->elem(3,0)=$b4;

declare @gens_matrices=($A,$B);
declare @gens_labels=("A","B");
declare @Cells=();

# apply valuation map coefficientwise
sub val_matrix($) {
    my $X=shift(@_);
    my $d=$X->rows();
    my $e=$X->cols();
    my $val_X=new Matrix<TropicalNumber<Min,Int>>($d,$e);
    for (my $i=0; $i<$d; ++$i) {
	for (my $j=0; $j<$d; ++$j) {
	    $val_X->elem($i,$j)=$X->elem($i,$j)->val();
	}
    }
    return $val_X;
}

# check if given (square) matrix is lower triangular;
# this is the stabilizer of a flag in the projective space
sub lowerTriangular($) {
    my $X=shift(@_);
    my $d=$X->rows();
    die "$X not a square matrix" unless $X->cols()==$d;
    for (my $i=0; $i<$d; ++$i) {
	for (my $j=$i+1; $j<$d; ++$j) {
	    return false if $X->elem($i,$j)!=0;
	}
    }
    return true;
}

# check if given (square) matrix is in the affine flag stabilizer;
# cf. Abramenko & Brown, section 6.9
sub affineFlagStab($) {
    my $X=shift(@_);
    my $d=$X->rows();
    die "$X not a square matrix" unless $X->cols()==$d;
    my $val_X=val_matrix($X);
    my $trop_unit=new TropicalNumber<Min,Int>(0);
    my $one_as_trop=new TropicalNumber<Min,Int>(1);
    for (my $i=0; $i<$d; ++$i) {
	return false if $val_X->elem($i,$i)!=$trop_unit;
    }
    for (my $i=0; $i<$d; ++$i) {
	for (my $j=$i+1; $j<$d; ++$j) {
	    return false if $X->elem($i,$j)<=$trop_unit;
	    return false if $X->elem($j,$i)<=$one_as_trop;
	}
    }
    return true;
}

# generate up to 2k+1 nodes of the Cayley graph (minus identifications)
sub CayleyGraph($) {
    my $k=shift(@_); my $n=2*$k+1; # safe upper bound for number of nodes
    my $cnt=0;
    my %M=();
    $M{"$I"}=$cnt++;
    @Cells=($I);
    my @NodeLabels=("");
    my @queue=(0);
    my $Cayley=new props::Graph<Directed>($n);
    for (my $i=0; $i<$k; ++$i) {
	my $c=shift(@queue);
	for (my $j=0; $j<=$#gens_matrices; ++$j) {
	    my $X=$gens_matrices[$j]*$Cells[$c];
	    my $strX="$X";
	    if (!defined($M{$strX})) {
		$M{$strX}=$cnt;
		$Cells[$cnt]=$X;
		$NodeLabels[$cnt]=$gens_labels[$j].$NodeLabels[$c];
		push @queue, $cnt;
		++$cnt;
	    }
	    my $cc=$M{$strX};
	    $Cayley->edge($c,$cc);
	}
    }
    $Cayley->squeeze_isolated();
    my $CG=new Graph(ADJACENCY=>$Cayley, N_NODES=>$cnt, NODE_LABELS=>\@NodeLabels);
    return $CG;
}

# jet scale node colors by degree; red signals degree=4
sub NodeColors($) {
    my $G=shift(@_);
    my $n=$G->N_NODES;
    my @colors=("","green","yellow","orange","red");
    my @NC = map { $colors[$_] } @{$G->NODE_DEGREES};
    return \@NC;
}

declare $CG=CayleyGraph(200);
$CG->VISUAL(NodeColor=>NodeColors($CG));
