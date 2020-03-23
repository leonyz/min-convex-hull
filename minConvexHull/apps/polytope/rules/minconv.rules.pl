#input: matrix X
#output: matrix val X obtained by taking entrywise valuations of X
user_function val_matrix($) {
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

# input: square matrix and lists of previous pivot rows and columns
# output: lists of pivot rows and columns with next pivot appended
user_function next_pivot( $ , @ , @ ){
   my($N, @R, @C) = @_;
   my $infinity = new TropicalNumber<Min,Int>("inf");
   my $val_N = val_matrix($N);
   my $num_pivots = scalar @R;
   for (my $i = 0; $i < $num_pivots; ++$i) {
   	$val_N -> elem(@R[$i], @C[$i]) = $infinity;
   }
   my $current_min = $infinity;
   my $current_pivot_row = 0;
   my $current_pivot_col = 0;
   my $d = $val_N->rows();
   for (my $i = 0; $i<$d; ++$i) {
   	for (my $j = 0; $j<$d; ++$j) {
   		if ($val_N->elem($i,$j) < $current_min) {
   			$current_pivot_row = $i;
   			$current_pivot_col = $j;
   			$current_min = $val_N->elem($i,$j);
   		}
   	}
   }
   push @R, $current_pivot_row;
   push @C, $current_pivot_col;
   return @R, @C;
}