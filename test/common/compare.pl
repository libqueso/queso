#!/usr/bin/perl
#
#------------------------------------------------------------------------------
# Simple solution comparison tool for QUESO generated output (Matlab format).
#------------------------------------------------------------------------------

use warnings;
$two_cols = 0;

my $DEFAULT_TOL = 1e-8;

if (@ARGV >= 2) {
    $file1 = shift @ARGV;
    $file2 = shift @ARGV;

} else {
    print "\nUsage: compare [solfile1] [solfile2]  <tolerance>\n\n";
    exit 0;
}

if (@ARGV == 1) {
    $TOL = shift @ARGV;
} else {
    $TOL = $DEFAULT_TOL;
}


my $count1 = 0;
my $count2 = 0;


#------------------------------------------
# Parsing routine
# -> support for 1 and 2 column formats
#------------------------------------------

sub parse_queso_matlab_output {

    my ($infile) = @_; 		# INPUT: filename to parse

    my @sola     = ();		# Column1 of solution input
    my @solb     = ();		# Column2 of solution input (optional)

    open($IN1, "<$infile") || die "Cannot open $file1\n";

    $line = <$IN1>;		# skip first line

    while ($line = <$IN1>)
    {

	# Terminate at the end of the array definition

	if($line =~ m/\];/) {
	    last;
	}

	# 2nd line of output is special: defines the array and has 1st value (2 col format)

	if($line =~ m/\S+ = \[(\S+) (\d+)\b/) {
	    push(@sola,$1);
	    push(@solb,$2);
	    $two_cols = 1;
	}

	# 2nd line of output is special: defines the array and has 1st value (1 col format)

	elsif($line =~ m/\S+ = \[(\S+)\b/) {
	    push(@sola,$1);
	}

	# Normal lines: two column format...

	elsif($line =~ m/(\S+) (\S+)\b/ ) {
	    push(@sola,$1);
	    push(@solb,$2);
	    $two_cols = 1;
	}

	# Normal lines: one column format...

	elsif($line =~ m/(\S+)\b/ ) {
	    push(@sola,$1);
	}

    }

    close $IN1;

#    my $num_lines = @solb;
#    print "local readb = $num_lines\n";

#    my $num_lines = @sola;
#    print "local reada = $num_lines\n";

    return (\@sola,\@solb);

}

#-------------------------------------------------------------
# Call routine to Store groovy contents of each solution file
#
# Note: slight bit of perl trickery here to get two array
#       references back from the parsing routine
#-------------------------------------------------------------

my ($sol1a_ref,$sol1b_ref) = parse_queso_matlab_output ( $file1 );
my ($sol2a_ref,$sol2b_ref) = parse_queso_matlab_output ( $file2 );

my @sol1a = @$sol1a_ref;
my @sol1b = @$sol1b_ref;

my @sol2a = @$sol2a_ref;
my @sol2b = @$sol2b_ref;

$num_sol1 = @sol1a;
$num_sol2 = @sol2a;

print "\nQUESO: Comparing $num_sol1 values using a tolerance of: $TOL\n";
print " --> ($file1 <-> $file2)\n";

my $success=1;

# Make sure we parsed some actual data.

if ($num_sol1 <= 0 ) {
    print "Error: No solution data available in $file1 (or format changed)\n";
    print "\nFAILED\n";
    exit -1;
}

# Make sure we have the same number of data points from both solutions.

if($num_sol1 != $num_sol2) {
    print "Error: # of solution variables do not agree ($num_sol1,$num_sol2)\n";
    print "\nFAILED\n";
    exit -1;
}

# Numerically diff the two solutions (each column individually).

for($count = 0;$count < $num_sol1; $count++) {

    if( $sol1a[$count] != $sol2a[$count] ) {
	my $delta = abs(($sol1a[$count] - $sol2a[$count]) / $sol2a[$count]);
	if ( $delta <= $TOL ) {
	    print " --> Within tol. ($sol1a[$count] - $sol2a[$count]) / $sol2a[$count] -> [delta = $delta] -> " .
		"Col 1 -> Index = $count\n";
	} else {
	    print " --> Different   ($sol1a[$count] - $sol2a[$count]) / $sol2a[$count] -> [delta = $delta] -> " .
		"Col 1 -> Index =$count\n";
	    $success=0;
	}
    }

    if($two_cols == 1) {
	if( $sol1b[$count] != $sol2b[$count] ) {

	    my $delta = abs(($sol1b[$count] - $sol2b[$count]) / $sol2b[$count]);
	    if ( $delta <= $TOL ) {
		print " --> Within tol. ($sol1b[$count] - $sol2b[$count]) / $sol2b[$count] -> [delta = $delta] -> " .
		    "Col 2 -> Index = $count\n";
	    } else {
		print " --> Different   ($sol1b[$count] - $sol2b[$count]) / $sol2b[$count] -> [delta = $delta] -> " .
		    "Col 2 -> Index = $count\n";
		$success=0;
	    }
	}
    }

}

if ($success == 1) {
    print " --> PASSED\n";
    exit 0;
} else {
    print " --> FAILED\n";
    exit -1;
}
