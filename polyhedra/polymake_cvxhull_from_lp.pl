#!/usr/bin/perl

# Polymake script to print the facet description of the convex hull of integer 
# points corresponding to feasible solutions to the given .lp file

use application "polytope";

my $argv_len = @ARGV;

if ($argv_len != 1)
{
    print("\nusage: polymake --script polymake_cvxhull_from_lp.pl [input_.lp_file_path]\n");
    print("\noptional: redirect std output with >[output_file_path]\n\n");
    exit;
}

my $f=lp2poly($ARGV[0]);
my $p = new Polytope<Rational>($f);

$p->LATTICE_POINTS;

my $s=new Polytope(POINTS=>$p->LATTICE_POINTS, COORDINATE_LABELS=>$p->COORDINATE_LABELS);

print_constraints($s);
