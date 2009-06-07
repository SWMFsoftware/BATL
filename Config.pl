#!/usr/bin/perl -i
use strict;
our @Arguments       = @ARGV;
our $Code            = "BATL";
our $MakefileDefOrig = 'src/Makefile.def';
require "share/Scripts/Config.pl";

# Variables inherited from share/Scripts/Config.pl
our $ERROR;
our $WARNING;
our $Help;
our $Show;
our $ShowGridSize;
our $NewGridSize;

# Source code directory
my $Src         = 'src';

# Grid size variables
my $NameGridFile = "$Src/BATL_size.f90";
my $GridSize;
my ($nDim, $nDimAmr, $nI, $nJ, $nK);

&print_help if $Help;

# Read previous grid size
&get_settings;

# Read previous grid size, equation and user module
&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

if($ShowGridSize or $Show){
    print "Config.pl -g=$nDim,$nDimAmr,$nI";
    print ",$nJ" if $nDim > 1;
    print ",$nK" if $nDim > 2;
    print "\n";
}

exit;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameGridFile
    open(MODSIZE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<MODSIZE>){
        next if /^\s*!/; # skip commented out lines
	$nDim=$1         if /\bnDim\s*=\s*(\d+)/i;
	$nDimAmr=$1      if /\bnDimAmr\s*=\s*(\d+)/i;
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
        $nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
        $nK=$1           if /\bnK\s*=\s*(\d+)/i;
    }
    close MODSIZE;

    die "$ERROR could not read nDimAmr from $NameGridFile\n" 
	unless length($nDimAmr);
    die "$ERROR could not read nDim from $NameGridFile\n" 
	unless length($nDim);
    die "$ERROR could not read nI from $NameGridFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameGridFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameGridFile\n" unless length($nK);

    $GridSize  = "$nDim,$nDimAmr,$nI";
    $GridSize .= ",$nJ" if $nDim > 1;
    $GridSize .= ",$nK" if $nDim > 2;
}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize =~ /^\d+(,\d+){2,4}$/){
	($nDim,$nDimAmr,$nI,$nJ,$nK) = split(',', $GridSize);
	if($nDim == 1){
	    die "$ERROR for nDim=1 ".
		"-g=$GridSize should contain 3 positive integers\n"
		if $nJ > 0;
	    $nJ = 1; $nK = 1;
	}elsif($nDim == 2){
	    die "$ERROR for nDim=2 ".
		"-g=$GridSize should contain 4 positive integers\n"
		if $nJ < 1 or $nK > 0;
	    $nK = 1;
	}elsif($nDim == 3){
	    die "$ERROR for nDim=3 ".
		"-g=$GridSize should contain 5 integers\n"
		if $nJ < 1 or $nK < 1;
	}
    }elsif($GridSize){
	die "$ERROR ".
	    "-g=$GridSize should be 3 to 5 integers separated by commas\n";
    }

    # Check the grid size (to be set)
    die "$ERROR nDim=$nDim must be 1, 2 or 3\n" if $nDim < 1 or $nDim > 3;
    die "$ERROR nDimAmr=$nDimAmr must be 1 to nDim=$nDim\n" 
	if $nDimAmr < 1 or $nDimAmr > $nDim;
    die "$ERROR nI=$nI must be 2 or more\n" if $nI < 2;
    die "$ERROR nJ=$nJ must be 2 or more\n" if $nJ < 2 and $nDim > 1;
    die "$ERROR nK=$nK must be 2 or more\n" if $nK < 2 and $nDim > 2;

    die "$ERROR nI=$nI must be an even integer\n" if $nI%2!=0;
    die "$ERROR nJ=$nJ must be an even integer\n" if $nJ%2!=0 and $nDimAmr > 1;
    die "$ERROR nK=$nK must be an even integer\n" if $nK%2!=0 and $nDimAmr > 2;

    print "Writing new grid size $GridSize into $NameGridFile...\n";

    @ARGV = ($NameGridFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nDim\s*=[^0-9]*)(\d+)/$1$nDim/i;
	s/\b(nDimAmr\s*=[^0-9]*)(\d+)/$1$nDimAmr/i;
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	print;
    }

}

##############################################################################

sub print_help{

    print "
Additional options for BATL/Config.pl:

-g[=NDIM,NDIMAMR,NI[,NJ[,NK]]]
    If -g is used without a value, it shows grid size. 
    Otherwise set grid and AMR dimensionality and the AMR block size.
    NDIM=1,2 or 3 is the dimensionality of the problem. 
    NDIMAMR is the dimensionality of the AMR tree: 1..NDIM.
    NI, NJ and NK are the number of cells in a block in the I, J and K 
    directions, respectively. The number of cells is always 1 in the 
    ignored dimensions, and should not be set.

Examples for BATS/Config.pl:

Show grid size:

    Config.pl -g

Set 3D domain with 3D AMR and block size 8x8x8 cells:

    Config.pl -g=3,3,8,8,8

Set 2D domain with 1D AMR and block size 40x10 cells:

    Config.pl -g=2,1,40,10
\n";
    exit 0;
}

