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
my ($nDim, $nI, $nJ, $nK, $MaxBlock);

&print_help if $Help;

# Read previous grid size
&get_settings;

# Read previous grid size, equation and user module
&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

print "Config.pl -g=$nDim,$nI,$nJ,$nK,$MaxBlock\n" 
    if $ShowGridSize or $Show;

exit;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameGridFile
    open(MODSIZE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<MODSIZE>){
        next if /^\s*!/; # skip commented out lines
	$nDim=$1         if /\bnDimTree\s*=\s*(\d+)/i;
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
        $nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
        $nK=$1           if /\bnK\s*=\s*(\d+)/i;
        $MaxBlock=$1     if /\bMaxBlock\s*=\s*(\d+)/i;
    }
    close MODSIZE;

    die "$ERROR could not read nDimTree from $NameGridFile\n" 
	unless length($nDim);
    die "$ERROR could not read nI from $NameGridFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameGridFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameGridFile\n" unless length($nK);
    die "$ERROR could not read MaxBlock from $NameGridFile\n"
        unless length($MaxBlock);

    $GridSize = "$nDim,$nI,$nJ,$nK,$MaxBlock";
}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+,\d+,\d+,\d+$/){
	($nDim,$nI,$nJ,$nK,$MaxBlock)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be 5 integers separated with commas\n";
    }

    # Check the grid size (to be set)
    die "$ERROR nDim=$nDim must be 1, 2 or 3\n" if $nDim < 1 or $nDim > 3;
    die "$ERROR nI=$nI must be 2 or more\n" if $nI < 2;
    die "$ERROR nJ=$nJ must be 1 or more\n" if $nJ < 1 and $nDim < 2;
    die "$ERROR nK=$nK must be 1 or more\n" if $nK < 1 and $nDim < 3;
    die "$ERROR nJ=$nJ must be 2 or more\n" if $nJ < 2 and $nDim > 1;
    die "$ERROR nK=$nK must be 2 or more\n" if $nK < 2 and $nDim > 2;

    die "$ERROR nI=$nI must be an even integer\n" if $nI%2!=0;
    die "$ERROR nJ=$nJ must be an even integer\n" if $nJ%2!=0 and $nDim > 1;
    die "$ERROR nK=$nK must be an even integer\n" if $nK%2!=0 and $nDim > 2;

    die "$ERROR MaxBlock=$MaxBlock must be a positive integer\n" 
	if $MaxBlock<1 or $MaxBlock != int($MaxBlock);

    print "Writing new grid size $GridSize into $NameGridFile...\n";

    @ARGV = ($NameGridFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nDimTree\s*=[^0-9]*)(\d+)/$1$nDim/i;
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	s/\b(MaxBlock\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	print;
    }

}

##############################################################################

sub print_help{

    print "
Additional options for BATL/Config.pl:

-g=NDIM,NI,NJ,NK,MAXBLK
                Set grid size. NDIM=1,2 or 3 is the dimensionality of the 
                AMR tree. NI, NJ and NK are the number of cells in a block 
                in the I, J and K directions, respectively. 
                MAXBLK is the maximum number of blocks per processor.

Examples for BATS/Config.pl:

Set the AMR tree to 3D, block size to 8x8x8, number of blocks to 400:

    Config.pl -g=3,8,8,8,400

Set the AMR tree to 2D, block size to 40x10x1, number of blocks to 1000:

    Config.pl -g=2,40,10,1,1000
\n";
    exit 0;
}

