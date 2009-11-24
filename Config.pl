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
our $AmrRatio;
our $NewAmrRatio;

our %Remaining; # Unprocessed arguments


# Source code directory
my $Src         = 'src';

# Grid size variables
my $NameGridFile = "$Src/BATL_size.f90";
my $GridSize;
my ($nI, $nJ, $nK, $iRatio, $jRatio, $kRatio);

&print_help if $Help;

# Read previous grid size
&get_settings;

foreach (@Arguments){
    if(/^-r=(.*)$/)     {$NewAmrRatio=$1; next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

# Set new grid size and AMR dimensions
&set_grid_size if ($NewGridSize and $NewGridSize ne $GridSize)
    or            ($NewAmrRatio and $NewAmrRatio ne $AmrRatio);

# Show grid size and AMR dimensions
if($ShowGridSize or $Show){
    print "Config.pl -g=$nI,$nJ,$nK -r=$iRatio,$jRatio,$kRatio\n";
}

exit;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameGridFile
    open(MODSIZE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<MODSIZE>){
        next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
        $nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
        $nK=$1           if /\bnK\s*=\s*(\d+)/i;
	$iRatio=$1       if /\biRatio\s*=\s*min\(\s*(\d)/;
	$jRatio=$1       if /\bjRatio\s*=\s*min\(\s*(\d)/;
	$kRatio=$1       if /\bkRatio\s*=\s*min\(\s*(\d)/;
    }
    close MODSIZE;

    die "$ERROR could not read nI from $NameGridFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameGridFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameGridFile\n" unless length($nK);

    die "$ERROR could not read iRatio from $NameGridFile\n" 
	unless length($iRatio);
    die "$ERROR could not read jRatio from $NameGridFile\n" 
	unless length($jRatio);
    die "$ERROR could not read kRatio from $NameGridFile\n" 
	unless length($kRatio);

    $GridSize  = "$nI,$nJ,$nK";
    $AmrRatio  = "$iRatio,$jRatio,$kRatio";
}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

    if($GridSize =~ /^[1-9]\d*,[1-9]\d*,[1-9]\d*$/){
	($nI,$nJ,$nK) = split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR ".
	    "-g=$GridSize should be 3 positive integers separated by commas\n";
    }

    $AmrRatio = $NewAmrRatio if $NewAmrRatio;

    if($AmrRatio =~ /^[12],[12],[12]$/){
	($iRatio,$jRatio,$kRatio) = split(',', $AmrRatio);
    }elsif($GridSize){
	die "$ERROR ".
	    "-r=$AmrRatio should be 3 integers = 1 or 2 separated by commas\n";
    }

    # Check the grid size (to be set)
    die "$ERROR nK=$nK must be 1 if nJ is 1\n" 
	if $nJ == 1 and $nK > 1;
    die "$ERROR nI=$nI must be an even integer >= 4 if iRatio=2\n" 
	if $iRatio==2 and ($nI<4 or $nI%2!=0);
    die "$ERROR nJ=$nJ must be 1 or an even integer >= 4 if jRatio=2\n" 
	if $jRatio==2 and $nJ>1 and ($nJ==2 or $nJ%2!=0);
    die "$ERROR nK=$nK must be 1 or an even integer >= 4 if kRatio=2\n" 
	if $kRatio==2 and $nK>1 and ($nK==2 or $nK%2!=0);;

    print "Writing new grid size $GridSize and AMR ratio $AmrRatio ".
	"into $NameGridFile...\n";

    @ARGV = ($NameGridFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines

	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	s/\b(iRatio\s*=[^0-9]*)(\d+)/$1$iRatio/i;
	s/\b(jRatio\s*=[^0-9]*)(\d+)/$1$jRatio/i;
	s/\b(kRatio\s*=[^0-9]*)(\d+)/$1$kRatio/i;
	print;
    }

}

##############################################################################

sub print_help{

    print "
Additional options for BATL/Config.pl:

-g=NI,NJ,NK
    If -g is used without a value, it shows grid size. 
    Otherwise set grid dimensionality and the grid block size.
    NI, NJ and NK are the number of cells in a block in the I, J and K 
    directions, respectively. If nK=1 the last dimension is ignored: 2D grid.
    If nJ=1 and nK=1 then the last two dimensions are ignored: 1D grid.

-r=IRATIO,JRATIO,KRATIO
    Set the AMR ratio for each (non-ignored) dimensions. The value can be
    1 (for no adaptation), or 2 for adaptation in the given direction.

Examples for BATS/Config.pl:

Show grid size:

    Config.pl -g

Set 3D domain with 3D AMR and block size 8x8x8 cells:

    Config.pl -g=8,8,8 -r=2,2,2

Set block size 40x10x1 cells (2D grid) with AMR in the first dimension only:

    Config.pl -g=40,10,1 -r=2,1,1
\n";
    exit 0;
}

