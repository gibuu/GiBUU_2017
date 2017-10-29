#!/usr/bin/perl
use strict;

my %ProvidedBy;
my %FileNeeds;
my %FilesObj;

my $fNameObj;

my $verbose = 0; # Overall flag for verbosity

###
### STEP 1: Scan all input files, build up the lists/hashes
###

foreach my $fName (@ARGV)
{
  $fNameObj = $fName;
  $fNameObj =~ s/\.f90$/\.o/;
  $fNameObj =~ s/\.F90$/\.o/;
  $fNameObj =~ s/\.f$/\.o/;
  $fNameObj =~ s/\.F$/\.o/;

  $FilesObj{$fName} = $fNameObj;

  if ($verbose) {print "\nREADING FILE '$fName' ...  ($fNameObj)\n";};

  open(file, $fName) || die "ERROR: can not open file '$fName'! STOP\n";

  my @ThisFileProvides;
  my @ThisFileNeeds;

  while (<file>)
  {
    if (/^\s*module\s+(\S+)\s*/i)
    {
#      print "<-- module '$1' is provided in file '$fName'\n";
      push @ThisFileProvides, lc($1);
    };

    if ( /^\s*use\s+(\S+)\s*,\s*only\s*:\s*/i|| /^\s*use\s+(\S+)\s*/i)
    {
#      print "--> module '$1' is needed by file '$fName'\n";
      push @ThisFileNeeds, lc($1);
    }

  }
  close(file);

### create a list of unique entries of needed modules:

  my %hash   = map { $_, 1 } @ThisFileNeeds;
  my @ThisFileNeedsUnique = keys %hash;

### just print the stuff:

  if ($verbose)
  {
    print "Provides:\n";
    foreach my $field (@ThisFileProvides) {print "   $field\n";}

    print "Needs:\n";
    foreach my $field (@ThisFileNeedsUnique) {print "   $field\n";}
  }
### store it in the overall lists/hashes:

  foreach my $field (@ThisFileProvides) {$ProvidedBy{$field} = $fName;}
  $FileNeeds{$fName} = [ @ThisFileNeedsUnique ];

}


###
### STEP 2: Produce Output
###

if ($verbose) {print "\n\n FINAL OUTPUT:\n\n";};

foreach my $field ( keys %FileNeeds )
{
  print "$FilesObj{$field} : $field";
  foreach my $i ( 0 .. $#{ $FileNeeds{$field} } ) 
  {
    if ($field ne $ProvidedBy{$FileNeeds{$field}[$i]})
    {
      print " $FilesObj{$ProvidedBy{$FileNeeds{$field}[$i]}}";
    }
  }
  print "\n";

}





