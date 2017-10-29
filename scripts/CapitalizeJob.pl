#!/usr/bin/perl
use strict;

my $RoboPath = "../Documentation/code/";
my $fNameIndex = $RoboPath."robo_namelist.html";
my $fNameHTML;
my $sNamelist;
my $fileHTML;
my $isOpen = 0;

my $verbose = 0; # Overall flag for verbosity

while (<>)
{
    my $line = $_;

# find line starting a new namelist:

    if (/^\s*&(\S+)\s*/)
    {
        if ($verbose) {print "Namelist: $1\n";};
        my $NamelistName = $1;

        open(file, $fNameIndex) || die "ERROR: can not open '$fNameIndex': $!.\n";

        while (<file>)
        {
            if ( />$NamelistName<\/a>/i )
            {
                if ( /<a href=.(\S+)\#\S+\s+\S+\s+>(\S+)<\/a>/ )
                {
                    $fNameHTML = $1;
                    $sNamelist = $2;
                }
                else
                {
                    die "ERROR: HTML file or Name not recognized.\n";
                }
                last;
            }
        }
        close(file);
        if ($verbose) {print "HTML: $fNameHTML\nName: $sNamelist\n";};

        my $fNameF = $fNameHTML;
        $fNameF =~ s/_f90.html/.f90/;
        $fNameF =~ s|\./|code/|;
        print "! file: $fNameF\n";

        $line =~ s/$sNamelist/$sNamelist/i; # Does the replacement!

        if ($isOpen)
        { 
            close($fileHTML);
            $isOpen = 0;
        }
        open($fileHTML, $RoboPath.$fNameHTML) || die "ERROR: can not open '$RoboPath.$fNameHTML': $!.\n";
        $isOpen = 1;
        
    }
    elsif (/^[!\s]\s*(\S+)\s*=/) 
    {
        if ($verbose) {print "Variable: |$1|\n";};
        my $VariableName = $1;

        seek($fileHTML, 0, 0); # goto beginning of file
        while (<$fileHTML>)
        {
            if ( /<li> <a \S+>$VariableName<\/a>/i )
            {
                if (/<li> <a href=\S+>(\S+)<\/a>/)
                {
                    my $NewName = $1;
                    if ($verbose) {print "$1\n";};
                    $line =~ s/$VariableName/$NewName/i; # Does the replacement!
                    last;
                }
            }
#            else
#            {
#                die "ERROR: Variable '$VariableName' not found.\n"
#            }
        }
    }


    print $line;
}
