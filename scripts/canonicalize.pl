#!/usr/bin/perl

# This script intends to move a source file into some canonical form.
#
# Some stuff it does:
# * for lines like ' !*** ' the stars extend up to column 79
# * Spelling of keywords standardized: 'if', 'end if'
# * a whitespace after if: 'if ('
# * ... many other stuff
#
# Some stuff it does not yet:
# * capitalized 'intent(in,out,inout)', 'integer', 'real', ...,
#   'pointer'
# * whitespace between function name and '('
# * no space in front of ',' [i.e. switch ' ,' to ', ']
# * only one space after 'call'
# * ... many other stuff
#


use strict;
use warnings;
#use re "debug";
use File::Copy;

my $verbose = 1; # Overall flag for verbosity
my $doReplace = 0; # Overall flag to replace the file by the new version

foreach my $fName (@ARGV)
{
    my $fNameNew = $fName.".canonical";

    if ($verbose) {print "\nREADING FILE '$fName' ...  ($fNameNew)\n";};
    open(fileI, $fName) || die "ERROR: can not open file '$fName'! STOP\n";

    open(fileO, '>', $fNameNew) || die "ERROR: can not open file '$fNameNew'! STOP\n";

    while(<fileI>)
    {
        my $L = $_;
        chomp($L);

        $L = CorrectStars($L);
        ($L,my @S) = RemoveStrings($L);
        ($L,my $C) = StripComment($L);

        # HERE YOU ADD YOUR STUFF:

        $L = Keywords($L);
        $L = IfSpace($L);
        $L = ThenSpace($L);
        $L = CaseSpace($L);
        $L = PrintWrite($L);
        $L = WriteSpace($L);
        $L = OpenClose($L);
        $L = UseOnly($L);

        $L = RenameVars($L);

        # Now the line is printed to the file:

        $L = RecreateLine($L,$C,\@S);
        $L =~ s/\s+$//; # remove trailing whitespaces
        print fileO "$L\n";
    }

    close(fileI);
    close(fileO);

    if ($doReplace)
    {
        move($fNameNew,$fName)  or die "Move failed: $!";
    };
}

#
# HERE YOU ADD YOUR STUFF:
#

# This guarantes, that for all lines like ' !*** ' the stars extend up
# to column 79 (80-1 for safety)
#
# Do this also for: '!===','!---','!###','!+++',
sub CorrectStars
{
    my $L = shift;

    if ($L =~ /^([\s]*)![*]+[\s]*$/ )
    {
        my $length = 78 - length($1);
        $L = $1."!".('*' x $length);
    }

    if ($L =~ /^([\s]*)![=]+[\s]*$/ )
    {
        my $length = 78 - length($1);
        $L = $1."!".('=' x $length);
    }

    if ($L =~ /^([\s]*)![#]+[\s]*$/ )
    {
        my $length = 78 - length($1);
        $L = $1."!".('#' x $length);
    }

    if ($L =~ /^([\s]*)![+]+[\s]*$/ )
    {
        my $length = 78 - length($1);
        $L = $1."!".('+' x $length);
    }

    if ($L =~ /^([\s]*)![-]+[\s]*$/ )
    {
        my $length = 78 - length($1);
        $L = $1."!".('-' x $length);
    }


    return($L);
}


# Find unique capitalization and spelling of keywords
#
# * keywords are lower case
# * space in 'end if' etc
#
sub Keywords
{
    my $L = shift;

    $L =~ s/([\s]+)else[\s]*if([\s]*)\(/$1else if$2\(/i; # "else if("

    $L =~ s/([\s]+)select[\s]*case([\s]*)\(/$1select case$2\(/i; # "select case"
    $L =~ s/([\s]+)case[\s]+default([\s]*)$/$1case default$2/i; # "case default"

    #
    # do not include "if" and "write" here, since we check these explicitely
    # against "if(" etc.
    #
    my @reserved=
        ("do", "while", "else", "then",
         "optional", "public", "private" );

    foreach (@reserved)
    {
        my $w = $_;

        $L =~ s/([\s]+)$w$/$1$w/i;
        $L =~ s/([\s]+)$w([\s]+)/$1$w$2/i;
    }

    #
    # words in combination with "("
    #
    my @reservedBr=
        ("if", "case", "write", "allocate", "deallocate",
         "open", "close" );

    foreach (@reservedBr)
    {
        my $w = $_;

        $L =~ s/([\s]+)$w([\s]*)\(/$1$w$2\(/i;
    }

    #
    # words, which can have ".(=+-," left of it and have "(" right of it:
    #
    my @reservedBr2=
        ("present", "allocated", "associated",
         "abs", "lbound", "ubound");

    foreach (@reservedBr2)
    {
        my $w = $_;

        $L =~ s/([\s\.\(,=+-]+)$w([\s]*)\(/$1$w$2\(/i;
    }

    #
    # words in combination with "end"
    #
    foreach ("do", "if", "select", "subroutine", "function", "module")
    {
        my $w = $_;

        $L =~ s/([\s]+)end[\s]*$w$/$1end $w/i;
        $L =~ s/([\s]+)end[\s]*$w([\s]+)/$1end $w$2/i;
    }

    return($L);
}

# Space between 'if' and '('
sub IfSpace
{
    my $L = shift;

    $L =~ s/([\s]+)(if)([\s]*)\(/$1$2 \(/i; # this 'eats' $3

    return($L);
}

# Space between ')' and 'then'
sub ThenSpace
{
    my $L = shift;

    $L =~ s/\)([\s]*)(then)/\) $2/i; # this 'eats' $1

    return($L);
}

# Space between 'case' and '('
sub CaseSpace
{
    my $L = shift;

    $L =~ s/([\s]+)(case)([\s]*)\(/$1$2 \(/i; # this 'eats' $3

    return($L);
}

# No space between 'write' and '('
### A space between 'write(...)' and further arguments # does not work yet!
# A space between 'write(*,*)' and further arguments
sub WriteSpace
{
    my $L = shift;

    $L =~ s/([\s]+)(write)([\s]*)\(/$1$2\(/i;

    # greediness is not the same as getting matching '(' and ')' !!!
    #    $L =~ s/([\s]+)(write)[\s]*\((.*)\)/$1$2\($3\) /i;

    # simple solution for '(*,*)':
    # (use [ \t] instead of [\s], since this would also eat a \n)
    $L =~ s/([\s]+)(write)[\s]*\(\*,\*\)[ \t]*/$1$2\(*,*\) /i;

    return($L);
}

# use lines have the format:
#  use XXX, only: ...
sub UseOnly
{
    my $L = shift;

    $L =~ s/([\s]+)(use)[\s]+([\S]+)[\s]*,[\s]*only[\s]*:/$1$2 $3, only:/i;

    return($L);
}

# open and close as : 'open(...', 'close(...'
sub OpenClose
{
    my $L = shift;

    # do not change the capitalization:
    #    $L =~ s/([\s]+)(open)([\s]*)\(/$1$2\(/i; # this 'eats' $3
    #    $L =~ s/([\s]+)(close)([\s]*)\(/$1$2\(/i; # this 'eats' $3

    # change the capitalization:
    $L =~ s/([\s]+)(open)([\s]*)\(/$1open\(/i; # this 'eats' $3
    $L =~ s/([\s]+)(close)([\s]*)\(/$1close\(/i; # this 'eats' $3

    return($L);
}

# replace "print *," by "write(*,*)"
sub PrintWrite
{
    my $L = shift;

    $L =~ s/([\s]+)(print)([\s]*)\*,/$1write(*,*)/i;

    return($L);
}

# replace a string with 'XXX', return a list of all extracted strings
sub RemoveStrings
{
    my $L = shift;
    my @strings;

    while ($L =~ /([\"'])(?:\\\1|.)*?\1/)
    {
        $L =~ s/([\"'])(?:\\\1|.)*?\1/XXX/;
        push @strings, $&;
    }

    #    use Data::Dumper;
    #    print Dumper(\@strings);

    return($L,@strings);
}

# split code and comment
sub StripComment
{
    my $L = shift;

    if ($L =~/^(.*?)(!.*)$/)
    {
#        print(">$1< >$2<\n");
        return( $1, $2 );
    }
    return( $L, "" );
}

# combine code and comment, reinsert strings
sub RecreateLine
{
    my $L = $_[0];
    my $C = $_[1];
    my @S = @{$_[2]};

    $L .= $C;
    foreach (@S)
    {
        my $str = $_;
        $L =~ s/XXX/$str/;
    }

    return($L);
}

# some replacements of variables
sub RenameVars
{
    my $L = shift;

    $L =~ s/teilchenin/partIn/ig;
    $L =~ s/teilchenout/partOut/ig;
    $L =~ s/mediumatcollision/mediumAtColl/ig;
    $L =~ s/momentumlrf/momLRF/ig;

    $L =~ s/nucleon_particle/partNucl/ig;
    $L =~ s/pion_particle/partPion/ig;
    $L =~ s/phi_particle/partPhi/ig;
    $L =~ s/omega_particle/partOmega/ig;



    return($L);
}




#$statement = '';
#$_ = $parexp;
#while (/( *)\.(gt|lt|eq|ne|le|ge|and|or|not|eqv|neqv)\.( *)/i) {
#    $begin = $`;
#    $middle = '.'.$2.'.';
#    $end = $';
#    if ($convert) {
#      MATCH:      {
#          if ($middle =~ /\.gt\./) { $middle = '>'; last MATCH; }
#          if ($middle =~ /\.ge\./) { $middle = '>='; last MATCH; }
#          if ($middle =~ /\.lt\./) { $middle = '<'; last MATCH; }
#          if ($middle =~ /\.le\./) { $middle = '<='; last MATCH; }
#          if ($middle =~ /\.eq\./) { $middle = '=='; last MATCH; }
#          if ($middle =~ /\.ne\./) { $middle = '/='; last MATCH; }
#        } }
#    $statement .= $begin.' '.$middle.' ';
#    $_ = $end;
