#!/usr/bin/perl -w

# This script is a modified version of egypt 1.6
# it does not produce output suitable to be processed via dot, but
# instead it generates a text output, where for every function its calling
# functions are listed
#
# modifications in respect to the original egypt script are marked with #!#


#!/usr/bin/perl -w

=head1 NAME

egypt - create call graph from gcc RTL dump

=head1 SYNOPISIS

 egypt [--omit function,function,...] [--include-external] <rtl-file>... | dotty -
 egypt [--omit function,function,...] [--include-external] <rtl-file>... | dot <dot-options>

=head1 DESCRIPTION

Egypt is a simple tool for creating call graphs of C programs.

=head1 OPTIONS

=over 8

=item omit

Omit the given functions from the call graph.  Multiple function names
may be given separated by commas.

=item include-external

Include calls to external functions in the call graph.  A function is
considered external if it is not defined in any of the input files.
For example, functions in the standard C library are external.  By
default, such calls are not shown.  Even with this option, only direct
calls will be shown; there is no way to show indirect calls to
external functions.

=back

=head1 HOW IT WORKS

The two major tasks in creating a call graph are analyzing the syntax
of the source code to find the function calls and laying out the
graph, but Egypt actually does neither.  Instead, it delegates the
source code analysis to GCC and the graph layout to Graphviz, both of
which are better at their respective jobs than egypt could ever hope
to be itself.  Egypt itself is just a small Perl script that acts as
glue between these existing tools.

Egypt takes advantage of GCC's capability to dump an intermediate
representation of the program being compiled into a file (a I<RTL
file>); this file is much easier to extract information from than a C
source file.  Egypt extracts information about function calls from the
RTL file and massages it into the format used by Graphviz.

=head1 GENERATING THE RTL FILE

Compile the program or source file you want to create a call graph for
with gcc, adding the option "-dr" to CFLAGS.  This option causes gcc
to dump its intermediate code representation of each file it compiles
into a a file.

For example, the following works for many programs:

   make clean
   make CFLAGS=-dr

Depending on the GCC version, the RTL file for a source file F<foo.c>
may be called something like F<foo.c.rtl>, F<foo.c.00.rtl>, or
F<foo.c.00.expand>.

=head1 VIEWING THE CALL GRAPH

To view the call graph in an X11 window, run egypt with one or
more .rtl files as command line arguments and pipe its output to the
B<dotty> program from the Graphviz package.  For example, if you
compiled F<foo.c> with C<gcc -dr> to generate F<foo.c.00.rtl>, use

    egypt foo.c.00.rtl | dotty -

=head1 PRINTING THE CALL GRAPH

To generate a PostScript version of the call graph for printing, use
the B<dot> program from the Graphviz package.  For example, to generate
a callgraph in the file F<callgraph.ps> fitting everything on a US
letter size page in landscape mode, try

   egypt foo.c.00.rtl | dot -Grotate=90 -Gsize=11,8.5 -Tps -o callgraph.ps

Sometimes, the graph will fit better if function calls go from left to
right instead of top to bottom.  The B<dot> option B<-Grankdir=LR>
will do that:

   egypt foo.c.00.rtl | dot -Gsize=8.5,11 -Grankdir=LR -Tps -o callgraph.ps

For nontrivial programs, the graph may end up too small
to comfortably read.  If that happens, try N-up printing:

   egypt foo.c.00.rtl | dot -Gpage=8.5,11 -Tps -o callgraph.ps

You can also try playing with other B<dot> options such as B<-Gratio>,
or for a different style of graph, try using B<neato> instead of
B<dot>.  See the Graphwiz documentation for more information about the
various options available for customizing the style of the graph.

=head1 READING THE CALL GRAPH

Function calls are displayed as solid arrows.  A dotted arrow means
that the function the arrow points from takes the address of the
function the arrow points to; this typically indicates that the latter
function is being used as a callback.

=head1 WHY IS IT CALLED EGYPT?

Egypt was going to be called B<rtlcg>, short for I<RTL Call Graph>,
but it turned out to be one of those rare cases where ROT13'ing the
name made it easier to remember and pronounce.

=head1 SEE ALSO

L<gcc>, L<dotty>, L<dot>, L<neato>

=head1 COPYRIGHT

Copyright 1994-2006 Andreas Gustafsson

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 AUTHOR

Andreas Gustafsson

=cut

use strict;
use Getopt::Long;

use vars qw($VERSION);

$VERSION = "1.6";

# A data structure containing information about potential function
# calls.  This is a reference to a hash table where the key is a
# the name of a function (the caller) and the value is a reference
# to another hash table indexed by the name of a symbol referenced
# by the caller (the potential callee) and a value of "direct"
# (if the reference is a direct function call) or "indirect"
# (if the reference is a non-function-call symbol reference;
# if the referenced symbol itself turns out to be a function,
# this will be considered an indirect function call).

my $calls = { };
my $called = { }; #!#

# A map from mangled C++ names to the corresponding demangled ones
my $demangle = { };

# The current function
my $curfunc;

# Functions to omit
my @omit = ();
my $include_external = 0;

# Mapping from symbol reference types to dot styles
my $styles = {
    direct => 'solid',
    indirect => 'dotted'
};

sub demangle {
    my ($name) = @_;
    return $demangle->{$name} || $name;
}

sub demangle2 {  #!#
  my ($name) = @_;#!#
  my ($name2);#!#
#!#
  if ($name =~ /__(\S+)_MOD_(\S+)/) {#!#
    $name2 = sprintf "(%s) %s",$1,$2;#!#
  }#!#
  else {#!#
#    $name2 = $demangle->{$name} || $name;#!#
    $name2 = $name;#!#
  }#!#
#!#
  return $name2;#!#
}#!#

sub demangle3 {#!#
    my ($name) = @_;#!#
    return $name;#!#
}#!#

sub CheckOK {#!#
  my ($name) = @_;#!#
  my ($val) = 0;#!#

  if ($name =~ /__(\S+)_MOD_(\S+)/) {#!#
    $val = 1;#!#
  } elsif ($name =~ /.+_$/) {#!#
    $val = 2;#!#
  } elsif ($name =~ /.+\.\d+$/) {#!#
    $val = 3;#!#
  }#!#
#!#
#  print " $name ==> $val\n";#!#
#!#
  return $val;#!#
}#!#



GetOptions('omit=s' => \@omit,
	   'include-external' => \$include_external);

@omit = split(/,/, join(',', @omit));

while (<>) {
    chomp;
    if (/^;; Function (\S+)\s*$/) {
	# pre-gcc4 style
	$curfunc = $1;
	$calls->{$curfunc} = { } if ! exists($calls->{$curfunc});
	$called->{$curfunc} = { } if ! exists($called->{$curfunc});#!#
    } elsif (/^;; Function (.*)\s+\((\S+)(,.*)?\).*$/) { #!# from egypt 1.10
#    } elsif (/^;; Function (.*)\s+\((\S+)\)$/) { #!#
	# gcc4 style
	$curfunc = $2;
	$demangle->{$curfunc} = $1;
	$calls->{$curfunc} = { } if ! exists($calls->{$curfunc});
	$called->{$curfunc} = { } if ! exists($called->{$curfunc});#!#
    }
    if (/^.*\(call.*"(.*)".*$/) {
	$calls->{$curfunc}->{$1} = "direct";
    } elsif (/^.*\(symbol_ref.*"(.*)".*$/) {
#	$calls->{$curfunc}->{$1} = "indirect";
    }
}

delete @$calls{@omit};

my %omit_map;
@omit_map{@omit} = ();

#!# from here til end of file: new output routine
foreach my $caller (keys %{$calls}) {
  foreach my $callee (keys %{$calls->{$caller}}) {
    my $reftype = $calls->{$caller}->{$callee};
    $called->{$callee}->{$caller} = $reftype;
  }
}

foreach my $callee (sort keys %{$called}) {
  next if exists $omit_map{$callee};
  next if ! CheckOK($callee);
#  my $hhh = CheckOK($callee);


  my $callee_d = demangle2($callee);
  print "\n$callee_d:\n";
  foreach my $caller (keys %{$called->{$callee}}) {

    my $reftype = $called->{$callee}->{$caller};
    my $style = $styles->{$reftype};
    my $caller_d = demangle2($caller);

    print "    <--- $caller_d\n";
  }
}

#!# The following is the original output routine:
#!#
# print "digraph callgraph {\n";
#
# foreach my $caller (keys %{$calls}) {
#     foreach my $callee (keys %{$calls->{$caller}}) {
# 	my $reftype = $calls->{$caller}->{$callee};
# 	# If the referenced symbol is not a defined function
# 	# or a direct call to an external function, ignore it.
# 	next unless exists($calls->{$callee}) or
# 	    $include_external and $reftype eq 'direct'
# 	        and ! exists $omit_map{$callee};
# 	my $style = $styles->{$reftype};
# 	my $caller_d = demangle3($caller);
# 	my $callee_d = demangle3($callee);	
# 	print "\"$caller_d\" -> \"$callee_d\" [style=$style];\n";
#     }
# }
