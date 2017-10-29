#!/usr/bin/perl
#
# This script generates as output a LaTeX file with the information
# about all namelists, the used variables and their default values
# and a description.
#
# In order to use this:
# [~/GiBUU/workingCode]$ make doku
# [~/GiBUU/workingCode]$ ./scripts/GenerateNamelistDoku.pl > namelists.tex
# [~/GiBUU/workingCode]$ pdflatex namelists.tex
#
# please repeat the call to pdflatex, until no warning messages about
# 'changed width' occurs anymore.

use strict;

print "\\documentclass{article}\n";
print "\\usepackage[a4paper,margin=0.7in,landscape]{geometry}\n";
print "\\usepackage{booktabs}\n";
print "\\usepackage{longtable}\n";

print "\\setlength{\\parindent}{0pt}\n";
print "\\makeatletter\n";
print "\\newcommand\\hlinenobreak{%\n";
print "  \\multispan\\LT\@cols\n";
print "  \\unskip\\leaders\\hrule\\\@height\\arrayrulewidth\\hfill\\\\*}";
print "\\makeatother\n";
print "\\begin{document}\n";

my $RoboPath = "../Documentation/code/";
my $fNameIndex = $RoboPath."robo_namelist.html";

my $fNameHTML;
my $sNamelist;
my $TAG;

my $fileHTML;
my $isOpen = 0;
my $isPre = 0;
my $isNotFirst = 0;

my $NamelistsHTML = { };
my $NamelistsTAG  = { };

my $VarsTAG = { };
my $VarsTyp = { };
my $VarsValue = { };
my $VarsText = { };

### Scan index, build up the list of namelists

open(file, $fNameIndex) || die "ERROR: can not open '$fNameIndex': $!.\n";
while (<file>)
{
    if ( /<a href=.(\S+)\#(\S+)\s+\S+\s+>(\S+)<\/a>/ )
    {
	$fNameHTML = $1;
	$TAG = $2;
	$sNamelist = $3;

	$NamelistsHTML->{$sNamelist} = $fNameHTML;
	$NamelistsTAG->{$sNamelist} = $TAG;

#	print "$fNameHTML   $sNamelist\n";
    }
}
close(file);

### Loop over all namelists
foreach my $sNamelist (sort  {lc $a cmp lc $b} keys %{$NamelistsHTML})
{
    $fNameHTML = $NamelistsHTML->{$sNamelist};
    $TAG = $NamelistsTAG->{$sNamelist};

    my $fNameF = $fNameHTML;
    $fNameF =~ s/_f90.html/.f90/;
    $fNameF =~ s|\./|code/|;

    my @Vars;

    print "% $sNamelist     $fNameF\n\n";
#    print "$fNameHTML   $sNamelist  $TAG\n";

### LaTeX output:
    my $sNamelistTeX = FormatTeX($sNamelist);
    my $fNameFTeX = FormatTeX($fNameF);

#    print "\\begin{table}{h}\n";

#    print "\\begin{longtable}{|l|l|l|l|}\n";
#    print "\\hlinenobreak\n";
#    print "\\textbf{{\\HUGE $sNamelistTeX}} & \\multicolumn{3}{|l|}{\\footnotesize{$fNameFTeX}}\\\\*\n";
#    print "\\hlinenobreak\n";

    print "\\begin{longtable}{llll}\n";
    print "\\toprule\n";
    print "\\textbf{\\large{$sNamelistTeX}} & \\multicolumn{3}{l}{\\footnotesize{$fNameFTeX}}\\\\*\n";
    print "\\midrule\n";
    print "\\endfirsthead\n";
    print "\\midrule\n";
    print "\\endhead\n";

    $isNotFirst = 0;
    $isOpen = 0;
    $VarsTAG = { };
    open(file, $RoboPath.$fNameHTML)  || die "ERROR: can not open '$fNameHTML': $!.\n";
    while (<file>)
    {
	# find the corresponding part...
	if ($isOpen > 0)
	{
	    if (/<\/ul>/)
	    {
#		print "closed\n";
		$isOpen = 0;
	    }
	}
	else
	{
	    if (/<a name="$TAG>/)
	    {
#		print "found. \n";
		$isOpen = 1;
	    }
	}

	# scan the corresponding part...
	if ($isOpen > 0)
	{
#	    print $_;

	    if (/<li> <a href="(\S+)">(\S+)<\/a>/)
	    {
#		print "  $2   $1\n";
		$VarsTAG->{$2} = $1;
		push @Vars, $2;

	    }
	    elsif (/<li>\s+(\S+)/)
	    {
#		print "  $1\n";
		$VarsTAG->{$1} = 'undocumented';
		push @Vars, $1;
	    }
	}


    }
    close(file);

### check, whether we can correct a link

    foreach my $var (sort {lc $a cmp lc $b} keys %{$VarsTAG})
    {
	$TAG = $VarsTAG->{$var};

	if ($TAG =~ /\S+#/)
	{
#	    print "$sNamelist     $fNameHTML: $var\n";

	    open(file, $RoboPath.$fNameHTML)  || die "ERROR: can not open '$fNameHTML': $!.\n";
	    while (<file>)
	    {
		if (/<li>\S+ <a href="(\S+)">\S+\/$var<\/a><\/li>/)
		{
#		    print "... $1\n";
		    $VarsTAG->{$var} = $1;
		}
	    }
	    close(file)

	}
    }

### Now try to extract additional informations

    foreach my $var (sort {lc $a cmp lc $b} keys %{$VarsTAG})
    {
	$TAG = $VarsTAG->{$var};
	$TAG =~ s/#//;

#	print "searching for $TAG ($var)\n";
	my $text = "";

	$isOpen = 0;
	$isPre = 0;
	open(file, $RoboPath.$fNameHTML)  || die "ERROR: can not open '$fNameHTML': $!.\n";
	while (<file>)
	{
	    my $line = $_;
	    my $VarTyp = "";
	    my $VarValue = "";

	    if (/<a name="$TAG">/)
	    {
		$isOpen = 1;
		$line = "";
	    }
	    elsif (/<a name="robo\d+">/)
	    {
		$isOpen = 0;
	    }
	    elsif (/<a name="\S+">/)
	    {
		$isOpen = 0;
	    }

	    if ($isOpen)
	    {
		$line =~ s/<a href="\S+">(\S+)<\/a>/$1/sig;
		if ($line =~ /^<p>\[ Top \]/) {$line = "";}
		if ($line =~ /^<p class="item_name">SOURCE<\/p>/) {$line = "";}
		if ($line =~ /^<hr \/>/) {$line = "";}

		$line =~ s/<span class="\S+">//sig;
		$line =~ s/<\/span>//sig;
		$line =~ s/<strong>//sig;
		$line =~ s/<\/strong>//sig;


		if ($line =~ /^<\/pre>/)
		{
		    $line = "";
		    $isPre = 0;
		}

		if ($line =~ /^<pre class="source">/)
		{
		    $isPre = 1;
		    $line =~ s/^<pre class="source">//sig;

		    $line =~ s/^\s*//sig;
		    $line =~ s/\!.*$//sig; # remove comment at end of line

		    $line =~ s/::/§/sig; # replace '::' by something else

#		    if ($line =~ /\s+(\S+)\s+::\s+(\S+)\s+=\s+(\S+)/)
		    if ($line =~ m/([^§]*)§\s*$var(.*)/)
		    {
			$VarTyp = $1;
			$VarValue = $2;

#			print "\n---- >$VarTyp< --- >$VarValue<\n";
			$VarTyp =~ s/SAVE//sig;
			$VarTyp =~ s/PUBLIC//sig;
			$VarTyp =~ s/,\s*,/, /sig;
			$VarTyp =~ s/\s*,\s*/, /sig;
			$VarTyp =~ s/^\s*,\s*//sig;
			$VarTyp =~ s/\s*,\s*$//sig;



			$VarValue =~ s/=//sig;
			$VarValue = trim($VarValue);

			# ensure same capitalization
			$VarTyp = EnsureCap($VarTyp);
			$VarValue = EnsureCap($VarValue);

			# remove some HTML stuff
#			$VarValue =~ s/&amp;/&/sig;
			$VarValue =~ s/&amp;/\\dots/sig;

			$VarsTyp->{$var} = $VarTyp;
			$VarsValue->{$var} = $VarValue;
#			print "---- >$VarTyp< --- >$VarValue<\n";

			$line = "";
		    }


		}

### Delete lines which are multiple 'source' lines
		if ($isPre)
		{
#		    $line = $line." ISPRE";
		    $line = "";
		}

		$line =~ s/\n//sig;
#		print ">|$line|<\n";
		$text = $text.$line;
	    }

	}
	$text =~ s/<p class="item_name">PURPOSE<\/p>//sig;
#	print ">|$text|<\n";
	$VarsText->{$var} = $text;
    }

#    print "\n";


### print the list of variables

#    foreach my $var (sort {lc $a cmp lc $b} keys %{$VarsTAG})
#    foreach my $var (keys %{$VarsTAG})
    foreach my $var (@Vars)
    {
	$TAG = $VarsTAG->{$var};

	my $mark = "";
	if ($TAG =~ /^#/)
	{
	    $mark = " ...";
	}
	elsif ($TAG =~ /#/)
	{
	    $mark = " ###";
	}
	else
	{
	    $mark = " 000";
	}

	my $varN = FormatTeX($var);

	my $typ = FormatTeX($VarsTyp->{$var});
	my $value = FormatTeX($VarsValue->{$var});
#	my $text = $VarsText->{$var};
	my $text = FormatTeX($VarsText->{$var});

#	print "$mark $var   $TAG  $VarsTyp->{$var}  $VarsValue->{$var}\n";
#	print "$mark $var   | $VarsTyp->{$var} | $VarsValue->{$var}\n";
#	print " $varN   | $typ | $value | $text\n";

### LaTeX output:
        if ($isNotFirst > 0)
        {
            print "\\midrule\n";
        }
        $isNotFirst = 1;
#	print "$varN & $typ & $value & \\begin{minipage}[t]{8cm}$text\\end{minipage}\\\\\n";
	print "$varN & \\begin{minipage}[t]{2cm}$typ\\end{minipage} & \\begin{minipage}[t]{2cm}$value\\end{minipage} & \\begin{minipage}[t]{12cm}$text\\end{minipage}\\\\*\n";


    }

### LaTeX output

#    print "\\hline\n";
    print "\\bottomrule\n";
    print "\\end{longtable}\n";
    print "{ }\n";

    print "\n";

    print "\n\n";


}


print "\\end{document}\n";

sub EnsureCap
{
    (my $text) = @_;

    # ensure same capitalization
    $text =~ s/logical/logical/sig;
    $text =~ s/integer/integer/sig;
    $text =~ s/real/real/sig;

    $text =~ s/.true./.true./sig;
    $text =~ s/.false./.false./sig;

    return $text;
}

sub trim
{
    (my $text) = @_;
    $text =~ s/^\s+//sig;
    $text =~ s/\s+$//sig;
    return $text;
}

sub FormatTeX
{
    (my $text) = @_;

    $text =~ s/<p>//g;
    $text =~ s/<\/p>/\\\\/g;

    $text =~ s/<ul>/\\begin{itemize}\\leftmargin0em\\itemsep0pt\\itemindent0pt/sig;
    $text =~ s/<\/ul>/\\end{itemize}/sig;

    $text =~ s/<li>/\\item/sig;
    $text =~ s/<\/li>//sig;

    $text =~ s/&gt;/\$>\$/g;
    $text =~ s/&lt;/\$<\$/g;
    $text =~ s/&amp;/ and /sig;
    $text =~ s/_/\\_/g;

    $text =~ s/~/\\~{}/g;
    $text =~ s/\^/\\^{}/g;
    $text =~ s/\#/\\#/g;
    $text =~ s/\%/\\%/g;
#    $text =~ s/\&/\\&/g;

    $text =~ s/\\rho/\$\\rho\$/g;

    $text =~ s/<p class="item\\_name">//g;

    $text =~ s/\\end\{itemize\}\\\\/\\end\{itemize\}/g;
    $text =~ s/\\\\\\begin\{itemize\}/\\begin\{itemize\}/g;


    $text = trim($text);

    $text =~ s/\\\\$//g; # remove trailing newline

    $text =~ s/\$<\$-\$>\$/\$\\leftrightarrow\$/g;
    $text =~ s/-\$>\$/\$\\rightarrow\$/g;

    return $text;
}
