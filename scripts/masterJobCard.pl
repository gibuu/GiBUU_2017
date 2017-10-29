#!/usr/bin/perl
use warnings;
use diagnostics;



# READ from standard input:
@Zeile=<STDIN>;



foreach $wort(@Zeile) {
  $_=$wort;
  # Setting up the current module name
  if (!m#\s*!.*#i) {
    if (m#.*program.*#i) {
      $moduleName="testProgramm";
    }
    if (!m#.*end.*module.*#i) {
      if (!m#.*use.*module.*#i) {
	if (!m#.*module.*procedure#i) {
	  if (m#.*module\s*(\w+)\s*#i) {
	    $moduleName=$1;
	  } 
	}
      }
    }
  }
  # Extracting namelists of the module
  if (!m#\s*!.*#i) {
    if (m#.*namelist.*/\s*(\w+)\s*/\s*#i) {
      $rest=$';
      chomp($rest);
      $name=$1;
      $flagNameList=1;
      $flagAnd=0;
      if (!($moduleName=~ m#.*testProgramm.*#i)) {
	print "\n";
      }
      $firstLine=1;
      if (m#.*&.*#i) {
	$flagAnd=1;
      }
    }
  }
  # Taking care of namelist stretching over several lines
  if (!$firstLine) {
    if ($flagAnd) {
      $_=$wort;
      # & at start and end
      if (m#&.*&#i) {
	$rest=$&;
	if ($rest=~m#\s*#i) {
	  $rest=$';
	  chomp($rest);
	}
      }
      if (!(m#&.*&#i)) {
	# Only & at the end
	if (m#\s*&\s*\n$#i) {
	  $rest=$`;
	  if ($rest=~m#\s*#i) {
	    $rest=$';
	    chomp($rest);
	  }
	}
	# Only & at start
	if (m#^\s*&\s*#i) {
	  $rest=$';
	  chomp($rest);
	}
	# No & at start or end
	if (!m#^\s*&#i) {
	  if (!m#&\s*\n$#i) {
	    $rest=$_;
	    if ($rest=~m#\s*#i) {
	      $rest=$';
	      chomp($rest);
	    }
	  }
	}
      }
      if (!m#&\s*$#i) {
	$flagAnd=0;
	#	    $flagNameList=0;
      }
    }
  }
  # Writing output
  if ($flagNameList) {
    if (!($moduleName=~ m#.*testProgramm.*#i)) {
      #	  print ($flagAnd.$flagNameList.$wort);
      if ($firstLine) {
	print "\&".$name."   ! Module: ".$moduleName."\n";
	$printEnd=1;
      }
#      print ($wort);
      @inputs=split/\s*,\s*/, $rest;
      foreach $input(@inputs) {
	# Take care that there is no "&" in the array
	$no_and=1; #Default value
	if ($input=~m#.*&.*#i) {
	  $no_and=0; # There is a & inside
	}
	if (!($input=~m#.*&.*#i)) {
	  $variablenName=$input;
	  $no_and=1; # there is no & inside
 	}	  
	# Take care of "& somethingUseful" cases
	if ($input=~m#.*&\s#i) {
	  $variablenName=$';
	  $no_and=1; # there is a & and something useful inside $input
	}
	if ($no_and) {
	  $_=$variablenName;
	  s/ //g;
	  $variablenName=$_;
	  # Now go and find the variable definition to deduce the type and the default value
	  $foundVariable=0;
	  $isArray=0;
	  $moduleName2="";
	  foreach $neueZeile(@Zeile) {
	    if (!$foundVariable) {
	      # Setting up the second module name
	      $_=$neueZeile;
	      if (!(m#\s*!.*#i)) {
		if (m#.*program.*#i) {
		  $moduleName2="testProgramm";
		}
		if (!m#.*end.*module.*#i) {
		  if (!m#.*use.*module.*#i) {
		    if (!m#.*module.*procedure#i) {
		      if (m#.*module\s*(\w+)\s*#i) {
			$moduleName2=$1;
		      }
		    }
		  }
		}
	      }
	      if ($moduleName eq $moduleName2) {
		if (m#.*save.*:+\s*($variablenName)\s*=#i) {
		  $foundVariable=1;
		}
#		print $foundVariable.$_;
		if (m#.*save.*\s+($variablenName)\s*=#i) {
		  $foundVariable=1;
		}
#		print $foundVariable.$_;
		if ($foundVariable) {
		  $value_save=$';
		  if (m#.*\(.*\).*$variablenName#i) {
		    $isArray=1;
		  }
		  if (m#$variablenName\s*\(.*\)#i) {
		    $isArray=1;
		  }
		  # check that it's no character array
		  if($isArray){
		    if (m#character#i) {
		      $isArray=0;
		    }
		  }
		  #print $neueZeile;
		  if (!(m#!.*($variablenName).*=.*#i)) {
		    $value=$value_save;
		  }
		  # strip off comments
		  if ($value=~m#\s*!.*#i) {
		    $value=$`;
		  }
		  chomp($value);
		  # get rid off blanks
		  $_=$value;
		  s/ //g;
		  $value=$_;
#		  print $value."\n";

		  if (m#.*logical.*#i) {
		    $type='logical';
		  } else  {
			if (m#.*real.*#i) {
		         $type='real   ';
                        } else {
		         if  (m#.*integer.*#i) {
		         $type='integer';
		         }
		       }
	}
		}
	      }
	    }
	  }
	  if (!$foundVariable) {
	    print "! WARNING. Variable not found: ".$variablenName."\n";
	  } else {
	    #	print $variablenName."  ,  ".$type."  ,  ".$value."\n";
	    if ($isArray) {
              # get rid of array brackets
              $_=$value;
              s/\(\///g;
              s/\/\)//g;
              $value=$_;
	      printf "! %40s \= %20s \n", $variablenName, $value;
	    } else {
	      printf "! %40s \= %20s \n", $variablenName, $value;
	    }
	  }
	}
      }
    }
  }
  
  if (!$flagAnd) {
    if ($printEnd) {
      print "\/\n" ;
      $printEnd=0;
    }
    $flagNameList=0;
  }
  $firstLine=0
}
#print "\n \n\n\n";


