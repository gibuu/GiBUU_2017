#!/bin/bash

# This script looks through all fortran files and searches for 
#    read(5,nml=...
# and replaces this with
#    write(107,nml=...
#    read(5,nml=...
#    write(108,nml=...
#
# Running the code, in file fort.107 and fort.108, all values of all used 
# namelists are reported, once before reading the jobcard, once after reading 
# the values.
# Thus, by diffing the two files, you may extract a 'minimal jobcard' doing
# the same as the given jobcard.
#
# You may undo all changes of this script by 'svn revert -R .'
#


for f in `find -path '*/.svn' -prune -o -type f -print0 | xargs -0 grep -I -l -i "read(5,nml="`; 
do
awk '{ if (/read\(5,nml(.)*/) { print gensub(/read\(5,nml/,"write(107,nml",1);  print $0; print gensub(/read\(5,nml/,"write(108,nml",1); } else { print $0 } }' $f > /tmp/aaa
ls -l $f
mv /tmp/aaa $f
done
