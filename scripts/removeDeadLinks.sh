#!/bin/bash

path=`pwd`

cd objects

for i in `find . -type l ! -exec test -r {} \; -print`;
do
  file $i
  rm $i

  j=`echo $i | sed 's/.f90$/.o/g'`
  if (test -r $j); then 
    file $j 
    rm $j
  fi;

  j=`echo $i | sed 's/.F90$/.o/g'`
  if (test -r $j); then
    file $j 
    rm $j
  fi;

  j=`echo $i | sed 's/.f$/.o/g'`
  if (test -r $j); then
    file $j 
    rm $j
  fi;

done;


cd $path