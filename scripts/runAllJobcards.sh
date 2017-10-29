#!/bin/bash
#*******************************************************************************
#****e* /runAllJobcards.sh
# NAME
# runAllJobcards.sh
# PURPOSE
# This script generates a subdirectory in testRun/ and then runs all jobards
# it can find in the directory testRun/jobCards
#*******************************************************************************

path=`pwd`
rundir=`date +"RunAll_%Y-%m-%d_%H-%M-%S"`
pathrundir="${path}/testRun/${rundir}"
echo "Generate: ${pathrundir}"
mkdir $pathrundir
cd $pathrundir



for job in `ls ${path}/testRun/jobCards/*.job`; do
    jobbase=`basename ${job} .job`
    echo "Running: ${jobbase} ..."
    mkdir ${jobbase}
    cd ${jobbase}
    ln -s ${path}/testRun/jobCards/022_ExternalSource.inp ExternalSource.inp
    /usr/bin/time ../../GiBUU.x < ${job} >& ../${jobbase}.rep
    cd ..
done
