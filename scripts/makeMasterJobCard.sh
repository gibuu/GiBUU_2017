echo ""
echo "*******************************************************************************************"
echo " This script generates a job card including all namelists which are available in the code."
echo "*******************************************************************************************"
echo ""

# cd to the GiBUU base directory, i.e. make the script callable from anywhere
s=`which "$0"`
d=`dirname "$s"`
cd $d/..

fileName=testRun/jobCards/masterJobCard

rm $fileName

echo "!*******************************************************************************">> $fileName
echo "! Master job card for GiBUU">> $fileName
echo "!" `date`>> $fileName
echo "!*******************************************************************************">> $fileName


echo "Generating the job card. This may take several seconds..."


for file in `find code/ -name '*.f90' -o -name '*.f'`
	do
      		echo "Processing " $file
             	scripts/masterJobCard.pl < $file >> $fileName
	done

echo "...done."

echo ""
echo 'You can find the job card in '`pwd`'/'$fileName
echo ""
