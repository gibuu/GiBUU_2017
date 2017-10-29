
rm -f File21.txt

for i in `ls jobGlauber000_*`
do
  cp $i jobGlauber

  sed s/'PartE = 1.50'/'PartE = 1.00'/ jobGlauber > jobGlauber.2
  mv jobGlauber.2 jobGlauber

  ./GlauberAbs.x < jobGlauber 
  cat fort.21 >> File21.txt
done

echo " " >>  File21.txt
echo " " >>  File21.txt

for i in `ls jobGlauber000_*`
do
  cp $i jobGlauber

  sed s/'PartE = 1.50'/'PartE = 2.00'/ jobGlauber > jobGlauber.2
  mv jobGlauber.2 jobGlauber

  ./GlauberAbs.x < jobGlauber 
  cat fort.21 >> File21.txt
done

echo " " >>  File21.txt
echo " " >>  File21.txt

for i in `ls jobGlauber000_*`
do
  cp $i jobGlauber

  sed s/'PartE = 1.50'/'PartE = 3.00'/ jobGlauber > jobGlauber.2
  mv jobGlauber.2 jobGlauber

  ./GlauberAbs.x < jobGlauber 
  cat fort.21 >> File21.txt
done

echo " " >>  File21.txt
echo " " >>  File21.txt

for i in `ls jobGlauber000_*`
do
  cp $i jobGlauber

  sed s/'PartE = 1.50'/'PartE = 5.00'/ jobGlauber > jobGlauber.2
  mv jobGlauber.2 jobGlauber

  ./GlauberAbs.x < jobGlauber 
  cat fort.21 >> File21.txt
done

echo " " >>  File21.txt
echo " " >>  File21.txt

for i in `ls jobGlauber000_*`
do
  cp $i jobGlauber

  sed s/'PartE = 1.50'/'PartE = 10.00'/ jobGlauber > jobGlauber.2
  mv jobGlauber.2 jobGlauber

  ./GlauberAbs.x < jobGlauber 
  cat fort.21 >> File21.txt
done

echo " " >>  File21.txt
echo " " >>  File21.txt