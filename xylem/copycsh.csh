#!/bin/csh
# Copies .csh files from current difectory to bin 
# removes .csh extension
# sets world read and execute premissions
cd csh
foreach file (*.csh)
  set name = $file:r
  cp $file ../bin/$name
  chmod a+rx ../bin/$name
end
