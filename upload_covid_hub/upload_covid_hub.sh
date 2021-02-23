#!/bin/tcsh

foreach method (CMU-TimeSeries Covid19Sim-Simulator COVIDhub-baseline COVIDhub-ensemble GT-DeepCOVID IHME-SEIR Karlen-pypm MOBS-GLEAM_COVID OliverWyman-Navigator UA-EpiCovDA UMass-MechBayes YYG-ParamSearch)

cp -r ../covid19-forecast-hub/data-processed/${method} .

end



