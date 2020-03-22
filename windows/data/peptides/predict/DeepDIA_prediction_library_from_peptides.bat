CALL C:\Users\OmicsolutionF\Anaconda3\Scripts\activate.bat C:\Users\OmicsolutionF\Anaconda3\envs\tensorflow

set pwd=%cd%
md charge2 charge3 irt

R CMD BATCH filter_undetected_peptides_to_library.r

#Predict MS/MS Spectra

copy charge2\*.csv ..\..\models\charge2\
copy charge3\*.csv ..\..\models\charge3\
copy irt\*.csv ..\..\models\irt\
md library

cd ..\..\models\charge2\ 
python ..\..\..\code\deepms2\py\predict.py
copy *.json %pwd%\library

cd ..\charge3\
python ..\..\..\code\deepms2\py\predict.py
copy *.json %pwd%\library

cd ..\irt\
python ..\..\..\code\deeprt\py\predict.py
copy *.csv %pwd%\library

cd %pwd%\library

#RScript.exe ..\..\..\..\code\init.R
#RScript.exe ..\..\..\..\code\generate_spectral_library_for_Spectronaut.R

#R
#source("../../../../code/init.R")
#source("../../../../code/generate_spectral_library_for_Spectronaut.R")
#quit(save = "yes")

copy ..\library.r .\
R CMD BATCH library.r

exit

