CALL C:\Users\OmicsolutionF\Anaconda3\Scripts\activate.bat C:\Users\OmicsolutionF\Anaconda3\envs\tensorflow

#The file name of fragment report should end with .FragmentReport.csv, e.g. HeLa.FragmentReport.csv.
R CMD BATCH "train ms2.r"
R CMD BATCH ..\..\..\..\code\deepms2\R\remove_redundant_ions.R
md charge2
move *charge2.ions.json charge2
md charge3
move *charge3.ions.json charge3
cd charge2
python ..\..\..\..\..\code\deepms2\py\train.py
cd ..\charge3
python ..\..\..\..\..\code\deepms2\py\train.py