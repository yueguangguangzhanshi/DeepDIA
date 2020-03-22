CALL C:\Users\OmicsolutionF\Anaconda3\Scripts\activate.bat C:\Users\OmicsolutionF\Anaconda3\envs\tensorflow
R CMD BATCH "get negative peptides.r"
python ..\..\..\..\code\deepdetect\py\train_hard_negative.py