#get coverage over 25 protein
Rscript.exe "get detectability from SpectroMine.r"

#get coverage over 25 protein fasta
#open powershell policy
##set-executionpolicy remotesigned
#PowerShell
$fasta=Resolve-Path $PWD\*.fasta
$protein=Resolve-Path $PWD\*proteinAccession.txt
..\..\..\..\code\deepdetect\Filter-Fasta.ps1 $fasta $protein
$name1=dir -name *.filtered.fasta
$test=Dir -name *.detectability.csv
$name2=($test -split '\.')[0]
rename-Item $name1 -NewName "$name2.fasta"
exit