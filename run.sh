# 
#$ -cwd 
#$ -j y 
#$ -pe openmp 24
#$ -S /bin/bash 
#
MATLABPATH=/usr/local/bin
FILEPATH=/cis/home/rtang/SBM3_Parallel_n150_p1/Part1
M_FILENAME=main.m
# Name the job #$ -N matlabScript #
cd $FILEPATH; 
$MATLABPATH/matlab -nosplash -nodisplay -r "main; quit;"
echo "" 
echo "Done at " `date` 