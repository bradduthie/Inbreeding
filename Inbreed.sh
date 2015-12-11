# SGE SUBMISSION SCRIPT 
# This is an embarrassingly parallel simulation
# Submit to 1 core -- should print to one output file
# I estimate that this program will take ca 48 hrs to run
#
#$ -V -cwd -o out.$JOB_ID
./InMult
