#model oprobit

#data Doprobit.txt

#verbose 1
#delta 0.01

#reltol 0.00001
#abstol 0.00001

#se none
#nboot 5000
regression X1 X2 X3 X4
Ycont Y
