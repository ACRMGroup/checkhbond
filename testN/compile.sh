INC=/home/bsm/martin/include
LIB=/home/bsm/martin/lib
OPT="-ansi -Wall -strict -D DEBUG"
cc -g -o mytest $OPT -I$INC -I$INC/bioplib -L$LIB mytest.c ReadPDB.c -lbiop -lgen -lm