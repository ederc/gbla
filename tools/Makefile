#  CC=gcc
#  CC=clang
#  CXX=g++
#  CXX=clang++
DBG=-O0 -g
WAR=-Wall -Wextra -ansi
OPT=-march=native -O3 -DNDEBUG -UDEBUG  -flto
#  STA=-static -O2 -DNDEBUG -UDEBUG  -flto
#  OPT=
LIB=-lm
#  PRF=-pg
PRF=

par=_par
PAR=-fopenmp

DMP=dump_matrix
CVT=converter
RUN=runit
opt=_opt
dbg=_dbg

all:
	make $(DMP)
	make $(CVT)


### dump matrix (to sms...)

$(DMP)$(dbg): $(DMP).c $(DMP).h
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)$(dbg)

$(DMP): $(DMP).c $(DMP).h
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)


### convert (old binary to new binary)
$(CVT)$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(CVT).c -o $(CVT)$(dbg)

$(CVT): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(CVT).c -o $(CVT)


clean:
	-rm *.o
	-rm `find . -type f -executable`


