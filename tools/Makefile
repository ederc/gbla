#  CC=gcc
#  CC=clang
#  CXX=g++
#  CXX=clang++
DBG=-O0 -g
WAR=-Wall -Wextra -ansi -Wno-unused-result
#  STA=
#  STA=-static
#  OPT=-O2 $(STA) -DNDEBUG -UDEBUG  -flto
OPT=-march=native -O3 $(STA) -DNDEBUG -UDEBUG  -flto
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

$(DMP)$(dbg): $(DMP).c $(DMP).h types.h tools.h selecter.h ouvrir.h print_helper.h macros.h
	$(CC) $(DBG) $(PRF) $(WAR)  $(DMP).c -o $(DMP)$(dbg)

$(DMP): $(DMP).c $(DMP).h types.h tools.h selecter.h ouvrir.h print_helper.h macros.h
	$(CC) $(OPT) $(PRF) $(WAR)  $(DMP).c -o $(DMP)


### convert (old binary to new binary)
$(CVT)$(dbg): $(CVT).c types.h tools.h selecter.h ouvrir.h macros.h

	$(CC) $(DBG) $(PRF) $(WAR)  $(CVT).c -o $(CVT)$(dbg)

$(CVT): $(CVT).c types.h tools.h selecter.h ouvrir.h macros.h
	$(CC) $(OPT) $(PRF) $(WAR)  $(CVT).c -o $(CVT)


clean:
	-rm *.o
	-rm `find . -type f -executable`


