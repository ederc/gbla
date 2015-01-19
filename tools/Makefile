CC=gcc
#  CC=clang
CXX=g++
#  CXX=clang++
DBG=-O0 -g
WAR=-Wall -Wextra -ansi
STA=
#  STA=-static
OPT=-march=native -O3 -DNDEBUG -UDEBUG  -flto
NEW=-D_PROPOSED_FORMAT
LIB=-lm
#  PRF=-pg
PRF=

par=_par
PAR=-fopenmp

TST=dumb_write
DMP=dump_matrix
CVT=converter
RUN=runit
opt=_opt
dbg=_dbg

all:
	make $(DMP)_new$(opt)
	make $(CVT)_rev_sor$(opt)


### dump matrix (to sms...)
$(DMP)$(dbg): $(DMP).c $(DMP).h 
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)$(dbg)

$(DMP)$(opt): $(DMP).c $(DMP).h 
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)$(opt)

$(DMP)_new$(dbg): $(DMP).c $(DMP).h
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(DMP).c -o $(DMP)_new$(dbg)

$(DMP)_new$(opt): $(DMP).c $(DMP).h 
	$(CC) $(OPT) $(PRF) $(NEW) $(WAR) $(STA) $(DMP).c -o $(DMP)_new$(opt)

### convert (old binary to new binary)
$(CVT)_rev$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(dbg)

$(CVT)$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -UREVERT -USORT $(CVT).c -o $(CVT)$(dbg)

$(CVT)_rev$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(opt)

$(CVT)$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -UREVERT -USORT $(CVT).c -o $(CVT)$(opt)

$(CVT)_rev_sor$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(dbg)

$(CVT)_sor$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(dbg)

$(CVT)_rev_sor$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(opt)

$(CVT)_sor$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(opt)



clean:
	-rm *.o
	-rm `find . -type f -executable`

check:
	-diff types.h ../draft/types.h
	-diff macros.h ../draft/macros.h
	-diff selecter.h ../draft/selecter.h
	-diff ouvrir.h ../draft/ouvrir.h

