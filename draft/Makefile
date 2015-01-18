CC=gcc -std=c11
#  CC=clang
CXX=g++ 
#  CXX=clang++
DBG=-O0 -g -mavx
WAR=-Wall -Wextra -ansi
STA=
#  STA=-static
OPT=-march=native -mtune=native -O3 -DNDEBUG -UDEBUG  -flto
NEW=-D_PROPOSED_FORMAT
LIB=-lm
#  PRF=-pg
PRF=

PAR=-fopenmp
#  PAR=

TST=dumb_write
DMP=dump_matrix
CVT=converter
RUN=runit
opt=_opt
dbg=_dbg

### dumb write (test matrix)
$(TST)_new$(dbg): $(TST).c
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(TST).c -o $(TST)_new$(dbg)

$(TST)_old$(dbg): $(TST).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(TST).c -o $(TST)_old$(dbg)

$(TST)_new$(opt): $(TST).c
	$(CC) $(OPT) $(PRF) $(NEW) $(WAR) $(STA) $(TST).c -o $(TST)_new$(opt)

$(TST)_old$(opt): $(TST).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(TST).c -o $(TST)_old$(opt)

### dump matrix (to sms...)
$(DMP)$(dbg): $(DMP).c $(DMP).h printer.h
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)$(dbg)

$(DMP)$(opt): $(DMP).c $(DMP).h printer.h
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)$(opt)

$(DMP)_new$(dbg): $(DMP).c $(DMP).h printer.h
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(DMP).c -o $(DMP)_new$(dbg)

$(DMP)_new$(opt): $(DMP).c $(DMP).h printer.h
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

### convert (old binary to new binary)
$(CVT)_rev_sor$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(dbg)

$(CVT)_sor$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(dbg)

$(CVT)_rev_sor$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(opt)

$(CVT)_sor$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(opt)


### run the algo
RED=fflas_reduce
FFLAS=/home/bboyer/usr/
FFLAS_CFLAGS=`sh $(FFLAS)/bin/fflas-ffpack-config --cflags-full`
FFLAS_LIBS=`sh $(FFLAS)/bin/fflas-ffpack-config --libs`


LIBS$(dbg): $(RED).cpp
	$(CXX) $(DBG) $(PRF) $(PAR) $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(dbg).o


OBJ$(dbg):
	$(CC) $(DBG) $(PRF) $(PAR) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(dbg).o



LIBS$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF)  $(STA) $(PAR) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(opt).o


OBJ$(opt):
	$(CC) $(OPT) $(PRF) $(WAR) $(PAR) $(STA) -c $(RUN).c -o $(RUN)$(opt).o




$(RUN)$(dbg): $(RUN).c
	make LIBS$(dbg) ;
	make OBJ$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(PAR) $(RED)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(dbg) -lm -lstdc++


$(RUN)$(opt): $(RUN).c
	make LIBS$(opt) ;
	make OBJ$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(opt).o $(PAR) $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(opt) -lm -lstdc++


clean:
	-rm *.o
	-rm `find . -type f -executable`

