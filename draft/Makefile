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

TST=dumb_write
DMP=dump_matrix
CVT=converter
RUN=runit

opt=_opt
dbg=_dbg
par=_par
seq=_seq

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

# SEQ

LIBS$(seq)$(dbg): $(RED).cpp
	$(CXX) $(DBG) $(PRF)  $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(seq)$(dbg).o


OBJ$(seq)$(dbg):
	$(CC) $(DBG) $(PRF)  $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(seq)$(dbg).o



LIBS$(seq)$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF)  $(STA)  $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(seq)$(opt).o


OBJ$(seq)$(opt):
	$(CC) $(OPT) $(PRF) $(WAR)  $(STA) -c $(RUN).c -o $(RUN)$(seq)$(opt).o




$(RUN)$(seq)$(dbg): $(RUN).c
	make LIBS$(seq)$(dbg) ;
	make OBJ$(seq)$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA)  $(RED)$(seq)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(seq)$(dbg) -lm -lstdc++


$(RUN)$(seq)$(opt): $(RUN).c
	make LIBS$(seq)$(opt) ;
	make OBJ$(seq)$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(seq)$(opt).o  $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(seq)$(opt) -lm -lstdc++



# PAR

LIBS$(par)$(dbg): $(RED).cpp
	$(CXX) $(DBG) $(PRF) $(PAR) $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(par)$(dbg).o


OBJ$(par)$(dbg):
	$(CC) $(DBG) $(PRF) $(PAR) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(par)$(dbg).o



LIBS$(par)$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF)  $(STA) $(PAR) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(par)$(opt).o


OBJ$(par)$(opt):
	$(CC) $(OPT) $(PRF) $(WAR) $(PAR) $(STA) -c $(RUN).c -o $(RUN)$(par)$(opt).o




$(RUN)$(par)$(dbg): $(RUN).c
	make LIBS$(par)$(dbg) ;
	make OBJ$(par)$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(PAR) $(RED)$(par)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(par)$(dbg) -lm -lstdc++


$(RUN)$(par)$(opt): $(RUN).c
	make LIBS$(par)$(opt) ;
	make OBJ$(par)$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(par)$(opt).o $(PAR) $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(par)$(opt) -lm -lstdc++


clean:
	-rm *.o
	-rm `find . -type f -executable`

