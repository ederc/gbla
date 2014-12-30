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
TGT=dump_matrix
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
$(TGT)$(dbg): $(TGT).c $(TGT).h printer.h
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(TGT).c -o $(TGT)$(dbg)

$(TGT)$(opt): $(TGT).c $(TGT).h printer.h
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(TGT).c -o $(TGT)$(opt)

$(TGT)_new$(dbg): $(TGT).c $(TGT).h printer.h
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(TGT).c -o $(TGT)_new$(dbg)

$(TGT)_new$(opt): $(TGT).c $(TGT).h printer.h
	$(CC) $(OPT) $(PRF) $(NEW) $(WAR) $(STA) $(TGT).c -o $(TGT)_new$(opt)

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
	$(CXX) $(DBG) $(PRF)  $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(dbg).o

LIBS$(dbg)$(par): $(RED).cpp
	$(CXX) $(DBG) $(PRF) $(PAR) $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(dbg)$(par).o


OBJ$(dbg):
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(dbg).o

OBJ$(dbg)$(par):
	$(CC) $(DBG) $(PRF) $(PAR) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(dbg)$(par).o


LIBS$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF)  $(STA) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(opt).o

LIBS$(opt)$(par): $(RED).cpp
	$(CXX) $(OPT) $(PRF)  $(STA) $(PAR) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(opt)$(par).o


OBJ$(opt):
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(opt).o

OBJ$(opt)$(par):
	$(CC) $(OPT) $(PRF) $(WAR) $(PAR) $(STA) -c $(RUN).c -o $(RUN)$(opt)$(par).o



$(RUN)$(dbg): $(RUN).c
	make LIBS$(dbg) ;
	make OBJ$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(RED)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(dbg) -lm -lstdc++

$(RUN)$(dbg)$(par): $(RUN).c
	make LIBS$(dbg)$(par) ;
	make OBJ$(dbg)$(par)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(PAR) $(RED)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(dbg)$(par) -lm -lstdc++


$(RUN)$(opt): $(RUN).c
	make LIBS$(opt) ;
	make OBJ$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(opt).o $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(opt) -lm -lstdc++

$(RUN)$(opt)$(par): $(RUN).c
	make LIBS$(opt)$(par) ;
	make OBJ$(opt)$(par)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(opt)$(par).o $(PAR) $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(opt)$(par) -lm -lstdc++


clean:
	-rm *.o
	-rm `find . -type f -executable`

