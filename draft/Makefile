DBG=-O0 -g
WAR=-Wall -Wextra -ansi
OPT=-march=native -O3 -DNDEBUG -UDEBUG
NEW=-D_PROPOSED_FORMAT
LIB=-lm
PRF=-pg
#  PRF=

TST=dumb_write
TGT=dump_matrix
CVT=converter
RUN=runit
opt=_opt
dbg=_dbg

### dumb write (test matrix)
$(TST)_new$(dbg): $(TST).c
	gcc $(DBG) $(PRF) $(NEW) $(WAR) $(TST).c -o $(TST)_new$(dbg)
$(TST)_old$(dbg): $(TST).c
	gcc $(DBG) $(PRF) $(WAR) $(TST).c -o $(TST)_old$(dbg)

$(TST)_new$(opt): $(TST).c
	gcc $(OPT) $(PRF) $(NEW) $(WAR) $(TST).c -o $(TST)_new$(opt)

$(TST)_old$(opt): $(TST).c
	gcc $(OPT) $(PRF) $(WAR) $(TST).c -o $(TST)_old$(opt)

### dump matrix (to sms...)
$(TGT)$(dbg): $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(PRF) $(WAR) $(TGT).c -o $(TGT)$(dbg)

$(TGT)$(opt): $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(PRF) $(WAR) $(TGT).c -o $(TGT)$(opt)

$(TGT)_new$(dbg): $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(PRF) $(NEW) $(WAR) $(TGT).c -o $(TGT)_new$(dbg)

$(TGT)_new$(opt): $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(PRF) $(NEW) $(WAR) $(TGT).c -o $(TGT)_new$(opt)

### convert (old binary to new binary)
$(CVT)_rev$(dbg): $(CVT).c
	gcc $(DBG) $(PRF) $(WAR) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(dbg)

$(CVT)$(dbg): $(CVT).c
	gcc $(DBG) $(PRF) $(WAR) -UREVERT -USORT $(CVT).c -o $(CVT)$(dbg)

$(CVT)_rev$(opt): $(CVT).c
	gcc $(OPT) $(PRF) $(WAR) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(opt)

$(CVT)$(opt): $(CVT).c
	gcc $(OPT) $(PRF) $(WAR) -UREVERT -USORT $(CVT).c -o $(CVT)$(opt)

### convert (old binary to new binary)
$(CVT)_rev_sor$(dbg): $(CVT).c
	gcc $(DBG) $(PRF) $(WAR) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(dbg)

$(CVT)_sor$(dbg): $(CVT).c
	gcc $(DBG) $(PRF) $(WAR) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(dbg)

$(CVT)_rev_sor$(opt): $(CVT).c
	gcc $(OPT) $(PRF) $(WAR) -DREVERT -DSORT $(CVT).c -o $(CVT)_rev_sor$(opt)

$(CVT)_sor$(opt): $(CVT).c
	gcc $(OPT) $(PRF) $(WAR) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(opt)


### run the algo
RED=fflas_reduce
FFLAS=/home/bboyer/usr/
FFLAS_CFLAGS=`sh $(FFLAS)/bin/fflas-ffpack-config --cflags-full`
FFLAS_LIBS=`sh $(FFLAS)/bin/fflas-ffpack-config --libs`
LIBS$(dbg): $(RED).cpp
	g++ $(DBG) $(PRF) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(dbg).o
OBJ$(dbg):
	gcc $(DBG) $(PRF) -c $(RUN).c -o $(RUN)$(dbg).o


LIBS$(opt): $(RED).cpp
	g++ $(OPT) $(PRF) $(FFLAS_CFLAGS)  -c $(RED).cpp  -o $(RED)$(opt).o
OBJ$(opt):
	gcc $(OPT) $(PRF) -c $(RUN).c -o $(RUN)$(opt).o



$(RUN)$(dbg): $(RUN).c
	make LIBS$(dbg) ;
	make OBJ$(dbg)  ;
	gcc $(DBG) $(PRF) $(WAR) $(FFLAS_LIBS) $(RED)$(dbg).o $(RUN).c -o $(RUN)$(dbg) -lm -lstdc++

$(RUN)$(opt): $(RUN).c
	make LIBS$(opt) ;
	make OBJ$(opt)  ;
	gcc $(OPT) $(PRF) $(WAR) $(FFLAS_LIBS) $(RED)$(opt).o $(RUN).c -o $(RUN)$(opt) -lm -lstdc++

clean:
	-rm *.o
	-rm `find . -type f -executable`

