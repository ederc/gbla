DBG=-O0 -g
WAR=-Wall -Wextra -ansi
OPT=-march=native -O3
NEW=-D_PROPOSED_FORMAT


TGT=dump_matrix
CVT=converter
RUN=runit
opt=_opt
dbg=_dbg

$(TGT)$(dbg): $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(WAR) $(TGT).c -o $(TGT)$(dbg)

$(TGT)$(opt): $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(WAR) $(TGT).c -o $(TGT)$(opt)

$(TGT)_new$(dbg): $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(NEW) $(WAR) $(TGT).c -o $(TGT)_new$(dbg)

$(TGT)_new$(opt): $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(NEW) $(WAR) $(TGT).c -o $(TGT)_new$(opt)

$(CVT)_revert$(dbg): $(CVT).c
	gcc $(DBG) $(WAR) -DREVERT $(CVT).c -o $(CVT)_revert$(dbg)

$(CVT)$(dbg): $(CVT).c
	gcc $(DBG) $(WAR) -UREVERT $(CVT).c -o $(CVT)$(dbg)

$(CVT)_revert$(opt): $(CVT).c
	gcc $(OPT) $(WAR) -DREVERT $(CVT).c -o $(CVT)_revert$(opt)

$(CVT)$(opt): $(CVT).c
	gcc $(OPT) $(WAR) -UREVERT $(CVT).c -o $(CVT)$(opt)

RED=fflas_reduce
FFLAS=/home/bboyer/usr/
LIBS: $(RED).cpp
	g++ `sh $(FFLAS)/bin/fflas-ffpack-config --cflags` `sh $(FFLAS)/bin/fflas-ffpack-config --libs` -c $(RED).cpp

$(RUN)_dgb: $(RUN).c
	gcc $(DBG) $(WAR) $(RUN).c -o $(RUN)$(dbg)

$(RUN)$(opt): $(RUN).c
	gcc $(OPT) $(AR) $(RUN).c -o $(RUN)$(opt)


clean:
	-rm *.o
	-rm `find . -type f -executable`

