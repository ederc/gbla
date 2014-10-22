DBG=-O0 -g
WAR=-Wall -Wextra -ansi
OPT=-O3
TGT=dump_matrix
NEW=-D_PROPOSED_FORMAT

dump-dbg: $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(WAR) $(TGT).c -o $(TGT)_dbg

dump-opt: $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(WAR) $(TGT).c -o $(TGT)_opt

new-dump-dbg: $(TGT).c $(TGT).h printer.h
	gcc $(DBG) $(NEW) $(WAR) $(TGT).c -o $(TGT)_dbg_new

new-dump-opt: $(TGT).c $(TGT).h printer.h
	gcc $(OPT) $(NEW) $(WAR) $(TGT).c -o $(TGT)_opt_new
