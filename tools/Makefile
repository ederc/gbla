#  CC=gcc
#  CC=clang
#  CXX=g++
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
	make $(DMP)_new
	make $(DMP)_old
	make $(CVT)_new


### dump matrix (to sms...)
$(DMP)_old$(dbg): $(DMP).c $(DMP).h 
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)_old$(dbg)

$(DMP)_old: $(DMP).c $(DMP).h 
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(DMP).c -o $(DMP)_old

$(DMP)_new$(dbg): $(DMP).c $(DMP).h
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(DMP).c -o $(DMP)_new$(dbg)

$(DMP)_new: $(DMP).c $(DMP).h 
	$(CC) $(OPT) $(PRF) $(NEW) $(WAR) $(STA) $(DMP).c -o $(DMP)_new

### convert (old binary to new binary)
$(CVT)_rev$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(dbg)

$(CVT)$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -UREVERT -USORT $(CVT).c -o $(CVT)$(dbg)

$(CVT)_rev$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -DREVERT -USORT $(CVT).c -o $(CVT)_rev$(opt)

$(CVT)$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -UREVERT -USORT $(CVT).c -o $(CVT)$(opt)

$(CVT)_new$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_new$(dbg)

$(CVT)_sor$(dbg): $(CVT).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(dbg)

$(CVT)_new: $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -DREVERT -DSORT $(CVT).c -o $(CVT)_new

$(CVT)_sor$(opt): $(CVT).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) -UREVERT -DSORT $(CVT).c -o $(CVT)_sor$(opt)



clean:
	-rm *.o
	-rm `find . -type f -executable`


