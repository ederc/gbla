CC=gcc -std=c11
#  CC=clang
CXX=g++ 
#  CXX=clang++
DBG=-march=native -mtune=native -O0 -g
WAR=-Wall -Wextra -ansi
STA=
#  STA=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive 
OPT=-march=native -mtune=native -O3 -DNDEBUG -UDEBUG  -flto
NEW=-D_PROPOSED_FORMAT
LIB=-lm
#  PRF=-pg
PRF=

PAR=-fopenmp

RUN=runit

opt=_opt
dbg=_dbg
par=_par
seq=_seq
ffl=_ffl

all:
	make $(RUN)$(seq)$(opt)
	make $(RUN)$(par)$(opt)

TETES=matrix.h field_ops.h selecter.h tools.h config.h types.h macros.h printer.h sparseops.h sparseops_block.h io.h reduce_B.h reduce_C.h echelon.h ouvrir.h


### run the algo
RED=fflas_reduce
FFLAS=~/usr/
FFLAS_CFLAGS=`sh $(FFLAS)/bin/fflas-ffpack-config --cflags-full`
FFLAS_LIBS=`sh $(FFLAS)/bin/fflas-ffpack-config --libs`

# SEQ

LIBS$(seq)$(dbg): $(RED).cpp
	$(CXX) $(DBG) $(PRF) $(FFLAS_CFLAGS)  $(STA)  -c $(RED).cpp  -o $(RED)$(seq)$(dbg).o

LIBS$(seq)$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF) $(FFLAS_CFLAGS)  $(STA)  -c $(RED).cpp  -o $(RED)$(seq)$(opt).o

OBJ$(seq)$(dbg): 
	$(CC)  $(DBG) $(PRF) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(seq)$(dbg).o

OBJ$(seq)$(opt): 
	$(CC)  $(OPT) $(PRF) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(seq)$(opt).o

$(RUN)$(seq)$(dbg): $(RUN).c $(TETES)
	make LIBS$(seq)$(dbg) ;
	make OBJ$(seq)$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(RED)$(seq)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(seq)$(dbg) -lm -lstdc++

$(RUN)$(seq)$(opt): $(RUN).c $(TETES)
	make LIBS$(seq)$(opt) ;
	make OBJ$(seq)$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(RED)$(seq)$(opt).o  $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(seq)$(opt) -lm -lstdc++



# PAR

LIBS$(par)$(dbg): $(RED).cpp
	$(CXX) $(DBG) $(PRF) $(PAR) $(FFLAS_CFLAGS)  $(STA) -c $(RED).cpp  -o $(RED)$(par)$(dbg).o

LIBS$(par)$(opt): $(RED).cpp
	$(CXX) $(OPT) $(PRF) $(PAR) $(FFLAS_CFLAGS)  $(STA) -c $(RED).cpp  -o $(RED)$(par)$(opt).o

OBJ$(par)$(dbg):
	$(CC)  $(DBG) $(PRF) $(PAR) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(par)$(dbg).o

OBJ$(par)$(opt):
	$(CC)  $(OPT) $(PRF) $(PAR) $(WAR) $(STA) -c $(RUN).c -o $(RUN)$(par)$(opt).o

$(RUN)$(par)$(dbg):  $(RUN).c $(TETES)
	make LIBS$(par)$(dbg) ;
	make OBJ$(par)$(dbg)  ;
	$(CC) $(DBG) $(PRF) $(WAR) $(PAR) $(STA)  $(RED)$(par)$(dbg).o $(FFLAS_LIBS)  $(RUN).c -o $(RUN)$(par)$(dbg) -lm -lstdc++


$(RUN)$(par)$(opt):  $(RUN).c $(TETES)
	make LIBS$(par)$(opt) ;
	make OBJ$(par)$(opt)  ;
	$(CC) $(OPT) $(PRF) $(WAR) $(PAR) $(STA) $(RED)$(par)$(opt).o $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(par)$(opt) -lm -lstdc++



# FFLAS only

OBJ$(par)$(ffl):
	$(CC) $(OPT) $(PRF) -DONLY_FFLAS $(WAR) $(PAR) $(STA) -c $(RUN).c -o $(RUN)$(par)$(ffl).o

$(RUN)$(par)$(ffl): 
	make LIBS$(par)$(opt) ;
	make OBJ$(par)$(ffl)  ;
	$(CC) $(OPT) $(PRF)  -DONLY_FFLAS $(WAR) $(PAR) $(STA) $(RED)$(par)$(opt).o  $(FFLAS_LIBS) $(RUN).c -o $(RUN)$(par)$(ffl) -lm -lstdc++


clean:
	-rm *.o
	-rm `find . -type f -executable`

check:
	-diff types.h    ../tools/types.h
	-diff macros.h   ../tools/macros.h
	-diff selecter.h ../tools/selecter.h
	-diff ouvrir.h   ../tools//ouvrir.h

