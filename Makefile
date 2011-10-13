#
#ILAを使用している場合、全エネルギーの計算がうまくいかないので、-DDEBUGは正しく動作しない。
#
#CC=icc
#CFLAGS  = -I/usr/include/mpi -Dreal=float -DMPI -O3 -Zp16 -tpp7 -xW -cxxlib-icc -static-libcxa # -Wall -Werror -g -pg# -DDEBUG# -DDEBUG
CFLAGS  = -O6 -Wall -DDOUBLEPRECISION #-Werror #-I/usr/include/mpi -DMPI #-g -pg -DSINGLEPRECISION # -DDEBUG# -DILA -DDEBUG 
LDFLAGS = -lm
MPICC   = mpicc.mpich
MPILDFLAGS = -lmpi
all: visualize search2 # protonRXMC pdice

search2: LattIce2.o SearchLib.o search.o PairList.o IntHash.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

visualize: search2
	ln -f search2 visualize

%.v.o: %.c
	$(CC) -c $(CFLAGS) -DVISUALIZE $< -o $@
#cluster
%.cl.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.d.o: %.c
	$(CC) -c $(CFLAGS) -DDEBUG $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) -DILA $< -o $@

ilacheck.%: search
	./search --lattice=iceIh50-30-30.ar3a --check=1 $* > $@


depend:
	makedepend *.c *.h
clean:
	-rm *.o
	cp Makefile Makefile.tmp
	#The following line must be the last line.
	sed -n -e '1,/--8x--END OF MAKEFILE--8x--/p' Makefile.tmp > Makefile
