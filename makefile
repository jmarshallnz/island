objects = main.o cluster.o cluster4.o cluster5.o cluster6.o
myutils = cmatrix.h DNA.h lotri_matrix.h \
	matrix.h mydouble.h myerror.h myutils.h \
	pause.h random.h sort.h tsv.h utils.h vector.h
headers = cluster.h $(myutils)
idir = "./"

## INSTRUCT MAKE TO MAKE MULTIPLE GOALS
.PHONY : all
all : isource

## LINK
isource : $(objects)
	gcc -w -O3 -o isource $(objects) -lstdc++

## COMPILE
main.o : main.cpp $(headers)
	gcc -w -O3 -c -o main.o main.cpp -I$(idir) -D _MODEL6

cluster.o : cluster.cpp $(headers)
	gcc -w -O3 -c -o cluster.o cluster.cpp -I$(idir) -D _MODEL6

cluster4.o : cluster4.cpp $(headers)
	gcc -w -O3 -c -o cluster4.o cluster4.cpp -I$(idir) -D _MODEL6

cluster5.o : cluster5.cpp $(headers)
	gcc -w -O3 -c -o cluster5.o cluster5.cpp -I$(idir) -D _MODEL6

cluster6.o : cluster6.cpp $(headers)
	gcc -w -O3 -c -o cluster6.o cluster6.cpp -I$(idir) -D _MODEL6

## MAKE CLEAN
clean :
	-rm $(objects)
