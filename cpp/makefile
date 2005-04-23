#LIBDIR = /usr/lib/gcc-lib/i686-pc-linux-gnu/pgcc-2.91.66:${HOME}/lib/cc
#MATLABLIBDIR=/usr/local/pkg/matlab/bin/glnx86
#MATLABBASE=/cygdrive/c/MATLAB701/
#MATLABLIB=${MATLABBASE}/extern/lib/win32/gcc
#MATLABRUNLIB=${MATLABBASE}/bin/win32
MATLABBASE=/usr/local/pkg/matlab/
MATLABLIB=${MATLABBASE}/bin/glnx86
MATLABRUNLIB=${MATLABBASE}/bin/glnx86
MATLABINCLUDE=${MATLABBASE}/extern/include
LDFLAGS = -L${MATLABLIB} -lmat -lmx -llapack -lblas -lg2c -lc  -Xlinker -rpath -Xlinker ${MATLABRUNLIB} 
LD=g++
CC=g++
FC=g77
INCLUDE = ${MATLABINCLUDE}
CFLAGS = -O3 -I$(INCLUDE) -DLINUX
#CFLAGS = -ggdb -DDBG -I$(INCLUDE) -DLINUX
G77FLAGS = -funroll-all-loops -O3
all: ivm_learn testCode testKern testDist testTransform testMatrix testNoise testIvm

ivm_learn: ivm_learn.o CIvm.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CKern.o CDist.o
	$(LD) ${XLINKERFLAGS} -o ivm_learn  ivm_learn.o CIvm.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CKern.o CDist.o $(LDFLAGS)

ivm_learn.o: ivm_learn.cpp CIvm.h CKern.h CMatrix.h
	$(CC) -c ivm_learn.cpp -o ivm_learn.o $(CFLAGS)

testCode: testCode.o CMatrix.o dtransr.o dgrand.o derf.o CKern.o CDist.o CTransform.o COptimisable.o 
	$(LD) ${XLINKERFLAGS} -o testCode testCode.o CMatrix.o dtransr.o dgrand.o derf.o CKern.o CDist.o CTransform.o COptimisable.o $(LDFLAGS) 

testCode.o: testCode.cpp CKern.h CDist.h CTransform.h CMatrix.h
	$(CC) -c testCode.cpp -o testCode.o $(CFLAGS)

testIvm: testIvm.o CIvm.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CKern.o CDist.o
	$(LD) ${XLINKERFLAGS} -o testIvm  testIvm.o CIvm.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CKern.o CDist.o $(LDFLAGS)

testIvm.o: testIvm.cpp CIvm.h CKern.h CMatrix.h
	$(CC) -c testIvm.cpp -o testIvm.o $(CFLAGS)


testNoise: testNoise.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CDist.o
	$(LD) ${XLINKERFLAGS} -o testNoise  testNoise.o CMatrix.o dtransr.o dgrand.o derf.o CNoise.o ndlutil.o CTransform.o COptimisable.o CDist.o $(LDFLAGS)

testNoise.o: testNoise.cpp CNoise.h CMatrix.h
	$(CC) -c testNoise.cpp -o testNoise.o $(CFLAGS)

testKern: testKern.o CMatrix.o dtransr.o dgrand.o derf.o CKern.o CTransform.o COptimisable.o CDist.o
	$(LD) ${XLINKERFLAGS} -o testKern testKern.o CMatrix.o dtransr.o dgrand.o derf.o CKern.o CTransform.o COptimisable.o CDist.o $(LDFLAGS) 

testKern.o: testKern.cpp CKern.h CDist.h CTransform.h CMatrix.h
	$(CC) -c testKern.cpp -o testKern.o $(CFLAGS)

testDist: testDist.o CMatrix.o dtransr.o dgrand.o derf.o CTransform.o COptimisable.o CDist.o
	$(LD) ${XLINKERFLAGS} -o testDist testDist.o CMatrix.o dtransr.o dgrand.o derf.o CTransform.o COptimisable.o CDist.o $(LDFLAGS) 

testDist.o: testDist.cpp CDist.h CTransform.h CMatrix.h
	$(CC) -c testDist.cpp -o testDist.o $(CFLAGS)

testTransform: testTransform.o CMatrix.o dtransr.o dgrand.o derf.o  CTransform.o 
	$(LD) ${XLINKERFLAGS} -o testTransform testTransform.o CMatrix.o dtransr.o dgrand.o derf.o CTransform.o $(LDFLAGS) 

testTransform.o: testTransform.cpp CTransform.h CMatrix.h
	$(CC) -c testTransform.cpp -o testTransform.o $(CFLAGS)

testMatrix: testMatrix.o CMatrix.o dtransr.o dgrand.o derf.o
	$(LD) ${XLINKERFLAGS} -o testMatrix testMatrix.o CMatrix.o dtransr.o dgrand.o derf.o $(LDFLAGS) 

testMatrix.o: testMatrix.cpp CMatrix.h
	$(CC) -c testMatrix.cpp  -o testMatrix.o $(CFLAGS)

CIvm.o: CIvm.cpp CIvm.h CMatrix.h CMatlab.h COptimisable.h
	$(CC) -c CIvm.cpp -o CIvm.o $(CFLAGS)

CNoise.o: CNoise.cpp CNoise.h COptimisable.h CTransform.h
	$(CC) -c CNoise.cpp -o CNoise.o $(CFLAGS)

CKern.o: CKern.cpp CKern.h CMatrix.h CTransform.h CMatlab.h
	$(CC) -c CKern.cpp -o CKern.o $(CFLAGS)

CTransform.o: CTransform.cpp CTransform.h CMatrix.h
	$(CC) -c CTransform.cpp -o CTransform.o $(CFLAGS)

COptimisable.o: COptimisable.cpp COptimisable.h
	$(CC) -c COptimisable.cpp -o COptimisable.o $(CFLAGS)
CDist.o: CDist.cpp CTransform.h CMatrix.h
	$(CC) -c CDist.cpp -o CDist.o $(CFLAGS)

CMatrix.o: CMatrix.cpp CMatrix.h CMatlab.h ndlutil.h CMatlab.h lapack.h xpose.h
	$(CC) -c CMatrix.cpp -o CMatrix.o $(CFLAGS)

ndlutil.o: ndlutil.cpp ndlutil.h 
	$(CC) -c ndlutil.cpp -o ndlutil.o $(CFLAGS)


# FORTRAN sampling algorithms
dgrand.o: dgrand.f
	$(FC) -c dgrand.f -o dgrand.o $(G77FLAGS)

# FORTRAN in place transpose algorithms
# algorithm 380
dtrans.o: dtrans.f 
	$(FC) -c dtrans.f -o dtrans.o $(G77FLAGS)

# algorithm 467
dxpose.o: dxpose.f 
	$(FC) -c dxpose.f -o dxpose.o $(G77FLAGS)

# algorithm 513
dtransr.o: dtransr.f 
	$(FC) -c dtransr.f -o dtransr.o $(G77FLAGS)

# William Cody's FORTRAN ERF, ERFC and ERFCX function
derf.o: derf.f
	$(FC) -c derf.f -o derf.o $(G77FLAGS) 

clean:
	rm *.o
