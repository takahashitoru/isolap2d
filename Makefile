#OBJ = test_isolap2d.o
OBJ = main.o

OBJ += envs.o opts.o double2.o gaussone.o closedcurve.o \
	input.o ld.o isolap2d.o uinc.o msolve.o output.o \
	quadtree.o mknf.o mktree.o chtree.o fmm.o bincoef.o \
	dxf.o \
	circle.o hexa.o cshape.o sshape.o \
	options.o incal.o

OBJ += nbie.o # normal boundary integral equation
##OBJ += hnbie.o
##OBJ += cbie.o # conventional boundary integral equation

OBJ += myintde2.o newtoncotes.o

CC = gcc # a.nuem.nagoya-u.ac.jp
#CC = /home/ttaka/gcc-4.6.1/bin/gcc # a.nuem.nagoya-u.ac.jp
#GCC_LIB = -L/home/ttaka/gcc-4.6.1/lib64
#GCC_INCLUDE = -L/home/ttaka/gcc-4.6.1/include

#CC = gcc44 # z.nuem.nagoya-u.ac.jp

OBJ += timer.o
ELAPSED_OBJ = elapsed.o
ELAPSED_FLAGS = -O0 -std=gnu99 # GCC
#ELAPSED_FLAGS = -O0 -std=c99 # ICC

#ELAPSED_FLAGS += -DCPU_CLOCK_GHZ=2.20 # Sakabe's machine (Dell Vostro 200)
ELAPSED_FLAGS += -DCPU_CLOCK_GHZ=3.07 # hpctech (192.168.1.197)

GMRES_HOME = .
GMRES =
#OBJ += dgmres110609.o
OBJ += dgmres120116.o

LIB = ${GMRES} ${GCC_LIB} -lm
INCLUDE = -I${GMRES_HOME} ${GCC_INCLUDE}

OPTS = -std=c99 -Wall # always
OPTS += -g -DMYDEBUG # for debug
#OPTS += -O3 -msse2 # for optimization
#OPTS += -O3 -mtune=core2 -mssse3 ##### for optimization my PC
#OPTS += -O3 -msse3 # for optimization
#OPTS += -O3 -march=corei7-avx -mtune=corei7-avx # for optimization
OPTS += -O3 -mavx # for optimization
#OPTS += -O3 -mtune=native -march=native # for optimization
#OPTS += -O3 -msse4_2

OPTS += -fopenmp # OpenMP (GCC)
LIB += -lgomp -lpthread #OpenMP

OPTS += -DCIRCLE # Solve a test problem for a circular boundary

OPTS += -DCIRCLE_NUM_CLOSED_CURVES=1 # this must be one if CIRCLE is enabled

#OPTS += -DCIRCLE_SPLINE_DEGREE=2
#OPTS += -DCIRCLE_NUM_CONTROL_POINTS=100
P = 2 # default
CP = 100 # default
OPTS += -DCIRCLE_SPLINE_DEGREE=${P} ##### P is given by work/go.sh #####
OPTS += -DCIRCLE_NUM_CONTROL_POINTS=${CP} ##### CP is given by work/go.sh #####


#OPTS += -DHEXA # hexagonal boundaries

NG = 9 # default
OPTS += -DNGAUSS=${NG} ##### NG is given by work/go.sh #####

ITMAX=1000 #default
OPTS += -DITMAX=${ITMAX}

TOL=1e-5 #default
OPTS += -DTOL=${TOL}

#TOL=1e-2 #default
#OPTS += -DTOL=${TOL}
#OPTS += -DDGMRES_RATE_CONVERGENCE
#OPTS += -DDGMRES_KEEP_RESIDUALS=3

METHOD=CONV # default
#METHOD=FMM
OPTS += -D${METHOD}
OPTS += -DFUKUI # default for FMM

MAXEPC=30 #default
NBETA=1 #default
NTERM=20 #default
OPTS += -DMAXEPC=${MAXEPC} # only for FMM
OPTS += -DNBETA=${NBETA} # only for FMM
OPTS += -DNTERM=${NTERM} # only for FMM

KIND=_DUMMY_ # default, i.e., NBIE
#KIND=ORIGINAL
#KIND=CBIE
OPTS += -D${KIND}

FLAG=-D_DUMMY_ # default
OPTS += ${FLAG} ##### FLAG is given by work/go.sh #####


#OPTS += -DSET_FREE_TERM # N21, p12: Sakabe did not use this option in his thesis


#OPTS += -DMYINTDE # enables DE formula for nbie
#OPTS += -DMYINTDE_EPS=1e-8 # 1e-15 is default
#OPTS += -DMYINTDE_LENAW=1000 # 8000 is default

#OPTS += -DADAPTIVE_GAUSS
#OPTS += -DADAPTIVE_GAUSS_NMIN=4
#OPTS += -DADAPTIVE_GAUSS_NMAX=24
#OPTS += -DADAPTIVE_GAUSS_EPS=1e-5

#OPTS += -DNEWTONCOTES # N21, p16

CFLAGS = ${INCLUDE} ${OPTS}

TAR = tmp.tar
EXE = a.out

${EXE}: ${OBJ} ${ELAPSED_OBJ}
	${CC} ${CFLAGS} -o ${EXE} ${OBJ} ${ELAPSED_OBJ} ${LIB}
${OBJ}:
${ELAPSED_OBJ}:
	${CC} ${ELAPSED_FLAGS} -c elapsed.c
clean : 
	rm -f *.o core
tar backup:
	tar cvfh ${TAR} Makefile *.c *.h *.sh work/go.sh
	gzip -f ${TAR}
diff:
	(./diff.sh ../120201/ ${OBJ:.o=.c} *.h)
