all: gsl_bia gsl_bia_inc gsl_dim gsl_dim_out hmt hmt_fort ccc_maker

LIB=-lgsl -lgslcblas -lm
GCC=gcc -O3 -std=c99 -Wall

gsl_bia: gsl_bia.c gsl_bia.h bia.h io.h nbody_bia.h nbody_chk.h
	$(GCC) -o gsl_bia gsl_bia.c $(LIB)

gsl_bia_inc: gsl_bia.c gsl_bia.h bia.h io.h nbody_bia.h nbody_chk.h
	$(GCC) -o gsl_bia_inc gsl_bia.c -D__BIA_INC__ $(LIB)

hmt: gsl_hmt.c io.h nbody_chk.h gsl_hmt.h
	$(GCC) -o hmt gsl_hmt.c $(LIB) -D__HMT__

hmt_fort: gsl_hmt.c io.h nbody_chk.h HermiteIntegrator2.o
	$(GCC) -o hmt_fort gsl_hmt.c HermiteIntegrator2.o $(LIB) -lgfortran

HermiteIntegrator2.o: HermiteIntegrator2.f90
	gfortran -O3 -c HermiteIntegrator2.f90

gsl_dim: gsl_dim.c gsl_dim.h io.h nbody_chk.h
	$(GCC) -o gsl_dim gsl_dim.c $(LIB)

gsl_dim_out: gsl_dim.c gsl_dim.h io.h nbody_chk.h
	$(GCC) -o gsl_dim_out gsl_dim.c -DINTERMEDIATE_VALUE $(LIB)

ccc_maker: ccc_maker.c
	$(GCC) -o ccc_maker ccc_maker.c $(LIB)

clean:
	rm -f gsl_bia gsl_bia_inc gsl_dim gsl_dim_out ccc_maker hmt hmt_fort HermiteIntegrator2.o
