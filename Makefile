all: gsl_bia gsl_bia_inc gsl_dim gsl_dim_out hmt ccc_maker HMT hmt

LIB=-lgsl -lgslcblas -lm
GCC=gcc -O3 -std=c99 -Wall

gsl_bia: gsl_bia.c gsl_bia.h bia.h io.h nbody_bia.h nbody_chk.h
	$(GCC) -o gsl_bia gsl_bia.c $(LIB)

gsl_bia_inc: gsl_bia.c gsl_bia.h bia.h io.h nbody_bia.h nbody_chk.h
	$(GCC) -o gsl_bia_inc gsl_bia.c -D__BIA_INC__ $(LIB)

HMT: gsl_hmt.c io.h nbody_chk.h gsl_hmt.h
	$(GCC) -o HMT gsl_hmt.c $(LIB) -D__HMT__

hmt: gsl_hmt.c io.h nbody_chk.h HermiteIntegrator2.o
	$(GCC) -o hmt gsl_hmt.c HermiteIntegrator2.o $(LIB) -lgfortran

HermiteIntegrator2.o: HermiteIntegrator2.f90
	f95 -O3 -c HermiteIntegrator2.f90

gsl_dim: gsl_dim.c gsl_dim.h io.h nbody_chk.h
	$(GCC) -o gsl_dim gsl_dim.c $(LIB)

gsl_dim_out: gsl_dim.c gsl_dim.h io.h nbody_chk.h
	$(GCC) -o gsl_dim_out gsl_dim.c -DINTERMEDIATE_VALUE $(LIB)

ccc_maker: ccc_maker.c
	$(GCC) -o ccc_maker ccc_maker.c $(LIB)

clean: 
	rm gsl_bia gsl_bia_inc gsl_dim gsl_dim_out ccc_maker hmt HermiteIntegrator2.o HMT
