#ifndef NEGPOTENTIAL_H
#define NEGPOTENTIAL_H

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#define MATHLIB_STANDALONE
//#include <Rmath.h>
/**
 * @file negpotential.h
 * @author Naihui Zhou (nzhou@iastate.edu)
 * Function to calculate negpotential function Q
 *
 **/
double negpotential(double *y, double *ystar, int size_y, double **neighbor,
		double alpha, double eta, double tau2);
double H_i(double *y, int i,  double *ystar, int size_y, double **neighbor,
		double alpha, double eta, double tau2);
double H_ij(double *y, int i, int j, double *ystar, int size_y,
		double **neighbor, double alpha, double eta, double tau2);
double vector_dot_product(int n, double *vector1, double *vector2) ;
#endif
