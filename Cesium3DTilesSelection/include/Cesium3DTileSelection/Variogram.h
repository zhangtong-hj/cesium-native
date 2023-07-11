// Original File: Variogram.h
// Author: Rodolfo Mora Z.
// Date: Jun 8, 2016
// Purpose: Multiple semivariogram models for Kriging Implementation
// https://github.com/CNCA-CeNAT/kriging/blob/master/src/Variogram.h
#pragma once
#ifndef VARIOGRAM_H
#define VARIOGRAM_H

//Standard and Math
#include <stdlib.h>
#include <math.h>
#include "Cesium3DTilesSelection/Vector3.h"
#include "Cesium3DTilesSelection/Matrix.h"
#include <corecrt_math_defines.h>
#define SPHERICAL 0
#define EXPONENTIAL 1
#define GAUSSIAN 2
#define WAVE 3
#define RATIONAL_Q 4
#define CIRCULAR 5

struct VariogramModel {
	double C0; //NUGGET
	double CX; //SILL (C0+C)
	double A;  //RANGE (Max distance to consider v(a)=SILL)
	int VAR;   //Type of variogram to use
};

/**
 * @brief Spherical Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double sphericalVariogram(double H, double C0, double CX, double A);


/**
 * @brief Exponential Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double exponentialVariogram(double H, double C0, double CX, double A);

/**
 * @brief Gaussian Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double gaussianVariogram(double H, double C0, double CX, double A);

/**
 * @brief Wave Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double waveVariogram(double H, double C0, double CX, double A);

/**
 * @brief Rational Quadratic Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double rationalQuadraticVariogram(double H, double C0, double CX, double A);

/**
 * @brief Circular Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double circularVariogram(double H, double C0, double CX, double A);

/** 
  * A, B arrays of points
  *	A/B[i] is the form of {x,y,z}
  *	z might be unknown
  * startA/startB is the first row to calculate, it is assumed that start >= 0
  * endA/endB is the last row to calculate, it is assumed that endA/endB < sizeA/sizeB,
  */
  
/**
 * @brief Calculate the mutual impact between the points of B and the points of A based on their planar xy coordinates
 * @param A sample points in the form of {x,y,z}
 * @param B predictant points in the form of {x,y,z}
 * @param M size of A
 * @param N size of B
 * @param startA first row to calculate, it is assumed that startA >= 0
 * @param endA last row to calculate, it is assumed that endA < M
 * @param startB first row to calculate, it is assumed that startB >= 0
 * @param endB last row to calculate, it is assumed that endB < N
 * @param v semivariogram model to apply
 * @return matrix D with the semivariogram outputs for each pair of points (A,B)
 */
matrix calculateVariogram(
    vector3 A[],
    vector3 B[],
    int M,
    int N,
    int startA,
    int endA,
    int startB,
    int endB,
    VariogramModel v);

/**
 * @brief Calculate the impact impact between a single point A and a block of points B based on their planar xy coordinates  
 * @param A the predictant point
 * @param B the sample points
 * @param M size of B
 * @param D output matrix
 * @param v variogram model
 * @return a matrix of size 1xM with the variogram outputs of the model between A and every point in B.
 */
matrix calculateVariogram(
    vector3 A,
    const vector3 B[],
    int M,
    matrix& D,
    const VariogramModel v);


#endif
