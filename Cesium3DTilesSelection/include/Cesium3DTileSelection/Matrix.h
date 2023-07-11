// Original File: Matrix.cxx
// Author: Jharrod LaFon
// Date: Spring 2011
// Purpose: Simple Matrix Class
// https://github.com/hpc/MPI-Examples/blob/master/newton-type-optimization/matrix.cxx

// Modified by: Rodolfo Mora
// Date: February 2015
// Purpose: Parallel implementation of Matrix Multiplication and Matrix Invertion algorithms
// https://github.com/CNCA-CeNAT/kriging/blob/master/src/Matrix.h
#pragma once
#ifndef MATRIX_CXX
#define MATRIX_CXX

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

using std::cout;
using std::endl;

#define EPSILON 0.00001

// Matrix class definition
class matrix
{
public:
	// Constructors
	matrix(int Row = 10,  int Col = 10);
	// Copies another
	matrix(matrix * m);
	// Constructs a square matrix 
	matrix(int Size = 10, bool One = false);

	// Destructor
	~matrix();

	// Accessors
	int rowCount();
	int colCount();

	// Print
	void print();

	// Print A:I
	void printAugmented();

	// Returns the identity value for this pivoted matrix on the required row and column
	float i(int r, int c);

	// Swaps two rows of the matrix
	void swap(int r1, int r2);

	// Clears permutation of a matrix, returning it to it's original row order
	void clearPermutation();

	// Applies LU Factorization to matrix
	void luFactorize();

	// Solves
	void solveLU(int c);

	// Replaces row r1 with r1-r2*s
	void setRow(int r1, int r2, float s);

	// Replaces row r with r*s
	void setRow(int r, float s);

	// Inverts a matrix using the gauss jordan algorithm
	void gaussInversion();

	// Returns the inverse of a matrix
	void invert();

	// Returns true if this is an identity matrix
	bool testIdentity();

	// Returns this matrix multiplied by another matrix
	matrix multiply(matrix m);

	// Swaps two matrix rows
	void swap_rows(const  int r1, const  int r2);

	// Multiplies the matrix by a vector
	matrix vector_multiply(const std::vector<float> v);

	// Reference Operator
	float  get(const  int row, const  int col);
	float & operator() (const  int row, const  int col);
	// Assignation Operator
	float & set(const  int row, const  int col, float value);

	matrix * L;
	matrix * U;
	matrix * I;
private:
	int row, col;
	int * r;
	float * array;
};



#endif
