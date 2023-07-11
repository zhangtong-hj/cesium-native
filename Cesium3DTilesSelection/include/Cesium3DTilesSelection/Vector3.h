// Original File: Vector3.h
// Author: Rodolfo Mora Z.
// Date: May 19, 2015
// Purpose: Simple TriDimentional Vector Class
// https://github.com/CNCA-CeNAT/kriging/blob/master/src/Vector3.h
#pragma once
#ifndef VECTOR3_H_
#define VECTOR3_H_

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>

using std::endl;
using namespace std;

class vector3
{
public:
	double x;
	double y;
	double z;

	// Constructors
	//Default
	vector3();

	vector3(double x, double y, double z);
	// Copies another
	vector3(vector3 * v);

	// Default values
	static vector3 zero(){return vector3(0,0,0);};
	static vector3 forward(){return vector3(0,0,1);};
	static vector3 up(){return vector3(0,1,0);};
	static vector3 right(){return vector3(1,0,0);};
	static vector3 one(){return vector3(1,1,1);};

	// Printing
	string toString();
	void print();

	vector3 operator+(vector3 other);
	vector3 operator-();
	vector3 operator-(vector3 other);
	vector3 operator*(double scalar);
	double operator*(vector3 other);

	double magnitude();
	static double distance(vector3 a, vector3 b){return (a-b).magnitude();};
};


#endif /* VECTOR3_H_ */
