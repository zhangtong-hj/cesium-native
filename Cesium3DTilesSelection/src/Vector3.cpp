#include "Cesium3DTilesSelection/Vector3.h"


vector3::vector3() {
  x = 0;
  y = 0;
  z = 0;
}

vector3::vector3(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

vector3::vector3(vector3* v) {
  this->x = v->x;
  this->y = v->y;
  this->z = v->z;
}

vector3 vector3::operator+(vector3 other) {
  return vector3(x + other.x, y + other.y, z + other.z);
}

vector3 vector3::operator-() { return vector3(-x, -y, -z); }

vector3 vector3::operator-(vector3 other) {
  return vector3(x - other.x, y - other.y, z - other.z);
}

vector3 vector3::operator*(double scalar) {
  return vector3(x * scalar, y * scalar, z * scalar);
}

double vector3::operator*(vector3 other) {
  return x * other.x + y * other.y + z * other.z;
}

string vector3::toString() {
  std::stringstream s;
  s.precision(2);
  s << "(" << x << "," << y << "," << z << ")" << endl;
  return s.str();
}

void vector3::print() {}

double vector3::magnitude() {
  double x2 = x * x;
  double y2 = y * y;
  double z2 = z * z;
  double p = x2 + y2 + z2;
  return sqrt(p);
}
