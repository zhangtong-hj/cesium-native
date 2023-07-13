#include "Cesium3DTilesSelection/Matrix.h"

// Constructors
matrix::matrix(int Row, int Col) {
  this->row = Row;
  this->col = Col;
  array = new float[row * col];

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; ++j) {
      array[i * col + j] = 0;
    }
  }

  r = new int[row];
  for (int i = 0; i < row; ++i) {
    r[i] = i;
  }
  L = NULL;
  U = NULL;
  I = NULL;
}

matrix::matrix(matrix* m) {
  cout << "COPY CONSTRUCTOR" << endl;
  row = m->row;
  col = m->col;
  array = new float[row * col];
  for (int i = 0; i < row * col; ++i) {
    array[i] = m->array[i];
  }

  r = new int[row];
  for (int i = 0; i < row; ++i) {
    r[i] = m->r[i];
  }
  if (m->L != NULL)
    L = new matrix(m->L);
  else
    L = NULL;
  if (m->U != NULL)
    U = new matrix(m->U);
  else
    U = NULL;
  if (m->I != NULL)
    I = new matrix(m->I);
  else
    I = NULL;
}

matrix::matrix(int Size, bool One) {
  this->row = Size;
  this->col = Size;

  array = new float[row * col];
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; ++j) {
      array[i * row + j] = 0;
    }
    if (One)
      array[i * row + i] = 1;
  }

  r = new int[row];
  for (int i = 0; i < row; ++i) {
    r[i] = i;
  }
  L = NULL;
  U = NULL;
  I = NULL;
}

// Destructor
matrix::~matrix() {
  // THIS LINE YIELDS A DOUBLE FREE OR CORRUPTION ERROR
  // delete [] array;
  // delete r;
  if (L != NULL)
    delete (L);
  if (U != NULL)
    delete (U);
  if (I != NULL)
    delete (I);
}

// Accessors
int matrix::rowCount() { return this->row; }
int matrix::colCount() { return this->col; }

float matrix::get(const int irow, const int icol) {
  assert(irow >= 0 && irow < this->row && icol >= 0 && icol < this->col);
  return array[r[irow] * this->col + icol];
}

float& matrix::set(const int irow, const int icol, float value) {
  assert(irow >= 0 && irow < this->row && icol >= 0 && icol < this->col);
  float rnd = round(value);
  float res = value - rnd;
  if (res < EPSILON && res > -EPSILON) {
    value = rnd;
  }
  array[r[irow] * this->col + icol] = value;
  return array[r[irow] * this->col + icol];
}

// Reference Operator
float& matrix::operator()(int irow, int icol) {
  assert(irow >= 0 && irow < this->row && icol >= 0 && icol < this->col);
  return array[r[irow] * this->col + icol];
}

//
void matrix::swap(int r1, int r2) {
  int t = r[r1];
  r[r1] = r[r2];
  r[r2] = t;
}

//
void matrix::clearPermutation() {
  for (int i = 0; i < row; ++i) {
    this->r[i] = i;
  }
}

//
float matrix::i(int rr, int c) {
  assert(rr >= 0 && rr < this->row && c >= 0 && c < this->col);
  if (c == this->r[rr]) {
    return 1;
  } else {
    return 0;
  }
}

// Multiplies the matrix by a vector
matrix matrix::vector_multiply(const std::vector<float> v) {
  assert(this->col == (int)v.size());
  matrix result(this->row, 1);
  int i, j;
  float temp;
  for (i = 0; i < this->row; i++) {
    temp = 0.0;
    for (j = 0; j < this->col; j++)
      temp += get(i, j) * v[j];
    result(i, 0) = temp;
  }
  return result;
}

// Multiply this matrix by another matrix m
matrix matrix::multiply(matrix m) {

  /*
  int n = row;
  int q = n/100;
  int p = q;*/

  matrix K = matrix(row, m.col);
  for (int i = 0; i < row; ++i) {
    /*
    if(q != 0 && i > p)
    {
            p += q;
            int prog = p/q;
            std::cout << "\r" << prog - 1 << "% completed: ";
            std::cout.flush();
    }*/

    for (int j = 0; j < m.col; ++j) {
      float v = 0;
      for (int k = 0; k < col; ++k) {
        v += get(i, k) * m.get(k, j);
      }
      K.set(i, j, v);
    }
  }
  return K;
}

//
void matrix::luFactorize() {
  int n = row;
  float val;

  // Pivoting
  for (int i = 1; i < n; ++i) {
    if (get(i, i) == 0) {
      swap(i - 1, i);
    }
  }

  L = new matrix(row, true);  // Unit Matrix
  U = new matrix(row, false); // Zero Matrix

  for (int p = 0; p < n - 1; ++p) {
    for (int k = 0; k < n; ++k) {
      for (int j = k; j < n; ++j) {
        val = 0;
        for (int m = 0; m < k; ++m) {
          val += L->get(k, m) * U->get(m, j);
        }
        val = get(k, j) - val;
        U->set(k, j, val);
      }
      for (int i = k + 1; i < n; ++i) {
        val = 0;
        for (int m = 0; m < k; ++m) {
          val += L->get(i, m) * U->get(m, k);
        }
        val = (get(i, k) - val) / U->get(k, k);
        L->set(i, k, val);
      }
    }
  }
}

// replaces row r1 with r1-(r2*s)
void matrix::setRow(int r1, int r2, float s) {
  // cout << "O: [" << r1 << "]" << " = [" << r1 << "] - [" << r2 << "] * " << s
  // << endl;
  for (int c = 0; c < col; ++c) {
    set(r1, c, get(r1, c) - (get(r2, c) * s));
    I->set(r1, c, I->get(r1, c) - (I->get(r2, c) * s));
  }
  // printAugmented();
}

// replaces row r1 with r1*s
void matrix::setRow(int rr, float s) {
  // cout << "O: [" << r << "] = [" << r << "] * " << s << endl;
  for (int c = 0; c < col; ++c) {
    set(rr, c, get(rr, c) * s);
    I->set(rr, c, I->get(rr, c) * s);
  }
  // printAugmented();
}

// Inverts a matrix using the gauss-jordan reduction method
void matrix::gaussInversion() {
  I = new matrix(this->row, true);

  int n = this->row;
  float q = float(n) / 100.0;
  float p = q;

  // Pivoting
  for (int i = 1; i < n; ++i) {
    if (get(i, i) == 0) {
      swap(i - 1, i);
      I->swap(i - 1, i);
    }
  }

  for (int i = 0; i < this->row; ++i) {

    if (q != 0 && i > p) {
      p += q;
      int prog = p / q;
      std::cout << "\r" << (prog - 1) << "% completed: ";
      std::cout.flush();
    }

    try {
      setRow(i, 1.0 / get(i, i));
    } catch (std::runtime_error& e) {
      cout << "AN ERROR HAS OCURRED: " << e.what();
      exit(0);
    } catch (int e) {
      std::cout << "AN UNIDENTIFIED ERROR OCURRED: " << e << endl;
      exit(0);
    }

    for (int rr = 0; rr < this->row; ++rr) {
      if (rr != i) {
        setRow(rr, i, get(rr, i));
      }
    }
  }

  cout << endl;
}

// Inverts a matrix. The original matrix is destroyed
// The matrix must be invertible.
void matrix::invert() {
  gaussInversion();
  I->clearPermutation();
}

bool matrix::testIdentity() {
  bool test = true;
  int rr = 0;
  int c = 0;
  while ((rr < this->row) && test) {
    c = 0;
    while ((c < this->col) && test) {
      if (c == rr)
        test = get(rr, c) == 1;
      else
        test = get(rr, c) == 0;
      ++c;
    }
    ++r;
  }
  if (!test) {
    cout << "TEST FAILED ON ROW: " << rr - 1 << ", COL: " << c - 1
         << " WITH VALUE: " << get(rr - 1, c - 1);
  }
  return test;
}

// Prints a matrix
void matrix::print() {
  for (int i = 0; i < this->row; i++) {
    std::cout << "[" << r[i] << "] ";
    for (int j = 0; j < this->col; j++) {
      if (get(i, j) >= 0)
        cout << " ";
      if (get(i, j) == -0.0)
        set(i, j, 0);
      std::cout << get(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

void matrix::printAugmented() {
  std::cout << std::fixed;
  std::cout << std::setprecision(2);

  for (int i = 0; i < this->row; i++) {
    std::cout << "[" << r[i] << "] ";
    for (int j = 0; j < this->col; j++) {
      if (get(i, j) >= 0)
        cout << " ";
      if (get(i, j) == -0.0)
        set(i, j, 0);
      std::cout << get(i, j) << " ";
    }

    cout << "\t : ";

    for (int j = 0; j < this->col; j++) {
      if (I->get(i, j) >= 0)
        cout << " ";
      if (I->get(i, j) == -0.0)
        I->set(i, j, 0);
      std::cout << I->get(i, j) << " ";
    }

    std::cout << std::endl;
  }
}
