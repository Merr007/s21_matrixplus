#ifndef SRC_S21_MATRIX_OOP_H
#define SRC_S21_MATRIX_OOP_H

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  double** matrix_;
  int rows_, cols_;

  // service methods
  void CreateMatrix(int rows, int cols);
  void RemoveMatrix();
  bool IsEqualDimensions(const S21Matrix& other) const;
  double GetDeterminant();
  void GetMatrix(int rows, int cols, const S21Matrix& res);

 public:
  // constructors and destructor
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  // getters and setters
  int getRows() const noexcept;
  int getCols() const noexcept;
  void setRows(int rows);
  void setCols(int cols);

  // operations
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double number);
  void MulMatrix(const S21Matrix& other);
  double Determinant();
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  // operators
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double number);
  double& operator()(int i, int j);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double number);
  bool operator==(const S21Matrix& other) const;
};

#endif