#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : matrix_(nullptr), rows_(), cols_() {}

S21Matrix::S21Matrix(int rows, int cols) : S21Matrix() {
  if (rows < 1 || cols < 1) {
    throw std::length_error("Matrix size can't be less than 1x1!");
  }
  this->CreateMatrix(rows, cols);
}

// Copy constructor
S21Matrix::S21Matrix(const S21Matrix& other) : S21Matrix() { *this = other; }

// Move constructor
S21Matrix::S21Matrix(S21Matrix&& other) noexcept : S21Matrix() {
  *this = std::move(other);
}

// Destructor
S21Matrix::~S21Matrix() { this->RemoveMatrix(); }

int S21Matrix::getRows() const noexcept { return rows_; }

int S21Matrix::getCols() const noexcept { return cols_; }

void S21Matrix::setRows(int rows) {
  int new_rows = rows_ < rows ? rows_ : rows;
  if (rows < 1) {
    throw std::length_error("Invalid size! Rows can't be less than 1!");
  }

  S21Matrix tmp(rows, cols_);
  for (int i{}; i < new_rows; i++) {
    for (int j{}; j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(tmp);
}

void S21Matrix::setCols(int cols) {
  int new_cols = cols_ < cols ? cols_ : cols;
  if (cols < 1) {
    throw std::length_error("Invalid size! Columns can't be less than 1!");
  }

  S21Matrix tmp(rows_, cols);
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < new_cols; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(tmp);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    RemoveMatrix();
    CreateMatrix(other.rows_, other.cols_);
    for (int i{}; i < rows_; i++) {
      for (int j{}; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    RemoveMatrix();
    CreateMatrix(other.rows_, other.cols_);
    for (int i{}; i < rows_; i++) {
      for (int j{}; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
    other.RemoveMatrix();
  }
  return *this;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool is_equal = true;

  if (!IsEqualDimensions(other)) {
    is_equal = false;
  } else {
    for (int i{}; i < rows_; i++) {
      for (int j{}; j < cols_; j++) {
        if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-07)
          is_equal = false;
      }
    }
  }

  return is_equal;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!IsEqualDimensions(other)) {
    throw std::logic_error("Can't sum matrices with different dimensions!");
  }

  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result{*this};
  return result += other;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!IsEqualDimensions(other)) {
    throw std::logic_error("Can't sub matrices with different dimensions!");
  }

  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result{*this};
  return result -= other;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Incorrect input! Index is out of range!");
  }

  return matrix_[i][j];
}

void S21Matrix::MulNumber(const double number) {
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      matrix_[i][j] *= number;
    }
  }
}

S21Matrix S21Matrix::operator*(const double number) {
  S21Matrix result{*this};
  return result *= number;
}

S21Matrix& S21Matrix::operator*=(const double number) {
  MulNumber(number);
  return *this;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (rows_ != other.cols_ || cols_ != other.rows_) {
    throw std::logic_error("Wrong dimensions for multiplication!");
  }

  S21Matrix result{rows_, other.cols_};
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < other.cols_; j++) {
      for (int k{}; k < other.rows_; k++)
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
    }
  }
  this->RemoveMatrix();
  *this = std::move(result);
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result{*this};
  return result *= other;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

double S21Matrix::Determinant() {
  double result = 0.0;
  if (rows_ != cols_) {
    throw std::logic_error("Wrong dimensions! Matrix must be square!");
  }

  result = GetDeterminant();
  return result;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::logic_error("Wrong dimensions! Matrix must be square!");
  }
  double det = 0.0;
  S21Matrix result(rows_, cols_);
  S21Matrix minor_matrix(rows_ - 1, cols_ - 1);
  for (int i{}; i < rows_; i++) {
    for (int j{}; j < cols_; j++) {
      GetMatrix(i, j, minor_matrix);
      det = minor_matrix.Determinant();
      result.matrix_[i][j] = pow(-1, (i + j)) * det;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) {
    throw std::logic_error("Wrong dimensions! Matrix must be square!");
  }
  double det = Determinant();
  if (std::fabs(det) < 1e-07) {
    throw std::logic_error("Determinant is zero!");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    S21Matrix temp = CalcComplements();
    result = temp.Transpose();
    result.MulNumber(1 / det);
  }
  return result;
}

double S21Matrix::GetDeterminant() {
  double result = 0.0;
  int sign = 1;
  if (rows_ == 1) {
    result = this->matrix_[0][0];
  } else {
    S21Matrix temp(rows_ - 1, cols_ - 1);
    for (int i{}; i < cols_; i++) {
      GetMatrix(0, i, temp);
      result += sign * matrix_[0][i] * temp.GetDeterminant();
      sign = -sign;
    }
  }
  return result;
}

void S21Matrix::GetMatrix(int rows, int cols, const S21Matrix& res) {
  int new_row = 0, new_col = 0;
  for (int i{}; i < rows_; i++) {
    if (i == rows) continue;
    for (int j{}; j < cols_; j++) {
      if (j == cols) continue;
      res.matrix_[new_row][new_col] = matrix_[i][j];
      new_col++;
    }
    new_col = 0;
    new_row++;
  }
}

bool S21Matrix::IsEqualDimensions(const S21Matrix& other) const {
  bool is_equal_dimensions = rows_ == other.rows_ && cols_ == other.cols_;
  return is_equal_dimensions;
}

void S21Matrix::CreateMatrix(int rows, int cols) {
  this->rows_ = rows;
  this->cols_ = cols;
  this->matrix_ = new double* [rows_] {};
  for (int i{}; i < rows_; i++) {
    this->matrix_[i] = new double[cols_]{};
  }
}

void S21Matrix::RemoveMatrix() {
  if (this->matrix_ != nullptr) {
    for (int i{}; i < this->rows_; i++) {
      delete[] this->matrix_[i];
    }
    delete[] this->matrix_;
    this->matrix_ = nullptr;
    this->rows_ = 0;
    this->cols_ = 0;
  }
}
