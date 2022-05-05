#pragma once

#include <algorithm>
#include <array>
#include <vector>
#include "Rational.h"

template<class ValueType = Rational<int64_t>>
class Matrix {
public:
    Matrix() = default;

    Matrix(const std::vector<std::vector<ValueType>> &new_matrix) {
        for (size_t i = 1; i < _matrix.size(); ++i) {
            if (_matrix[i - 1].size() != _matrix[i].size()) {
                throw std::invalid_argument("Matrix rows' sizes are not equal");
            }
        }
        _matrix = new_matrix;
    }

    Matrix(const Matrix &B) {
        _matrix = B._matrix;
    }

    Matrix &operator=(const Matrix &B) {
        _matrix = B._matrix;
        return *this;
    }

    Matrix(Matrix &&B) {
        _matrix = B._matrix;
        B._matrix = std::vector<std::vector<ValueType>>();
    }

    Matrix &operator=(Matrix &&B) {
        _matrix = B._matrix;
        B._matrix = std::vector<std::vector<ValueType>>();
        return *this;
    }

    std::vector<ValueType> &operator[](const size_t i) {
        return _matrix[i];
    }

    const std::vector<ValueType> &operator[](const size_t i) const {
        return _matrix[i];
    }

    Matrix operator*(const Matrix &B) const {
        if (_matrix.empty()) {
            throw std::invalid_argument("Empty _matrix multiplication exception");
        }
        if (B._matrix.empty()) {
            throw std::invalid_argument("Empty _matrix multiplication exception");
        }
        if (this->_matrix[0].size() != B._matrix.size()) {
            throw std::invalid_argument("Wrong matrices sizes A(n1, m1) * B(n2, m2) m1 != n2");
        }
        Matrix C;
        C._matrix.resize(_matrix.size());
        for (size_t i = 0; i < _matrix.size(); ++i) {
            C[i].resize(B[0].size());
        }
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < B._matrix.size(); ++j) {
                for (size_t k = 0; k < B[0].size(); ++k) {
                    C[i][k] += _matrix[i][j] * B[j][k];
                }
            }
        }
        return C;
    }

    Matrix &operator*=(const Matrix &B) {
        Matrix C = *this * B;
        *this = C;
        return *this;
    }

    Matrix &operator*=(const ValueType &number) {
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[i].size(); ++j) {
                _matrix[i][j] *= number;
            }
        }
        return *this;
    }

    Matrix operator*(const ValueType &number) const {
        Matrix C = *this;
        C *= number;
        return C;
    }

    friend Matrix operator*(const ValueType &number, const Matrix &B) {
        return B * number;
    }

    Matrix &operator+=(const Matrix &B) {
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[i].size(); ++j) {
                _matrix[i][j] += B[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix &B) const {
        Matrix C = *this;
        C += B;
        return C;
    }

    Matrix &operator-=(const Matrix &B) {
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[i].size(); ++j) {
                _matrix[i][j] += B[i][j];
            }
        }
        return *this;
    }

    Matrix operator-(const Matrix &B) const {
        Matrix C = *this;
        C -= B;
        return C;
    }

    Matrix operator-() const {
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[i].size(); ++j) {
                _matrix[i][j] = -_matrix[i][j];
            }
        }
        return *this;
    }

    bool operator==(const Matrix &B) const {
        return _matrix == B._matrix;
    }

    std::array<size_t, 2> Sizes() {
        if (_matrix.empty()) {
            return {0, 0};
        } else {
            return {_matrix.size(), _matrix[0].size()};
        }
    }

    void Print() const {
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[0].size(); ++j) {
                _matrix[i][j].Print();
                std::cout << ' ';
            }
            std::cout << '\n';
        }
    }

    /*
     * Adding j-th row to i-th row with lambda coefficient
     */
    void RowAddition(const size_t i, const size_t j, const ValueType &lambda) {
        for (size_t k = 0; k < _matrix[0].size(); ++k) {
            _matrix[i][k] += _matrix[j][k] * lambda;
        }
    }

    /*
     * Adding j-th column to i-th column with lambda coefficient
     */
    void ColumnAddition(const size_t i, const size_t j, const ValueType &lambda) {
        for (size_t k = 0; k < _matrix.size(); ++k) {
            _matrix[k][i] += _matrix[k][j] * lambda;
        }
    }

    void RowSwap(const size_t i, const size_t j) {
        std::swap(_matrix[i], _matrix[j]);
    }

    void ColumnSwap(const size_t i, const size_t j) {
        for (size_t k = 0; k < _matrix.size(); ++k) {
            std::swap(_matrix[k][i], _matrix[k][j]);
        }
    }

    void RowMultiplication(const size_t i, const ValueType &lambda) {
        for (size_t k = 0; k < _matrix[0].size(); ++k) {
            _matrix[i][k] *= lambda;
        }
    }

    void ColumnMultiplication(const size_t i, const ValueType &lambda) {
        for (size_t k = 0; k < _matrix.size(); ++k) {
            _matrix[k][i] *= lambda;
        }
    }

    Matrix Transpose() {
        if (_matrix.empty()) {
            return *this;
        }
        Matrix B;
        B._matrix.resize(_matrix[0].size());
        for (size_t i = 0; i < _matrix[0].size(); ++i) {
            B[i].resize(_matrix.size());
        }
        for (size_t i = 0; i < _matrix.size(); ++i) {
            for (size_t j = 0; j < _matrix[0].size(); ++j) {
                B[j][i] = _matrix[i][j];
            }
        }
        return B;
    }

    /*
     * Output: reduced row echelon form of the _matrix
     */
    Matrix Gauss() const {
        if (_matrix.empty()) {
            return *this;
        }
        Matrix A = *this;
        std::vector<std::pair<size_t, size_t>> main_variables;
        for (size_t first_row = 0, column = 0; first_row < A._matrix.size() && column < A[0].size();) {
            for (size_t row = first_row; row < A._matrix.size(); ++row) {
                if (A[row][column] != 0) {
                    A.RowSwap(row, first_row);
                    break;
                }
            }
            ValueType first_element = A[first_row][column];
            if (first_element == 0) {
                ++column;
                continue;
            }
            main_variables.push_back({first_row, column});
            A.RowMultiplication(first_row, ValueType(1) / first_element);
            for (size_t row = first_row + 1; row < A._matrix.size(); ++row) {
                ValueType lambda = -A[row][column];
                A.RowAddition(row, first_row, lambda);
            }
            ++first_row;
            ++column;
        }
        std::reverse(main_variables.begin(), main_variables.end());
        for (const auto[last_row, column]: main_variables) {
            for (size_t row = 0; row < last_row; ++row) {
                ValueType lambda = -A[row][column];
                A.RowAddition(row, last_row, lambda);
            }
        }
        return A;
    }

    ValueType Determinant() const {
        if (_matrix.empty()) {
            return 0;
        }
        if (_matrix.size() != _matrix[0].size()) {
            throw std::invalid_argument("Determinant of non-square _matrix exception");
        }
        Matrix A = *this;
        ValueType det = 1;
        for (size_t first_row = 0, column = 0; first_row < A._matrix.size() && column < A[0].size();) {
            for (size_t row = first_row; row < A._matrix.size(); ++row) {
                if (A[row][column] != 0) {
                    A.RowSwap(row, first_row);
                    if (row != first_row) {
                        det *= -1;
                    }
                    break;
                }
            }
            ValueType first_element = A[first_row][column];
            if (first_element == 0) {
                ++column;
                continue;
            }
            A.RowMultiplication(first_row, ValueType(1) / first_element);
            det *= first_element;
            for (size_t row = first_row + 1; row < A._matrix.size(); ++row) {
                ValueType lambda = -A[row][column];
                A.RowAddition(row, first_row, lambda);
            }
            ++first_row;
            ++column;
        }
        for (size_t i = 0; i < A._matrix.size(); ++i) {
            det *= A[i][i];
        }
        return det;
    }

    size_t Rank() const {
        if (_matrix.empty()) {
            return 0;
        }
        Matrix A = *this;
        for (size_t first_row = 0, column = 0; first_row < A._matrix.size() && column < A[0].size();) {
            for (size_t row = first_row; row < A._matrix.size(); ++row) {
                if (A[row][column] != 0) {
                    A.RowSwap(row, first_row);
                    break;
                }
            }
            ValueType first_element = A[first_row][column];
            if (first_element == 0) {
                ++column;
                continue;
            }
            A.RowMultiplication(first_row, ValueType(1) / first_element);
            for (size_t row = first_row + 1; row < A._matrix.size(); ++row) {
                ValueType lambda = -A[row][column];
                A.RowAddition(row, first_row, lambda);
            }
            ++first_row;
            ++column;
        }
        size_t rank = 0;
        for (size_t i = 0; i < A._matrix.size(); ++i) {
            if (A[i][i] != 0) {
                rank += 1;
            }
        }
        return rank;
    }

    /*
     * Jacobi method of diagonalizing symmetrical bilinear form
     * Input: symmetrical bilinear form's _matrix
     * Output: pair(normalized _matrix, transition _matrix)
     */
    std::pair<Matrix, Matrix> Jacobi() const {
        Matrix E;
        Matrix A = *this;
        if (A._matrix.empty()) {
            return {A, E};
        }
        if (A._matrix.size() != A[0].size()) {
            throw std::invalid_argument("Jacobi diagonalizing non-square _matrix exception");
        }
        for (size_t i = 0; i < A._matrix.size(); ++i) {
            for (size_t j = i + 1; j < A._matrix.size(); ++j) {
                if (A[i][j] != A[j][i]) {
                    throw std::invalid_argument("Jacobi diagonalizing asymmetrical _matrix exception");
                }
            }
        }
        E._matrix.resize(A._matrix.size());
        for (size_t i = 0; i < A._matrix.size(); ++i) {
            E[i].resize(A._matrix.size(), ValueType(0));
            E[i][i] = 1;
        }
        for (size_t iteration = 0; iteration < A._matrix.size(); ++iteration) {
            size_t first_row = iteration;
            size_t column = iteration;
            ValueType first_element = A[first_row][column];
            if (first_element == 0) {
                for (size_t i = first_row + 1; i < A._matrix.size(); ++i) {
                    if (A[i][i] != 0) {
                        A.RowSwap(first_row, i);
                        E.RowSwap(first_row, i);
                        A.ColumnSwap(first_row, i);
                        break;
                    }
                }
                first_element = A[first_row][column];
                if (first_element == 0) {
                    for (size_t row = first_row + 1; row < A._matrix.size(); ++row) {
                        if (A[row][column] != 0) {
                            A.RowAddition(first_row, row, ValueType(1));
                            E.RowAddition(first_row, row, ValueType(1));
                            A.ColumnAddition(first_row, row, ValueType(1));
                            break;
                        }
                    }
                }
            }
            first_element = A[first_row][column];
            if (first_element == 0) {
                continue;
            }
            for (size_t row = first_row + 1; row < A._matrix.size(); ++row) {
                ValueType lambda = -A[row][column] / first_element;
                A.RowAddition(row, first_row, lambda);
                E.RowAddition(row, first_row, lambda);
                A.ColumnAddition(row, first_row, lambda);
            }
        }
        return {A, E.Transpose()};
    }

    /*
     * Orthogonalization of column _matrix space using Gram-Schmidt process
     * Input: vectors in columns of the input _matrix
     * Output: orthogonal vectors in columns of the output _matrix
     */
    Matrix Orthogonalization() const {
        Matrix A = *this;
        if (_matrix.empty()) {
            return A;
        }
        for (size_t i = 1; i < A[0].size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                ValueType a_norm = 0;
                ValueType dot_product = 0;
                for (size_t k = 0; k < _matrix.size(); ++k) {
                    a_norm += A[k][j] * A[k][j];
                    dot_product += A[k][j] * _matrix[k][i];
                }
                if (a_norm == 0) {
                    continue;
                }
                for (size_t k = 0; k < A._matrix.size(); ++k) {
                    A[k][i] -= dot_product / a_norm * A[k][j];
                }
            }
        }
        return A;
    }

private:
    std::vector<std::vector<ValueType>> _matrix;
};
