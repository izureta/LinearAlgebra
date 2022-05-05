#include "Matrix.h"

namespace LinearAlgebra {
    template<class ValueType = Rational<int64_t>>
    Matrix<ValueType> Eye(size_t n) {
        std::vector<std::vector<ValueType>> e;
        e.resize(n);
        for (size_t i = 0; i < n; ++i) {
            e[i].resize(n, ValueType(0));
            e[i][i] = 1;
        }
        return Matrix(e);
    }

    template<class ValueType = Rational<int64_t>>
    Matrix<ValueType> Zero(size_t n) {
        std::vector<std::vector<ValueType>> zero;
        zero.resize(n);
        for (size_t i = 0; i < n; ++i) {
            zero[i].resize(n, ValueType(0));
        }
        return Matrix(zero);
    }
}
