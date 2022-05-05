#pragma once

#include <iostream>

template<class ValueType = int64_t>
class Rational {

public:
    Rational() {}

    Rational(const ValueType &numerator, const ValueType &denominator) : _numerator(numerator),
                                                                         _denominator(denominator) {
        Reduce();
    }

    Rational(const ValueType &number) {
        _numerator = number;
        _denominator = 1;
    }

    Rational(const Rational &number) {
        _numerator = number._numerator;
        _denominator = number._denominator;
    }

    Rational &operator=(const Rational &number) {
        _numerator = number._numerator;
        _denominator = number._denominator;
        return *this;
    }

    Rational(Rational &&number) {
        _numerator = number._numerator;
        _denominator = number._denominator;
    }

    Rational &operator=(Rational &&number) {
        _numerator = number._numerator;
        _denominator = number._denominator;
        return *this;
    }

    Rational &operator*=(const Rational &number) {
        _numerator *= number._numerator;
        _denominator *= number._denominator;
        Reduce();
        return *this;
    }

    Rational operator*(const Rational &number) const {
        Rational product = *this;
        product *= number;
        return product;
    }

    Rational &operator/=(const Rational &number) {
        if (number._numerator == 0) {
            throw std::overflow_error("Division by zero exception");
        }
        _numerator *= number._denominator;
        _denominator *= number._numerator;
        Reduce();
        return *this;
    }

    Rational operator/(const Rational &number) const {
        Rational product = *this;
        product /= number;
        return product;
    }

    Rational &operator+=(const Rational &number) {
        _numerator = _numerator * number._denominator + _denominator * number._numerator;
        _denominator *= number._denominator;
        Reduce();
        return *this;
    }

    Rational operator+(const Rational &number) const {
        Rational product = *this;
        product += number;
        return product;
    }

    Rational &operator-=(const Rational &number) {
        _numerator = _numerator * number._denominator - _denominator * number._numerator;
        _denominator *= number._denominator;
        Reduce();
        return *this;
    }

    Rational operator-(const Rational &number) const {
        Rational product = *this;
        product -= number;
        return product;
    }

    Rational operator-() const {
        Rational number = *this;
        number._numerator *= -1;
        return number;
    }

    bool operator==(const Rational &number) const {
        return (_numerator == number._numerator) && (_denominator == number._denominator);
    }

    bool operator!=(const Rational &number) const {
        return !((*this) == number);
    }

    void Print() const {
        if (_denominator == 1 || _numerator == 0) {
            std::cout << _numerator;
        } else {
            std::cout << _numerator << '/' << _denominator;
        }
    }

private:
    ValueType GCD(ValueType a, ValueType b) {
        if (b == 0) {
            return a;
        }
        return GCD(b, a % b);
    }

    void Reduce() {
        ValueType g = GCD(_numerator, _denominator);
        if (g == 0) {
            throw std::overflow_error("Division by zero exception");
        }
        _numerator /= g;
        _denominator /= g;
        if (_denominator < 0) {
            _numerator *= -1;
            _denominator *= -1;
        }
    }

    ValueType _numerator = 0;
    ValueType _denominator = 1;
};
