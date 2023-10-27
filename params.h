#ifndef PARAMS_H
#define PARAMS_H

#include <cmath>
#include <stdexcept>
#include <type_traits>

template <typename T>
class Safe {
    // Solo admite 'float', 'double', o 'long double'
    static_assert(std::is_floating_point<T>::value, "Tipo de template debe ser de punto flotante");

private:
    T value;
    static constexpr T EPSILON = static_cast<T>(1e-5);

public:
    Safe() : value(static_cast<T>(0)) {}
    Safe(T value) : value(value) {}

    // Sobrecarga de operadores de comparación
    bool operator==(const Safe& other) const {
        return std::abs(value - other.value) < EPSILON;
    }
    bool operator==(const T& scalar) const {
        return std::abs(value - scalar) < EPSILON;
    }

    bool operator!=(const Safe& other) const {
        return !(*this == other);
    }
    bool operator!=(const T& scalar) const {
        return !(*this == scalar);
    }

    bool operator<(const Safe& other) const {
        return value < other.value - EPSILON;
    }
    bool operator<(const T& scalar) const {
        return value < scalar - EPSILON;
    }

    bool operator<=(const Safe& other) const {
        return value <= other.value + EPSILON;
    }
    bool operator<=(const T& scalar) const {
        return value <= scalar + EPSILON;
    }

    bool operator>(const Safe& other) const {
        return value > other.value + EPSILON;
    }
    bool operator>(const T& scalar) const {
        return value > scalar + EPSILON;
    }

    bool operator>=(const Safe& other) const {
        return value >= other.value - EPSILON;
    }
    bool operator>=(const T& scalar) const {
        return value >= scalar - EPSILON;
    }

    // Sobrecarga de operadores aritméticos
    Safe operator+(const Safe& other) const {
        return Safe(value + other.value);
    }
    Safe operator-(const Safe& other) const {
        return Safe(value - other.value);
    }
    Safe operator*(const Safe& other) const {
        return Safe(value * other.value);
    }
    Safe operator*(const double& other) const {
        return Safe(value * other);
    }
    Safe operator*(const float& other) const {
        return Safe(value * other);
    }
    Safe operator/(const Safe& other) const {
        if (std::abs(other.value) < EPSILON) {
            throw std::runtime_error("División por cero");
        }
        return Safe(value / other.value);
    }
    Safe operator-() const {
        return Safe(-value);
    }

    // Operadores de asignación
    Safe& operator+=(const Safe& other) {
        value += other.value;
        return *this;
    }
    Safe& operator-=(const Safe& other) {
        value -= other.value;
        return *this;
    }
    Safe& operator*=(const Safe& other) {
        value *= other.value;
        return *this;
    }
    Safe& operator/=(const Safe& other) {
        if (std::abs(other.value) < EPSILON) {
            throw std::runtime_error("División por cero");
        }
        value /= other.value;
        return *this;
    }


    // Métodos para manipulación directa de valores
    T getValue() const {
        return value;
    }
    void setValue(T value) {
        this->value = value;
    }
    
    // Funciones matemáticas
    static Safe abs(const Safe& other) {
        return Safe(std::abs(other.value));
    }
    static Safe sqrt(const Safe& other) {
        if (other.value < 0) {
            throw std::runtime_error("Intento de calcular la raíz cuadrada de un número negativo");
        }
        return Safe(std::sqrt(other.value));
    }
    static Safe pow(const Safe& base, const int& exponent) {
        return Safe(std::pow(base.value, exponent));
    }
    static Safe min(const Safe& a, const Safe& b) {
        return a < b ? a : b;
    }
    static Safe max(const Safe& a, const Safe& b) {
        return a > b ? a : b;
    }

    // Print
    friend std::ostream& operator<<(std::ostream& os, const Safe& other) {
        os << other.value;
        return os;
    }
    friend std::stringstream& operator<<(std::stringstream& os, const Safe& other) {
        os << other.value;
        return os;
    }
    friend std::stringstream& operator>>(std::stringstream& is, Safe& other) {
        is >> other.value;
        return is;
    }

    // Valores especiales
    static Safe max_value() {
        return Safe(std::numeric_limits<T>::max());
    }
    static Safe min_value() {
        return Safe(std::numeric_limits<T>::min());
    }
};

template <typename T>
Safe<T> abs(const Safe<T>& x) {
    return Safe<T>::abs(x);
}
template <typename T>
Safe<T> sqrt(const Safe<T>& x) {
    return Safe<T>::sqrt(x);
}
template <typename T>
Safe<T> pow(const Safe<T>& base, const int& exponent) {
    return Safe<T>::pow(base, exponent);
}
template <typename T>
Safe<T> min(const Safe<T>& a, const Safe<T>& b) {
    return Safe<T>::min(a, b);
}
template <typename T>
Safe<T> max(const Safe<T>& a, const Safe<T>& b) {
    return Safe<T>::max(a, b);
}

using NType = Safe<float>;

namespace Settings {
    inline size_t M = 8;
    inline size_t m = 4;
}


#endif // PARAMS_H