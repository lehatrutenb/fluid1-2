#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include <functional>
#include <limits>

using namespace std;
#define TYPES FIXED(10, 10),FIXED(20, 20),FLOAT
/*
class Fixed_w {
public:
    virtual void setFType(auto func) = 0;
};

template<std::size_t N, std::size_t K, bool F>
class Fixed_wc : Fixed_w {
    Fixed_wc(){};
    setFType(auto func) override {
        f.template operator()<Fixed<N,K,F>>();
    }
};


std::array<, 2>


constexpr std::vector<Fixed_t> getTypes() {
    std::vector<Fixed_t> res;
    std::string s = "TYPES";
    for (int i = 0; i < s.size(); i++) {
        if (s[i] == '(') {
            std::string cur = "";
            std::size_t n, k;
            for (int j = i + 1; j < s.size(); j++) {
                if (s[j] == ' ') {
                    n = stoi(cur);
                    cur = "";
                }
                if (s[j] == ')') {
                    k = stoi(cur);
                    res.push_back({n, k, true});
                    res.push_back({n, k, false});
                    Fixed<n,k, false> res;
                    break;
                }
            }
        }
    }
    return res;
}*/

template<std::size_t N, std::size_t K, bool F>
struct Fixed {
    using T = std::conditional<N <= 8,  int8_t, typename std::conditional<N <= 16,  int16_t, typename std::conditional<N <= 32,  int32_t, int64_t>::type>::type>::type;
    using T2 = std::conditional<N <= 8,  int16_t, typename std::conditional<N <= 16,  int32_t, int64_t>::type>::type;
    constexpr Fixed(int v_): v(v_ << K) {}
    constexpr Fixed(float f_): v(f_ * (1 << K)) {}
    constexpr Fixed(double f_): v(f_ * (1 << K)) {}
    template<std::size_t No, std::size_t Ko, bool Fo>
    constexpr Fixed(Fixed<No, Ko, Fo> o): v(o.to_double()*(1<<K)) {}
    constexpr Fixed(): v(0) {}

    constexpr double to_double() {
        return v / (double) (1 << K);
    }

    static constexpr Fixed from_raw(T x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }
    T v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    //static constexpr Fixed inf = Fixed::from_raw(std::numeric_limits<T>::max());
    //static constexpr Fixed eps = Fixed::from_raw(deltas.size());

    Fixed operator+(Fixed o) {
        return from_raw(v + o.v);
    }
    Fixed operator-(Fixed o) {
        return from_raw(v - o.v);
    }
    Fixed operator*(Fixed o) {
        return Fixed::from_raw(((T2) v * o.v) >> K);
    }
    Fixed operator/(Fixed o) {
        return Fixed::from_raw(((T2) v << K) / o.v);
    }
    Fixed& operator=(Fixed o) {
        this->v = o.v;
        return *this;
    }
    Fixed &operator+=(Fixed o) {
        return *this = *this + o;
    }

    Fixed &operator-=(Fixed o) {
        return *this = *this - o;
    }

    Fixed &operator*=(Fixed o) {
        return *this = *this * o;
    }

    Fixed &operator/=(Fixed o) {
        return *this = *this / o;
    }

    Fixed operator-() {
        return Fixed::from_raw(-v);
    }
};

template<std::size_t N, std::size_t K, bool F>
Fixed<N,K,F> abs(Fixed<N,K,F> f) {
    if (f.v < 0) {
        f.v = -f.v;
    }
    return f;
}

template<std::size_t N, std::size_t K, bool F>
double operator+(double o, Fixed<N,K,F> t) {
    return t.to_double() + o;
}

template<std::size_t N, std::size_t K, bool F>
double operator-(double o, Fixed<N,K,F> t) {
    return o - t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double operator*(double o, Fixed<N,K,F> t) {
    return o * t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double operator/(double o, Fixed<N,K,F> t) {
    return o / t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double& operator+=(double& o, Fixed<N,K,F> t) {
    return o = t.to_double() + o;
}

template<std::size_t N, std::size_t K, bool F>
double& operator-=(double o, Fixed<N,K,F> t) {
    return o = o - t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double& operator*=(double o, Fixed<N,K,F> t) {
    return o = o * t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double operator/(double& o, Fixed<N,K,F> t) {
    return o = o / t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
std::ostream& operator<<(std::ostream &out, const Fixed<N,K,F>& f) {
    return out << f.v / (double) (1 << K);
}

//typename T, 
/*
template<std::size_t N, std::size_t K, bool F>
struct Fixed {
    using T = std::conditional<N <= 8,  int8_t, typename std::conditional<N <= 16,  int16_t, typename std::conditional<N <= 32,  int32_t, int64_t>::type>::type>::type;
    constexpr Fixed(int v_): v(v_ << K) {}
    constexpr Fixed(float f_): v(f_ * (1 << K)) {}
    constexpr Fixed(double f_): v(f_ * (1 << K)) {}
    template<std::size_t No, std::size_t Ko, bool Fo>
    constexpr Fixed(Fixed<No, Ko, Fo> o): v(o.to_double()*(1<<K)) {}
    constexpr Fixed(): v(0) {}

    constexpr double to_double() {
        return v / (double) (1 << K);
    }

    static constexpr Fixed from_raw(T x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    T v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    //static constexpr Fixed inf = Fixed::from_raw(std::numeric_limits<int32_t>::max());
    Fixed operator+(Fixed o) {
        return from_raw(o.v + v);
    }
    Fixed& operator=(Fixed o) {
        this->v = o.v;
        return *this;
    }
};


template<std::size_t N, std::size_t K, bool F>
double operator+(double o, Fixed<N,K,F> t) {
    return t.to_double() + o;
}

template<std::size_t N, std::size_t K, bool F>
std::ostream& operator<<(std::ostream &out, const Fixed<N,K,F>& f) {
    return out << f.v / (double) (1 << K);
}*/


/*
template <std::size_t N, std::size_t K, bool F>
constexpr Fixed<N,K,F>::Fixed(double f) {//: v(64 * (1 << 32))
}

template<>
template <std::size_t N=32, std::size_t K=32, bool F=false>
Fixed<N,K,F>::Fixed(double f) {//: v(64 * (1 << 32))
}*/

//template<>
//constexpr Fixed<64, 32, false>::Fixed<64, 32, false>(double f): v(f * (1 << K)) {}


//constexpr const size_t possibleTypes[1][2] = {{10, 10}};

//const std::vector<Fixed_t> types = getTypes();

template<std::size_t N, std::size_t K, bool F>
void setTypeFixed(auto& func, int n, int k, bool f) {
    if (n == N && k == K && f == F) {
        func.template operator()<Fixed<N,K,F>>();
    }
}

void setTypeFloat(auto& func, std::string_view type) {
    if (type == "FLOAT") {
        func.template operator()<float>();
    }
}

void setTypeDouble(auto& func, std::string_view type) {
    if (type == "DOUBLE") {
        func.template operator()<double>();
    }
}

void setType(char* type_b, std::string type, auto func) {
    std::size_t n = std::numeric_limits<size_t>::max(), k = n;
    bool f;
    if (type != "FLOAT" && type != "DOUBLE") {
        int res = sscanf(type_b, "FIXED(%zu, %zu)", &n, &k);
        if (res == 0) {
            throw std::runtime_error("error");
            return;
        }
    }
    #define FIXED(N, K) setTypeFixed<N,K,false>(func, n, k, f)
    #define FLOAT setTypeFloat(func, type)
    #define DOUBLE setTypeDouble(func, type)
    TYPES;
    #undef FIXED
    #undef FLOAT
    #undef DOUBLE
}

template<typename T1, typename T2, typename T3>
void test_func(T1 a, T2 b, T3 c) {
    a /= a + b;
    b -= b + abs(b);
    c = c + c;
    std::cout << a <<  ' ' << b << ' ' << c << std::endl;
}

int main() {
    //Fixed<10, 10, false> f = 1;
    //f = f + 2.3;
    auto first_type = "FIXED(10, 10)";
    char first_type_b[] = "FIXED(10, 10)";
    auto second_type = "FIXED(20, 20)";
    char second_type_b[] = "FIXED(20, 20)";
    auto third_type = "FLOAT";
    char third_type_b[] = "FLOAT";
    setType(first_type_b, first_type, [&]<typename T1>() {
        setType(second_type_b, second_type, [&]<typename T2>() {
            setType(third_type_b, third_type, [&]<typename T3>() {
                T1 a = T1(1);
                T2 b = T2(2.2);
                T3 c = 3.5;
                test_func(a, b, c);
            });
        });
    });
}