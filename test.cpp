//#include <map>
#include <string>
#include <iostream>
#include <cstdio>
#define TYPES FIXED(10, 10),FIXED(20, 20)
#define COMMA ,
#define TCOMMA ;
#define COMMA TCOMMA
TYPES
#define FIXED(N, K) 

using namespace std;

template<typename T, std::size_t N, std::size_t K>
struct Fixed {
    constexpr Fixed(int v): v(v << K) {}
    constexpr Fixed(float f): v(f * (1 << K)) {}
    constexpr Fixed(double f): v(f * (1 << K)) {}
    constexpr Fixed(): v(0) {}

    static constexpr Fixed from_raw(int32_t x) {
        Fixed ret;
        ret.v = x;
        return ret;
    } 

    T v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    Fixed operator+(Fixed b) {
        return Fixed::from_raw(v + b.v);
    }

    Fixed& operator=(Fixed b) {
        v = b.v;
        return *this;
    }

    Fixed operator-(Fixed b) {
        return Fixed::from_raw(v - b.v);
    }
};

template<typename T, std::size_t N, std::size_t K>
std::ostream& operator<<(std::ostream &out, const Fixed<T,N,K>& f) {
    return out << f.v / (double) (1 << K);
}

template<typename T1, typename T2, typename T3>
void test_func(T1 a, T2 b, T3 c) {
    a = a + a;
    b = b + b;
    c = c + c;
    std::cout << a <<  ' ' << b << ' ' << c << std::endl;
}


//constexpr const size_t possibleTypes[1][2] = {{10, 10}};

void FixedSetType(size_t N, size_t K, auto f) {
    f.template operator()<Fixed<int, 10, 10>>();
    /*for (int ind = 0; ind < 1; ind++) {
        if (N == possibleTypes[ind][0] && K == possibleTypes[ind][1]) {
            f.template operator()<Fixed<int8_t, N, K>>();
            return;
        }
        if (N == 10 && K == 10) {
            f.template operator()<Fixed<int, 10, 10>>();
        }
    }*/
}

void setType(char* type_b, std::string_view type, auto f) {
    if (type == "FLOAT") {
        f.template operator()<float>();
    } else if (type == "DOUBLE") {
        f.template operator()<double>();
    } else {
        size_t N, K;
        int res = sscanf(type_b, "FIXED(%zu, %zu)", &N, &K);
        if (res == 0) {
            throw std::runtime_error("error");
            //scanf(type, "FAST_FIXED(%u, %u)", &N, &K);
            //f.template operator()<FastFixed<N, K>>();
            return;
        }
        FixedSetType(N, K, f);
        ////(type == "fixed") {
        //f.template operator()<Fixed<10, 10>>();
    }
}

int main() {
    auto first_type = "FIXED(10, 10)";
    char first_type_b[] = "FIXED(10, 10)";
    auto second_type = "FIXED(10, 10)";
    char second_type_b[] = "FIXED(10, 10)";
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
    return 0;
}