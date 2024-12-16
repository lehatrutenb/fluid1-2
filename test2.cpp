#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include <functional>

using namespace std;
//#define TYPES FIXED(10, 10),FIXED(20, 20)
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

//typename T, 
template<std::size_t N, std::size_t K, bool F>
struct Fixed {
    using T = std::conditional<N <= 8,  int8_t, std::conditional<N <= 16,  int16_t, int16_t>>::type;
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

template<std::size_t N, std::size_t K, bool F>
std::ostream& operator<<(std::ostream &out, const Fixed<N,K,F>& f) {
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

//const std::vector<Fixed_t> types = getTypes();

template<std::size_t N, std::size_t K, bool F>
void setType(auto& func, int n, int k, bool f) {
    if (n == N && k == K && f == F) {
        func.template operator()<Fixed<N,K,F>>();
    }
}

void FixedSetType(size_t n, size_t k, bool f, auto& func) {
    //f.template operator()<Fixed<10, 10, false>>();
    /*for (auto& t : types) {
        if (N == t.N && K == t.K && F == t.F) {
            f.template operator()<Fixed<t.N, t.K, t.F>>();
            return;
        }
        if (N == 10 && K == 10) {
            f.template operator()<Fixed<int, 10, 10>>();
        }
    }*/
   //#define FIXED(N, K) if (n == N && k == K && !f) { func.template operator()<Fixed<N,K,false>>();}
   //TYPES;
   //#undef FIXED
   //#define FIXED(N, K) if (n == N && k == K && !f) { func.template operator()<Fixed<N,K,false>>();}
   //#undef FIXED,
    #define FIXED(N, K) setType<N,K,false>(f, n, k, f)
    TYPES;
    #undef FIXED
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
        FixedSetType(N, K, false, f);
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
}