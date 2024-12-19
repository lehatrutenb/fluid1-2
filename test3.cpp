#include <array>
#include <vector>
#include <iostream>

template<bool t>
struct A{
void test() {
    if constexpr (t) {
        std::cout << 1;
    } else {
        std::cout << 2;
    }
}
};
int main() {
    A<false>{}.test();
    std::cout << std::endl;
}