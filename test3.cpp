#include <array>
#include <vector>

constexpr int f(std::vector<std::vector<int>>& v) {
    v[1][1] = 2;
    return v[1][1];
}

int main() {
    constexpr int N = 4, M = 3;
    constexpr std::vector<std::vector<int>> v(N, constexpr std::vector<int>(M));
    constexpr const int res = v[1][1];

    //std::array<std::array<int, N>, M> a2;
    //a2[0][0] = 1;
    //a2[2][3] = 1;
}