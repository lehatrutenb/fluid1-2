#pragma GCC optimize("Ofast,unroll-loops")
#include <cstdint>
#include <limits>
#include <utility>
#include <cassert>
#include <random>
#include <iostream>
#include <ranges>
#include <algorithm>
#include <cstring>
#define TYPES FAST_FIXED(20, 10), FIXED(20, 10),FIXED(20, 20),FLOAT

using namespace std;

//constexpr size_t N = 36, M = 84; // размеры поля ? TODO
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000; // кол-во тиков которые будут сделаны
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}}; // куда из каждой точки попытаемся потечь?
constexpr long double eps = 1e-6;

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

mt19937 rnd(1303);

template<std::size_t N, std::size_t K, bool F>
struct Fixed {
    using T_ = std::conditional<N <= 8, int8_t, typename std::conditional<N <= 16, int16_t, typename std::conditional<N <= 32, int32_t, int64_t>::type>::type>::type;
    using T2_ = std::conditional<N <= 8, int16_t, typename std::conditional<N <= 16, int32_t, int64_t>::type>::type;

    using TF_ = std::conditional<N <= 8, int_fast8_t, typename std::conditional<N <= 16, int_fast16_t, typename std::conditional<N <= 32, int_fast32_t, int_fast64_t>::type>::type>::type;
    using TF2_ = std::conditional<N <= 8, int_fast16_t, typename std::conditional<N <= 16, int_fast32_t, int_fast64_t>::type>::type;

    using T = std::conditional<F, TF_, T_>::type;
    using T2 = std::conditional<F, TF2_, T2_>::type;

    constexpr Fixed(int v_): v(v_ << K) {}
    constexpr Fixed(float f_): v(f_ * (1 << K)) {}
    constexpr Fixed(double f_): v(f_ * (1 << K)) {}
    template<std::size_t No, std::size_t Ko, bool Fo>
    constexpr Fixed(Fixed<No, Ko, Fo> o): v(o.to_double()*(1<<K)) {}
    constexpr Fixed(): v(0) {}
    //operator double() const { return to_double(); }
    //operator float() const { return to_double(); }

    constexpr double to_double() const {
        return v / (double) (1 << K);
    }

    static constexpr Fixed from_raw(T x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }
    T v;

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
auto operator<=>(const Fixed<N,K,F>& a, const double& b) {
    if (a == b) {
        return 0;
    } else {
        return a.to_double() > b ? 1 : -1;
    }
}

template<std::size_t N, std::size_t K, bool F>
bool operator==(const Fixed<N,K,F>& a, const double& b) {
    return abs(a.to_double() - b) < eps;
}

template<std::size_t N, std::size_t K, bool F, std::size_t No, std::size_t Ko, bool Fo>
auto operator<=>(const Fixed<N,K,F>& a, const Fixed<No,Ko,Fo>& b) {
    if (a == b) {
        return 0;
    } else {
        return a.to_double() > b.to_double() ? 1 : -1;
    }
}

template<std::size_t N, std::size_t K, bool F, std::size_t No, std::size_t Ko, bool Fo>
bool operator==(const Fixed<N,K,F>& a, const Fixed<No,Ko,Fo>& b) {
    return abs(a.to_double() - b.to_double()) < eps;
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N,K,F> abs(Fixed<N,K,F> f) {
    if (f.v < 0) {
        f.v = -f.v;
    }
    return f;
}

template<std::size_t N, std::size_t K, bool F>
constexpr void tryConv(Fixed<N,K,F> x, double& res) {
    res = x.to_double();
}

template<std::size_t N, std::size_t K, bool F>
constexpr void tryConv(Fixed<N,K,F> x, float& res) {
    res = x.to_double();
}

template<std::size_t N, std::size_t K, bool F>
constexpr void tryConv(double x, Fixed<N,K,F>& res) {
    res = Fixed<N,K,F>(x);
}

constexpr void tryConv(double x, float& res) {
    res = x;
}

constexpr void tryConv(double x, double& res) {
    res = x;
}

template<std::size_t N, std::size_t K, bool F, std::size_t No, std::size_t Ko, bool Fo>
constexpr void tryConv(Fixed<No,Ko,Fo> x, Fixed<N,K,F>& res) {
    res = Fixed<N,K,F>(x);
}


template<std::size_t N, std::size_t K, bool F>
constexpr double operator+(double o, Fixed<N,K,F> t) {
    return t.to_double() + o;
}

template<std::size_t N, std::size_t K, bool F>
constexpr double operator-(double o, Fixed<N,K,F> t) {
    return o - t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
constexpr double operator*(double o, Fixed<N,K,F> t) {
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
float& operator+=(float& o, Fixed<N,K,F> t) {
    return o = t.to_double() + o;
}

template<std::size_t N, std::size_t K, bool F>
double& operator-=(double& o, Fixed<N,K,F> t) {
    return o = o - t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
float& operator-=(float& o, Fixed<N,K,F> t) {
    return o = o - t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double& operator*=(double& o, Fixed<N,K,F> t) {
    return o = o * t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
double operator/=(double& o, Fixed<N,K,F> t) {
    return o = o / t.to_double();
}

template<std::size_t N, std::size_t K, bool F>
std::ostream& operator<<(std::ostream &out, const Fixed<N,K,F>& f) {
    return out << f.v / (double) (1 << K);
}


template<std::size_t N, std::size_t M, typename PT>
struct FieldStatic {
    FieldStatic(std::size_t _, std::size_t __){}
    constexpr static char f[N][M];
    static PT point[N][M]{};
    static PT old_point[N][M];
    constexpr const static std::size_t N_ = N;
    constexpr const static std::size_t M_ = M;
    constexpr const static bool isStatic = true;
};


template<typename PT>
struct FieldDinamic {
    FieldDinamic(std::size_t N__, std::size_t M__) : N_(N__), M_(M__) {

        f.resize(0);
        f.resize(N_, std::vector<char> (M_));
        point.resize(N_, std::vector<PT> (M_));
        // assign to point zeroes???? TODO CHECK
        old_point.resize(N_, std::vector<PT> (M_));
    }
    FieldDinamic(){}
    static std::vector<std::vector<char>> f;
    static std::vector<std::vector<PT>> point;
    static std::vector<std::vector<PT>> old_point;
    const std::size_t N_;
    const std::size_t M_;
    constexpr const static bool isStatic = false;
};

template<typename PT, typename VT, typename VFT, typename FT, std::size_t NStatic, std::size_t MStatic>
struct fluidEmulator {

//std::size_t N, std::size_t M;
//const static auto& field = FT::f;
#define isStatic FT::isStatic
#define field FT::f
#define point_ FT::point
#define old_point FT::old_point
#define N FT{}.N_
#define M FT{}.M_

/*
char field[N][M + 1] = { // само поле M + 1 хз почему
    "####################################################################################",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                       .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............################                     #                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "##################################################################                 #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "####################################################################################",
};
*/
VT rho[256]; // какие-то константы? TODO

template<typename T>
struct VectorField { // changed Fixed->T
    //array<T, deltas.size()> v[N][M]; // TODO вероятно для каждой клетки какой-то поток в разные стороны
    //using TStatic = std::array<std::array<std::array<T, deltas.size()>, M>, N>;
    using TDinamic = std::vector<std::vector<array<T, deltas.size()>>>;
    using TArr = std::conditional<isStatic, std::array<std::array<std::array<T, deltas.size()>, MStatic>, NStatic>, TDinamic>::type;

    template<bool isStatic_, enable_if<!isStatic_>>
    VectorField() {
        v.resize(N, std::vector<array<T, deltas.size()>>(M));
    }

    template<bool isStatic_, enable_if<isStatic_>>
    VectorField() {

    }

    TArr v;

    T &add(int x, int y, int dx, int dy, T dv) {
        return get(x, y, dx, dy) += dv;
    }

    T &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

//VectorField<VT> testt = VectorField<VT> ();
VectorField<VT> velocity<isStatic>{};
VectorField<VFT> velocity_flow = VectorField<VFT>(); // TODO ?

using LUDinamic = std::vector<std::vector<int>>;
using LUArr = std::conditional<isStatic,std::array<std::array<int, MStatic>, NStatic>, LUDinamic>::type;

//int last_use[N][M]{}; // базовый used для dfs
LUArr last_use;

template<enable_if<!isStatic>>
void initLastUse() {
    last_use.resize(N, std::vector<int>(M));
}

template<enable_if<isStatic>>
constexpr void initLastUse() {
    for (int i = 0; i < NStatic; i++) {
        for (int j = 0; j < MStatic; j++) {
            last_use[i][j] = 0;
        }
    }
}

int UT = 0;  // видимо переменная для dfs (или не проверял что там) чтобы не запускаться из одной вершины много раз

tuple<VFT, bool, pair<int, int>> propagate_flow(int x, int y, VFT lim) { // координаты рассм клетки и lim TODO?
    last_use[x][y] = UT - 1; // помечаем клетку посещённой
    VFT ret = 0;
    for (auto [dx, dy] : deltas) { // перебираем соседние клетки к нашей
        int nx = x + dx, ny = y + dy; // коорд соседних клеток
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT) { // если рассм клетка не препятствие и не была помечена раньше то идём внутрь
            VT cap = velocity.get(x, y, dx, dy);
            VFT flow = velocity_flow.get(x, y, dx, dy);
            if (flow == cap) { // TODO?
                continue;
            }
            // assert(v >= velocity_flow.get(x, y, dx, dy));
            VFT vp;
            if (lim <= (cap - flow)) {
                vp = lim;
            } else {
                tryConv(cap-flow, vp);
            }

            if (last_use[nx][ny] == UT - 1) { // видимо если была подсчитама в пред разы то пересчёт простой TODO?
                velocity_flow.add(x, y, dx, dy, vp); // добавляем к вытеканию сколько вытечет в напр dx dy TODDO?
                last_use[x][y] = UT; // помечаем изначальную клетку рассмотренной
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                return {vp, 1, {nx, ny}}; // TODO?
            }
            auto [t, prop, end] = propagate_flow(nx, ny, vp); // идём в соседнюю клетку рескурсивно 
            ret += t; // TODO?
            if (prop) {
                velocity_flow.add(x, y, dx, dy, t);
                last_use[x][y] = UT; // помечаем изначальную клетку рассмотренной
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                return {t, prop && end != pair(x, y), end}; // TODO?
            }
        }
    }
    last_use[x][y] = UT; // помечаем изначальную клетку рассмотренной
    return {ret, 0, {0, 0}}; // TODO?
}

double random01() {
    return rnd()/(double)std::numeric_limits<uint_fast32_t>::max();
}

// прикол что помечает какие-то клетки рассмотренным, но какие TODO ? 
void propagate_stop(int x, int y, bool force = false) { // коорд текущей клетки и force TODO?
    if (!force) { // TODO ?
        bool stop = true;
        for (auto [dx, dy] : deltas) { // перебираем соседей
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > VT(0)) { // если есть хотя бы одна клетка из которой ещё не всё вылилось (пересчиталось) что должно TODO ?
                stop = false; // типо все рассмотрели - говорим что можно остановиться
                break;
            }
        }
        if (!stop) { // если в цикле решили остановиться - останавливаемся
            return;
        }
    }
    last_use[x][y] = UT; // помечаем клетку рассмотренной если <= velocity.get(x, y, dx, dy) или это бордюрчик или уже рассмотренна на последнем тике
    for (auto [dx, dy] : deltas) { // итерируемся по соседним клеткам
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > VT(0)) { // если бордюрчик/уже помечена/ещё не всё вытекло ?TODO - не рассматриваем
            continue;
        }
        propagate_stop(nx, ny); // рекурсивно смотрим соседнюю
    }
}

VT move_prob(int x, int y) { // координаты клетки
    VT sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) { // если бордюр или уже рассмотренна - не трогаем
            continue;
        }
        auto v = velocity.get(x, y, dx, dy); // смотрим сколько чего-то а ней TODO?
        if (v < VT(0)) {
            continue;
        }
        sum += v; // суммируем чего-то в клетке, если положительное число TODO ?
    }
    return sum;
}

struct ParticleParams { // класс чтобы инициализировать/менять клетки поля TODO?
    char type;
    PT cur_p;
    array<VT, deltas.size()> v;

    void swap_with(auto point__, VectorField<VT>& velocity, int x, int y) {
        swap(field[x][y], type);
        swap(point__[x][y], cur_p);
        swap(velocity.v[x][y], v);
    }
};

bool propagate_move(int x, int y, bool is_first) { // клетки поля и is_first - вероятно корневая ли вершина TODO?
    last_use[x][y] = UT - is_first; // TODO?
    bool ret = false;
    int nx = -1, ny = -1;
    do {
        std::array<VT, deltas.size()> tres;
        VT sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) { // если бордюр или уже помеченная то TODO?
                tres[i] = sum;
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < VT(0)) {
                tres[i] = sum;
                continue;
            }
            sum += v;
            tres[i] = sum;
        }

        if (sum == VT(0)) {
            break;
        }

        double p = random01() * sum;
        size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

        auto [dx, dy] = deltas[d];
        nx = x + dx;
        ny = y + dy;
        assert(velocity.get(x, y, dx, dy) > VT(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

        ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < VT(0)) { // если не бордюр, не была ещё рассмотренная и что-то в ней есть TODO?
            propagate_stop(nx, ny);
        }
    }
    if (ret) { // TODO?
        if (!is_first) { // TODO?
            ParticleParams pp{}; // TODO?
            pp.swap_with(point_, velocity, x, y);
            pp.swap_with(point_, velocity, nx, ny);
            pp.swap_with(point_, velocity, x, y);
        }
    }
    return ret;
}

//double dirs[N][M]{}; // TODO?

using DDinamic = std::vector<std::vector<double>>;
using DArr = std::conditional<isStatic,std::array<std::array<double, MStatic>, NStatic>, LUDinamic>::type;

DArr dirs;

template<enable_if<!isStatic>>
void initDirs() {
    dirs.resize(N, std::vector<double>(M));
}

template<enable_if<isStatic>>
constexpr void initDirs() {
    for (int i = 0; i < NStatic; i++) {
        for (int j = 0; j < MStatic; j++) {
            dirs[i][j] = 0;
        }
    }
}

void run() {
    initLastUse();
    initDirs();
    rho[' '] = 0.01; // задаём константы
    rho['.'] = 1000;
    //VT g = 0.1; // типо g физическая
    VT g = 1;

    for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
            if (field[x][y] == '#')
                continue;
            for (auto [dx, dy] : deltas) {
                dirs[x][y] += (field[x + dx][y + dy] != '#');
            }
        }
    }

    for (size_t i = 0; i < T; ++i) { // итерируемся по тикам
        
        PT total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                if (field[x + 1][y] != '#')
                    velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(old_point, point_, sizeof(point_));
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && old_point[nx][ny] < old_point[x][y]) {
                        PT delta_p = old_point[x][y] - old_point[nx][ny];
                        PT force = delta_p;
                        VT &contr = velocity.get(nx, ny, -dx, -dy);
                        if (contr * rho[(int) field[nx][ny]] >= force) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        
                        VT in_;
                        tryConv((force) / rho[(int) field[x][y]], in_);
                        velocity.add(x, y, dx, dy, in_);
                        point_[x][y] -= force / VT(dirs[x][y]);
                        total_delta_p -= force / VT(dirs[x][y]);
                    }
                }
            }
        }

        // Make flow from velocities
        velocity_flow = {};
        bool prop = false;
        do {
            UT += 2;
            prop = 0;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        auto [t, local_prop, _] = propagate_flow(x, y, 1);
                        if (t > VFT(0)) {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    auto old_v = velocity.get(x, y, dx, dy);
                    auto new_v = velocity_flow.get(x, y, dx, dy);
                    if (old_v > VT(0)) {
                        assert(new_v <= old_v);
                        VT in_;
                        tryConv(new_v, in_);
                        velocity.get(x, y, dx, dy) = in_;
                        auto force = (old_v - new_v) * rho[(int) field[x][y]];
                        if (field[x][y] == '.')
                            force *= 0.8;
                        if (field[x + dx][y + dy] == '#') {
                            point_[x][y] += force / PT(dirs[x][y]);
                            total_delta_p += force / PT(dirs[x][y]);
                        } else {
                            point_[x + dx][y + dy] += force / PT(dirs[x + dx][y + dy]);
                            total_delta_p += force / PT(dirs[x + dx][y + dy]);
                        }
                    }
                }
            }
        }

        UT += 2;
        prop = false;
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] != '#' && last_use[x][y] != UT) {
                    //std::cout << random01() << ' ' << move_prob(x, y) << std::endl;
                    if (random01() < move_prob(x, y)) {
                        prop = true;
                        propagate_move(x, y, true);
                    } else {
                        propagate_stop(x, y, true);
                    }
                }
            }
        }

        if (prop) {
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < N; ++x) {
                cout << field[x] << "\n";
            }
        }
    }
}

#undef point_
#undef old_point
#undef field
#undef N
#undef M

};

template<typename FT>
void parseField(const char* fieldFile) {
    freopen(fieldFile, "r", stdin);
    
    std::size_t N, M;
    std::cin >> N >> M;
    FT{}(N, M);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            std::cin >> FT::f[i][j];
        }
    }

    fclose(stdin);
}

namespace typeSetter {

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

void setType(char* type_b, std::string_view type, auto func) {
    std::size_t n = std::numeric_limits<size_t>::max(), k = n;
    bool f;
    if (type != "FLOAT" && type != "DOUBLE") {
        int res = sscanf(type_b, "FIXED(%zu, %zu)", &n, &k);
        if (res == 0) {
            int res = sscanf(type_b, "FAST_FIXED(%zu, %zu)", &n, &k);
            if (res == 0) {
                throw std::runtime_error("error");
                return;
            }
        }
    }
    #define FIXED(N, K) setTypeFixed<N,K,false>(func, n, k, f)
    #define FAST_FIXED(N, K) setTypeFixed<N,K,true>(func, n, k, f)
    #define FLOAT setTypeFloat(func, type)
    #define DOUBLE setTypeDouble(func, type)
    TYPES;
    #undef FIXED
    #undef FLOAT
    #undef DOUBLE
}


template<std::size_t N, std::size_t M, typename T1, typename T2, typename T3>
void setNumber(std::size_t CN, std::size_t CM, std::string fieldFile) {
    if (CN == N && CM == M) {
        parseField<FieldStatic<N, M, T1>>(fieldFile);
        fluidEmulator<T1, T2, T3, FieldStatic<N, M, T1>, N, M>{}.run();
    }
}

}

int main() {
    auto first_type = "FAST_FIXED(20, 10)";
    char first_type_b[] = "FAST_FIXED(20, 10)";
    auto second_type = "FAST_FIXED(20, 10)";
    char second_type_b[] = "FAST_FIXED(20, 10)";
    auto third_type = "FAST_FIXED(20, 10)";
    char third_type_b[] = "FAST_FIXED(20, 10)";

    std::size_t CN = 36;
    std::size_t CM = 84;

    const char* fieldFile = "field";

    typeSetter::setType(first_type_b, first_type, [&]<typename T1>() {
        typeSetter::setType(second_type_b, second_type, [&]<typename T2>() {
            typeSetter::setType(third_type_b, third_type, [&]<typename T3>() {
                #define S(N, M) typeSetter::setNumber<N,M>(CN, CM, field)
                #undef S

                parseField<FieldDinamic<T1>>(fieldFile);
                fluidEmulator<T1, T2, T3, FieldDinamic<T1>, 0, 0>{}.run();
            });
        });
    });
}