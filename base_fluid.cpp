#include <cstdint>
#include <limits>
#include <utility>
#include <cassert>
#include <random>
#include <iostream>
#include <ranges>
#include <algorithm>
#include <cstring>
#include <chrono>

#define NOT_TO_PRINT

using namespace std;

constexpr size_t N = 36, M = 84; // размеры поля ? TODO
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1000;//1'000'000; // кол-во тиков которые будут сделаны
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}}; // куда из каждой точки попытаемся потечь?

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

struct Fixed {
    constexpr Fixed(int v): v(v << 16) {}
    constexpr Fixed(float f): v(f * (1 << 16)) {}
    constexpr Fixed(double f): v(f * (1 << 16)) {}
    constexpr Fixed(): v(0) {}

    static constexpr Fixed from_raw(int32_t x) {
        Fixed ret;
        ret.v = x;
        return ret;
    } 

    int32_t v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
};

static constexpr Fixed inf = Fixed::from_raw(std::numeric_limits<int32_t>::max());
static constexpr Fixed eps = Fixed::from_raw(deltas.size());

Fixed operator+(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v + b.v);
}

Fixed operator-(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v - b.v);
}

Fixed operator*(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v * b.v) >> 16);
}

Fixed operator/(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v << 16) / b.v);
}

Fixed &operator+=(Fixed &a, Fixed b) {
    return a = a + b;
}

Fixed &operator-=(Fixed &a, Fixed b) {
    return a = a - b;
}

Fixed &operator*=(Fixed &a, Fixed b) {
    return a = a * b;
}

Fixed &operator/=(Fixed &a, Fixed b) {
    return a = a / b;
}

Fixed operator-(Fixed x) {
    return Fixed::from_raw(-x.v);
}

Fixed abs(Fixed x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}

ostream &operator<<(ostream &out, Fixed x) {
    return out << x.v / (double) (1 << 16);
}

Fixed rho[256]; // какие-то константы? TODO

Fixed p[N][M]{}, old_p[N][M]; // ? TODO

struct VectorField {
    array<Fixed, deltas.size()> v[N][M]; // TODO вероятно для каждой клетки какой-то поток в разные стороны
    Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
        return get(x, y, dx, dy) += dv;
    }

    Fixed &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

VectorField velocity{}, velocity_flow{}; // TODO ?
int last_use[N][M]{}; // базовый used для dfs
int UT = 0;  // видимо переменная для dfs (или не проверял что там) чтобы не запускаться из одной вершины много раз


mt19937 rnd(1337);

tuple<Fixed, bool, pair<int, int>> propagate_flow(int x, int y, Fixed lim) { // координаты рассм клетки и lim TODO?
    last_use[x][y] = UT - 1; // помечаем клетку посещённой
    Fixed ret = 0;
    for (auto [dx, dy] : deltas) { // перебираем соседние клетки к нашей
        int nx = x + dx, ny = y + dy; // коорд соседних клеток
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT) { // если рассм клетка не препятствие и не была помечена раньше то идём внутрь
            auto cap = velocity.get(x, y, dx, dy);
            auto flow = velocity_flow.get(x, y, dx, dy);
            if (flow == cap) { // TODO?
                continue;
            }
            // assert(v >= velocity_flow.get(x, y, dx, dy));
            auto vp = min(lim, cap - flow); // TODO?
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

Fixed random01() {
    return Fixed::from_raw((rnd() & ((1 << 16) - 1)));
}

// прикол что помечает какие-то клетки рассмотренным, но какие TODO ? 
void propagate_stop(int x, int y, bool force = false) { // коорд текущей клетки и force TODO?
    if (!force) { // TODO ?
        bool stop = true;
        for (auto [dx, dy] : deltas) { // перебираем соседей
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) { // если есть хотя бы одна клетка из которой ещё не всё вылилось (пересчиталось) что должно TODO ?
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
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) { // если бордюрчик/уже помечена/ещё не всё вытекло ?TODO - не рассматриваем
            continue;
        }
        propagate_stop(nx, ny); // рекурсивно смотрим соседнюю
    }
}

Fixed move_prob(int x, int y) { // координаты клетки
    Fixed sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) { // если бордюр или уже рассмотренна - не трогаем
            continue;
        }
        auto v = velocity.get(x, y, dx, dy); // смотрим сколько чего-то а ней TODO?
        if (v < 0) {
            continue;
        }
        sum += v; // суммируем чего-то в клетке, если положительное число TODO ?
    }
    return sum;
}

struct ParticleParams { // класс чтобы инициализировать/менять клетки поля TODO?
    char type;
    Fixed cur_p;
    array<Fixed, deltas.size()> v;

    void swap_with(int x, int y) {
        swap(field[x][y], type);
        swap(p[x][y], cur_p);
        swap(velocity.v[x][y], v);
    }
};

bool propagate_move(int x, int y, bool is_first) { // клетки поля и is_first - вероятно корневая ли вершина TODO?
    last_use[x][y] = UT - is_first; // TODO?
    bool ret = false;
    int nx = -1, ny = -1;
    do {
        std::array<Fixed, deltas.size()> tres;
        Fixed sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) { // если бордюр или уже помеченная то TODO?
                tres[i] = sum;
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                tres[i] = sum;
                continue;
            }
            sum += v;
            tres[i] = sum;
        }

        if (sum == 0) {
            break;
        }

        Fixed p = random01() * sum;
        size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

        auto [dx, dy] = deltas[d];
        nx = x + dx;
        ny = y + dy;
        assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

        ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) { // перебираем соседей
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) { // если не бордюр, не была ещё рассмотренная и что-то в ней есть TODO?
            propagate_stop(nx, ny);
        }
    }
    if (ret) { // TODO?
        if (!is_first) { // TODO?
            ParticleParams pp{}; // TODO?
            pp.swap_with(x, y);
            pp.swap_with(nx, ny);
            pp.swap_with(x, y);
        }
    }
    return ret;
}

int dirs[N][M]{}; // TODO?

int main() {
    rho[' '] = 0.01; // задаём константы
    rho['.'] = 1000;
    Fixed g = 0.1; // типо g физическая

    for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
            if (field[x][y] == '#')
                continue;
            for (auto [dx, dy] : deltas) {
                dirs[x][y] += (field[x + dx][y + dy] != '#');
            }
        }
    }

    const auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < T; ++i) {
        
        Fixed total_delta_p = 0;
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
        memcpy(old_p, p, sizeof(p));
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                        auto delta_p = old_p[x][y] - old_p[nx][ny];
                        auto force = delta_p;
                        auto &contr = velocity.get(nx, ny, -dx, -dy);
                        if (contr * rho[(int) field[nx][ny]] >= force) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                        p[x][y] -= force / dirs[x][y];
                        total_delta_p -= force / dirs[x][y];
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
                        if (t > 0) {
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
                    if (old_v > 0) {
                        assert(new_v <= old_v);
                        velocity.get(x, y, dx, dy) = new_v;
                        auto force = (old_v - new_v) * rho[(int) field[x][y]];
                        if (field[x][y] == '.')
                            force *= 0.8;
                        if (field[x + dx][y + dy] == '#') {
                            p[x][y] += force / dirs[x][y];
                            total_delta_p += force / dirs[x][y];
                        } else {
                            p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                            total_delta_p += force / dirs[x + dx][y + dy];
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
            #ifndef NOT_TO_PRINT
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < N; ++x) {
                cout << field[x] << "\n";
            }
            #endif
        }
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff = end - start;
    std::cout << "Time spent: " << diff << '\n';
}