#ifndef COMMON_HPP
#define COMMON_HPP

#include <array>
#include <vector>
#include <optional>
#include <cmath>

typedef std::array<std::array<double, 2>, 3> VolRegion;
typedef std::array<double, 3> Coord;
typedef std::vector<Coord> Coords;
typedef std::array<int, 3> CoordIndex;
typedef std::tuple<std::optional<int>, double, std::optional<int> > IQHb;

struct GridValue {
    double vdw_or_clash;
    double elec;
    std::optional<double> hd_donor;
    std::optional<double> hd_acceptor;

    GridValue(double vdw, double e, std::optional<double> donor, std::optional<double> acceptor)
        : vdw_or_clash(vdw), elec(e), hd_donor(donor), hd_acceptor(acceptor) {}
};

struct TerAtom {
    Coord coord;
    double charge;
};

// 添加一个角度规范化函数
/*
inline double normalize_angle(double angle) {
    while (angle > 180.0) angle -= 360.0;
    while (angle < -180.0) angle += 360.0;
    return angle;
}*/
inline double normalize_angle(double angle) {
    angle = std::fmod(angle + 180.0, 360.0);
    return angle >= 0 ? angle - 180.0 : angle + 180.0;
}
#endif // COMMON_HPP
