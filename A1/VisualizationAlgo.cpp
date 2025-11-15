//
// Created by arsen on 12.11.2025.
//
#include <bits/stdc++.h>
#include <cstdint>
using namespace std;

struct Circle {
    double x, y, r;
};
static bool point_in_all(const Circle &c0, const Circle &c1, const Circle &c2, double px, double py) {
    double dx = px - c0.x;
    double dy = py - c0.y;
    if (dx * dx + dy * dy > c0.r * c0.r)
        return false;
    dx = px - c1.x;
    dy = py - c1.y;
    if (dx * dx + dy * dy > c1.r * c1.r)
        return false;
    dx = px - c2.x;
    dy = py - c2.y;
    if (dx * dx + dy * dy > c2.r * c2.r)
        return false;
    return true;
}
static void compute_wide_box(const array<Circle, 3> &C, double &xmin, double &xmax, double &ymin, double &ymax) {
    xmin = C[0].x - C[0].r;
    xmax = C[0].x + C[0].r;
    ymin = C[0].y - C[0].r;
    ymax = C[0].y + C[0].r;
    for (int i = 1; i < 3; ++i) {
        xmin = min(xmin, C[i].x - C[i].r);
        xmax = max(xmax, C[i].x + C[i].r);
        ymin = min(ymin, C[i].y - C[i].r);
        ymax = max(ymax, C[i].y + C[i].r);
    }
}
static void compute_tight_box(const array<Circle, 3> &C, double wide_xmin, double wide_xmax, double wide_ymin,
                              double wide_ymax, double &tight_xmin, double &tight_xmax, double &tight_ymin,
                              double &tight_ymax, int gridN = 801) {
    bool any = false;
    tight_xmin = wide_xmin;
    tight_xmax = wide_xmax;
    tight_ymin = wide_ymin;
    tight_ymax = wide_ymax;
    for (int iy = 0; iy < gridN; ++iy) {
        double y = wide_ymin + (wide_ymax - wide_ymin) * (double(iy) / (gridN - 1));
        for (int ix = 0; ix < gridN; ++ix) {
            double x = wide_xmin + (wide_xmax - wide_xmin) * (double(ix) / (gridN - 1));
            if (point_in_all(C[0], C[1], C[2], x, y)) {
                if (!any) {
                    tight_xmin = tight_xmax = x;
                    tight_ymin = tight_ymax = y;
                    any = true;
                } else {
                    tight_xmin = min(tight_xmin, x);
                    tight_xmax = max(tight_xmax, x);
                    tight_ymin = min(tight_ymin, y);
                    tight_ymax = max(tight_ymax, y);
                }
            }
        }
    }
    if (!any) {
        tight_xmin = wide_xmin;
        tight_xmax = wide_xmax;
        tight_ymin = wide_ymin;
        tight_ymax = wide_ymax;
        return;
    }
    double dx = (wide_xmax - wide_xmin) / double(gridN - 1);
    double dy = (wide_ymax - wide_ymin) / double(gridN - 1);
    tight_xmin = max(wide_xmin, tight_xmin - 1.5 * dx);
    tight_xmax = min(wide_xmax, tight_xmax + 1.5 * dx);
    tight_ymin = max(wide_ymin, tight_ymin - 1.5 * dy);
    tight_ymax = min(wide_ymax, tight_ymax + 1.5 * dy);
}
int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    bool has_exact = false;
    double exactS = 0.0;
    for (int i = 1; i + 1 < argc; ++i) {
        string a = argv[i];
        if (a == "--exact") {
            exactS = stod(argv[i + 1]);
            has_exact = true;
            ++i;
        }
    }
    array<Circle, 3> C;
    for (int i = 0; i < 3; ++i) {
        if (!(cin >> C[i].x >> C[i].y >> C[i].r)) {
            cerr << "Ошибка: ожидаются 3 строки x y r\n";
            return 1;
        }
    }
    double wide_xmin, wide_xmax, wide_ymin, wide_ymax;
    compute_wide_box(C, wide_xmin, wide_xmax, wide_ymin, wide_ymax);
    if (!(wide_xmax > wide_xmin && wide_ymax > wide_ymin)) {
        cerr << "Плохая геометрия: ширина/высота рамки <= 0\n";
        cout << "N,S_wide,relerr_wide,S_tight,relerr_tight\n";
        return 0;
    }
    double tight_xmin, tight_xmax, tight_ymin, tight_ymax;
    compute_tight_box(C, wide_xmin, wide_xmax, wide_ymin, wide_ymax, tight_xmin, tight_xmax, tight_ymin, tight_ymax,
                      801);
    vector<uint64_t> Ns;
    for (uint64_t n = 100; n <= 100000; n += 500)
        Ns.push_back(n);
    if (Ns.back() != 100000)
        Ns.push_back(100000);
    uint64_t Nmax = Ns.back();
    auto run_once = [&](double xmin, double xmax, double ymin, double ymax, uint64_t seed) -> vector<uint32_t> {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> dx(xmin, xmax);
        std::uniform_real_distribution<double> dy(ymin, ymax);
        vector<uint32_t> cumsum;
        cumsum.reserve(Nmax);
        uint32_t count = 0;
        for (uint64_t i = 0; i < Nmax; ++i) {
            double x = dx(rng);
            double y = dy(rng);
            if (point_in_all(C[0], C[1], C[2], x, y))
                ++count;
            cumsum.push_back(count);
        }
        return cumsum;
    };
    const uint64_t SEED_WIDE = 42;
    const uint64_t SEED_TIGHT = 4242;
    auto csum_wide = run_once(wide_xmin, wide_xmax, wide_ymin, wide_ymax, SEED_WIDE);
    auto csum_tight = run_once(tight_xmin, tight_xmax, tight_ymin, tight_ymax, SEED_TIGHT);
    double area_wide = (wide_xmax - wide_xmin) * (wide_ymax - wide_ymin);
    double area_tight = (tight_xmax - tight_xmin) * (tight_ymax - tight_ymin);
    cout.setf(std::ios::fixed);
    cout << setprecision(12);
    cout << "N,S_широкая,отн_ошибка_широкая,S_узкая,отн_ошибка_узкая\n";
    for (size_t idx = 0; idx < Ns.size(); ++idx) {
        uint64_t N = Ns[idx];
        uint32_t hits_wide = csum_wide[N - 1];
        uint32_t hits_tight = csum_tight[N - 1];
        double est_wide  = (double)hits_wide  / (double)N * area_wide;
        double est_tight = (double)hits_tight / (double)N * area_tight;
        if (has_exact) {
            double rel_wide = (exactS != 0.0)
                              ? fabs(est_wide - exactS) / exactS
                              : fabs(est_wide - exactS);

            double rel_tight = (exactS != 0.0)
                               ? fabs(est_tight - exactS) / exactS
                               : fabs(est_tight - exactS);

            cout << N << "," << est_wide << "," << rel_wide
                 << "," << est_tight << "," << rel_tight << "\n";

        } else {
            cout << N << "," << est_wide << ",,"
                 << est_tight << ",\n";
        }
    }
    const int visNx = 500, visNy = 500;
    ofstream fout("visualization.csv");
    if (fout) {
        fout << "x,y,inside\n";
        for (int iy = 0; iy < visNy; ++iy) {
            double y = tight_ymin + (tight_ymax - tight_ymin) * ((iy + 0.5) / double(visNy));
            for (int ix = 0; ix < visNx; ++ix) {
                double x = tight_xmin + (tight_xmax - tight_xmin) * ((ix + 0.5) / double(visNx));
                int inside = point_in_all(C[0], C[1], C[2], x, y) ? 1 : 0;
                fout << x << "," << y << "," << inside << "\n";
            }
        }
        fout.close();
    } else {
        cerr << "Не удалось создать visualization.csv\n";
    }
    cerr << "SEEDs: wide=" << SEED_WIDE << " tight=" << SEED_TIGHT << "\n";
    cerr << "wide box: x in [" << wide_xmin << ", " << wide_xmax << "], y in [" << wide_ymin << ", " << wide_ymax
         << "], area=" << area_wide << "\n";
    cerr << "tight box: x in [" << tight_xmin << ", " << tight_xmax << "], y in [" << tight_ymin << ", " << tight_ymax
         << "], area=" << area_tight << "\n";
    cerr << "Nmax used = " << Nmax << "\n";

    return 0;
}
