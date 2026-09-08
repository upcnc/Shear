#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
namespace fs = std::filesystem;

/*
 * Discrete-unit IBM detachment in a straight microchannel,
 * following Wang, Fu, Wu, Li, Wang,
 * Int. Biodeterior. Biodegrad. 207 (2026) 106245.
 *
 * Each biofilm unit is a circular ring of Lagrangian IB points.
 * Units connect when centre distance < Tc and break when
 * |L - L0| > Tf = 0.6 r (their experimental calibration).
 *
 * Cases (paper Fig. 6):
 *   1  semicircle, layered units
 *   2  semicircle, random units
 *   3  rectangle, layered / staggered units
 *
 *   g++ -O3 -std=c++17 -fopenmp wang_ibm_detach.cpp -o wang_detach
 *   ./wang_detach out_case1 1 1.0e6
 */

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

constexpr int Q = 9;
const int ex[Q] = {0,1,0,-1,0,1,-1,-1,1};
const int ey[Q] = {0,0,1,0,-1,1,1,-1,-1};
const double w[Q] = {4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};
const int opp[Q] = {0,3,4,1,2,7,8,5,6};

constexpr double DX = 2.0e-6;
constexpr int NX = 220;
constexpr int NY = 60;
constexpr double NU_PHYS = 1.0e-6;
constexpr double RHO_PHYS = 1000.0;
constexpr double MU_PHYS = NU_PHYS * RHO_PHYS;
constexpr double NU_L = 0.10;
constexpr double DT = NU_L * DX * DX / NU_PHYS;
constexpr double TAU_F = 3.0 * NU_L + 0.5;

constexpr double UMAX_PHYS = 5.0e-4;
constexpr double UMAX_LU = UMAX_PHYS * DT / DX;

constexpr int N_RING = 8;
constexpr double R_UNIT = 4.0e-6;
constexpr double TC = 2.35 * R_UNIT;
constexpr double TF = 0.60 * R_UNIT;
constexpr double KS_RING = 2.0e-3;
constexpr double KBEND = 4.0e-4;
constexpr double A0 = 4.0e-12;
constexpr double KUW = 8.0e-3;
constexpr double KT_WALL = 0.80;

constexpr int LBM_STEPS = 25000;
constexpr int OUTPUT_EVERY = 500;
constexpr unsigned SEED = 20260908u;

struct Marker {
    double x = 0, y = 0;
    double fx = 0, fy = 0;
    int unit = -1;
    bool wall = false;
    bool tether = false;
    double x0 = 0, y0 = 0;
};

struct Bond {
    int a = 0, b = 0;
    double L0 = 0;
    bool inter = false;
    bool wallb = false;
    bool alive = true;
};

struct Unit {
    vector<int> pts;
    double cx = 0, cy = 0;
    bool anchored = false;
};

vector<double> fin, fout;
vector<double> rho, ux, uy;
vector<char> solid;
vector<Marker> mk;
vector<Bond> bonds;
vector<Unit> units;
mt19937 rng(SEED);

inline int id(int i, int j) { return i * NY + j; }
inline bool inside(int i, int j) { return i >= 0 && i < NX && j >= 0 && j < NY; }

double feq(int k, double r, double u, double v) {
    const double eu = ex[k] * u + ey[k] * v;
    const double uu = u * u + v * v;
    return w[k] * r * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uu);
}

void alloc_fluid() {
    const int n = NX * NY;
    fin.assign(n * Q, 0.0);
    fout.assign(n * Q, 0.0);
    rho.assign(n, 1.0);
    ux.assign(n, 0.0);
    uy.assign(n, 0.0);
    solid.assign(n, 0);
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            solid[id(i, j)] = (j == 0 || j == NY - 1);
}

void init_poiseuille() {
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j) {
            const int p = id(i, j);
            rho[p] = 1.0;
            if (solid[p]) { ux[p] = uy[p] = 0.0; }
            else {
                const double eta = (j + 0.5) / NY;
                ux[p] = 4.0 * UMAX_LU * eta * (1.0 - eta);
                uy[p] = 0.0;
            }
            for (int k = 0; k < Q; ++k) fin[p * Q + k] = feq(k, rho[p], ux[p], uy[p]);
        }
}

int add_marker(double x, double y, int uid, bool wall, bool tether) {
    Marker m;
    m.x = x; m.y = y; m.x0 = x; m.y0 = y;
    m.unit = uid; m.wall = wall; m.tether = tether;
    mk.push_back(m);
    return (int)mk.size() - 1;
}

void add_bond(int a, int b, bool inter, bool wallb) {
    Bond bd;
    bd.a = a; bd.b = b; bd.inter = inter; bd.wallb = wallb; bd.alive = true;
    bd.L0 = hypot(mk[a].x - mk[b].x, mk[a].y - mk[b].y);
    bonds.push_back(bd);
}

int make_unit(double cx, double cy, bool anchor) {
    Unit u;
    u.cx = cx; u.cy = cy; u.anchored = anchor;
    const int uid = (int)units.size();
    for (int k = 0; k < N_RING; ++k) {
        const double th = 2.0 * M_PI * k / N_RING;
        const int p = add_marker(cx + R_UNIT * cos(th), cy + R_UNIT * sin(th), uid, false, false);
        u.pts.push_back(p);
    }
    units.push_back(u);
    for (int k = 0; k < N_RING; ++k) {
        add_bond(u.pts[k], u.pts[(k + 1) % N_RING], false, false);
        add_bond(u.pts[k], u.pts[(k + N_RING / 2) % N_RING], false, false);
    }
    if (anchor) {
        for (int k = 0; k < 3; ++k) {
            const int p = u.pts[(N_RING / 2 + k - 1 + N_RING) % N_RING];
            const int w = add_marker(mk[p].x, 1.02 * DX, -1, true, true);
            add_bond(p, w, false, true);
        }
    }
    return uid;
}

void build_case1() {
    const double cx = 0.38 * NX * DX;
    const double R = 28.0e-6;
    const double dy = 2.15 * R_UNIT;
    int row = 0;
    for (double y = 1.6 * R_UNIT; y <= R; y += dy, ++row) {
        const double half = sqrt(max(0.0, R * R - y * y));
        const double dxu = 2.15 * R_UNIT;
        for (double x = cx - half + R_UNIT; x <= cx + half - 0.5 * R_UNIT; x += dxu)
            make_unit(x, y, row == 0);
    }
}

void build_case2() {
    const double cx = 0.38 * NX * DX;
    const double R = 28.0e-6;
    uniform_real_distribution<double> ur(-1.0, 1.0);
    vector<pair<double, double>> pos;
    for (int t = 0; t < 400 && (int)pos.size() < 28; ++t) {
        const double rr = 0.15 + 0.80 * fabs(ur(rng));
        const double th = M_PI * 0.5 * (1.0 + ur(rng));
        const double x = cx + rr * R * cos(th);
        const double y = 1.4 * R_UNIT + rr * R * sin(th);
        if (y < 1.3 * R_UNIT || y > R) continue;
        bool ok = true;
        for (auto &q : pos)
            if (hypot(q.first - x, q.second - y) < 2.05 * R_UNIT) ok = false;
        if (ok) pos.push_back({x, y});
    }
    sort(pos.begin(), pos.end(), [](auto &a, auto &b) { return a.second < b.second; });
    for (size_t i = 0; i < pos.size(); ++i)
        make_unit(pos[i].first, pos[i].second, pos[i].second < 3.2 * R_UNIT);
}

void build_case3() {
    const double x0 = 0.28 * NX * DX;
    const double L = 72.0e-6, H = 28.0e-6;
    const double dxu = 2.15 * R_UNIT, dy = 1.90 * R_UNIT;
    int row = 0;
    for (double y = 1.6 * R_UNIT; y <= H; y += dy, ++row) {
        const double shift = (row % 2) ? 0.5 * dxu : 0.0;
        for (double x = x0 + shift; x <= x0 + L; x += dxu)
            make_unit(x, y, row == 0);
    }
}

void connect_nearby() {
    for (size_t i = 0; i < units.size(); ++i)
        for (size_t j = i + 1; j < units.size(); ++j) {
            double d = hypot(units[i].cx - units[j].cx, units[i].cy - units[j].cy);
            if (d > TC) continue;
            int ia = units[i].pts[0], ib = units[j].pts[0];
            double best = 1e300;
            for (int a : units[i].pts)
                for (int b : units[j].pts) {
                    double dd = hypot(mk[a].x - mk[b].x, mk[a].y - mk[b].y);
                    if (dd < best) { best = dd; ia = a; ib = b; }
                }
            bool exists = false;
            for (auto &bd : bonds)
                if (bd.inter && bd.alive &&
                    ((bd.a == ia && bd.b == ib) || (bd.a == ib && bd.b == ia)))
                    exists = true;
            if (!exists) add_bond(ia, ib, true, false);
        }
}

void update_unit_centres() {
    for (auto &u : units) {
        u.cx = u.cy = 0.0;
        for (int p : u.pts) { u.cx += mk[p].x; u.cy += mk[p].y; }
        u.cx /= max(1, (int)u.pts.size());
        u.cy /= max(1, (int)u.pts.size());
    }
}

double spring_ks_inter(double E_pa) {
    const double scale = DT * DT / (RHO_PHYS * DX * DX * DX);
    return (E_pa * A0) * scale;
}

void accumulate_forces(double E_pa) {
    for (auto &m : mk) m.fx = m.fy = 0.0;
    const double ksi = spring_ks_inter(E_pa);
    for (auto &bd : bonds) {
        if (!bd.alive) continue;
        Marker &A = mk[bd.a], &B = mk[bd.b];
        const double dx = B.x - A.x, dy = B.y - A.y;
        const double L = hypot(dx, dy);
        if (L < 1e-16) continue;
        if (bd.inter && fabs(L - bd.L0) > TF) { bd.alive = false; continue; }
        const double k = bd.wallb ? KUW : (bd.inter ? ksi : KS_RING);
        const double f = k * (L - bd.L0);
        const double fx = f * dx / L, fy = f * dy / L;
        A.fx += fx; A.fy += fy;
        B.fx -= fx; B.fy -= fy;
    }
    for (auto &u : units) {
        const int n = (int)u.pts.size();
        for (int k = 0; k < n; ++k) {
            const int im = u.pts[(k - 1 + n) % n];
            const int i = u.pts[k];
            const int ip = u.pts[(k + 1) % n];
            const double ax = mk[i].x - mk[im].x, ay = mk[i].y - mk[im].y;
            const double bx = mk[ip].x - mk[i].x, by = mk[ip].y - mk[i].y;
            const double cross = ax * by - ay * bx;
            const double target = R_UNIT * R_UNIT * sin(2.0 * M_PI / n);
            const double c = KBEND * (cross - target);
            mk[i].fx += c * (-ay + by);
            mk[i].fy += c * (ax - bx);
        }
    }
    for (auto &m : mk) if (m.tether) {
        m.fx += KT_WALL * (m.x0 - m.x);
        m.fy += KT_WALL * (m.y0 - m.y);
    }
}

vector<double> Fx, Fy;

void spread_to_euler() {
    const int n = NX * NY;
    Fx.assign(n, 0.0); Fy.assign(n, 0.0);
    auto phi = [](double r) {
        r = fabs(r);
        if (r >= 2.0) return 0.0;
        if (r <= 1.0) return 0.125 * (3.0 - 2.0 * r + sqrt(max(0.0, 1.0 + 4.0 * r - 4.0 * r * r)));
        return 0.125 * (5.0 - 2.0 * r - sqrt(max(0.0, -7.0 + 12.0 * r - 4.0 * r * r)));
    };
    const double inv = 1.0 / (DX * DX);
    for (const auto &m : mk) {
        const double gx = m.x / DX - 0.5, gy = m.y / DX - 0.5;
        const int i0 = (int)floor(gx) - 1, j0 = (int)floor(gy) - 1;
        for (int di = 0; di <= 3; ++di)
            for (int dj = 0; dj <= 3; ++dj) {
                int i = i0 + di, j = j0 + dj;
                if (!inside(i, j) || solid[id(i, j)]) continue;
                const double s = phi(gx - i) * phi(gy - j) * inv;
                Fx[id(i, j)] += s * m.fx;
                Fy[id(i, j)] += s * m.fy;
            }
    }
}

void interpolate_move() {
    auto phi = [](double r) {
        r = fabs(r);
        if (r >= 2.0) return 0.0;
        if (r <= 1.0) return 0.125 * (3.0 - 2.0 * r + sqrt(max(0.0, 1.0 + 4.0 * r - 4.0 * r * r)));
        return 0.125 * (5.0 - 2.0 * r - sqrt(max(0.0, -7.0 + 12.0 * r - 4.0 * r * r)));
    };
    for (auto &m : mk) {
        if (m.tether) { m.x = m.x0; m.y = m.y0; continue; }
        const double gx = m.x / DX - 0.5, gy = m.y / DX - 0.5;
        const int i0 = (int)floor(gx) - 1, j0 = (int)floor(gy) - 1;
        double u = 0, v = 0, ssum = 0;
        for (int di = 0; di <= 3; ++di)
            for (int dj = 0; dj <= 3; ++dj) {
                int i = i0 + di, j = j0 + dj;
                if (!inside(i, j) || solid[id(i, j)]) continue;
                const double s = phi(gx - i) * phi(gy - j);
                u += s * ux[id(i, j)];
                v += s * uy[id(i, j)];
                ssum += s;
            }
        if (ssum > 1e-12) { u /= ssum; v /= ssum; }
        m.x += u * DX;
        m.y += v * DX;
        m.x = min(max(m.x, DX), (NX - 1.01) * DX);
        m.y = min(max(m.y, 1.05 * DX), (NY - 1.05) * DX);
    }
}

void collide_stream() {
    const double g = 8.0 * NU_L * UMAX_LU / ((NY - 2) * (NY - 2));
    vector<double> nxt(NX * NY * Q, 0.0);
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j) {
            const int p = id(i, j);
            if (solid[p]) {
                for (int k = 0; k < Q; ++k) fout[p * Q + k] = fin[p * Q + opp[k]];
                continue;
            }
            double r = 0, mx = 0, my = 0;
            for (int k = 0; k < Q; ++k) {
                r += fin[p * Q + k];
                mx += ex[k] * fin[p * Q + k];
                my += ey[k] * fin[p * Q + k];
            }
            r = max(r, 1e-12);
            const double fx = Fx[p] * (DT * DT / (RHO_PHYS * DX)) + g;
            const double fy = Fy[p] * (DT * DT / (RHO_PHYS * DX));
            ux[p] = (mx + 0.5 * fx) / r;
            uy[p] = (my + 0.5 * fy) / r;
            rho[p] = r;
            for (int k = 0; k < Q; ++k) {
                const double eu = ex[k] * ux[p] + ey[k] * uy[p];
                const double fsrc = (1.0 - 0.5 / TAU_F) * w[k] *
                    (3.0 * ((ex[k] - ux[p]) * fx + (ey[k] - uy[p]) * fy) + 9.0 * eu * (ex[k] * fx + ey[k] * fy));
                fout[p * Q + k] = fin[p * Q + k] - (fin[p * Q + k] - feq(k, r, ux[p], uy[p])) / TAU_F + fsrc;
            }
        }
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < Q; ++k) {
                int ni = i + ex[k], nj = j + ey[k];
                if (ni < 0) ni += NX;
                if (ni >= NX) ni -= NX;
                if (nj < 0 || nj >= NY) continue;
                nxt[id(ni, nj) * Q + k] += fout[id(i, j) * Q + k];
            }
    fin.swap(nxt);
}

int count_alive_inter() {
    int n = 0;
    for (auto &b : bonds) if (b.inter && b.alive) ++n;
    return n;
}

int count_anchored_units() {
    int n = 0;
    for (auto &u : units) {
        bool hold = false;
        for (auto &b : bonds)
            if (b.wallb && b.alive)
                for (int p : u.pts) if (p == b.a || p == b.b) hold = true;
        if (hold) ++n;
    }
    return n;
}

void write_frame(const fs::path &dir, int step, double t, double E) {
    ostringstream name;
    name << "state_" << setw(6) << setfill('0') << step << ".dat";
    ofstream out(dir / name.str());
    out << "TITLE = \"WangIBM\"\nVARIABLES = \"X_um\",\"Y_um\",\"U\",\"V\",\"Phase\"\n";
    out << "ZONE T=\"Field\", I=" << NX << ", J=" << NY << ", F=POINT\n";
    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            const int p = id(i, j);
            int ph = solid[p] ? 1 : 0;
            out << i * DX * 1e6 << ' ' << j * DX * 1e6 << ' '
                << ux[p] * DX / DT << ' ' << uy[p] * DX / DT << ' ' << ph << '\n';
        }
    out << "ZONE T=\"Units\", I=" << max(1, (int)units.size()) << ", F=POINT\n";
    if (units.empty()) out << "0 0 0 0 2\n";
    else for (auto &u : units)
        out << u.cx * 1e6 << ' ' << u.cy * 1e6 << " 0 0 2\n";
    int nb = 0;
    for (auto &b : bonds) if (b.alive) ++nb;
    out << "ZONE T=\"Bonds\", I=" << max(1, 2 * nb) << ", F=POINT\n";
    if (!nb) out << "0 0 0 0 3\n";
    else for (auto &b : bonds) if (b.alive) {
        out << mk[b.a].x * 1e6 << ' ' << mk[b.a].y * 1e6 << " 0 0 3\n";
        out << mk[b.b].x * 1e6 << ' ' << mk[b.b].y * 1e6 << " 0 0 3\n";
    }
    (void)E; (void)t;
}

int main(int argc, char **argv) {
    const fs::path outdir = argc > 1 ? fs::path(argv[1]) : fs::path("WangDetach");
    const int cse = argc > 2 ? atoi(argv[2]) : 1;
    const double E_pa = argc > 3 ? atof(argv[3]) : 1.0e6;
    fs::create_directories(outdir);
    alloc_fluid();
    init_poiseuille();
    if (cse == 2) build_case2();
    else if (cse == 3) build_case3();
    else build_case1();
    update_unit_centres();
    connect_nearby();

    ofstream hist(outdir / "history.csv");
    hist << "step,time_s,n_units,n_inter_bonds,n_anchored,mean_cx_um,mean_cy_um,frac_detached\n";

    cout << "Wang-style IBM detach  case=" << cse << "  E=" << E_pa
         << " Pa  units=" << units.size() << "  bonds=" << bonds.size() << "\n";

    for (int n = 0; n <= LBM_STEPS; ++n) {
        update_unit_centres();
        if (n % 40 == 0) connect_nearby();
        accumulate_forces(E_pa);
        spread_to_euler();
        collide_stream();
        interpolate_move();

        if (n % OUTPUT_EVERY == 0 || n == LBM_STEPS) {
            double mx = 0, my = 0;
            for (auto &u : units) { mx += u.cx; my += u.cy; }
            mx /= max(1, (int)units.size());
            my /= max(1, (int)units.size());
            const int anc = count_anchored_units();
            const double frac = 1.0 - double(anc) / max(1, (int)units.size());
            hist << n << ',' << n * DT << ',' << units.size() << ','
                 << count_alive_inter() << ',' << anc << ','
                 << mx * 1e6 << ',' << my * 1e6 << ',' << frac << '\n';
            write_frame(outdir, n, n * DT, E_pa);
            cout << "n=" << n << " t=" << scientific << n * DT
                 << " bonds=" << count_alive_inter()
                 << " anchored=" << anc
                 << " detached=" << frac << "\n";
        }
    }
    ofstream sum(outdir / "summary.txt");
    sum << "case = " << cse << "\nE_Pa = " << E_pa
        << "\nTf_over_r = 0.6\nunits = " << units.size()
        << "\nfinal_anchored = " << count_anchored_units()
        << "\nfinal_inter_bonds = " << count_alive_inter() << "\n";
    cout << "done -> " << outdir << "\n";
    return 0;
}
