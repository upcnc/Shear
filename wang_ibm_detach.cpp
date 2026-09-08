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

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

constexpr int Q = 9;
const int ex[Q] = {0,1,0,-1,0,1,-1,-1,1};
const int ey[Q] = {0,0,1,0,-1,1,1,-1,-1};
const double w[Q] = {4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};
const int opp[Q] = {0,3,4,1,2,7,8,5,6};

constexpr double DX = 2.0e-6;
constexpr int NX = 240;
constexpr int NY = 56;
constexpr double NU_PHYS = 1.0e-6;
constexpr double RHO_PHYS = 1000.0;
constexpr double NU_L = 0.08;
constexpr double DT = NU_L * DX * DX / NU_PHYS;
constexpr double TAU_F = 3.0 * NU_L + 0.5;
constexpr double UMAX_PHYS = 2.0e-3;
constexpr double UMAX_LU = UMAX_PHYS * DT / DX;

constexpr int N_RING = 8;
constexpr double R_UNIT = 2.6e-6;
constexpr double SPACING = 5.4e-6;
constexpr double TC = 2.20 * R_UNIT;
constexpr double TF = 0.60 * R_UNIT;
constexpr double TF_WALL = 0.45 * R_UNIT;

constexpr int LBM_STEPS = 40000;
constexpr int OUTPUT_EVERY = 1000;
constexpr unsigned SEED = 20260908u;

struct Marker {
    double x=0,y=0,fx=0,fy=0,x0=0,y0=0;
    int unit=-1;
    bool tether=false;
};
struct Bond {
    int a=0,b=0,ua=-1,ub=-1;
    double L0=0;
    bool inter=false,wallb=false,alive=true;
};
struct Unit {
    vector<int> pts;
    double cx=0,cy=0,cx0=0,cy0=0;
    bool started_anchored=false;
};

vector<double> fin,fout,rho,ux,uy,Fx,Fy,alpha;
vector<char> solid;
vector<Marker> mk;
vector<Bond> bonds;
vector<Unit> units;
mt19937 rng(SEED);

inline int id(int i,int j){return i*NY+j;}
inline bool inside(int i,int j){return i>=0&&i<NX&&j>=0&&j<NY;}

double feq(int k,double r,double u,double v){
    const double eu=ex[k]*u+ey[k]*v, uu=u*u+v*v;
    return w[k]*r*(1.0+3.0*eu+4.5*eu*eu-1.5*uu);
}
double peskin(double r){
    r=fabs(r);
    if(r>=2.0) return 0.0;
    if(r<=1.0) return 0.125*(3.0-2.0*r+sqrt(max(0.0,1.0+4.0*r-4.0*r*r)));
    return 0.125*(5.0-2.0*r-sqrt(max(0.0,-7.0+12.0*r-4.0*r*r)));
}

void alloc_fluid(){
    const int n=NX*NY;
    fin.assign(n*Q,0); fout.assign(n*Q,0);
    rho.assign(n,1); ux.assign(n,0); uy.assign(n,0);
    Fx.assign(n,0); Fy.assign(n,0); alpha.assign(n,0); solid.assign(n,0);
    for(int i=0;i<NX;++i) for(int j=0;j<NY;++j) solid[id(i,j)]=(j==0||j==NY-1);
}
void init_poiseuille(){
    for(int i=0;i<NX;++i) for(int j=0;j<NY;++j){
        const int p=id(i,j); rho[p]=1.0;
        if(solid[p]) ux[p]=uy[p]=0;
        else { const double eta=(j+0.5)/NY; ux[p]=4.0*UMAX_LU*eta*(1.0-eta); uy[p]=0; }
        for(int k=0;k<Q;++k) fin[p*Q+k]=feq(k,rho[p],ux[p],uy[p]);
    }
}
int add_marker(double x,double y,int uid,bool tether){
    Marker m; m.x=x; m.y=y; m.x0=x; m.y0=y; m.unit=uid; m.tether=tether;
    mk.push_back(m); return (int)mk.size()-1;
}
void add_bond(int a,int b,bool inter,bool wallb,int ua=-1,int ub=-1){
    Bond bd; bd.a=a; bd.b=b; bd.inter=inter; bd.wallb=wallb; bd.ua=ua; bd.ub=ub; bd.alive=true;
    bd.L0=hypot(mk[a].x-mk[b].x,mk[a].y-mk[b].y);
    if(bd.L0<1e-12) bd.L0=R_UNIT;
    bonds.push_back(bd);
}
int make_unit(double cx,double cy,bool anchor){
    Unit u; u.cx=u.cx0=cx; u.cy=u.cy0=cy; u.started_anchored=anchor;
    const int uid=(int)units.size();
    const int pc=add_marker(cx,cy,uid,false);
    u.pts.push_back(pc);
    for(int k=0;k<N_RING;++k){
        const double th=2.0*M_PI*k/N_RING;
        const int p=add_marker(cx+R_UNIT*cos(th),cy+R_UNIT*sin(th),uid,false);
        u.pts.push_back(p); add_bond(pc,p,false,false,uid,uid);
    }
    for(int k=0;k<N_RING;++k)
        add_bond(u.pts[1+k],u.pts[1+(k+1)%N_RING],false,false,uid,uid);
    units.push_back(u);
    if(anchor){
        const int w=add_marker(cx,1.05*DX,-1,true);
        add_bond(pc,w,false,true,uid,-1);
    }
    return uid;
}
void build_case1(){
    const double cx=0.40*NX*DX, R=50.0e-6; int row=0;
    for(double y=1.35*R_UNIT; y<=R-0.4*R_UNIT; y+=SPACING*0.87, ++row){
        const double half=sqrt(max(0.0,R*R-y*y));
        const double shift=(row%2)?0.5*SPACING:0.0;
        for(double x=cx-half+R_UNIT+shift; x<=cx+half-R_UNIT; x+=SPACING)
            make_unit(x,y,row==0);
    }
}
void build_case2(){
    const double cx=0.40*NX*DX, R=50.0e-6;
    uniform_real_distribution<double> ur(0.0,1.0);
    vector<pair<double,double>> pos;
    for(int t=0;t<8000 && (int)pos.size()<90;++t){
        const double th=M_PI*ur(rng), rr=R*sqrt(ur(rng));
        const double x=cx+rr*cos(th), y=1.2*R_UNIT+rr*sin(th);
        if(y<1.15*R_UNIT||y>R) continue;
        bool ok=true;
        for(auto &q:pos) if(hypot(q.first-x,q.second-y)<0.96*SPACING){ok=false;break;}
        if(ok) pos.push_back({x,y});
    }
    sort(pos.begin(),pos.end(),[](auto&a,auto&b){return a.second<b.second;});
    for(auto &q:pos) make_unit(q.first,q.second,q.second<3.0*R_UNIT);
}
void build_case3(){
    const double x0=0.28*NX*DX, L=100.0e-6, H=50.0e-6; int row=0;
    for(double y=1.35*R_UNIT; y<=H; y+=SPACING*0.87, ++row){
        const double shift=(row%2)?0.5*SPACING:0.0;
        for(double x=x0+shift; x<=x0+L; x+=SPACING) make_unit(x,y,row==0);
    }
}
void update_unit_centres(){
    for(auto &u:units){
        u.cx=u.cy=0;
        for(int p:u.pts){u.cx+=mk[p].x; u.cy+=mk[p].y;}
        const double inv=1.0/max(1,(int)u.pts.size());
        u.cx*=inv; u.cy*=inv;
    }
}
void connect_nearby(){
    for(size_t i=0;i<units.size();++i) for(size_t j=i+1;j<units.size();++j){
        if(hypot(units[i].cx-units[j].cx,units[i].cy-units[j].cy)>TC+R_UNIT) continue;
        bool exists=false;
        for(auto &bd:bonds)
            if(bd.inter&&bd.alive&&((bd.ua==(int)i&&bd.ub==(int)j)||(bd.ua==(int)j&&bd.ub==(int)i)))
                exists=true;
        if(exists) continue;
        int ia=units[i].pts[0], ib=units[j].pts[0]; double best=1e300;
        for(int a:units[i].pts) for(int b:units[j].pts){
            const double dd=hypot(mk[a].x-mk[b].x,mk[a].y-mk[b].y);
            if(dd<best){best=dd; ia=a; ib=b;}
        }
        if(best<TC) add_bond(ia,ib,true,false,(int)i,(int)j);
    }
}
double k_inter_lu(double E_pa){
    const double A0=DX*DX;
    const double k_phys=E_pa*A0/max(SPACING,DX);
    return k_phys*DT*DT/(RHO_PHYS*DX*DX);
}
double k_ring_lu(){ return 8.0*k_inter_lu(4.0e6); }
double k_wall_lu(){ return 6.0*k_inter_lu(4.0e6); }

void accumulate_forces(double E_pa){
    for(auto &m:mk) m.fx=m.fy=0;
    const double ki=k_inter_lu(E_pa), kr=k_ring_lu(), kw=k_wall_lu();
    for(auto &bd:bonds){
        if(!bd.alive) continue;
        Marker &A=mk[bd.a], &B=mk[bd.b];
        const double dx=B.x-A.x, dy=B.y-A.y, L=hypot(dx,dy);
        if(L<1e-16) continue;
        const double lim=bd.wallb?TF_WALL:TF;
        if((bd.inter||bd.wallb) && fabs(L-bd.L0)>lim){ bd.alive=false; continue; }
        const double k=bd.wallb?kw:(bd.inter?ki:kr);
        const double f=k*(L-bd.L0);
        A.fx+=f*dx/L; A.fy+=f*dy/L; B.fx-=f*dx/L; B.fy-=f*dy/L;
    }
    for(auto &m:mk) if(m.tether){ m.fx+=kw*(m.x0-m.x); m.fy+=kw*(m.y0-m.y); }
}
void spread_to_euler(){
    fill(Fx.begin(),Fx.end(),0); fill(Fy.begin(),Fy.end(),0); fill(alpha.begin(),alpha.end(),0);
    for(const auto &m:mk){
        const double gx=m.x/DX-0.5, gy=m.y/DX-0.5;
        const int i0=(int)floor(gx)-1, j0=(int)floor(gy)-1;
        for(int di=0;di<=3;++di) for(int dj=0;dj<=3;++dj){
            const int i=i0+di, j=j0+dj;
            if(!inside(i,j)||solid[id(i,j)]) continue;
            const double s=peskin(gx-i)*peskin(gy-j);
            Fx[id(i,j)]+=s*m.fx; Fy[id(i,j)]+=s*m.fy;
        }
    }
    const double a0=0.35;
    for(const auto &u:units){
        const int i0=max(1,(int)floor((u.cx-1.2*R_UNIT)/DX));
        const int i1=min(NX-2,(int)ceil((u.cx+1.2*R_UNIT)/DX));
        const int j0=max(1,(int)floor((u.cy-1.2*R_UNIT)/DX));
        const int j1=min(NY-2,(int)ceil((u.cy+1.2*R_UNIT)/DX));
        for(int i=i0;i<=i1;++i) for(int j=j0;j<=j1;++j){
            const double x=(i+0.5)*DX, y=(j+0.5)*DX;
            if(hypot(x-u.cx,y-u.cy)<=0.95*R_UNIT) alpha[id(i,j)]=max(alpha[id(i,j)],a0);
        }
    }
}
void interpolate_move(){
    for(auto &m:mk){
        if(m.tether){ m.x=m.x0; m.y=m.y0; continue; }
        const double gx=m.x/DX-0.5, gy=m.y/DX-0.5;
        const int i0=(int)floor(gx)-1, j0=(int)floor(gy)-1;
        double u=0,v=0,ssum=0;
        for(int di=0;di<=3;++di) for(int dj=0;dj<=3;++dj){
            const int i=i0+di, j=j0+dj;
            if(!inside(i,j)||solid[id(i,j)]) continue;
            const double s=peskin(gx-i)*peskin(gy-j);
            u+=s*ux[id(i,j)]; v+=s*uy[id(i,j)]; ssum+=s;
        }
        if(ssum>1e-12){u/=ssum; v/=ssum;}
        m.x+=u*DX; m.y+=v*DX;
        m.x=min(max(m.x,1.2*DX),(NX-1.2)*DX);
        m.y=min(max(m.y,1.15*DX),(NY-1.2)*DX);
    }
}
void collide_stream(){
    const double g=8.0*NU_L*UMAX_LU/((NY-2.0)*(NY-2.0));
    vector<double> nxt(NX*NY*Q,0.0);
    for(int i=0;i<NX;++i) for(int j=0;j<NY;++j){
        const int p=id(i,j);
        if(solid[p]){ for(int k=0;k<Q;++k) fout[p*Q+k]=fin[p*Q+opp[k]]; continue; }
        double r=0,mx=0,my=0;
        for(int k=0;k<Q;++k){ r+=fin[p*Q+k]; mx+=ex[k]*fin[p*Q+k]; my+=ey[k]*fin[p*Q+k]; }
        r=max(r,1e-12);
        const double fx=Fx[p]+g-alpha[p]*mx/r;
        const double fy=Fy[p]-alpha[p]*my/r;
        ux[p]=(mx+0.5*fx)/r; uy[p]=(my+0.5*fy)/r; rho[p]=r;
        for(int k=0;k<Q;++k){
            const double eu=ex[k]*ux[p]+ey[k]*uy[p];
            const double fsrc=(1.0-0.5/TAU_F)*w[k]*
                (3.0*((ex[k]-ux[p])*fx+(ey[k]-uy[p])*fy)+9.0*eu*(ex[k]*fx+ey[k]*fy));
            fout[p*Q+k]=fin[p*Q+k]-(fin[p*Q+k]-feq(k,r,ux[p],uy[p]))/TAU_F+fsrc;
        }
    }
    for(int i=0;i<NX;++i) for(int j=0;j<NY;++j) for(int k=0;k<Q;++k){
        int ni=i+ex[k], nj=j+ey[k];
        if(ni<0) ni+=NX; if(ni>=NX) ni-=NX;
        if(nj<0||nj>=NY) continue;
        nxt[id(ni,nj)*Q+k]+=fout[id(i,j)*Q+k];
    }
    fin.swap(nxt);
}
int count_alive_inter(){int n=0; for(auto &b:bonds) if(b.inter&&b.alive) ++n; return n;}
int count_still_anchored(){
    int n=0;
    for(size_t i=0;i<units.size();++i){
        if(!units[i].started_anchored) continue;
        bool hold=false;
        for(auto &b:bonds) if(b.wallb&&b.alive&&b.ua==(int)i) hold=true;
        if(hold) ++n;
    }
    return n;
}
int count_started_anchored(){int n=0; for(auto &u:units) if(u.started_anchored) ++n; return n;}
double mean_speed_reduction(){
    double u0=0,u1=0; int n0=0,n1=0;
    const double xc=0.40*NX*DX;
    for(int i=0;i<NX;++i) for(int j=1;j<NY-1;++j){
        const double x=(i+0.5)*DX, u=ux[id(i,j)]*DX/DT;
        if(x<xc-80e-6){u0+=u;++n0;}
        if(x>xc-20e-6&&x<xc+20e-6){u1+=u;++n1;}
    }
    if(n0<10||n1<10) return 0.0;
    return 1.0-(u1/n1)/max(1e-12,u0/n0);
}
void write_frame(const fs::path &dir,int step){
    ostringstream name; name<<"state_"<<setw(6)<<setfill('0')<<step<<".dat";
    ofstream out(dir/name.str());
    out<<"TITLE = \"WangIBM\"\nVARIABLES = \"X_um\",\"Y_um\",\"U\",\"V\",\"Phase\"\n";
    out<<"ZONE T=\"Field\", I="<<NX<<", J="<<NY<<", F=POINT\n";
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i){
        const int p=id(i,j); int ph=solid[p]?1:0;
        if(alpha[p]>0.05) ph=2;
        out<<i*DX*1e6<<' '<<j*DX*1e6<<' '<<ux[p]*DX/DT<<' '<<uy[p]*DX/DT<<' '<<ph<<'\n';
    }
    out<<"ZONE T=\"Units\", I="<<max(1,(int)units.size())<<", F=POINT\n";
    if(units.empty()) out<<"0 0 0 0 2\n";
    else for(auto &u:units) out<<u.cx*1e6<<' '<<u.cy*1e6<<" 0 0 2\n";
    int nb=0; for(auto &b:bonds) if(b.alive&&(b.inter||b.wallb)) ++nb;
    out<<"ZONE T=\"Bonds\", I="<<max(1,2*nb)<<", F=POINT\n";
    if(!nb) out<<"0 0 0 0 3\n";
    else for(auto &b:bonds) if(b.alive&&(b.inter||b.wallb)){
        out<<mk[b.a].x*1e6<<' '<<mk[b.a].y*1e6<<" 0 0 3\n";
        out<<mk[b.b].x*1e6<<' '<<mk[b.b].y*1e6<<" 0 0 3\n";
    }
}
int main(int argc,char **argv){
    const fs::path outdir=argc>1?fs::path(argv[1]):fs::path("WangDetach");
    const int cse=argc>2?atoi(argv[2]):1;
    const double E_pa=argc>3?atof(argv[3]):1.0e6;
    fs::create_directories(outdir);
    alloc_fluid(); init_poiseuille();
    if(cse==2) build_case2(); else if(cse==3) build_case3(); else build_case1();
    update_unit_centres(); connect_nearby();
    ofstream hist(outdir/"history.csv");
    hist<<"step,time_s,n_units,n_inter_bonds,n_anchored,n_started_anchored,mean_cx_um,mean_cy_um,frac_detached,u_reduction\n";
    cout<<"Wang IBM two-way  case="<<cse<<"  E="<<E_pa<<" Pa  units="<<units.size()
        <<"  markers="<<mk.size()<<"  bonds="<<bonds.size()<<"\n";
    for(int n=0;n<=LBM_STEPS;++n){
        update_unit_centres();
        if(n%20==0) connect_nearby();
        accumulate_forces(E_pa); spread_to_euler(); collide_stream(); interpolate_move();
        if(n%OUTPUT_EVERY==0 || n==LBM_STEPS){
            double mx=0,my=0; for(auto &u:units){mx+=u.cx; my+=u.cy;}
            mx/=max(1,(int)units.size()); my/=max(1,(int)units.size());
            const int anc0=count_started_anchored(), anc=count_still_anchored();
            const double frac=anc0?1.0-double(anc)/anc0:0.0;
            const double red=mean_speed_reduction();
            hist<<n<<','<<n*DT<<','<<units.size()<<','<<count_alive_inter()<<','
                <<anc<<','<<anc0<<','<<mx*1e6<<','<<my*1e6<<','<<frac<<','<<red<<'\n';
            hist.flush(); write_frame(outdir,n);
            cout<<"n="<<n<<" t="<<scientific<<n*DT<<" units="<<units.size()
                <<" bonds="<<count_alive_inter()<<" anchored="<<anc<<"/"<<anc0
                <<" detach="<<fixed<<setprecision(3)<<frac<<" u_red="<<red<<"\n";
        }
    }
    ofstream sum(outdir/"summary.txt");
    sum<<"case = "<<cse<<"\nE_Pa = "<<E_pa<<"\nunits = "<<units.size()
        <<"\nmarkers = "<<mk.size()<<"\nfinal_anchored = "<<count_still_anchored()
        <<"\nstarted_anchored = "<<count_started_anchored()
        <<"\nfinal_inter_bonds = "<<count_alive_inter()
        <<"\nu_reduction = "<<mean_speed_reduction()<<"\n";
    cout<<"done -> "<<outdir<<"\n";
    return 0;
}
