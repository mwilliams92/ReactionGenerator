// BeamMC.cpp
// Standalone Monte Carlo for beam energy loss, straggling, multiple scattering,
// and optional two-body reaction recoils.
// Build: g++ -O2 -std=c++11 BeamMC.cpp -o BeamMC

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static const double PI = 3.14159265358979323846;
static const double AMU_MEV = 931.49410242;
static const double C_MM_PER_NS = 299.792458;

struct Element {
    int Z = 0, A = 0;
    double molarMass = 0.0;
    double fraction = 0.0;              // mass fraction
    std::string srimFile;
};

struct Material {
    double thicknessMm = 0.0;
    double density = 0.0;               // g/cm^3
    std::vector<Element> elements;
    std::string energyUnit = "keV";
    std::string stoppingUnit = "MeV_cm2_g";
};

struct BeamConfig {
    int Z = 0, A = 0;
    double energyMeV = 0.0;
    double energySpreadMeV = 0.0;       // Gaussian sigma
    double emittanceX = 0.0;            // geometric rms, mm mrad
    double emittanceY = 0.0;            // geometric rms, mm mrad
    double sigmaXmm = 1.0, sigmaYmm = 1.0;
    int particles = 10000;
    double stepMm = 0.01;
    unsigned long long seed = 0;
    std::string output = "beam_output.csv";
    double stragglingCoefficient = 0.02;
};

struct DegraderConfig : Material {
    bool enabled = false;
    std::string output;
};

struct ReactionConfig {
    int beamZ = 0, beamA = 0;
    int targetZ = 0, targetA = 0;
    int ejectileZ = 0, ejectileA = 0;
    int recoilZ = 0, recoilA = 0;
    double excitationMeV = 0.0;
    double qValueMeV = 0.0;              // optional; default 0 if omitted
    std::string output = "reaction_output.csv";
};

struct Point { double x=0, y=0, z=0; };
struct Vec3 { double x=0, y=0, z=1; };

static std::string trim(const std::string& s) {
    const std::string ws = " \t\r\n";
    size_t a = s.find_first_not_of(ws), b = s.find_last_not_of(ws);
    return a == std::string::npos ? "" : s.substr(a, b-a+1);
}

static std::string stripComment(std::string s) {
    size_t p = s.find('#');
    if (p != std::string::npos) s.resize(p);
    p = s.find(';');
    if (p != std::string::npos) s.resize(p);
    return trim(s);
}

static std::map<std::string,std::string> readKeyValues(const std::string& file,
                                                        std::vector<std::string>& elementLines) {
    std::ifstream in(file.c_str());
    if (!in) throw std::runtime_error("Cannot open " + file);
    std::map<std::string,std::string> kv;
    std::string line;
    while (std::getline(in,line)) {
        line = stripComment(line);
        if (line.empty()) continue;
        size_t p = line.find('=');
        if (p == std::string::npos) p = line.find(':');
        if (p != std::string::npos) {
            kv[trim(line.substr(0,p))] = trim(line.substr(p+1));
        } else if (line.find("element") == 0 || line.find("Element") == 0) {
            elementLines.push_back(line);
        }
    }
    return kv;
}

static double getD(const std::map<std::string,std::string>& kv, const std::string& k,
                   double def=0.0) {
    std::map<std::string,std::string>::const_iterator it=kv.find(k);
    return it==kv.end() ? def : std::stod(it->second);
}
static int getI(const std::map<std::string,std::string>& kv, const std::string& k, int def=0) {
    return static_cast<int>(getD(kv,k,def));
}
static std::string getS(const std::map<std::string,std::string>& kv, const std::string& k,
                        const std::string& def="") {
    std::map<std::string,std::string>::const_iterator it=kv.find(k);
    return it==kv.end() ? def : it->second;
}

static std::vector<std::string> words(const std::string& s) {
    std::istringstream is(s); std::vector<std::string> v; std::string x;
    while (is >> x) v.push_back(x);
    return v;
}

static std::vector<Element> parseElements(const std::vector<std::string>& lines) {
    std::vector<Element> out;
    for (size_t i=0;i<lines.size();++i) {
        std::string s=lines[i]; size_t p=s.find('=');
        if (p!=std::string::npos) s=s.substr(p+1);
        else { p=s.find_first_of(" \t"); if (p!=std::string::npos) s=s.substr(p+1); }
        std::vector<std::string> w=words(s);
        if (w.size()<5) throw std::runtime_error("Element line needs: Z A molar_mass mass_fraction srim_file");
        Element e; e.Z=std::stoi(w[0]); e.A=std::stoi(w[1]); e.molarMass=std::stod(w[2]);
        e.fraction=std::stod(w[3]); e.srimFile=w[4]; out.push_back(e);
    }
    if (out.empty()) throw std::runtime_error("No element lines found");
    return out;
}

static Material readMaterial(const std::string& file) {
    std::vector<std::string> elines; std::map<std::string,std::string> kv=readKeyValues(file,elines);
    Material m; m.thicknessMm=getD(kv,"thickness_mm"); m.density=getD(kv,"density_g_cm3");
    m.energyUnit=getS(kv,"srim_energy_unit","keV");
    m.stoppingUnit=getS(kv,"srim_stopping_unit","MeV_cm2_g");
    m.elements=parseElements(elines);
    if (m.thicknessMm<=0 || m.density<=0) throw std::runtime_error("Invalid material thickness or density in "+file);
    double sum=0; for (size_t i=0;i<m.elements.size();++i) sum+=m.elements[i].fraction;
    if (sum<=0) throw std::runtime_error("Element fractions sum to zero in "+file);
    for (size_t i=0;i<m.elements.size();++i) m.elements[i].fraction/=sum;
    return m;
}

static BeamConfig readBeam(const std::string& file) {
    std::vector<std::string> e; std::map<std::string,std::string> k=readKeyValues(file,e);
    BeamConfig b; b.Z=getI(k,"Z"); b.A=getI(k,"A"); b.energyMeV=getD(k,"energy_MeV");
    b.energySpreadMeV=getD(k,"energy_spread_MeV"); b.emittanceX=getD(k,"emittance_x_mm_mrad",getD(k,"emittance_mm_mrad"));
    b.emittanceY=getD(k,"emittance_y_mm_mrad",b.emittanceX); b.sigmaXmm=getD(k,"beam_sigma_x_mm",1.0);
    b.sigmaYmm=getD(k,"beam_sigma_y_mm",b.sigmaXmm); b.particles=getI(k,"particles",10000);
    b.stepMm=getD(k,"step_mm",0.01); b.seed=static_cast<unsigned long long>(getD(k,"seed",0));
    b.output=getS(k,"output",b.output); b.stragglingCoefficient=getD(k,"straggling_coefficient",0.02);
    if (b.Z<=0 || b.A<=0 || b.energyMeV<=0) throw std::runtime_error("Invalid beam Z, A, or energy");
    if (b.seed==0) b.seed=std::random_device()();
    return b;
}

static DegraderConfig readDegrader(const std::string& file) {
    std::vector<std::string> e; std::map<std::string,std::string> k=readKeyValues(file,e);
    DegraderConfig d; d.enabled=getI(k,"enabled",0)!=0; d.thicknessMm=getD(k,"thickness_mm");
    d.density=getD(k,"density_g_cm3"); d.energyUnit=getS(k,"srim_energy_unit","keV");
    d.stoppingUnit=getS(k,"srim_stopping_unit","MeV_cm2_g"); d.elements=parseElements(e);
    double sum=0; for(size_t i=0;i<d.elements.size();++i) sum+=d.elements[i].fraction;
    for(size_t i=0;i<d.elements.size();++i) d.elements[i].fraction/=sum;
    if (d.enabled && (d.thicknessMm<=0 || d.density<=0)) throw std::runtime_error("Invalid enabled degrader");
    return d;
}

static ReactionConfig readReaction(const std::string& file) {
    std::vector<std::string> e; std::map<std::string,std::string> k=readKeyValues(file,e);
    ReactionConfig r; r.beamZ=getI(k,"beam_Z"); r.beamA=getI(k,"beam_A"); r.targetZ=getI(k,"target_Z"); r.targetA=getI(k,"target_A");
    r.ejectileZ=getI(k,"ejectile_Z"); r.ejectileA=getI(k,"ejectile_A"); r.recoilZ=getI(k,"recoil_Z"); r.recoilA=getI(k,"recoil_A");
    r.excitationMeV=getD(k,"recoil_excitation_MeV",0); r.qValueMeV=getD(k,"q_value_MeV",0); r.output=getS(k,"output",r.output);
    if(r.beamA<=0||r.targetA<=0||r.ejectileA<=0||r.recoilA<=0) throw std::runtime_error("Incomplete reaction input");
    return r;
}

struct Table { std::vector<double> E, S; };
static Table loadSRIM(const std::string& file, const std::string& eunit) {
    std::ifstream in(file.c_str()); if(!in) throw std::runtime_error("Cannot open SRIM file "+file);
    Table t; std::string line;
    while(std::getline(in,line)) {
        for(size_t i=0;i<line.size();++i) if(line[i]==',' || line[i]==';' || line[i]=='\t') line[i]=' ';
        std::istringstream is(line); std::vector<double> nums; double x;
        while(is>>x) nums.push_back(x);
        if(nums.size()>=2 && nums[0]>0 && nums.back()>0) {
            double E=nums[0]; if(eunit=="keV"||eunit=="kev") E/=1000.0;
            else if(eunit=="eV"||eunit=="ev") E/=1e6;
            t.E.push_back(E); t.S.push_back(nums.back());
        }
    }
    if(t.E.size()<2) throw std::runtime_error("Could not find two numeric columns in "+file);
    std::vector<size_t> idx(t.E.size()); for(size_t i=0;i<idx.size();++i)idx[i]=i;
    std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return t.E[a]<t.E[b];});
    Table q; for(size_t i:idx){ if(q.E.empty()||t.E[i]>q.E.back()){q.E.push_back(t.E[i]);q.S.push_back(t.S[i]);} } return q;
}

static double interp(const Table& t,double E) {
    if(E<=t.E.front()) return t.S.front(); if(E>=t.E.back()) return t.S.back();
    std::vector<double>::const_iterator it=std::upper_bound(t.E.begin(),t.E.end(),E); size_t i=it-t.E.begin()-1;
    double f=(E-t.E[i])/(t.E[i+1]-t.E[i]); return t.S[i]+f*(t.S[i+1]-t.S[i]);
}

struct MaterialModel { Material m; std::vector<Table> tables; double x0gcm2=1e30; };
static double atomicX0(int Z,double A) { return 716.4*A/(Z*(Z+1)*std::log(287.0/std::sqrt((double)Z))); }
static MaterialModel buildModel(const Material& m) {
    MaterialModel mm; mm.m=m; double inv=0;
    for(size_t i=0;i<m.elements.size();++i){ mm.tables.push_back(loadSRIM(m.elements[i].srimFile,m.energyUnit)); inv+=m.elements[i].fraction/atomicX0(m.elements[i].Z,m.elements[i].A); }
    mm.x0gcm2=1.0/inv; return mm;
}

static double stoppingMeVmm(const MaterialModel& mm,double E,double rho) {
    double s=0; for(size_t i=0;i<mm.m.elements.size();++i) {
        double v=interp(mm.tables[i],std::max(E,1e-9));
        if(mm.m.stoppingUnit=="MeV_cm2_mg") v*=1000.0;
        else if(mm.m.stoppingUnit=="keV_um") v*=1.0; // 1 keV/um = 1 MeV/mm
        s+=mm.m.elements[i].fraction*v;
    }
    return s*rho/10.0; // MeV cm^2/g * g/cm^3 -> MeV/cm -> MeV/mm
}

static Vec3 unit(Vec3 v){double n=std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z);return {v.x/n,v.y/n,v.z/n};}
static Vec3 add(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static Vec3 scale(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static double betaFromE(double E,double mass){return std::sqrt(std::max(0.0,2*E/mass));}

struct Particle { Point p; Vec3 dir; double E=0; double ax=0, ay=0; };
static void propagate(Particle& p,const MaterialModel& mm,double length,double step,double zcharge,double straggle,std::mt19937_64& rng) {
    if(length<=0||p.E<=0) return; std::normal_distribution<double> gauss(0,1); double left=length;
    double mass=mm.m.elements[0].A*AMU_MEV;
    while(left>1e-12 && p.E>1e-9){ double ds=std::min(step,left); double dEdx=stoppingMeVmm(mm,p.E,mm.m.density); double dE=std::min(p.E,std::max(0.0,dEdx*ds));
        double sigE=std::sqrt(std::max(0.0,straggle*dE*std::max(p.E,1e-6))); p.E=std::max(0.0,p.E-dE+sigE*gauss(rng));
        p.p.x+=p.dir.x*ds; p.p.y+=p.dir.y*ds; p.p.z+=p.dir.z*ds;
        double beta=std::max(1e-5,betaFromE(p.E,mass)); double pMeV=std::max(1e-5,std::sqrt(2*mass*p.E));
        double xg=mm.m.density*ds/10.0; double theta=13.6*zcharge/(beta*pMeV)*std::sqrt(std::max(0.0,xg/mm.x0gcm2));
        if(theta>0){p.ax+=theta*gauss(rng)/std::sqrt(2.0);p.ay+=theta*gauss(rng)/std::sqrt(2.0);p.dir=unit({p.dir.x+p.ax*ds/length,p.dir.y+p.ay*ds/length,p.dir.z});}
        left-=ds;
    }
}

static Particle initialBeam(const BeamConfig& b,std::mt19937_64& rng) {
    Particle p; std::normal_distribution<double> g(0,1); p.p={b.sigmaXmm*g(rng),b.sigmaYmm*g(rng),0};
    double sx=b.emittanceX/b.sigmaXmm/1000.0, sy=b.emittanceY/b.sigmaYmm/1000.0;
    p.dir=unit({sx*g(rng),sy*g(rng),1.0}); p.E=std::max(0.0,b.energyMeV+b.energySpreadMeV*g(rng)); return p;
}

static Particle reactionRecoil(const Particle& in,const ReactionConfig& r,std::mt19937_64& rng) {
    Particle out; double ma=r.beamA*AMU_MEV, mt=r.targetA*AMU_MEV, mr=r.recoilA*AMU_MEV, me=r.ejectileA*AMU_MEV;
    double Ecm=in.E*mt/(ma+mt), avail=Ecm+r.qValueMeV-r.excitationMeV; if(avail<=0) return out;
    double mu=me*mr/(me+mr), pf=std::sqrt(2*mu*avail), pcm=std::sqrt(2*ma*in.E);
    std::normal_distribution<double> g(0,1); double c=2.0*std::generate_canonical<double,64>(rng)-1, s=std::sqrt(1-c*c), ph=2*PI*std::generate_canonical<double,64>(rng);
    Vec3 w={s*std::cos(ph),s*std::sin(ph),c}; Vec3 z=unit(in.dir), ref=(std::fabs(z.z)<0.9?Vec3{0,0,1}:Vec3{0,1,0});
    Vec3 ex=unit({ref.y*z.z-ref.z*z.y,ref.z*z.x-ref.x*z.z,ref.x*z.y-ref.y*z.x}); Vec3 ey={z.y*ex.z-z.z*ex.y,z.z*ex.x-z.x*ex.z,z.x*ex.y-z.y*ex.x};
    Vec3 pfv=scale(add(add(scale(ex,w.x),scale(ey,w.y)),scale(z,w.z)),pf); Vec3 P=scale(z,pcm); Vec3 pre=add(P,scale(pfv,-1));
    double vcm=pcm/(ma+mt); out.p=in.p; out.dir=unit(add(scale(pre,1.0/mr),scale(z,0))); out.E=pre.x*pre.x+pre.y*pre.y+pre.z*pre.z; out.E/=2*mr; // momentum is in MeV/c; non-relativistic kinetic energy
    return out;
}

static void checkReactionBeam(const BeamConfig& b,const ReactionConfig& r){
    if(b.Z!=r.beamZ||b.A!=r.beamA) std::cerr<<"Warning: beam_input and input_reaction beam species differ\n";
}

int main(int argc,char**argv){
    if(argc!=2 || (std::string(argv[1])!="beam" && std::string(argv[1])!="reaction")){std::cerr<<"Usage: "<<argv[0]<<" beam|reaction\n";return 2;}
    try{
        BeamConfig b=readBeam("beam_input.txt"); Material target=readMaterial("target_input.dat"); MaterialModel tm=buildModel(target);
        DegraderConfig dc=readDegrader("degrader_input.txt"); MaterialModel dm; if(dc.enabled) dm=buildModel(dc);
        std::mt19937_64 rng(b.seed); std::ofstream out; std::string mode=argv[1];
        ReactionConfig rc; if(mode=="reaction"){rc=readReaction("input_reaction.txt");checkReactionBeam(b,rc);out.open("reaction_output.csv");out<<"E_recoil_MeV,theta_x_mrad,theta_y_mrad,x_mm,y_mm,z_mm\n";}else out.open(b.output.c_str()),out<<"E_out_MeV,theta_x_mrad,theta_y_mrad,x_mm,y_mm,z_mm\n";
        for(int n=0;n<b.particles;++n){ Particle p=initialBeam(b,rng); double depth=(mode=="reaction"?target.thicknessMm*std::generate_canonical<double,64>(rng):target.thicknessMm); depth=std::max(std::min(depth,target.thicknessMm),std::min(b.stepMm,target.thicknessMm)); propagate(p,tm,depth,b.stepMm,b.Z,b.stragglingCoefficient,rng);
            if(mode=="beam"){if(dc.enabled)propagate(p,dm,dc.thicknessMm,b.stepMm,b.Z,b.stragglingCoefficient,rng);out<<std::setprecision(10)<<p.E<<","<<std::atan2(p.dir.x,p.dir.z)*1000<<","<<std::atan2(p.dir.y,p.dir.z)*1000<<","<<p.p.x<<","<<p.p.y<<","<<p.p.z<<"\n";}
            else {Particle recoil=reactionRecoil(p,rc,rng); double remain=target.thicknessMm-depth; if(recoil.E>0){double dz=std::max(1e-9,recoil.dir.z);propagate(recoil,tm,remain/dz,b.stepMm,rc.recoilZ,b.stragglingCoefficient,rng);if(dc.enabled)propagate(recoil,dm,dc.thicknessMm/std::max(1e-9,recoil.dir.z),b.stepMm,rc.recoilZ,b.stragglingCoefficient,rng);}out<<std::setprecision(10)<<recoil.E<<","<<std::atan2(recoil.dir.x,recoil.dir.z)*1000<<","<<std::atan2(recoil.dir.y,recoil.dir.z)*1000<<","<<recoil.p.x<<","<<recoil.p.y<<","<<recoil.p.z<<"\n";}
        }
        std::cerr<<"Completed "<<b.particles<<" particles using seed "<<b.seed<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;} return 0;
}
