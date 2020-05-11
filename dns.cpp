 /* nsintegrator.cpp: time-integration classes for spectral Navier-Stokes DNS
  */

#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/chebyshev.h"
#include "channelflow/vector.h"
#include "channelflow/periodicfunc.h"
//#include "orrsommfunc.h"
//#include <fstream> // tmp debugging need

using namespace std;
using namespace channelflow;


//*************Functions for calc of stress moments by interpolation****************************************************************//
void reverseforce(int i,FlowField& force);
void third_order_stress_moments(int i,FlowField& force, FlowField& Dxx1, FlowField& Dyy1, FlowField& Dzz1);
void rotational_stress_tensor(int i, FlowField& force, FlowField& dx1, FlowField& dy1, FlowField& dz1);
void symmetric_stress_tensor(int i,FlowField& u,FlowField& utmp,FlowField& force,FlowField& sx1, const Vector& x_grid, const Vector&  z_grid);
//**************End of Functions for calc of stress moments by interpolation**********************************************************//

//**************Functions for calc of stress moments by spectral method***************************************************************//
//void reverseforce_spectral_method(int i, FlowField& force);
//void third_order_stress_moments_spectral(int i,FlowField& Dxx,FlowField& Dyy,FlowField& Dzz,FlowField& Dxx1,FlowField& Dyy1,FlowField& Dzz1);
//extern void rotational_stress_tensor_spectral(int i,FlowField& dx, FlowField& dy, FlowField& dz, FlowField& dx1, FlowField& dy1, FlowField& dz1);
//void symmetric_stress_tensor_spectral_method(int i,FlowField& u,FlowField& sx,FlowField& sx1,const Vector& x_grid, const Vector&  z_grid);*/
//**************End Functions for calc of stress moments by spectral method***************************************************************//

namespace channelflow {

const Real EPSILON = 1e-11;

TimeStep::TimeStep(Real dt, Real dtmin, Real dtmax, Real dT,
		   Real CFLmin, Real CFLmax) 
  :
  n_(int(dT/dt + EPSILON)),
  nmin_(int(dT/dtmax + EPSILON)),
  nmax_(int(dT/dtmin + EPSILON)),
  dT_(dT),
  CFLmin_(CFLmin),
  CFL_((CFLmax+CFLmin)/2), // will take on meaninful value after first adjust
  CFLmax_(CFLmax)
{
  assert(dt>0 && dt<=dT);
  assert(dt>=dtmin && dt<=dtmax);
  assert(n_>=nmin_ &&  n_<=nmax_);
  assert(nmin_ > 0);
}


// relations
// n*dt = dT
// 
bool TimeStep::adjust(Real CFL, bool verbose) {
  CFL_ = CFL;
  if (CFLmin_ <= CFL && CFL <= CFLmax_)
    return false;
  
  // Potential new value of n. Potential new value of CFL is kept in CFL.
  int n = n_;   

  // Decrease CFL by increasing n, decreasing dt = T/n
  // CFL' == CFL * dt'/dt == CFL * n/n' == CFL * n/(n+1);
  while (CFL > CFLmax_ && n < nmax_) {
    CFL *= Real(n)/Real(n+1);
    ++n;
  }

  // Increase CFL by decreasing n, increasing dt = T/n
  // CFL' == CFL * dt'/dt == CFL * n/n' == CFL * n/(n-1);
  while (CFL < CFLmin_ && n > nmin_) {
    CFL *= Real(n)/Real(n-1);
    --n;
  }

  // Check to see if above loops exited before getting CFL in proper range
  // or if somehow n got bigger than nmax (latter should never happen).
  if (CFL > CFLmax_ || n > nmax_) {
    CFL *= Real(n)/Real(nmax_);
    n = nmax_;
    cerr << "TimeStep::adjust(CFL) : dt bottomed out at\n" 
	 << " dt  == " << dT_/n << endl
	 << " CFL == " << CFL  << endl
	 << " n   == " << n <<endl;
  }
  else if (CFL < CFLmin_ || n < nmin_) {
    CFL *= Real(n)/Real(nmin_);
    n = nmin_;
    cerr << "TimeStep::adjust(CFL) : dt topped out at\n"
	 << " dt  == " << dT_/n << endl
	 << " CFL == " << CFL  << endl
	 << " n   == " << n << endl;
  }

  // If final choice for n differs from original n_, reset internal values
  bool adjustment = (n == n_) ? false : true;
  if (adjustment && verbose) {
    cerr << "TimeStep::adjust(CFL) { " << endl;
    cerr << "   n : " << n_ << " -> " << n << endl;
    cerr << "  dt : " << dT_/n_ << " -> " << dT_/n << endl;
    cerr << " CFL : " << CFL_ << " -> " << (n_*CFL_)/n << endl;
    cerr << "}" << endl;
    CFL_ *= Real(n_)/Real(n);
    n_ = n;
  }
  return adjustment;
}
    
Real TimeStep::CFL() const {return CFL_;}
Real TimeStep::dt() const {return dT_/n_;}
Real TimeStep::dT() const {return dT_;}
TimeStep::operator Real() const {return dT_/n_;}
int TimeStep::n() const {return n_;}

//====================================================================
DNS::DNS() 
  :
  main_algorithm_(0),
  init_algorithm_(0)
{}

DNS::DNS(FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt, 
	 const DNSFlags& flags, Real t)
  :
  main_algorithm_(0),
  init_algorithm_(0) 
{


 main_algorithm_ = newAlgorithm(u, Ubase, nu, dt, flags, t);

  if (!main_algorithm_->full() && 
      flags.initstepping != flags.timestepping) {
    DNSFlags initflags = flags;
    initflags.timestepping = flags.initstepping;
    init_algorithm_ = newAlgorithm(u,Ubase,nu,dt,initflags,t);

    // Safety check
    if (init_algorithm_->Ninitsteps() != 0) 
      cerr << "DNS::DNS(u, Ubase, nu, dt, flags, t) :\n" << flags.initstepping
	   << " can't initialize " << flags.timestepping
	   << " since it needs initialization itself.\n";
  }
}
                   

DNS::DNS(const DNS& dns) 
  :
  main_algorithm_(dns.main_algorithm_ ? dns.main_algorithm_->clone() : 0),
  init_algorithm_(dns.init_algorithm_ ? dns.init_algorithm_->clone() : 0)
{}
 
DNSAlgorithm* DNS::newAlgorithm(FlowField& u, const ChebyCoeff& Ubase, 
				Real nu, Real dt, const DNSFlags& flags,
				Real t) { 

  DNSAlgorithm* alg = 0;
  switch (flags.timestepping) {
  case CNFE1:
  case SBDF1:
  case SBDF2:
  case SBDF3:
  case SBDF4:
    alg = new MultistepDNS(u, Ubase, nu, dt, flags, t);
    break;
  case CNRK2:
    alg = new RungeKuttaDNS(u, Ubase, nu, dt, flags, t);
    break;
  case SMRK2:
  case CNAB2:  
    alg = new CNABstyleDNS(u, Ubase, nu, dt, flags, t);
    break;
  default: 
    cerr << "DNS::newAlgorithm : algorithm " << flags.timestepping 
	 << " is unimplemented" << endl;
    exit(1);
  }
// cout<< " in newAlgorithm "<<alg<<endl;
 return alg;


}

DNS::~DNS() {
  delete main_algorithm_;
  delete init_algorithm_;
}

DNS& DNS::operator=(const DNS& dns) {
  delete main_algorithm_;
  delete init_algorithm_;
  main_algorithm_ = dns.main_algorithm_ ? dns.main_algorithm_->clone() : 0;
  init_algorithm_ = dns.init_algorithm_ ? dns.init_algorithm_->clone() : 0;
  return *this;
}  

void DNS::advance(FlowField& u, FlowField& q, int Nsteps) {
  assert(main_algorithm_);

  // Error check
  if (!main_algorithm_->full() && !init_algorithm_) {
    cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
	 << "and the initialization algorithm is not set. This should not be\n"
	 << "possible. Please submit a bug report (see documentation)." 
	 << endl;
    exit(1);
  }
  int n=0;
  while (!main_algorithm_->full() && n<Nsteps) {
    init_algorithm_->advance(u,q,1);
    main_algorithm_->push(u);
    if (main_algorithm_->full()) {
      delete init_algorithm_;
      init_algorithm_ = 0;
    }
    ++n;
  }
  main_algorithm_->advance(u,q,Nsteps-n);
}

//void DNS::reset() {
//}

void DNS::reset_dt(Real dt) {
  assert(main_algorithm_);
  main_algorithm_->reset_dt(dt);

  DNSFlags mainflags = main_algorithm_->flags();
  if (!main_algorithm_->full() && 
      mainflags.initstepping != mainflags.timestepping) {

    DNSFlags initflags = mainflags;
    initflags.timestepping = mainflags.initstepping;

    // An initialization algorithm needs only u's parameters at construction,
    // not its data, so we can construct it from a zero-valued u of right size
    FlowField u(main_algorithm_->Nx(), main_algorithm_->Ny(), 
		main_algorithm_->Nz(), 3,
		main_algorithm_->Lx(), main_algorithm_->Lz(),
		main_algorithm_->a(),  main_algorithm_->b());
		

    init_algorithm_ = newAlgorithm(u, main_algorithm_->Ubase(), 
				   main_algorithm_->nu(), dt, initflags, 
				   main_algorithm_->time());

    // Safety check
    if (init_algorithm_->Ninitsteps() != 0) 
      cerr << "DNS::DNS(u, Ubase, nu, dt, flags, t) :\n"
	   << mainflags.initstepping  << " can't initialize " 
	   << mainflags.timestepping
	   << " since it needs initialization itself.\n";
  }
  // No need for safety check on init_algorithm, has already been ok'd in ctor
}

// The mindless hassle of wrapper classes in C++ follows  
void DNS::reset_time(Real t) {
  assert(main_algorithm_);
  main_algorithm_->reset_time(t);
  if (init_algorithm_)
    init_algorithm_->reset_time(t);
}
void DNS::reset_dPdx(Real dPdx) {
  assert(main_algorithm_);
  main_algorithm_->reset_dPdx(dPdx);
  if (init_algorithm_)
    init_algorithm_->reset_dPdx(dPdx);

}
void DNS::reset_Ubulk(Real Ubulk) {
  assert(main_algorithm_);
  main_algorithm_->reset_Ubulk(Ubulk);
  if (init_algorithm_)
    init_algorithm_->reset_Ubulk(Ubulk);
}
//void DNS::reset_uj(const FlowField& uj, int j) {
//assert(main_algorithm_);
//main_algorithm_->reset_uj(uj, j);
//}
bool DNS::push(const FlowField& u) {
  assert(main_algorithm_);
  return main_algorithm_->push(u);
}
bool DNS::full() const {
  assert(main_algorithm_);
  return main_algorithm_->full();
}
int DNS::order() const { 
  assert(main_algorithm_);
  return main_algorithm_->order();
}
int DNS::Ninitsteps() const {
  assert(main_algorithm_);
  return main_algorithm_->Ninitsteps();
}
Real DNS::dt() const {
  assert(main_algorithm_);
  return main_algorithm_->dt();
}
Real DNS::CFL() const {
  assert(main_algorithm_);
  return main_algorithm_->CFL();
}
Real DNS::time() const {
  assert(main_algorithm_);
  return main_algorithm_->time();
}
Real DNS::dPdx() const {
  assert(main_algorithm_);
  return main_algorithm_->dPdx();
}
Real DNS::Ubulk() const {
  assert(main_algorithm_);
  return main_algorithm_->Ubulk();
}
Real DNS::dPdxRef() const {
  assert(main_algorithm_);
  return main_algorithm_->dPdxRef();
}
Real DNS::UbulkRef() const {  // the bulk velocity enforced during integ.
  assert(main_algorithm_);
  return main_algorithm_->UbulkRef();
}
const DNSFlags& DNS::flags() const {
  assert(main_algorithm_);
  return main_algorithm_->flags();
}
TimeStepMethod DNS::timestepping() const {
  assert(main_algorithm_);
  return main_algorithm_->timestepping();
}

/*********************************
void DNS::printStack() const {
  assert(main_algorithm_);
  main_algorithm_->printStack();
}
***********************************/


//====================================================================

DNSAlgorithm::~DNSAlgorithm() {}

DNSAlgorithm::DNSAlgorithm() 
  :
  Nx_(0),
  Ny_(0),
  Nz_(0),
  Mx_(0),
  Mz_(0),
  Nyd_(0),
  kxd_max_(0),
  kzd_max_(0),
  Lx_(0),
  Lz_(0),
  a_(0),
  b_(0),
  flags_(), 
  order_(0),
  Ninitsteps_(0),
  nu_(0),
  dt_(0),
  t_(0),
  cfl_(0),
  dPdxRef_(0),
  dPdxAct_(0),
  UbulkRef_(0), 
  UbulkAct_(0), 
  UbulkBase_(0), 
  Ubase_(),
  Ubaseyy_(),
  tmp_(),
  uk_(),
  vk_(),
  wk_(),
  Pk_(),
  Pyk_(),
  Rxk_(),
  Ryk_(),
  Rzk_()
{}  

DNSAlgorithm::DNSAlgorithm(FlowField& u, const ChebyCoeff& Ubase,
			   Real nu, Real dt, const DNSFlags& flags, 
			   Real t0)
  :
  Nx_(u.Nx()),
  Ny_(u.numYmodes()),
  Nz_(u.Nz()),
  Mx_(u.numXmodes()),
  Mz_(u.numZmodes()),
  Nyd_(flags.dealias_y() ? 2*(u.numYmodes()-1)/3 + 1 : u.numYmodes()),
  //kxd_max_(flags.dealias_xz() ? 2*(u.kxmax()/3) : u.kxmax()),
  //kzd_max_(flags.dealias_xz() ? 2*(u.kzmax()/3) : u.kzmax()),
  kxd_max_(flags.dealias_xz() ? u.Nx()/3-1 : u.kxmax()),
  kzd_max_(flags.dealias_xz() ? u.Nz()/3-1 : u.kzmax()),
  Lx_(u.Lx()),
  Lz_(u.Lz()),
  a_(u.a()),
  b_(u.b()),
  flags_(flags),
  order_(0),
  Ninitsteps_(0),
  nu_(nu),
  dt_(dt),
  t_(t0),
  cfl_(0),
  dPdxRef_(0),
  dPdxAct_(0),
  UbulkRef_(0), 
  UbulkAct_(0), 
  UbulkBase_(0), 
  Ubase_(Ubase),
  Ubaseyy_(), 
  tmp_(),
  uk_(Nyd_,a_,b_,Spectral),
  vk_(Nyd_,a_,b_,Spectral),
  wk_(Nyd_,a_,b_,Spectral),
  Pk_(Nyd_,a_,b_,Spectral),
  Pyk_(Nyd_,a_,b_,Spectral),
  Rxk_(Nyd_,a_,b_,Spectral),
  Ryk_(Nyd_,a_,b_,Spectral),
  Rzk_(Nyd_,a_,b_,Spectral)
{
  //tmp_.setToZero();

  u.assertState(Spectral, Spectral);
  assert(u.vectorDim() == 3);
  // Set the aliased modes to zero
  Complex zero = 0.0 + 0.0*I;
  if (flags_.dealias_xz()) {
    for (int i=0; i<3; ++i)
      for (int mx=0; mx<Mx_; ++mx)
	for (int mz=0; mz<Mz_; ++mz) {
	  if (isAliasedMode(u.kx(mx), u.kz(mz)))
	    for (int ny=0; ny<Nyd_; ++ny)
	      u.cmplx(mx,ny,mz,i) = zero;
	    for (int ny=Nyd_; ny<Ny_; ++ny)
	    u.cmplx(mx,ny,mz,i) = zero;
	}
  }
	    
  // Calculate Ubasey_ and related quantities if nonzero. Require that 
  // base flow is quadratic (solves an equilibrium problem).
  ChebyTransform trans(Ny_);
  ChebyCoeff Ubasey;
  if (Ubase_.length() != 0) {
    Ubase_.makeSpectral(trans);
    UbulkBase_ = Ubase_.mean();
    Ubasey   = diff(Ubase_);
    Ubaseyy_ = diff(Ubasey);
    Ubase_.makePhysical(trans);
    Ubasey.makePhysical(trans);
    // keep Ubaseyy_ as spectral for calc of R in advance() method..
  }
  else
    Ubase_.setState(Physical);


  // These methods require a 9d (3x3) tmp flowfield
  if (flags_.nonlinearity ==  Alternating || 
      flags_.nonlinearity ==  Alternating_ || 
      flags_.nonlinearity ==  Convection || 
      flags_.nonlinearity ==  Divergence || 
      flags_.nonlinearity ==  SkewSymmetric) 
    tmp_.resize(u.Nx(), u.Ny(), u.Nz(), 9, u.Lx(), u.Lz(), u.a(), u.b());
  else 
    tmp_.resize(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());

 cfl_ = u.CFLfactor(Ubase_);                                    
 cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_; 

  // Determine actual Ubulk and dPdx from initial data Ubase + u. 
  ChebyCoeff u00(Ny_,a_,b_,Spectral);
  for (int ny=0; ny<Ny_; ++ny)
    u00[ny] = Re(u.cmplx(0,ny,0,0));
  ChebyCoeff du00dy = diff(u00);

  UbulkAct_ = UbulkBase_ + u00.mean();
  dPdxAct_  = nu*(du00dy.eval_b() - du00dy.eval_a())/(b_-a_);
  if (Ubase_.length() != 0)
    dPdxAct_  += nu*(Ubasey.eval_b() - Ubasey.eval_a())/(b_-a_);

  /****************************************************************
  // Old code
  // Set whichever reference value is being held const.
  if (flags_.constraint == BulkVelocity)
    UbulkRef_ = UbulkAct_;
  else {
    // But if dPdx is held const, and dPdxAct is very near zero,
    // set it to zero, to correct for numerical error in above calc.
    if (abs(dPdxAct_) < EPSILON)
      dPdxAct_ = 0.0;
    dPdxRef_ = dPdxAct_;
  }
  ********************************************************************/
  // Replaced with 
  if (flags_.constraint == BulkVelocity)
    UbulkRef_ = flags_.Ubulk;
  else {
    dPdxAct_ = flags_.dPdx;
    dPdxRef_ = flags_.dPdx;
  }
}

bool DNSAlgorithm::push(const FlowField& u) {return true;}
bool DNSAlgorithm::full() const {return true;}

//void DNSAlgorithm::reset_uj(const FlowField& uj, int j) {;}

void DNSAlgorithm::reset_time(Real t) {t_=t;}

void DNSAlgorithm::reset_dPdx(Real dPdx) {
  flags_.constraint = PressureGradient;
  flags_.dPdx = dPdx;
  flags_.Ubulk = 0.0;
  dPdxRef_ = dPdx;
  UbulkRef_ = 0.0;
}
void DNSAlgorithm::reset_Ubulk(Real Ubulk) {
  flags_.constraint = BulkVelocity;
  flags_.Ubulk = Ubulk;
  flags_.dPdx = 0.0;
  UbulkRef_ = Ubulk;
  dPdxRef_ = 0.0;
}

int DNSAlgorithm::Nx() const {return Nx_;}
int DNSAlgorithm::Ny() const {return Ny_;}
int DNSAlgorithm::Nz() const {return Nz_;}
Real DNSAlgorithm::Lx() const {return Lx_;}
Real DNSAlgorithm::Lz() const {return Lz_;}
Real DNSAlgorithm::a() const {return a_;}
Real DNSAlgorithm::b() const {return b_;}
Real DNSAlgorithm::dt() const {return dt_;}
Real DNSAlgorithm::nu() const {return nu_;}
Real DNSAlgorithm::CFL() const {return cfl_;}
Real DNSAlgorithm::time() const {return t_;}
Real DNSAlgorithm::dPdx() const {return dPdxAct_;}
Real DNSAlgorithm::dPdxRef() const {return dPdxRef_;}
Real DNSAlgorithm::Ubulk() const {return UbulkAct_;}
Real DNSAlgorithm::UbulkRef() const {return UbulkRef_;}
int DNSAlgorithm::order() const {return order_;}
int DNSAlgorithm::Ninitsteps() const {return Ninitsteps_;}
const DNSFlags& DNSAlgorithm::flags() const {return flags_;}
const ChebyCoeff&  DNSAlgorithm::Ubase() const {return Ubase_;}
TimeStepMethod DNSAlgorithm::timestepping() const {return flags_.timestepping;}
int DNSAlgorithm::kxmaxDealiased() const {return kxd_max_;}
int DNSAlgorithm::kzmaxDealiased() const {return kzd_max_;}
bool DNSAlgorithm::isAliasedMode(int kx, int kz) const {
  return (abs(kx) > kxd_max_ || (abs(kz) > kzd_max_)) ? true : false;
}
/***************************************************
void DNSAlgorithm::printStack() const {
  cout << "DNSAlgorithm::printStack()" << endl;
}
*************************************************/

// ====================================================================
// Multistep algorithms
MultistepDNS::MultistepDNS()
  :
  DNSAlgorithm()
{}

MultistepDNS::MultistepDNS(const MultistepDNS& dns) 
  :
  DNSAlgorithm(dns),
  eta_(dns.eta_),
  alpha_(dns.alpha_),
  beta_(dns.beta_),
  u_(dns.u_),
  f_(dns.f_)
{
  // Copy tausolvers
  tausolver_ = new TauSolver*[Mx_];       // new #1
  for (int mx=0; mx<Mx_; ++mx) {
    tausolver_[mx] = new TauSolver[Mz_];  // new #2
    for (int mz=0; mz<Mz_; ++mz) 
      tausolver_[mx][mz] = dns.tausolver_[mx][mz];
  }
}

MultistepDNS::MultistepDNS(FlowField& u, const ChebyCoeff& Ubase,
			   Real nu, Real dt, const DNSFlags& flags, Real t) 
			   
  :
  DNSAlgorithm(u,Ubase,nu,dt,flags, t)
{

  TimeStepMethod algorithm = flags.timestepping;
  switch (algorithm) {
  case CNFE1:
  case SBDF1:
    order_ = 1;
    eta_ = 1.0;
    alpha_.resize(order_);
    beta_.resize(order_);
    alpha_[0] = -1.0;
    beta_[0]  =  1.0; 
    break;
  case SBDF2: 
    order_ = 2;
    alpha_.resize(order_);
    beta_.resize(order_);
    eta_ = 1.5;
    alpha_[0] = -2.0; alpha_[1] =  0.5;
    beta_[0]  =  2.0;  beta_[1] = -1.0;
    break;
  case SBDF3: 
    order_ = 3;
    alpha_.resize(order_);
    beta_.resize(order_);
    eta_ = 11.0/6.0;
    alpha_[0] = -3.0;  alpha_[1] = 1.5; alpha_[2] = -1.0/3.0;
    beta_[0]  =  3.0;   beta_[1] = -3.0; beta_[2] = 1.0;
    break;
  case SBDF4: 
    order_ = 4;
    alpha_.resize(order_);
    beta_.resize(order_);
    eta_ = 25.0/12.0;
    alpha_[0] = -4.0; alpha_[1] =  3.0; alpha_[2] = -4.0/3.0; alpha_[3] = 0.25;
    beta_[0]  =  4.0;  beta_[1] = -6.0;  beta_[2] =  4.0;      beta_[3] = -1.0;
    break;
  default: 
    cerr << "MultistepDNS::MultistepDNS(un,Ubase,nu,dt,flags,t0)\n"
	 << "error: flags.timestepping == " << algorithm 
	 << "is a non-multistepping algorithm" << endl;
    exit(1);
  }

  // Configure tausolvers
  tausolver_ = new TauSolver*[Mx_];       // new #1
  for (int mx=0; mx<Mx_; ++mx) 
    tausolver_[mx] = new TauSolver[Mz_];  // new #2

  reset_dt(dt_);

  // Initialize arrays of previous u's and f's
  FlowField tmp(u);
  tmp.setToZero();
  
  u_.resize(order_);
  f_.resize(order_);
  for (int j=0; j<order_; ++j) {
    u_[j] = tmp;
    f_[j] = tmp;
  }
  if (order_ > 0) { // should always be true
    u_[0] = u;
    //navierstokesNL(u_[0], Ubase_, f_[0], tmp_, flags_.nonlinearity);
  }
  cfl_ = u.CFLfactor(Ubase_);
  cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;

  Ninitsteps_ = order_;
  countdown_ = Ninitsteps_;
}

MultistepDNS::~MultistepDNS() {
  if (tausolver_) {
    for (int mx=0; mx<Mx_; ++mx) 
      delete[] tausolver_[mx];  // undo new #2
    delete[] tausolver_;        // undo new #1
  }
  tausolver_ = 0;
}

DNSAlgorithm* MultistepDNS::clone() const {
  return new MultistepDNS(*this);
}

void MultistepDNS::reset_dt(Real dt) {
  cfl_ *= dt/dt_;
  //nu_ = nu;
  dt_ = dt;
  const Real c = 4.0*square(pi)*nu_;
  const int kxmax = tmp_.kxmax();
  const int kzmax = tmp_.kzmax();

  // This loop replaces the TauSolver objects at tausolver_[substep][mx][mz] 
  // with new TauSolver objects, with the given parameters.
  for (int mx=0; mx<Mx_; ++mx) {
    int kx = tmp_.kx(mx);
    for (int mz=0; mz<Mz_; ++mz) {
      int kz = tmp_.kz(mz);
      Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));

      // When using dealiasing, some modes get set to zero, rather than 
      // updated with momentum eqns. Don't initialize TauSolvers for these.
      if ((kx != kxmax && kz != kzmax) && 
	  (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
	
	tausolver_[mx][mz] = TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, 
				       nu_, Nyd_, flags_.taucorrection);
    }
  }
  // Start from beginning on initialization
  countdown_ = Ninitsteps_;
}

// This calculation follows Peyret section 4.5.1(b) pg 131.
void MultistepDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {

  const int kxmax = un.kxmax();
  const int kzmax = un.kzmax();
  for (int step=0; step<Nsteps; ++step) {

    //cout << "MultistepDNS::advance(u,q,N) stack : " << endl;

    if (order_ > 0) {
      u_[0] = un;
      navierstokesNL(u_[0], Ubase_, f_[0], tmp_, flags_.nonlinearity);
        }
    // Update each Fourier mode with time-stepping algorithm
    for (int mx=0; mx<Mx_; ++mx) {
      const int kx = un.kx(mx);

      for (int mz=0; mz<Mz_; ++mz) {
	const int kz = un.kz(mz);

	// Zero out the aliased modes and break to next kx,kz
	if ((kx == kxmax || kz == kzmax) || 
	    flags_.dealias_xz() && isAliasedMode(kx,kz)) {
	  for (int ny=0; ny<Nyd_; ++ny) {
	    u_[0].cmplx(mx,ny,mz,0) = 0.0;
	    u_[0].cmplx(mx,ny,mz,1) = 0.0;
	    u_[0].cmplx(mx,ny,mz,2) = 0.0;
	    qn.cmplx(mx,ny,mz,0) = 0.0;
	  }
	  break;
	}

	// For non-aliased modes
	Rxk_.setToZero();
	Ryk_.setToZero();
	Rzk_.setToZero();

	for (int j=0; j<order_; ++j) { 
	  const Real a = -alpha_[j]/dt_;
	  const Real b = -beta_[j];
	  for (int ny=0; ny<Nyd_; ++ny) {
	    Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
	    Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
	    Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
	  }
	}

	// Solve the tau solutions
	if (kx!=0 || kz!=0) 
	  tausolver_[mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);

	else { // kx,kz == 0,0

	  if (Ubaseyy_.length() > 0)
	    for (int ny=0; ny<Ny_; ++ny)
	      Rxk_.re[ny] += nu_*Ubaseyy_[ny]; // Rx has additional term

	  if (flags_.constraint == PressureGradient) {
	    // pressure is supplied, put on RHS of tau eqn
	    Rxk_.re[0] -= dPdxRef_;  

	    // Solve the tau equations
	    tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);

	    // Bulk vel is free variable determined from soln of tau eqn
	    UbulkAct_ = UbulkBase_ + uk_.re.mean(); 
	    dPdxAct_ = dPdxRef_;
	  }
	  else { // const bulk velocity
	    // bulk velocity is supplied, put on RHS of tau eqn
	    // dPdxAct (i.e. at next time step is solved for)
	    // constraint:    UbulkBase + mean(u) = UbulkRef.
	    tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, 
					Rxk_, Ryk_, Rzk_, 
					UbulkRef_ - UbulkBase_);
	    UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
	    //UbulkAct_ = UbulkRef_; 
	  }
	}

	// Load solutions back into the external 3d data arrays.
	// Because of FFTW complex symmetries
	// The 0,0 mode must be real.
	// For Nx even, the kxmax,0 mode must be real
	// For Nz even, the 0,kzmax mode must be real
	// For Nx,Nz even, the kxmax,kzmax mode must be real
	if ((kx == 0 && kz == 0) ||
	    (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
	    (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
	    (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {

	  for (int ny=0; ny<Nyd_; ++ny) {
	    un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
	    un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
	    un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
	    qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
	  }
	}
	// The normal case, for general kx,kz
	else
	  for (int ny=0; ny<Nyd_; ++ny) {
	    un.cmplx(mx,ny,mz,0) = uk_[ny];
	    un.cmplx(mx,ny,mz,1) = vk_[ny];
	    un.cmplx(mx,ny,mz,2) = wk_[ny];
	    qn.cmplx(mx,ny,mz,0) = Pk_[ny];
	  }

	// And now set the y-aliased modes to zero.
	for (int ny=Nyd_; ny<Ny_; ++ny) {
	  un.cmplx(mx,ny,mz,0) = 0.0;
	  un.cmplx(mx,ny,mz,1) = 0.0;
	  un.cmplx(mx,ny,mz,2) = 0.0;
	  qn.cmplx(mx,ny,mz,0) = 0.0;
	}
      }
    }
    
    // The solution is stored in un. Shift entire u and f arrays in time
    // and push un into u_[0]. Ie shift u_[K] <- u_[K-1] <- ... <- u_[0] <- un
    for (int j=order_-1; j>0; --j) {
      swap(f_[j], f_[j-1]);
      swap(u_[j], u_[j-1]);
    }

    t_ += dt_;

    if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll) 
      cout << t_ << ' ' << flush;
    else if (flags_.verbosity == PrintTicks) 
      cout << '.' << flush;
  }
  
  cfl_ = u_[0].CFLfactor(Ubase_);
  cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
    
  // If using dealiasing, set flag in FlowField that compactifies binary IO 
  un.setDealiased(flags_.dealias_xz()); 
  qn.setDealiased(flags_.dealias_xz()); 

  if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
      flags_.verbosity == PrintTicks) 
    cout << endl;
    
  return;
}


/*************************************************************
void MultistepDNS::printStack() const {
  cout << "Multistep::printStack() {" << endl;
  cout << "        t == " << t_ << endl;
  cout << "countdown == " << countdown_ << endl;
  cout << "     full == " << full() << endl;

  for (int j=order_-1; j>=0; --j) 
    printf("j=%2d t=%5.2f L2(uj)=%13.10f L2(fj)=%13.10f\n",
	   j, t_-j*dt_, L2Norm(u_[j]), L2Norm(f_[j]));    
  cout << endl;
  cout << "}" << endl;
}
*****************************************************/

bool MultistepDNS::push(const FlowField& un) {

  // Let K = order-1. Arrays are then u_[0:K], f_[0:K]
  // Shift u_[K] <- u_[K-1] <- ... <- u_[0] <- un  
  for (int j=order_-1; j>0; --j) {
    swap(u_[j], u_[j-1]);
    swap(f_[j], f_[j-1]);
  }

  if (order_ > 0) {
    u_[0] = un;
    navierstokesNL(u_[0], Ubase_, f_[0], tmp_, flags_.nonlinearity);
    cfl_ = u_[0].CFLfactor(Ubase_);
    cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
  }

  t_ += dt_;
  --countdown_;

  return full();
}

bool MultistepDNS::full() const {
  return (countdown_ == 0) ? true : false;
}

/**************************************
void MultistepDNS::reset_uj(const FlowField& uj, int j) {
  if (j<0 || j>=order_) {
    cerr << "error in MultistepDNS::reset_uj(uj, j) : j = " << j 
	 << " is out of bounds" << endl;
    exit(1);
  }
  u_[j] = uj;
  navierstokesNL(u_[j], Ubase_, f_[j], tmp_, flags_.nonlinearity);
}
**********************************************/
// ==============================================================
// Runge-Kutta algorithms

RungeKuttaDNS::RungeKuttaDNS()
  :
  DNSAlgorithm(),
  Qj1_(),
  Qj_()
{}

RungeKuttaDNS::RungeKuttaDNS(const RungeKuttaDNS& dns)
  :
  DNSAlgorithm(dns),
  Nsubsteps_(dns.Nsubsteps_),
  Qj1_(dns.Qj1_),
  Qj_(dns.Qj_),
  A_(dns.A_),
  B_(dns.B_),
  C_(dns.C_)
{
  // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver arrays
  // and copy tausolvers from dns argument
  tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
  for (int j=0; j<Nsubsteps_; ++j) {
    tausolver_[j] = new TauSolver*[Mx_];       // new #2
    for (int mx=0; mx<Mx_; ++mx) {
      tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
      for (int mz=0; mz<Mz_; ++mz)
	tausolver_[j][mx][mz] = dns.tausolver_[j][mx][mz];
	
    }
  }
}

// This algorithm is described in "Spectral Methods for Incompressible Viscous
// Flow", Roger Peyret, Springer-Verlag Appl Math Sci series vol 148, 2002.
// section 4.5.2.c.2 "Three-stage scheme (RK3/CN)". I use C for his B'
RungeKuttaDNS::RungeKuttaDNS(FlowField& u, const ChebyCoeff& Ubase,
			     Real nu, Real dt, const DNSFlags& flags, Real t) 
  :
  DNSAlgorithm(u,Ubase,nu,dt,flags,t),
  Qj1_(u),
  Qj_(u)
{
  Qj1_.setToZero();
  Qj_.setToZero();

  TimeStepMethod algorithm = flags.timestepping;
  switch (algorithm) {
  case CNRK2: 
    order_ = 2;   
    Nsubsteps_ = 3;
    Ninitsteps_ = 0;
    A_.resize(Nsubsteps_);
    B_.resize(Nsubsteps_);
    C_.resize(Nsubsteps_);
    A_[0] = 0.0;     A_[1] = -5.0/9.0;  A_[2] = -153.0/128.0; // Peyret A
    B_[0] = 1.0/3.0; B_[1] = 15.0/16.0; B_[2] = 8.0/15.0;     // Peyret B
    C_[0] = 1.0/6.0; C_[1] = 5.0/24.0;  C_[2] = 1.0/8.0;      // Peyret B'
    D_[0] = 0.0;     D_[1] = 1.0/3.0;   D_[2] = 3.0/4.0;
    break;
  default: 
    cerr << "RungeKuttaDNS::RungeKuttaDNS(un,Ubase,nu,dt,flags,t0)\n"
	 << "error: flags.timestepping == " << algorithm 
	 << " is a non-runge-kutta algorithm" << endl;
    exit(1);
  }

  // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
  tausolver_ = new TauSolver**[Nsubsteps_];           // new #1
  for (int j=0; j<Nsubsteps_; ++j) {
    tausolver_[j] = new TauSolver*[Mx_];       // new #2
    for (int mx=0; mx<Mx_; ++mx) 
      tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
  }
  reset_dt(dt_);
}

RungeKuttaDNS::~RungeKuttaDNS() {
  if (tausolver_) {
    for (int j=0; j<Nsubsteps_; ++j) {
      for (int mx=0; mx<Mx_; ++mx) 
	delete[] tausolver_[j][mx];  // undo new #3
      delete[] tausolver_[j];        // undo new #2
    } 
    delete[] tausolver_;                   // undo new #1
  }
  tausolver_ = 0;
}

DNSAlgorithm* RungeKuttaDNS::clone() const {
  return new RungeKuttaDNS(*this);
}

void RungeKuttaDNS::reset_dt(Real dt) {
  cfl_ *= dt/dt_;
  //nu_ = nu;
  dt_ = dt;
  const Real c = 4.0*square(pi)*nu_;
  const int kxmax = tmp_.kxmax();
  const int kzmax = tmp_.kzmax();

  // This loop replaces the TauSolver objects at tausolver_[i][mx][mz] 
  // with new TauSolver objects configured with appropriate parameters
  for (int j=0; j<Nsubsteps_; ++j) {
    for (int mx=0; mx<Mx_; ++mx) {
      int kx = tmp_.kx(mx);
      for (int mz=0; mz<Mz_; ++mz) {
	int kz = tmp_.kz(mz);
	Real lambda = 1.0/(C_[j]*dt_) + c*(square(kx/Lx_)+square(kz/Lz_));

	if ((kx != kxmax || kz != kzmax) && 
	    (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))

	  tausolver_[j][mx][mz] = 
	    TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, nu_, Nyd_, 
		      flags_.taucorrection);
      }
    }
  }
}

void RungeKuttaDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {

  const int kxmax = un.kxmax();
  const int kzmax = un.kzmax();

  for (int n=0; n<Nsteps; ++n) {
    for (int j=0; j<Nsubsteps_; ++j) {
      
      FlowField& uj(un); // Store uj in un during substeps, reflect in notation
          // Efficient implementation of 
      // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
      // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
      // Q_j = Q_{j+1}
      Qj_ *= A_[j];
      navierstokesNL(uj, Ubase_, Qj1_, tmp_, flags_.nonlinearity);
      Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
   
      // Update each Fourier mode with time-stepping algorithm
      for (int mx=0; mx<Mx_; ++mx) {
	const int kx = un.kx(mx);

	for (int mz=0; mz<Mz_; ++mz) {
	  const int kz = un.kz(mz);

	  // Zero out the aliased modes
	  if ((kx == kxmax || kz == kzmax) || 
	      flags_.dealias_xz() && isAliasedMode(kx,kz)) {
	    for (int ny=0; ny<Nyd_; ++ny) {
	      un.cmplx(mx,ny,mz,0) = 0.0;
	      un.cmplx(mx,ny,mz,1) = 0.0;
	      un.cmplx(mx,ny,mz,2) = 0.0;
	      qn.cmplx(mx,ny,mz,0) = 0.0;
	    }
	    break;
	  }

	  Rxk_.setToZero();
	  Ryk_.setToZero();
	  Rzk_.setToZero();
	  
	  // Make the following assignments in prep for computation of RHS

	  // nu (uk,vk,wk) = nu ujn(0,1,2)
	  //            Pk = qn 

	  // Goal is to compute
	  // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
	  //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj

	  // Extract relevant Fourier modes of uj and qj for computations
	  for (int ny=0; ny<Nyd_; ++ny) {
	    uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
	    vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
	    wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
	    Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
	  }

	  // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
	  diff2(uk_, Rxk_, Pyk_);
	  diff2(vk_, Ryk_, Pyk_);
	  diff2(wk_, Rzk_, Pyk_);

	  // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
	  diff(Pk_, Pyk_); 

	  // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
	  const Real c = 
	    1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
	  const Real B_C = B_[j]/C_[j];
	  const Complex Dx = un.Dx(mx);
	  const Complex Dz = un.Dz(mz);

	  for (int ny=0; ny<Nyd_; ++ny) {
	    Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
	    Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
	    Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
	  }
	  
	  // Do the tau solutions
	  if (kx!=0 || kz!=0) 
	    tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);

	  else { // kx,kz == 0,0
	    // Rx has additional terms, nu Uyy at both t=j and t=j+1
	    const Real c = 2*nu_;
	    if (Ubaseyy_.length() > 0)
	      for (int ny=0; ny<Ny_; ++ny)
		Rxk_.re[ny] += c*Ubaseyy_[ny];
	  
	    if (flags_.constraint == PressureGradient) {
	      // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
	      Rxk_.re[0] -= dPdxAct_ + dPdxRef_;

	      // Solve the tau equations
	      tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);

	      // Bulk velocity is free variable on LHS solved by tau eqn
	      UbulkAct_ = UbulkBase_ + uk_.re.mean();
	      dPdxAct_ = dPdxRef_;
	    }
	    else { // const bulk velocity
	      // Add the previous time-step's -dPdx to the RHS. The next 
	      // timestep's dPdx term appears on LHS as unknown.
	      Rxk_.re[0] -= dPdxAct_; 

	      // Use tausolver with additional variable and constraint:
	      // free variable: dPdxAct at next time-step, 
	      // constraint:    UbulkBase + mean(u) = UbulkRef.
	      tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, 
					  Rxk_, Ryk_, Rzk_, 
					  UbulkRef_ - UbulkBase_);
	      UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
	      //UbulkAct_ = UbulkRef_;
	    }
	    // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
	    // Pk_.set(0, Complex(0.0, 0.0));
	  }

	  // Load solutions back into the external 3d data arrays.
	  // Because of FFTW complex symmetries
	  // The 0,0 mode must be real.
	  // For Nx even, the kxmax,0 mode must be real
	  // For Nz even, the 0,kzmax mode must be real
	  // For Nx,Nz even, the kxmax,kzmax mode must be real
	  if ((kx == 0 && kz == 0) ||
	      (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
	      (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
	      (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {

	    for (int ny=0; ny<Nyd_; ++ny) {
	      un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
	      un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
	      un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
	      qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
	    }
	  }
	  // The normal case, for general kx,kz
	  else
	    for (int ny=0; ny<Nyd_; ++ny) {
	      un.cmplx(mx,ny,mz,0) = uk_[ny];
	      un.cmplx(mx,ny,mz,1) = vk_[ny];
	      un.cmplx(mx,ny,mz,2) = wk_[ny];
	      qn.cmplx(mx,ny,mz,0) = Pk_[ny];
	    }

	  // And now set the y-aliased modes to zero.
	  for (int ny=Nyd_; ny<Ny_; ++ny) {
	    un.cmplx(mx,ny,mz,0) = 0.0;
	    un.cmplx(mx,ny,mz,1) = 0.0;
	    un.cmplx(mx,ny,mz,2) = 0.0;
	    qn.cmplx(mx,ny,mz,0) = 0.0;
	  }
	}
      }
    }
    t_ += dt_;

    if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll) 
      cout << t_ << ' ' << flush;
    else if (flags_.verbosity == PrintTicks) 
      cout << '.' << flush;
  }
    
  cfl_ = un.CFLfactor(Ubase_);
  cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
  
  // If using dealiasing, set flag in FlowField that compactifies binary IO 
  un.setDealiased(flags_.dealias_xz()); 
  qn.setDealiased(flags_.dealias_xz()); 

  if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll) 
    cout << endl;
    
  return;
}


// ==============================================================
// CNAB-style algorithms

CNABstyleDNS::CNABstyleDNS()
  :
  DNSAlgorithm(),
  full_(false),
  fj1_(),
  fj_()
{}

CNABstyleDNS::CNABstyleDNS(const CNABstyleDNS& dns) 
  :
  DNSAlgorithm(dns),
  Nsubsteps_(dns.Nsubsteps_),
  full_(dns.full_),
  fj1_(dns.fj1_),
  fj_(dns.fj_),
  alpha_(dns.alpha_),
  beta_(dns.beta_),
  gamma_(dns.gamma_),
  zeta_(dns.zeta_)
{
  // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver arrays
  // and copy tausolvers from dns argument
  tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
  for (int j=0; j<Nsubsteps_; ++j) {
    tausolver_[j] = new TauSolver*[Mx_];       // new #2
    for (int mx=0; mx<Mx_; ++mx) {
      tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
      for (int mz=0; mz<Mz_; ++mz)
	tausolver_[j][mx][mz] = dns.tausolver_[j][mx][mz];
	
    }
  }
}

CNABstyleDNS::CNABstyleDNS(FlowField& u, const ChebyCoeff& Ubase,
			   Real nu, Real dt, const DNSFlags& flags, Real t) 
  :
  DNSAlgorithm(u,Ubase,nu,dt,flags, t),
  full_(false),
  fj1_(u),
  fj_(u)
{
  fj1_.setToZero();
  fj_.setToZero();


  TimeStepMethod algorithm = flags.timestepping;
  switch (algorithm) {
  case CNAB2: 
    order_ = 2;   
    Nsubsteps_ = 1;
    Ninitsteps_ = 1;
    full_ = false;
    alpha_.resize(Nsubsteps_);
    beta_.resize(Nsubsteps_);
    gamma_.resize(Nsubsteps_);
    zeta_.resize(Nsubsteps_);
    alpha_[0] = 0.5;
    beta_[0]  = 0.5;
    gamma_[0] = 1.5;
    zeta_[0]  = -0.5;
    
    break;
  case SMRK2:
    order_ = 2;
    Nsubsteps_ = 3;
    Ninitsteps_ = 0;
    full_ = true;
    alpha_.resize(Nsubsteps_);
    beta_.resize(Nsubsteps_);
    gamma_.resize(Nsubsteps_);
    zeta_.resize(Nsubsteps_);
    alpha_[0] = 29.0/96.0;  alpha_[1] = -3.0/40.0;  alpha_[2] = 1.0/6.0;
    beta_[0]  = 37.0/160.0;  beta_[1] =  5.0/24.0;  beta_[2]  = 1.0/6.0;
    gamma_[0] = 8.0/15.0;   gamma_[1] =  5.0/12.0;  gamma_[2] = 3.0/4.0;
    zeta_[0]  = 0.0;         zeta_[1] = -17.0/60.0;  zeta_[2] = -5.0/12.0;  
    break;
  default: 
    cerr << "CNABstyleDNS::CNABstyleDNS(un,Ubase,nu,dt,flags,t0)\n"
	 << "error: flags.timestepping == " << algorithm 
	 << " is not a CNAB-style algorithm." << endl;
    exit(1);
  }

  // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
  tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
  for (int j=0; j<Nsubsteps_; ++j) {
    tausolver_[j] = new TauSolver*[Mx_];       // new #2
    for (int mx=0; mx<Mx_; ++mx) 
      tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
  }
  reset_dt(dt_);
}

CNABstyleDNS::~CNABstyleDNS() {
  if (tausolver_) {
    for (int j=0; j<Nsubsteps_; ++j) {
      for (int mx=0; mx<Mx_; ++mx) 
	delete[] tausolver_[j][mx];  // undo new #3
      delete[] tausolver_[j];        // undo new #2
    } 
    delete[] tausolver_;             // undo new #1
  }
  tausolver_ = 0;
           
           
}

DNSAlgorithm* CNABstyleDNS::clone() const {
  return new CNABstyleDNS(*this);
}


void CNABstyleDNS::reset_dt(Real dt) {
  cfl_ *= dt/dt_;
  //nu_ = nu;
  dt_ = dt;
  const Real c = 4.0*square(pi)*nu_;
  const int kxmax = tmp_.kxmax();
  const int kzmax = tmp_.kzmax();
 

  // This loop replaces the TauSolver objects at tausolver_[i][mx][mz] 
  // with new TauSolver objects configured with appropriate parameters
  for (int j=0; j<Nsubsteps_; ++j) {
    for (int mx=0; mx<Mx_; ++mx) {
      int kx = tmp_.kx(mx);
      for (int mz=0; mz<Mz_; ++mz) {
	int kz = tmp_.kz(mz);

	Real lambda = 1.0/(beta_[j]*dt_)+c*(square(kx/Lx_)+square(kz/Lz_));

	if ((kx != kxmax) && (kz != kzmax) && 
	    (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))

	  tausolver_[j][mx][mz] = 
	    TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, nu_, Nyd_, 
		      flags_.taucorrection);
      }
    }
  }
  // For some forms of CNABstyle, need to reinitialize again
  switch (flags_.timestepping) {
  case CNAB2:
    full_ = false;
    break;
  default:
    ;
  }
}

bool CNABstyleDNS::push(const FlowField& un) {
  swap(fj_, fj1_);
   fj_.setToZero();
    fj1_.setToZero();
  navierstokesNL(un, Ubase_, fj_, tmp_, flags_.nonlinearity);
  t_ += dt_;
  full_ = true;
  return full_;
}

bool CNABstyleDNS::full() const {
  return full_;
}
// Incorporating Reverse Force by volume inverse rule method//



void CNABstyleDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {
  const int kxmax = tmp_.kxmax();
  const int kzmax = tmp_.kzmax();
 
  for (int n=0; n<Nsteps; ++n) 	{
    for (int j=0; j<Nsubsteps_; ++j) {

      // Store substeps uj,qj in un,qn; reflect this in notation
      FlowField& uj(un); 
      FlowField& qj(qn);
        
      swap(fj_, fj1_);
    
      navierstokesNL(un, Ubase_, fj_, tmp_, flags_.nonlinearity);
      
FlowField force(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                un.Lz(), un.a(), un.b(), Physical, Physical);

FlowField utmp(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                un.Lz(), un.a(), un.b(), Physical, Physical);

FlowField Dxx1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                     un.Lz(), un.a(), un.b(), Physical, Physical);
        
FlowField Dyy1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                     un.Lz(), un.a(), un.b(), Physical, Physical);            

FlowField Dzz1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                     un.Lz(), un.a(), un.b(), Physical, Physical);        

FlowField dx1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                     un.Lz(), un.a(), un.b(), Physical, Physical);        

FlowField dy1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                     un.Lz(), un.a(), un.b(), Physical, Physical);       

FlowField dz1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                    un.Lz(), un.a(), un.b(), Physical, Physical);       

FlowField sx1(un.Nx(), un.Ny(), un.Nz(), 3, un.Lx(),
                    un.Lz(), un.a(), un.b(), Physical, Physical);


//*************************************************************************//
//********End for Func calling for Stress moment calc by interpolation*****//
//*************************************************************************//

//*************************************************************************//
//*********Func calling for Stress moment calc by spectral method**********//
//*************************************************************************//
int i=0;
  Vector x = periodicpoints(un.Nx(), un.Lx());
  Vector z = periodicpoints(un.Nz(), un.Lz());
 
  reverseforce(i,force);

force.makePhysical();
               double totalforce_x = 0.0;
             for(int ny=0;ny<force.Ny();++ny)
               for(int nx=0;nx<force.Nx();++nx)
                  for(int nz=0;nz<force.Nz();++nz)
                         totalforce_x += force(nx,ny,nz,0);

               cout<<totalforce_x<<endl; 
                         
                  force.makeSpectral();

 third_order_stress_moments(i,force,Dxx1,Dyy1,Dzz1);
 rotational_stress_tensor(i,force,dx1,dy1,dz1);
 symmetric_stress_tensor(i,un,utmp,force,sx1,x,z);
//*************************************************************************//
//*****End for Func calling for Stress moment calc by spectral method*****//
//*************************************************************************//

                         // Set convenience variables 
      Real a_b = alpha_[j]/beta_[j]; 
      Real g_b = gamma_[j]/beta_[j]; 
      Real z_b = zeta_[j]/beta_[j];  
      Real anu_b = a_b*nu_;           
      Real anu = alpha_[j]*nu_;         
      
      // Update each Fourier mode with time-stepping algorithm
      for (int mx=0; mx<Mx_; ++mx) {
	const int kx = uj.kx(mx);

	for (int mz=0; mz<Mz_; ++mz) {
	  const int kz = uj.kz(mz);

	  // Zero out the aliased modes
	  if ((kx == kxmax || kz == kzmax) || 
	    flags_.dealias_xz() && isAliasedMode(kx,kz)) {
	    for (int ny=0; ny<Nyd_; ++ny) {
	      uj.cmplx(mx,ny,mz,0) = 0.0;
	      uj.cmplx(mx,ny,mz,1) = 0.0;
	      uj.cmplx(mx,ny,mz,2) = 0.0;
	      qj.cmplx(mx,ny,mz,0) = 0.0;
	    }
	    break;
	  }

	  Rxk_.setToZero();
	  Ryk_.setToZero();
	  Rzk_.setToZero();
	  
	  // Goal is to compute
	  // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
	  //     - g/b fj - z/b fj1
	  //
	  //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
	  //     - a/b grad qj - g/b fj - z/b fj1

	  // Extract relevant Fourier modes of uj and qj for computations
	  // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj

	  for (int ny=0; ny<Nyd_; ++ny) {
	    uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
	    vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
	    wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));             
	    Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
           
	  }

	  // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)

	  diff2(uk_, Rxk_, Pyk_);
	//  diff2(vk_, Ryk_, Pyk_);
	//  diff2(wk_, Rzk_, Pyk_);

	  // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
	  diff(Pk_, Pyk_); 
      
	  // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj 
	  //          - g/b fj - z/b fj1 to R, completing calculation of R
	  const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
	  const Real c = 1.0/(anu*dt_) - kappa2;
	  const Complex Dx = un.Dx(mx);
	  const Complex Dz = un.Dz(mz);
          double fluid_scale = 1/beta_[j];
	  for (int ny=0; ny<Nyd_; ++ny) {
	    Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]                     
                     + fluid_scale*(force.cmplx(mx,ny,mz,0)
                     + sx1.cmplx(mx,ny,mz,0)
                     + (dx1.cmplx(mx,ny,mz,1) + dx1.cmplx(mx,ny,mz,2))
                     + (Dxx1.cmplx(mx,ny,mz,0) + Dyy1.cmplx(mx,ny,mz,0) + Dzz1.cmplx(mx,ny,mz,0)))
		     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0)); 
                  
	   Ryk_.add(ny, c*vk_[ny] - Pyk_[ny] 
                    + fluid_scale*(force.cmplx(mx,ny,mz,1)
                    + sx1.cmplx(mx,ny,mz,1)
                   + (dy1.cmplx(mx,ny,mz,0) + dy1.cmplx(mx,ny,mz,2))
                   + (Dxx1.cmplx(mx,ny,mz,1) + Dyy1.cmplx(mx,ny,mz,1) + Dzz1.cmplx(mx,ny,mz,1)))
		     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));

	    Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny] 
                     + fluid_scale*(force.cmplx(mx,ny,mz,2)
                     + sx1.cmplx(mx,ny,mz,2)
                     + (dz1.cmplx(mx,ny,mz,0) + dz1.cmplx(mx,ny,mz,1))
                     + (Dxx1.cmplx(mx,ny,mz,2) + Dyy1.cmplx(mx,ny,mz,2) + Dzz1.cmplx(mx,ny,mz,2)))
                    - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
	 
          }  

  
		  // Solve the tau solutions
	  if (kx!=0 || kz!=0) 
	    tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
             //  tausolver_[j][mx][mz].solve(uk_,0,0,0, Rxk_,0,0);
	  else { // kx,kz == 0,0
	    // Rx has additional terms, nu Uyy at both t=j and t=j+1
	    Real a_b = alpha_[j]/beta_[j];
	    Real c = nu_*(a_b+1.0);
	    if (Ubaseyy_.length() > 0)
	      for (int ny=0; ny<Ny_; ++ny){
		Rxk_.re[ny] += c*Ubaseyy_[ny];
                
	  }
	    if (flags_.constraint == PressureGradient) {
	      // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
	      Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;

	      // Solve the tau equations
	      tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
              

	      // Bulk velocity is free variable on LHS solved by tau eqn
	      UbulkAct_ = UbulkBase_ + uk_.re.mean();                 
	      dPdxAct_ = dPdxRef_;
            
	    }
	    else { // const bulk velocity
	      // Add the previous time-step's -dPdx to the RHS. The next 
	      // timestep's dPdx term appears on LHS as unknown.
	      //Rxk_.re[0] -= a_b*dPdxAct_;
              
	      // Use tausolver with additional variable and constraint:
	      // free variable: dPdxAct at next time-step, 
	      // constraint:    UbulkBase + mean(u) = UbulkRef.
	      tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, 
					  Rxk_, Ryk_, Rzk_, 
					  UbulkRef_ - UbulkBase_);

	                                 UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
   
	      //UbulkAct_ = UbulkRef_;
	    }
	    // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
	    // Pk_.set(0, Complex(0.0, 0.0));
	  }

	  // Load solutions back into the external 3d data arrays.
	  // Because of FFTW complex symmetries
	  // The 0,0 mode must be real.
	  // For Nx even, the kxmax,0 mode must be real
	  // For Nz even, the 0,kzmax mode must be real
	  // For Nx,Nz even, the kxmax,kzmax mode must be real
	  if ((kx == 0 && kz == 0) ||
	      (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
	      (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
	      (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {

	    for (int ny=0; ny<Nyd_; ++ny) {
	      un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
	      un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
	      un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
	      qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
           
	    }
	  }
           
	  // The normal case, for general kx,kz
	  else
	    for (int ny=0; ny<Nyd_; ++ny) {
	      un.cmplx(mx,ny,mz,0) = uk_[ny];
	      un.cmplx(mx,ny,mz,1) = vk_[ny];
	      un.cmplx(mx,ny,mz,2) = wk_[ny];
	      qn.cmplx(mx,ny,mz,0) = Pk_[ny];
	    }
            
	  // And now set the y-aliased modes to zero.
	  for (int ny=Nyd_; ny<Ny_; ++ny) {
	    un.cmplx(mx,ny,mz,0) = 0.0;
	    un.cmplx(mx,ny,mz,1) = 0.0;
	    un.cmplx(mx,ny,mz,2) = 0.0;
	    qn.cmplx(mx,ny,mz,0) = 0.0;
	}
	}
      }

    }
    t_ += dt_;
       if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll) 
      cout << t_ << ' ' << flush;
    else if (flags_.verbosity == PrintTicks) 
      cout << '.' << flush;
  }


  cfl_ = un.CFLfactor(Ubase_);
  cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
 
  // If using dealiasing, set flag in FlowField that compactifies binary IO 
  un.setDealiased(flags_.dealias_xz()); 
  qn.setDealiased(flags_.dealias_xz()); 

  if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll) 
    cout << endl;
    
  return;

    
}



void changeBaseFlow(const ChebyCoeff& ubase0, const FlowField& ufluc0, 
		const FlowField& q0arg, 
		const ChebyCoeff& ubase1, FlowField& u1, FlowField& q1){
  ChebyCoeff& U0 = (ChebyCoeff&) ubase0;
  fieldstate U0state = U0.state();

  ChebyCoeff& U1 = (ChebyCoeff&) ubase1;
  fieldstate U1state = U1.state();

  FlowField& u0 = (FlowField&) ufluc0;
  fieldstate u0xzstate = u0.xzstate();
  fieldstate u0ystate = u0.ystate();
  
  FlowField& q0 = (FlowField&) q0arg;
  fieldstate q0xzstate = q0.xzstate();
  fieldstate q0ystate = q0.ystate();
  
  int Nx=u0.numXgridpts();
  int Ny=u0.numYgridpts();
  int Nz=u0.numZgridpts();

  u1 = u0; // want u1 FPF
  u1.makeState(Spectral, Physical);
  u0.makePhysical();
  q0.makePhysical();
  q1 = q0; // want q1 physical

  // At this point 
  // u1 == utot - U0
  // q1 == p + 1/2 u0 dot u0
  
  // Remove 1/2 u0 dot u0 from q1
  for (int ny=0; ny<Ny; ++ny)
    for (int nx=0; nx<Nx; ++nx)
      for (int nz=0; nz<Nz; ++nz) 
	q1(nx,ny,nz,0) -= 0.5*(square(u0(nx,ny,nz,0)) + 
			       square(u0(nx,ny,nz,1)) + 
			       square(u0(nx,ny,nz,2)));
  // At this point 
  // u1 == utot - U0
  // q1 == p 
  
  ChebyTransform t(U0.numModes());
  U0.makePhysical(t);
  U1.makePhysical(t);

  // Add U0-U1 to u1 
  ChebyCoeff delta_U(U0);
  delta_U -= U1;
  u1 += delta_U;
  u1.makePhysical();

  // At this point 
  // u1 == utot - U1
  // q1 == p 
  
  // Add 1/2 u1 dot u1 to q1
  for (int ny=0; ny<Ny; ++ny)
    for (int nx=0; nx<Nx; ++nx)
      for (int nz=0; nz<Nz; ++nz) 
	q1(nx,ny,nz,0) += 0.5*(square(u1(nx,ny,nz,0)) + 
			      square(u1(nx,ny,nz,1)) + 
			      square(u1(nx,ny,nz,2)));
  // At this point 
  // u1 == utot - U1
  // q1 == p + 1/2 u1 dot u1
  // et, voila

  U0.makeState(U0state,t);
  U1.makeState(U1state,t);
  u0.makeState(u0xzstate, u0ystate); 
  q0.makeState(q0xzstate, q0ystate); 
  u1.makeState(u0xzstate, u0ystate); 
  q1.makeState(q0xzstate, q0ystate); 
}
  
DNSFlags::DNSFlags(MeanConstraint  constraint_,
		   TimeStepMethod  timestepping_,
		   TimeStepMethod  initstepping_,
		   NonlinearMethod nonlinearity_,
		   Dealiasing      dealiasing_,
		   bool            taucorrection_,
		   Verbosity       verbosity_,
		   Real            dPdx_,
		   Real            Ubulk_)
  :
  constraint(constraint_),
  timestepping(timestepping_),
  initstepping(initstepping_),
  nonlinearity(nonlinearity_),
  dealiasing(dealiasing_),
  taucorrection(taucorrection_),
  verbosity(verbosity_),
  dPdx(dPdx_),
  Ubulk(Ubulk_)
{
  if (dealias_y() && (nonlinearity != Rotational)) {
    cerr << "DNSFlags::DNSFlags: DealiasY and DealiasXYZ work only with\n";
    cerr << "Rotational nonlinearity in the current version of channelflow.\n";
    cerr << "Setting nonlinearity to Rotational." << endl;
    nonlinearity = Rotational;
  }
  // maybe should print warnings about initstepping only mattering for SBDF and CNAB
}

bool DNSFlags::dealias_xz() const {
  return ((dealiasing == DealiasXZ || dealiasing == DealiasXYZ) ? true:false);
}

bool DNSFlags::dealias_y() const {
  return ((dealiasing == DealiasY || dealiasing == DealiasXYZ) ? true:false);
}

ostream& operator<<(ostream& os, Dealiasing d) {
  string s;
  switch(d) {
  case NoDealiasing: s="NoDealiasing"; break;
  case DealiasXZ: s="DealiasXZ"; break;
  case DealiasY: s="DealiasY"; break;
  case DealiasXYZ: s="DealiasXYZ"; break;
  default: s="Invalid Dealiasing value: please submit bug report";
  }
  os << s;
  return os;
}

ostream& operator<<(ostream& os, MeanConstraint m) {
  string s;
  switch(m) {
  case PressureGradient: s="PressureGradient"; break;
  case BulkVelocity: s="BulkVelocity"; break;
  default: s="Invalid MeanConstraint value: please submit bug report";
  }
  os << s;
  return os;
}

ostream& operator<<(ostream& os, TimeStepMethod t) {
  string s;
  switch(t) {
  case CNFE1: s="CNFE1"; break;
  case CNAB2: s="CNAB2"; break;
  case CNRK2: s="CNRK2"; break;
  case SMRK2: s="SMRK2"; break;
  case SBDF1: s="SBDF1"; break;
  case SBDF2: s="SBDF2"; break;
  case SBDF3: s="SBDF3"; break;
  case SBDF4: s="SBDF4"; break;
  default: s="Invalid TimeStepMethod value: please submit bug report";
  }
  os << s;
  return os;
}
ostream& operator<<(ostream& os, Verbosity v) {
  string s;
  switch(v) {
  case Silent: s="Silent"; break;
  case PrintTicks: s="PrintTicks"; break;
  case PrintTime: s="PrintTime"; break;
  case VerifyTauSolve: s="VerifyTauSolve"; break;
  case PrintAll: s="PrintAll"; break;
  default: s="Invalid Verbosity value: please submit bug report";
  }
  os << s;
  return os;
}

ostream& operator<<(ostream& os, DNSFlags& flags) {
  string s(", ");
  string tau = (flags.taucorrection) ? "TauCorrection" : "NoTauCorrection";
  os << flags.timestepping << s 
     << flags.initstepping << s 
     << flags.nonlinearity << s
     << flags.dealiasing << s 
     << flags.constraint << s 
     << tau << s 
     << flags.verbosity;
  return os;
}

} //namespace channelflow

