/* nsintegrator.h: time-integration class for spectral Navier-Stokes DNS
 
 */

// NSIntegrator is a class for integrating Navier-Stokes equation.
// DNSFlags is used to specify the integration parameters of NSIntegrator.
// TimeStep manages variable time-stepping, adjusting dt to keep CFL in range
// TimeStepParams holds the numerical parameters of timestepping algorithms


#ifndef CHANNELFLOW_DNS_H
#define CHANNELFLOW_DNS_H

#include "channelflow/mathdefs.h"
#include "channelflow/vector.h"
#include "channelflow/chebyshev.h"
#include "channelflow/array.h"
#include "channelflow/flowfield.h"
#include "channelflow/diffops.h"
#include "channelflow/tausolver.h"


namespace channelflow {

// Enum types for specifying the behavior of NSIntegrator, fields of DNSFlags.
enum Dealiasing      {NoDealiasing, DealiasXZ, DealiasY, DealiasXYZ};
enum TimeStepMethod  {CNFE1, CNAB2, CNRK2, SMRK2,
		      SBDF1, SBDF2, SBDF3, SBDF4}; 
enum MeanConstraint  {PressureGradient, BulkVelocity};
enum Verbosity       {Silent, PrintTicks, PrintTime, VerifyTauSolve, PrintAll};

// CNFE1  == Crank-Nicolson Forward-Euler order 1   (no init steps needed)
// CNAB2  == Crank-Nicolson Adams-Bashforth order 2 (needs 1 init steps)
// SMRK2  == Spalart,Moser,R?, Runge-Kutta order 2  (no init steps needed)
// SBDFn  == Semimplicit backwards-differentiation order n (needs n-1 init steps)
// See channelflow manual for detailed description of these algorithms
// Note: CNFE1 and SBDF1 are the same algorithm.

/*********************************************************************
// Now declared in diffops.h
enum NonlinearMethod {Rotational, SkewSymmetric, Convection, Divergence, 
                      Alternating, Linearized}; 
**********************************************************************/

// Specify the behavior of NSIntegrators by setting fields of DNSFlags.
class DNSFlags {
public: 
  //       Option type     Option name      Default value
  DNSFlags(MeanConstraint  constraint     = PressureGradient,
	   TimeStepMethod  timestepping   = SBDF3,
	   TimeStepMethod  initstepping   = SMRK2,
	   NonlinearMethod nonlinearity   = Rotational, 
	   Dealiasing      dealiasing     = NoDealiasing,
	   bool            taucorrection  = true,
	   Verbosity       verbosity      = PrintTicks,
	   Real            dPdx           = 0.0,
	   Real            Ubulk          = 0.0);
       //    bool            reverseforce   = true);

  MeanConstraint  constraint;   // Enforce const press grad or const bulk vel
  TimeStepMethod  timestepping; // Time-stepping algorithm
  TimeStepMethod  initstepping; // Algorithm for initializing multistep methods
  NonlinearMethod nonlinearity; // Method of calculating nonlinearity of NS eqn
  Dealiasing      dealiasing;   // Use 3/2 rule to eliminate aliasing
  bool            taucorrection;// Remove divergence caused by discretization
  Verbosity       verbosity;    // Print diagnostics, times, ticks, or nothing
  Real            dPdx;         // Constraint value
  Real            Ubulk;        // Constraint value
 // bool            reverseforce;

  bool dealias_xz() const;
  bool dealias_y() const;
};


// TimeStep keeps dt between dtmin and dtmax, and CFL between CFLminand CFLmax,
// in such a way that dt*n = dT for some integer n. That's useful if you
// want to plot/save data at regular dT intervals, but use a variable timestep 
// dt for efficiency. For example of use, see example codes. 

class TimeStep {
public:
  TimeStep(Real dt, Real dtmin, Real dtmax, Real dT, Real CFLmin, Real CFLmax);

  // adjust dt to keep CFLmin<=CFL<=CFLmax (soft), dtmin<=dt<=dtmax (hard)
  // returns true if dt changes, false if dt is unchanged.
  bool adjust(Real CFL, bool verbose=true); 

  Real CFL() const;
  int n() const;         // n*dt == dT
  Real dt() const;       // integration timestep
  Real dT() const;       // plot interval
  operator Real() const; // same as dt()

private:
  int n_;     // n_    == int(dT/dt)freverse_x2, number of dt steps per plot interval
  int nmin_;  // nmin_ == int(dT/dtmax)
  int nmax_;  // nmax_ == int(dT/dtmin)
  Real dT_;   // dT_ == n_*dt_, plot interval
  Real CFLmin_;
  Real CFL_;
  Real CFLmax_;
};


class DNSAlgorithm;

// DNS is a wrapper class for DNSAlgorithms. It's the main class for 
// integrating the Navier-Stokes equations in top-level programs. 
// Specify the integration algorithm and other parameters in the DNSFlags.
// If you like, you can construct and use specific DNS algorithms like 
// MultiStepDNS in top-level programs --any class derived from DNSAlgorithm. 
// Look in example codes for examples of initialization and use.

class DNS {
public:
  DNS();
  DNS(const DNS& dns);
  DNS(FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt, 
      const DNSFlags& flags, Real t=0.0);
  DNS(FlowField& u, FlowField& force, Real nu, Real dt,
      const DNSFlags& flags, Real t=0.0);
  ~DNS();

  DNS& operator=(const DNS& dns);
  
  void advance(FlowField& u, FlowField& q, int nSteps=1); 

  //void reset();                  // flush state, prepare for new integration
  void reset_dt(Real dt);
  void reset_time(Real t);
  void reset_dPdx(Real dPdx);    // change dPdx and enforce const dPdx
  void reset_Ubulk(Real Ubulk);  // change Ubulk and enforce const Ubulk

  //void reset_uj(const FlowField& uj, int j);  // set u[j]=u(t-j*dt)
  bool push(const FlowField& u); // push into u[j] stack, true when full, t+=dt
  bool full() const;             // pushed enough init data into u[j]?

  int order() const;             // err should scale as dt^order
  int Ninitsteps() const;        // number of steps needed to initialize

  Real dt() const;
  Real CFL() const;
  Real time() const;
  Real dPdx() const;      // the mean pressure gradient at the current time
  Real Ubulk() const;     // the actual bulk velocity at the current time

  Real dPdxRef() const;   // the mean press grad enforced during integration
  Real UbulkRef() const;  // the bulk velocity enforced during integ.

  const DNSFlags& flags() const;
  TimeStepMethod timestepping() const;

  //void printStack() const;

 private:
  DNSAlgorithm* main_algorithm_;        // same as 
  DNSAlgorithm* init_algorithm_;

  DNSAlgorithm* newAlgorithm(FlowField& u, const ChebyCoeff& Ubase, Real nu, 
			     Real dt, const DNSFlags& flags, Real t);
  //DNSAlgorithm* newAlgorithm(FlowField& u, FlowField& force, Real nu,
    //                        Real dt, const DNSFlags& flags, Real t);
};


// DNSAlgorithm is a base class for classes representing time-stepping 
// algorithms for the Navier-Stokes equations, using a Fourier x Chebyshev 
// x Fourier FlowField for spatial discretization and finite-differencing 
// and tau method for temporal discretization. 

class DNSAlgorithm {
public:
  DNSAlgorithm();
  DNSAlgorithm(FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt, 
	       const DNSFlags& flags, Real t=0);
  virtual ~DNSAlgorithm();

  DNSAlgorithm& operator=(const DNSAlgorithm& dns);
  
  virtual void advance(FlowField& u, FlowField& q, int nSteps=1) = 0; 

  //virtual void reset();             // flush state, prepare for new integration
  virtual void reset_dt(Real dt) = 0;     // somewhat expensive
  virtual bool push(const FlowField& u);  // push u onto u[j] stack, t += dt
  virtual bool full() const;              // have enough init data?
  
  void reset_time(Real t);
  void reset_dPdx(Real dPdx);    // change dPdx and enforce const dPdx
  void reset_Ubulk(Real Ubulk);  // change Ubulk and enforce const Ubulk

  int order() const;             // err should scale as dt^order
  int Ninitsteps() const;        // number of steps needed to initialize

  int Nx() const;
  int Ny() const;
  int Nz() const;
  int Mx() const;
  int My() const;
  int Mz() const;
  Real Lx() const;
  Real Lz() const;
  Real a() const;
  Real b() const;
  Real nu() const;
  Real dt() const;
  Real CFL() const;
  Real time() const;
  Real dPdx() const;      // the mean pressure gradient at the current time
  Real Ubulk() const;     // the actual bulk velocity at the current time
  Real dPdxRef() const;   // the mean press grad enforced during integration
  Real UbulkRef() const;  // the bulk velocity enforced during integ.

  const DNSFlags& flags() const;
  const ChebyCoeff& Ubase() const;
  TimeStepMethod timestepping() const;

  virtual DNSAlgorithm* clone() const = 0;  // new copy of *this

protected:
  // Spatial parameters 
  int Nx_;      // number of X gridpoints
  int Ny_;      // number of Chebyshev T(y) modes
  int Nz_;      // number of Z gridpoints
  int Mx_;      // number of X modes
  int Mz_;      // number of Z modes
  int Nyd_;     // number of dealiased Chebyshev T(y) modes 
  int kxd_max_; // maximum value of kx among dealiased modes
  int kzd_max_; // maximum value of kz among dealiased modes
  Real Lx_;
  Real Lz_;
  Real a_;
  Real b_;

  // Temporal integration parameters 
  DNSFlags flags_; // User-defined integration parameters
  int order_;
  int Ninitsteps_; // number of initialization steps required
  Real nu_;
  Real dt_;
  Real t_;         // time in convective units
  Real cfl_;       // CFL number
  Real dPdxRef_;   // Enforced mean pressure gradient (0.0 if unused).
  Real dPdxAct_;   // Actual mean pressure gradient at previous timestep.
  Real UbulkRef_;  // Enforced total bulk velocity (0.0 if unused).
  Real UbulkAct_;  // Actual total bulk velocity bulk obtained.
  Real UbulkBase_; // Bulk velocity of base flow

  ChebyCoeff Ubase_;   // baseflow physical
  ChebyCoeff Ubaseyy_; // baseflow'' physical

  FlowField tmp_;   // used in calculation of nonlinearity
  FlowField ubase_; // used when linearizing about genl field.

  // These variables are used as temp storage when solving indpt tau problems.
  ComplexChebyCoeff uk_;   // profile of u_{kx,kz} (y) at t = n dt
  ComplexChebyCoeff vk_;
  ComplexChebyCoeff wk_;
  ComplexChebyCoeff Pk_;   // profile of P_{kx,kz} (y)
  ComplexChebyCoeff Pyk_;  // profile of dP_{kx,kz}/dy (y)
  ComplexChebyCoeff Rxk_;
  ComplexChebyCoeff Ryk_;
  ComplexChebyCoeff Rzk_;
  
  int kxmaxDealiased() const;
  int kzmaxDealiased() const;
  bool isAliasedMode(int kx, int kz) const;
};



// Multistep algorithms, all using Backwards Differentiation: SBDFk
// Based on Peyret section. The order is set by flags.timestepping.
class MultistepDNS : public DNSAlgorithm {
public:
  MultistepDNS();
  MultistepDNS(const MultistepDNS& dns);
  MultistepDNS(FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt, 
	       const DNSFlags& flags, Real t=0);
  ~MultistepDNS();
  
  MultistepDNS& operator=(const MultistepDNS& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1); 

  //virtual void reset();             // flush state, prepare for new integration
  virtual void reset_dt(Real dt); 
  virtual bool push(const FlowField& u); // for initialization
  virtual bool full() const;             // have enough init data?

  virtual DNSAlgorithm* clone() const;  // new copy of *this
  //virtual void setreverseforce (FlowField& force);

protected:
  Real eta_;
  array<Real> alpha_;
  array<Real> beta_;
  array<FlowField> u_;  // u[j] == u at t-j*dt for multistep algorithms
  array<FlowField> f_;  // f[j] == f at t-j*dt for multistep algorithms
  
  TauSolver** tausolver_;  // 2d array of tausolvers, indexed by [mx][mz]

  int countdown_;
};

// CNRK2 and hopefully another. Based on algorithm in Peyret pg 149
class RungeKuttaDNS : public DNSAlgorithm {
public:
  RungeKuttaDNS();
  RungeKuttaDNS(const RungeKuttaDNS& dns);
  RungeKuttaDNS(FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt, 
		const DNSFlags& flags, Real t=0);
  ~RungeKuttaDNS();

  RungeKuttaDNS& operator=(const RungeKuttaDNS& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1); 

  //virtual void reset();             // flush state, prepare for new integration
  virtual void reset_dt(Real dt); 

  virtual DNSAlgorithm* clone() const;  // new copy of *this
  //void setreverseforce (FlowField& force);

protected:
  int Nsubsteps_;
  FlowField Qj1_;  // Q_{j-1} (Q at previous substep)
  FlowField Qj_;   // Q_j     (Q at current  substep)
  array<Real> A_;  // Q_{j+1} = A_j Q_j + N(u_j)
  array<Real> B_;  // u_{j+1} = u_j + dt B_j Q_j + dt C_j (L u_j + L u_{j+1})
  array<Real> C_;
  array<Real> D_;  

  TauSolver*** tausolver_; // 3d array indexed by [step][mx][mz]
};

// A generalization of CNAB2 with substeps. Implements CNAB2 and SMRK2
class CNABstyleDNS : public DNSAlgorithm{
public:
  CNABstyleDNS();
  CNABstyleDNS(const CNABstyleDNS& dns);
  CNABstyleDNS(FlowField& u, const ChebyCoeff& Ubase, Real nu, 
	       Real dt, const DNSFlags& flags, Real t=0);
  ~CNABstyleDNS();

  CNABstyleDNS& operator=(const CNABstyleDNS& dns);

  virtual void advance(FlowField& conc, FlowField& q, int nSteps=1); 
  virtual void reset_dt(Real dt); 
  virtual bool push(const FlowField& u);  // push u onto u[j] stack, t += dt
  virtual bool full() const;              // have enough init data?
  void setreverseforce (FlowField& force);
  void set_third_order_stress_moment_x(FlowField& Dxx, FlowField& Dxx1);
  void set_third_order_stress_moment_y(FlowField& Dyy, FlowField& Dyy1);
  void set_third_order_stress_moment_z(FlowField& Dzz, FlowField& Dzz1);
  void set_rotational_stress_tensor_x(FlowField& dx, FlowField& dx1);
  void set_rotational_stress_tensor_y(FlowField& dy, FlowField& dy1);
  void set_rotational_stress_tensor_z(FlowField& dz, FlowField& dz1);
  void set_symmetric_tensor(FlowField& u, FlowField& sx, FlowField& sx1);

private:
  int Nsubsteps_;
  bool full_;
  FlowField fj1_;      // f_{j-1} (f at previous substep)
  FlowField fj_;       // f_j     (f at current  substep)
  array<Real> alpha_;  // u_{j+1} = u_j + dt L (alpha_j u_j + beta_j u_{j+1})
  array<Real> beta_;    //           + dt gamma_j N(u_j) + dt zeta N(u_{j-1})
  array<Real> gamma_;
  array<Real> zeta_;
 

  TauSolver*** tausolver_; // 3d array indexed by [i][mx][mz]

  virtual DNSAlgorithm* clone() const;    // new copy of *this
};


 

std::ostream& operator<<(std::ostream& os, Dealiasing d);
std::ostream& operator<<(std::ostream& os, MeanConstraint m);
std::ostream& operator<<(std::ostream& os, TimeStepMethod t);
std::ostream& operator<<(std::ostream& os, Verbosity v);
std::ostream& operator<<(std::ostream& os, DNSFlags& flags);


// Given a baseflow, fluctation, modified pressure triple (U0,u0,q0) and
// a new baseflow U1, compute new fluctuation u1 and modified pressure q1.
//     utot == U0 + u0 == U1 + u1
// pressure == q0 - 1/2 || u0 ||^2 == q1 - 1/2 || u1 ||^2 
void changeBaseFlow(const ChebyCoeff& U0, const FlowField& u0, const FlowField& q0, 
		    const ChebyCoeff& U1, FlowField& u1, FlowField& q1);

} //namespace channelflow
#endif
