#ifndef _DHSHARED_H_
#define _DHSHARED_H_

#include <sstream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "particle.h"
#include "dh.h"

#if 1   /* set to 1 in order to use Symplectic Corrector */
#define USE_SYMPC
#endif

struct Orbit
{
  real inc;
  real sma;
  real ecc;
  real Omega;
  real lambda;

  Orbit(const vec3 &R, const vec3 &V, const real Mtot)
  {
    const vec3 L = R%V;
    const real h2 = L.norm2();
    const real h  = std::sqrt(h2); 

    const real inc = std::acos(L.z/h);
    const real fac = std::sqrt(L.x*L.x + L.y*L.y)/h;
    const real LSMALL = 1.0e-10;

    real capom, u;
    if (fac < LSMALL)
    {
      capom = 0.0;
      u     = std::atan2(R.y, R.x);
      if (std::abs(inc - M_PI) < 10*LSMALL) 
        u = -u;
    } 
    else
    {
      capom = std::atan2(L.x, -L.y);
      u = std::atan2( R.z/std::sin(inc) , R.x*std::cos(capom) + R.y*std::sin(capom));
    }

    if (capom < 0.0) capom += 2.0*M_PI;
    if (u     < 0.0) u     += 2.0*M_PI;

    const real r  = R.abs();
    const real v2 = V.norm2();
    const real vdotr = R*V;
    const real  energy = 0.5*v2 -  Mtot/r; 
    assert(energy < 0.0);

    real e, f, a, omega, capm;
    {

      a = -0.5*Mtot/energy;
      const real fac = 1.0 - h2/(Mtot*a);
      assert(a > 0.0);


      real cape; 
      if (fac > LSMALL)
      {
        e = std::sqrt(fac);
        const real face = (a-r)/(a*e);
        cape = face > 0.9999999 ? 0.0 : (face > -0.9999999 ? std::acos(face) : M_PI);

        if (vdotr < 0.0) cape = 2.0*M_PI - cape;

        const real cf = (std::cos(cape) - e)/(1.0 - e*std::cos(cape));
        const real sf = std::sqrt(1-e*e)*std::sin(cape)/(1.0 - e*std::cos(cape));
        f = std::atan2(sf, cf);
        if (f < 0.0) f += 2.0*M_PI;
      }
      else
      {
        e = 0.0;
        f = u;
        cape = u;
      }

      capm = cape - e*std::sin(cape);
      omega = u - f;
      if (omega < 0) omega += 2.0*M_PI;
      omega = omega - int(omega/(2.0*M_PI))*2.0*M_PI;  /* longitude of pericentre */
    }

    real wp = capom + omega;
    real lambda = capm + wp;
    lambda = lambda - int(lambda/(2.0*M_PI))*2.0*M_PI;

    this->inc    = inc;
    this->sma    = a;
    this->ecc    = e;
    this->Omega  = wp;
    this->lambda = lambda;
  }
};

struct Nbody
{
  /*****************/

  unsigned long long iteration;
  unsigned long long iteration1;
  real time;
  Particle::Vector ptcl;

  /*****************/

  real dt;
  DH dh;
  
  unsigned long long flops;
  real tbeg;
  real dt_force;
  real dt_step, dt_multi, dt_extra, dt_err;

  int    NCSTEP;
  double a[17], b[17];

  void reset_counters()
  {
    flops = 0.0;
    dt_force = 0.0;
    dt_step = dt_multi = dt_extra = dt_err = 0.0;
    tbeg = mytimer::get_wtime();
  }
  
  real get_gflops() const
  {
    return flops/(mytimer::get_wtime() - tbeg+1.0/HUGE)/1e9;
  }

  real get_gflops_force() const
  {
    return flops / (dt_force + 1.0/HUGE) / 1e9;
  }

  Nbody(const unsigned long long i, const real _tepoch, const real dt, const Particle::Vector &_ptcl) : 
    iteration(i), time(_tepoch), ptcl(_ptcl) 
  {
    iteration1 = 0;
    assert(ptcl.size() > 1);
    reset_counters();

    move_to_CM(ptcl);
    
    /* bring to heliocentric positions */

    const real Mcentre = ptcl[0].mass;
    fprintf(stderr, "Mcentre= %g \n", Mcentre);

    const Particle p0 = ptcl[0];
    assert(p0.dt == -1.0);

    const int n = ptcl.size();
    real dt_min = HUGE;
    for (int i = 0; i < n-1; i++)
    {
      ptcl[i]      = ptcl[i+1];
      ptcl[i].pos -= p0.pos;
      assert(ptcl[i].dt > 0.0);
      dt_min = std::min(dt_min, ptcl[i].dt);
    }
    ptcl.resize(n-1);
    for (int i = 0; i < n-1; i++)
      ptcl[i].dt = dt_min;

    /* setting up the DH solver */

    dh        = DH(Mcentre);

    /* setting symplectic corrector */

    const real alpha = std::sqrt(7.0/40.0);
    const real  beta = 1.0/(48.0*alpha);

#if 1
    NCSTEP = 17;
    {
      a[8] = 1.0*alpha;
      b[8] = 45815578591785473.0/24519298961757600.0*beta;

      a[7] = 2.0*alpha;
      b[7] = -104807478104929387.0/80063017017984000.0*beta;

      a[6] = 3.0*alpha;
      b[6] = 422297952838709.0/648658702692000.0*beta;

      a[5] = 4.0*alpha;
      b[5] = -27170077124018711.0/112088223825177600.0*beta;

      a[4] = 5.0*alpha;
      b[4] = 102433989269.0/1539673404192.0*beta;

      a[3] = 6.0*alpha;
      b[3] = -33737961615779.0/2641809989145600.0*beta;

      a[2] = 7.0*alpha;
      b[2] = 26880679644439.0/17513784972684000.0*beta;

      a[1] = 8.0*alpha;
      b[1] = +682938344463443.0/7846175667762432000.0*beta;

      for (int i = 9; i <= 16; i++)
      {
        a[i] = -a[17-i];
        b[i] = -b[17-i];
      }
    }
#endif

#if 0
    NCSTEP = 2;
    {
      a[2] = alpha;
      b[2] = 0.5*beta;
      a[1] = -a[2];
      b[1] = -b[2];
    }
#endif

#if 0
    NCSTEP = 4;
    {
      a[2] = alpha;
      b[2] = 5.0/6.0 * beta;
      a[1] = 2.0*alpha;
      b[1] = -1.0/6.0 * beta;
      for (int i = 3; i <= 4; i++)
      {
        a[i] = -a[5-i];
        b[i] = -b[5-i];
      }
    }
#endif

#if 0
    NCSTEP = 6;
    {
      a[3] =     alpha;
      a[2] = 2.0*alpha;
      a[1] = 3.0*alpha;
      b[3] = 53521.0/49392.0 * beta;
      b[2] = -22651.0/61740.0 * beta;
      b[1] = 12361.0/246960.0*beta;
      for (int i = 4; i <= 6; i++)
      {
        a[i] = -a[7-i];
        b[i] = -b[7-i];
      }
    }
#endif

#if 0
    NCSTEP = 10;
    {
      a[5] = 1.0*alpha;
      a[4] = 2.0*alpha;
      a[3] = 3.0*alpha;
      a[2] = 4.0*alpha;
      a[1] = 5.0*alpha;
      b[5] = 3394141.0/2328480.0 * beta;
      b[4] = -14556229.0/19015920.0 * beta;
      b[3] = 895249.0/3622080.0 * beta;
      b[2] = -329447.0/6985440.0*beta;
      b[1] = 2798927.0/684573120.0 * beta;
      for (int i = 6; i <= 10; i++)
      {
        a[i] = -a[11-i];
        b[i] = -b[11-i];
      }
    }
#endif

#if 1
    ptcl = cvt2symp(ptcl);
#endif

  }

  private:

  void move_to_CM(Particle::Vector&ptcl) const
  {
    vec3 cm_pos(0.0), cm_vel(0.0);
    real Mtot = 0.0;
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      assert(it->mass > 0.0);
      Mtot   += it->mass;
      cm_pos += it->mass*it->pos;
      cm_vel += it->momentum();
    }

    cm_pos *= 1.0/Mtot;
    cm_vel *= 1.0/Mtot;

    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      it->pos -= cm_pos;
      it->vel -= cm_vel;
    }
  }

  real Epot(const Particle::Vector &ptcl) const
  {
    const int nbody = ptcl.size();
    real gpot = 0.0;

    for (int i = 0; i < nbody-1; i++)
      for (int j = i+1; j < nbody; j++)
        gpot -= ptcl[i].mass*ptcl[j].mass/(ptcl[i].pos - ptcl[j].pos).abs();

    for (int i = 0; i < nbody; i++)
      gpot -= dh.Mcentre*ptcl[i].mass/ptcl[i].pos.abs();

    return gpot;
  }

  real Ekin(const Particle::Vector &ptcl) const
  {
    real Ekin = 0.0;
    vec3 cmom(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      Ekin += it->Ekin();
      cmom += it->momentum();
    }

    const real Eckin = cmom.norm2()*0.5/dh.Mcentre;

    return Ekin + Eckin;
  }

  public:

  real Etot() const
  {
#if 1
    const Particle::Vector ptcl = cvt2phys(this->ptcl);
#else
    const Particle::Vector &ptcl = this->ptcl;
#endif
    return Ekin(ptcl) + Epot(ptcl);
  }

  std::string print_orbit(const int i) const
  {
#if 1
    const Particle::Vector ptcl = cvt2phys(this->ptcl);
#else
    const Particle::Vector &ptcl = this->ptcl;
#endif

    vec3 cVel(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
      cVel += it->momentum();
    cVel *= -1.0/dh.Mcentre;

    const Particle &p = ptcl[i];
    const vec3 R = p.pos;
    const vec3 V = p.vel - cVel;
    const real Mtot = p.mass + dh.Mcentre;

    const Orbit orb(R, V, Mtot);
    std::stringstream oss;
    oss <<  "#" <<  i << ": I= " << orb.inc << " a= " << orb.sma << " e= " <<  orb.ecc << 
      " o= " << orb.Omega << " l= " << orb.lambda;
    return oss.str();
  }

  std::string print_output() const
  {
#if 1
    Particle::Vector ptcl = cvt2phys(this->ptcl);
#else
    Particle::Vector ptcl(this->ptcl);
#endif

    const int n = ptcl.size();
    for (int i = 1; i < n; i++)
      ptcl[i] = this->ptcl[i-1];

    ptcl[0].mass = dh.Mcentre;
    ptcl[0].pos = ptcl[0].vel = 0.0;


    move_to_CM(ptcl);

    std::stringstream oss;
    for (int i = 0; i < n; i++)
    {
      oss.precision(15);
      oss 
        << i+1 << " "
        << std::scientific << " "
        << ptcl[i].mass << " "
        << ptcl[i].pos.x << " "
        << ptcl[i].pos.y << " "
        << ptcl[i].pos.z << " "
        << ptcl[i].vel.x << " "
        << ptcl[i].vel.y << " "
        << ptcl[i].vel.z << " "
        << std::endl;
    }
    return oss.str();
  }

  void Xstep(const real a, const real b, Particle::Vector &ptcl) const
  {
    static std::vector<int> active_list;
    if (active_list.empty())
    {
      const int n = ptcl.size();
      active_list.resize(n);
      for (int i = 0; i < n; i++)
        active_list[i] = i;
    }
    dh.kepler_drift(ptcl, active_list,  a);
    dh.linear_drift(ptcl, active_list,  b);
    dh.kick        (ptcl, active_list,  b);
    dh.kepler_drift(ptcl, active_list,  -a);
  }
  void Zstep(const real a, const real b, Particle::Vector &ptcl) const
  {
    Xstep( a,  b, ptcl);
    Xstep(-a, -b, ptcl);
  }

  Particle::Vector cvt2symp(const Particle::Vector &ptcl_in) const
  {
    Particle::Vector ptcl(ptcl_in);

#ifndef USE_SYMPC
    return ptcl;
#endif

    for (int i = NCSTEP-1; i >= 1; i--)
      Zstep(a[i], -b[i], ptcl);

    return ptcl;
  }

  Particle::Vector cvt2phys(const Particle::Vector &ptcl_in) const
  {
    Particle::Vector ptcl(ptcl_in);

#ifndef USE_SYMPC
    return ptcl;
#endif

    for (int i = 1; i <= NCSTEP-1; i++)
      Zstep(a[i], b[i], ptcl);

    return ptcl;
  }



  void iterate(const int nsteps = 32)
  {
#if 0
    ptcl = cvt2symp(ptcl);
#endif

    const int n = ptcl.size();
    static std::vector<int> active_list;
    if (active_list.empty())
    {
      active_list.resize(n);
      for (int i = 0; i < n; i++)
        active_list[i] = i;
    }

    dh.drift1(ptcl, active_list);

    const real time0 = time;
    for (int i = 0; i < nsteps; i++)
    {
      const real dt = ptcl[0].dt;

      dh.kick(ptcl, active_list);

      if (i < nsteps-1) dh.staggered_drift (ptcl, active_list);
      else              dh.          drift2(ptcl, active_list);

      time += dt;
      iteration1++;
    }

    dt = time - time0;
    iteration++;

#if 0
    change_dt();
#endif

#if 0
    ptcl = cvt2phys(ptcl);
#endif
  }
  
  void change_dt(const real eta = 1.0/20, const real f1 = 0.01)
  {
    const int n = ptcl.size();
    
    vec3 cVel(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
      cVel += it->momentum();
    cVel *= -1.0/dh.Mcentre;

    real dt_new = HUGE;
    for (int i = 0; i < n; i++)
    {
      const Particle &p = ptcl[i];
      const vec3 R = p.pos;
      const vec3 V = p.vel - cVel;
      const real Mtot = p.mass + dh.Mcentre;

      const Orbit orb(R, V, Mtot);
      const real rp = orb.sma*(1.0 - orb.ecc);
      assert(rp > 0.0);
      
      const real tp = eta*2*M_PI*std::sqrt(rp*rp*rp/Mtot);
      dt_new = std::min(dt_new, tp);
    }

    const real dt = ptcl[0].dt;
    if (std::abs(dt_new - dt) > f1*dt)
    {
      fprintf(stderr, " -- converting : dt_old= %g  dt_new= %g-- \n", dt, dt_new);
      ptcl = cvt2phys(ptcl);
      for (int i = 0; i < n; i++)
        ptcl[i].dt = dt_new;
      ptcl = cvt2symp(ptcl);
    }
    
  }

};


#endif /* _DHSHARED_H_ */
