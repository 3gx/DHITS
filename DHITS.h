#ifndef _DHSHARED_H_
#define _DHSHARED_H_

#include <sstream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "particle.h"
#include "dh.h"
#include "Scheduler.h"

#if 1   /* set to 1 in order to use Symplectic Corrector */
#define USE_SYMPC
#endif

#define RUNGMAX 16

struct Nbody
{
  /*****************/

  unsigned long long iteration;
  unsigned long long iteration1;
  real time;
  Particle::Vector ptcl;

  /*****************/

  DH dh;
  real dt;
  
  unsigned long long flops;
  real tbeg;
  real dt_force;
  real dt_step, dt_multi, dt_extra, dt_err;

  Scheduler<RUNGMAX> scheduler;

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
    for (int i = 0; i < n-1; i++)
    {
      ptcl[i]      = ptcl[i+1];
      ptcl[i].pos -= p0.pos;
      assert(ptcl[i].dt > 0.0);
    }
    ptcl.resize(n-1);

    /* setting up the scheduler */

    real dt_max = 1.0;
    while (dt > dt_max) dt_max *= 2.0;
    while (dt < dt_max) dt_max *= 0.5;
    
    scheduler = Scheduler<RUNGMAX>(dt_max);
    dh        = DH(Mcentre);

    for (int i = 0; i < n-1; i++)
      scheduler.push_particle(i, ptcl[i].dt);

    const real alpha = std::sqrt(7.0/40.0);
    const real  beta = 1.0/(48.0*alpha);

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

#if 0
    for (int i = 1; i <= 8; i++)
    {
      fprintf(stderr, "a[%d]= %lg  b[%d]= %15.15lg\n", i, a[i]/alpha, i, b[i]/beta);
    }
#endif

    ptcl = cvt2symp(ptcl);
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
    const Particle::Vector ptcl = cvt2phys(this->ptcl);
    return Ekin(ptcl) + Epot(ptcl);
  }

  std::string print_orbit(const int i) const
  {
    const Particle::Vector ptcl = cvt2phys(this->ptcl);

    vec3 cVel(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
      cVel += it->momentum();
    cVel *= -1.0/dh.Mcentre;

    const Particle &p = ptcl[i];
    const vec3 R = p.pos;
    const vec3 V = p.vel - cVel;
    const real Mtot = p.mass + dh.Mcentre;

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

    std::stringstream oss;
    oss <<  "#" <<  i << ": I= " << inc << " a= " << a << " e= " <<  e << 
      " o= " << wp << " l= " << lambda;
    return oss.str();
  }

  std::string print_output() const
  {
    Particle::Vector ptcl = cvt2phys(this->ptcl);

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

#if 0
    for (int i = 1; i < (const int)ptcl.size(); i++)
      ptcl[i].dt = ptcl[0].dt;
    for (int i = 16; i >= 1; i--)
      Zstep(a[i], -b[i], ptcl);
    for (int i = 1; i < (const int)ptcl.size(); i++)
      ptcl[i].dt = ptcl_in[i].dt;
#else
    for (int i = 16; i >= 1; i--)
      Zstep(a[i], -b[i], ptcl);
#endif

    return ptcl;
  }

  Particle::Vector cvt2phys(const Particle::Vector &ptcl_in) const
  {
    Particle::Vector ptcl(ptcl_in);

#ifndef USE_SYMPC
    return ptcl;
#endif

#if 0
    for (int i = 1; i < (const int)ptcl.size(); i++)
      ptcl[i].dt = ptcl[0].dt;
    for (int i = 1; i <= 16; i++)
      Zstep(a[i], b[i], ptcl);
    for (int i = 1; i < (const int)ptcl.size(); i++)
      ptcl[i].dt = ptcl_in[i].dt;
#else
    for (int i = 1; i <= 16; i++)
      Zstep(a[i], b[i], ptcl);
#endif

    return ptcl;
  }

#if 0

  void iterate(const int nsteps = 16)
  {
    const int n = ptcl.size();
    std::vector<int> active_list(n);
    for (int i = 0; i < n; i++)
      active_list[i] = i;

    dh.drift1(ptcl, active_list);

    const real time0 = time;
    int step = 1;
    while (1)
    {
      step++;
      active_list.clear();
      const real dt = scheduler.pull_active_list(active_list);
      const bool full_step = active_list.size() == ptcl.size() && step > nsteps;
      assert(!active_list.empty());

      dh.kick(ptcl, active_list);

      if (!full_step) dh.staggered_drift (ptcl, active_list);
      else            dh.          drift2(ptcl, active_list);

      for (std::vector<int>::const_iterator it = active_list.begin(); it != active_list.end(); it++)
        scheduler.push_particle(*it, ptcl[*it].dt);

      time += dt;
      iteration1++;
      if (full_step)
        break;
    }

    dt = time - time0;
    iteration++;
  }

#else  /* experimental!!! */

  double compute_dt(const Particle::Vector &ptcl) const
  {
    vec3 cVel(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
      cVel += it->momentum();

    cVel *= -1.0/dh.Mcentre;

    double dt_min = HUGE;

    const int n = ptcl.size();

    for (int i = 0; i < n; i++)
    {
      /* compute curvature */
      const vec3 pos = ptcl[i].pos;
      const vec3 vel = ptcl[i].vel - cVel;

      const double r =  pos.abs();
      const double v =  vel.abs();
      const vec3 acc = -dh.Mcentre/(r*r*r)*pos;

      const vec3 VxA      = vel % acc;
      const double vxa      = VxA.abs();
      const double R        = v*v*v/vxa;
      const double dphi     = 2*M_PI/100;
      const double dt_curv  = dphi * R/v;

      const vec3 jrk      = -3.0*acc*(vel*pos)/(r*r) - dh.Mcentre/(r*r*r)*vel;
      const double Rdot     =  3.0*v*(vel*acc)/vxa - v*v*v/(vxa*vxa*vxa) * (VxA*(vel%jrk));
      const double dt_Rdot  =  0.01*R/(std::abs(Rdot) + 1.0/HUGE);

      dt_min = std::min(dt_min, 1.0/(1.0/dt_curv + 1.0/dt_Rdot));
    }

    //    return dt_max;
#if 1
    return dt_min; //std::min(dt_min, dt_max);
#else
    double dt1 = dt;
    while (    dt1 > dt_min) dt1 *= 0.5;
    while (2.0*dt1 < dt_min) dt1 *= 2.0;
    return std::min(dt1, dt_max);
#endif
  }

  void iterate()
  {
    static std::vector<int> active_list;
    if (active_list.empty())
    {
      const int n = ptcl.size();
      active_list.resize(n);
      for (int i = 0; i < n; i++)
        active_list[i] = i;

#if 1
      for (int i = 1; i < (const int)ptcl.size(); i++)
        ptcl[i].dt = ptcl[0].dt;
#endif

      scheduler.flush_list();
    }
    

#if 1
    ptcl = cvt2phys(ptcl);

    const real dt_step = compute_dt(ptcl);
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
      it->dt = dt_step;

    ptcl = cvt2symp(ptcl);
#endif


    const real time0 = time;

    dh.drift1(ptcl, active_list);
    const int Niter = 16;
    for (int n = 0; n < Niter; n++)
    {
      dh.kick(ptcl, active_list);
      if (n < Niter-1) dh.staggered_drift (ptcl, active_list);
      else             dh.          drift2(ptcl, active_list);

      time += ptcl[0].dt;
      iteration1++;
    }

    dt = time - time0;
    iteration++;
  }
#endif

};


#endif /* _DHSHARED_H_ */
