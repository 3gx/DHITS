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
  }

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

  real Epot() const
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

  real Ekin() const
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

  real Etot() const
  {

    return Ekin() + Epot();
  }

  std::string print_orbit(const int i) const
  {
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
    Particle::Vector ptcl(this->ptcl.size() + 1);
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

  void iterate()
  {
    const int n = ptcl.size();
    std::vector<int> active_list(n);
    for (int i = 0; i < n; i++)
      active_list[i] = i;

    dh.drift1(ptcl, active_list);

    real time0 = time;
    while (1)
    {
      active_list.clear();
      const real dt = scheduler.pull_active_list(active_list);
      const bool full_step = active_list.size() == ptcl.size();
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

};


#endif /* _DHSHARED_H_ */
