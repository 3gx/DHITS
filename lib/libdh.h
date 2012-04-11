#ifndef _LIBDH_H_
#define _LIBDH_H_

#include "particle.h"
#include "dh.h"
  
inline real Ekin(const Particle::Vector &ptcl, const real M) 
{
  real Ekin = 0.0;
  vec3 cmom(0.0);
  for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
  {
    Ekin += it->Ekin();
    cmom += it->momentum();
  }

  const real Eckin = cmom.norm2()*0.5/M;

  return Ekin + Eckin;
}

inline real Epot(const Particle::Vector &ptcl, const real M)
{
  const int nbody = ptcl.size();
  real gpot = 0.0;

  for (int i = 0; i < nbody-1; i++)
    for (int j = i+1; j < nbody; j++)
      gpot -= ptcl[i].mass*ptcl[j].mass/(ptcl[i].pos - ptcl[j].pos).abs();

  for (int i = 0; i < nbody; i++)
    gpot -= M*ptcl[i].mass/ptcl[i].pos.abs();

  return gpot;
}

inline real Etot(const Particle::Vector &ptcl, const real M) 
{
  return Ekin(ptcl, M) + Epot(ptcl, M);
}


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

struct libDH
{
  DH dh;
  int NCSTEP;
  double a[18], b[17];

  libDH(const real M, const int sc_mode = 17) : dh(DH(M)) 
  {
    set_sc(sc_mode);
  };

  void set_sc(const int sc_mode)
  {
    const real alpha = std::sqrt(7.0/40.0);
    const real  beta = 1.0/(48.0*alpha);
    switch (sc_mode)
    {
      case(2):
        fprintf(stderr, " -- libDH: 2nd order symplectic corrector ");
        NCSTEP = 2;
        a[2] = alpha;
        b[2] = 0.5*beta;
        a[1] = -a[2];
        b[1] = -b[2];
        break;

      case(4):
        fprintf(stderr, " -- libDH: 4th order symplectic corrector ");
        NCSTEP = 4;
        a[2] = alpha;
        b[2] = 5.0/6.0 * beta;
        a[1] = 2.0*alpha;
        b[1] = -1.0/6.0 * beta;
        for (int i = 3; i <= 4; i++)
        {
          a[i] = -a[5-i];
          b[i] = -b[5-i];
        }
        break;

      case(6):
        fprintf(stderr, " -- libDH: 6th order symplectic corrector ");
        NCSTEP = 6;
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
        break;

      case(10):
        fprintf(stderr, " -- libDH: 10th order symplectic corrector ");
        NCSTEP = 10;
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
        break;

      case(17):
        fprintf(stderr, " -- libDH: 17th order symplectic corrector ");
        NCSTEP = 17;
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
        break;

      default:
        fprintf(stderr, " -- libDH: disabling symplectic corrector ");
        NCSTEP = -1;
    }
  }

  void iterate(Particle::Vector &ptcl, const real dt, real &t, unsigned long long &i, int nstep = 32) const
  {
    dh.drift1(ptcl, dt);
    while (nstep > 0)
    {
      nstep--;
      dh.kick(ptcl, dt);

      if (nstep > 0) dh.staggered_drift (ptcl, dt);
      else           dh.          drift2(ptcl, dt);

      t += dt;
      i++;
    }
  }

  void Xstep(const real a, const real b, Particle::Vector &ptcl) const
  {
    dh.kepler_drift(ptcl,  a);
    dh.linear_drift(ptcl,  b);
    dh.kick        (ptcl,  b);
    dh.kepler_drift(ptcl, -a);
  }
  void Zstep(const real a, const real b, Particle::Vector &ptcl) const
  {
    Xstep( a,  b, ptcl);
    Xstep(-a, -b, ptcl);
  }
  Particle::Vector cvt2symp(const Particle::Vector &ptcl_in, const real dt) const
  {
    Particle::Vector ptcl(ptcl_in);
    if (-1 == NCSTEP) return ptcl;
    for (int i = NCSTEP - 1; i >= 1; i--)
      Zstep(a[i]*dt, -b[i]*dt, ptcl);
    return ptcl;
  }
  Particle::Vector cvt2phys(const Particle::Vector &ptcl_in, const real dt) const
  {
    Particle::Vector ptcl(ptcl_in);
    if (-1 == NCSTEP) return ptcl;
    for (int i = 1; i <= NCSTEP-1; i++)
      Zstep(a[i]*dt, b[i]*dt, ptcl);
    return ptcl;
  }
};

#endif /* _LIBDH_H_ */
