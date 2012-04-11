#ifndef __DH_H__
#define __DH_H__

#include <iostream>
#include "vector3.h"
#include "particle.h"
#include "kepler.h"
#include "mytimer.h"

struct DH
{
  real Mcentre;

  DH(const real M = 0.0) : Mcentre(M) {}
  
  void kepler_drift(Particle::Vector &ptcl, const real dt) const
  {
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      Particle &p = *it;
      const Kepler kep(p.pos, p.vel, Mcentre, dt);
      p.pos = kep.pos;
      p.vel = kep.vel;
    }
  }
  
  void linear_drift(Particle::Vector &ptcl, const real dt) const
  {
    vec3 cmom(0.0);
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
      cmom += it->momentum();

    const vec3 dr = cmom*(dt/Mcentre);

    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
      it->pos += dr;
  }


  void drift1(Particle::Vector &ptcl, const real dt) const
  {
    kepler_drift(ptcl, 0.5*dt);
    linear_drift(ptcl, 0.5*dt);
  }
  void drift2(Particle::Vector &ptcl, const real dt) const
  {
    linear_drift(ptcl, 0.5*dt);
    kepler_drift(ptcl, 0.5*dt);
  }
  
  void staggered_drift(Particle::Vector &ptcl, const real dt) const
  {
    linear_drift(ptcl, 0.5*dt);
    kepler_drift(ptcl, 1.0*dt);
    linear_drift(ptcl, 0.5*dt);
  }

#ifndef _SSE_
  void kick(Particle::Vector &ptcl, const real dt) const
  {
    const int n = ptcl.size();
    for (int i = 0; i < n-1; i++)
      for (int j = i+1; j < n; j++)
      {
        const vec3 dr   = ptcl[j].pos - ptcl[i].pos;
        const real ds2  = dr.norm2();
        const real ids2 = 1.0/ds2;
        const real ids  = std::sqrt(ids2);
        const real ids3 = ids2*ids; 
        const vec3 aij  = ids3*dr;
        ptcl[i].vel += aij*(+ptcl[j].mass*dt);
        ptcl[j].vel += aij*(-ptcl[i].mass*dt);
      }
  }
#else
  void kick(Particle::Vector &ptcl, const real sdt) const
  {
    const int n = ptcl.size();
    for (int i = 0; i < n-1; i++)
    {
      const v2df dt      = (v2df){          sdt,           sdt};
      const v2df massidt = (v2df){-ptcl[i].mass, -ptcl[i].mass}*dt;
      for (int j = i+1; j < n; j += 2)
      {
        const int jp     = std::min(j+1, n-1);

        const v2df massj = {ptcl[j].mass, (double)(jp-j)*ptcl[jp].mass};
        const v2df dx    = {ptcl[j].pos.x - ptcl[i].pos.x, ptcl[jp].pos.x - ptcl[i].pos.x};
        const v2df dy    = {ptcl[j].pos.y - ptcl[i].pos.y, ptcl[jp].pos.y - ptcl[i].pos.y};
        const v2df dz    = {ptcl[j].pos.z - ptcl[i].pos.z, ptcl[jp].pos.z - ptcl[i].pos.z};
        const v2df r2    = dx*dx + dy*dy + dz*dz;
        const v2df invr  = __rsqrtpd(r2);
        const v2df invr2 = invr*invr;
        const v2df invr3 = invr*invr2;

        const v2df ax    = invr3*dx;
        const v2df ay    = invr3*dy;
        const v2df az    = invr3*dz;

        ptcl[i].vel.x += reduce(ax*massj*dt);
        ptcl[i].vel.y += reduce(ay*massj*dt);
        ptcl[i].vel.z += reduce(az*massj*dt);

        const v2df axj = ax * massidt;
        const v2df ayj = ay * massidt;
        const v2df azj = az * massidt;
        ptcl[j].vel.x += __builtin_ia32_vec_ext_v2df(axj,0);
        ptcl[j].vel.y += __builtin_ia32_vec_ext_v2df(ayj,0);
        ptcl[j].vel.z += __builtin_ia32_vec_ext_v2df(azj,0);

        if (jp != j)
        {
          ptcl[jp].vel.x += __builtin_ia32_vec_ext_v2df(axj,1);
          ptcl[jp].vel.y += __builtin_ia32_vec_ext_v2df(ayj,1);
          ptcl[jp].vel.z += __builtin_ia32_vec_ext_v2df(azj,1);
        }
      }
    }
  }
#endif

};

#endif // __DH_H__
