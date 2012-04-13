#include "libdh.h"

namespace libDH_cpp
{
  libDH *dh = NULL;

  void open(const double M)
  {
    assert(NULL == dh);
    dh = new libDH(M);
  };

  void close()
  {
    assert (NULL != dh);
    delete dh;
    dh = NULL;
  };

  void iterate(
      const int n,
      const double dt,
      const int nstep,
      const double m[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[],
      double &time,
      int    &iteration)
  {
    assert(NULL != dh);
    static Particle::Vector ptcl;
    ptcl.clear();
    for (int i = 0; i < n; i++)
      ptcl.push_back(Particle(m[i], vec3(x[i], y[i], z[i]), vec3(vx[i], vy[i], vz[i])));

    unsigned long long i = 0;

    dh->iterate(ptcl, dt, time, i, nstep);
    iteration += i;

    for (int i = 0; i < n; i++)
    {
      x [i] = ptcl[i].pos.x;
      y [i] = ptcl[i].pos.y;
      z [i] = ptcl[i].pos.z;
      vx[i] = ptcl[i].vel.x;
      vy[i] = ptcl[i].vel.y;
      vz[i] = ptcl[i].vel.z;
    }
  }

  void cvt2phys(
      const int n,
      const double dt,
      const double m[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[])
  {
    assert(NULL != dh);
    static Particle::Vector ptcl;
    ptcl.clear();
    for (int i = 0; i < n; i++)
      ptcl.push_back(Particle(m[i], vec3(x[i], y[i], z[i]), vec3(vx[i], vy[i], vz[i])));

    ptcl = dh->cvt2phys(ptcl, dt);
    
    for (int i = 0; i < n; i++)
    {
      x [i] = ptcl[i].pos.x;
      y [i] = ptcl[i].pos.y;
      z [i] = ptcl[i].pos.z;
      vx[i] = ptcl[i].vel.x;
      vy[i] = ptcl[i].vel.y;
      vz[i] = ptcl[i].vel.z;
    }
  }
  
  void cvt2symp(
      const int n,
      const double dt,
      const double m[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[])
  {
    assert(NULL != dh);
    static Particle::Vector ptcl;
    ptcl.clear();
    for (int i = 0; i < n; i++)
      ptcl.push_back(Particle(m[i], vec3(x[i], y[i], z[i]), vec3(vx[i], vy[i], vz[i])));

    ptcl = dh->cvt2symp(ptcl, dt);
    
    for (int i = 0; i < n; i++)
    {
      x [i] = ptcl[i].pos.x;
      y [i] = ptcl[i].pos.y;
      z [i] = ptcl[i].pos.z;
      vx[i] = ptcl[i].vel.x;
      vy[i] = ptcl[i].vel.y;
      vz[i] = ptcl[i].vel.z;
    }
  }
  
  double Ekin(
      const int n,
      const double m[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[])
  {
    assert(NULL != dh);
    static Particle::Vector ptcl;
    ptcl.clear();
    for (int i = 0; i < n; i++)
      ptcl.push_back(Particle(m[i], vec3(x[i], y[i], z[i]), vec3(vx[i], vy[i], vz[i])));

    return Ekin(ptcl, dh->dh.Mcentre);
  }
  
  double Epot(
      const int n,
      const double m[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[])
  {
    assert(NULL != dh);
    static Particle::Vector ptcl;
    ptcl.clear();
    for (int i = 0; i < n; i++)
      ptcl.push_back(Particle(m[i], vec3(x[i], y[i], z[i]), vec3(vx[i], vy[i], vz[i])));

    return Epot(ptcl, dh->dh.Mcentre);
  }
}

extern "C"
{
  void dh_open_(double *M) {libDH_cpp::open(*M);}
  void dh_close_()         {libDH_cpp::close() ;}
  void dh_iterate_(
      int *n,
      double *dt,
      int *nstep,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz,
      double *time,
      int *iteration)
  {
    libDH_cpp::iterate(*n, *dt, *nstep, m, x,y,z, vx,vy,vz, *time, *iteration);
  }

  void dh_cvt2phys_(
      int *n,
      double *dt,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz)
  {
    libDH_cpp::cvt2phys(*n, *dt, m, x,y,z, vx, vy, vz);
  }
  
  void dh_cvt2symp_(
      int *n,
      double *dt,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz)
  {
    libDH_cpp::cvt2symp(*n, *dt, m, x,y,z, vx, vy, vz);
  }
  
  double dh_Ekin_(
      int *n,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz)
  {
    return libDH_cpp::Ekin(*n, m, x,y,z, vx,vy,vz);
  }
  
  double dh_Epot_(
      int *n,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz)
  {
    return libDH_cpp::Epot(*n, m, x,y,z, vx,vy,vz);
  }

}


