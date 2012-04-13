#include <cstdlib>
#include <fstream>

#if 1
#define _SSE_
#endif

#include <cassert>
#include <iostream>
#include <vector>
#include "vector3.h"
#include "kepler.h"
#include "libdh_inc.h"
#include "mytimer.h"

#ifndef __MACOSX_
#define __LINUX__
#endif

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
  _mm_setcsr( _MM_MASK_MASK &~
              (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) {
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
crap
void fpe_catch(void) {}
#endif


double read_xyz(
    int    *id,
    double *mass,
    double *x,
    double *y,
    double *z,
    double *vx,
    double *vy,
    double *vz,
    int nbody)
{
  assert(nbody > 0);
  double mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);
  
  double Mcentre;
  std::cin >> Mcentre;
  Mcentre *= mass_scale;
  fprintf(stderr, " Mcentre= %g \n", Mcentre);

  for (int i = 0; i < nbody; i++)
  {
    std::cin >> id[i] >> mass[i] >> 
       x[i] >>  y[i] >>  z[i] >> 
      vx[i] >> vy[i] >> vz[i];

    mass[i] *= mass_scale;
    x [i] *= pos_scale;
    y [i] *= pos_scale;
    z [i] *= pos_scale;
    vx[i] *= vel_scale;
    vy[i] *= vel_scale;
    vz[i] *= vel_scale;
  }

  fprintf(stderr, "read= %d  particles \n", nbody);

  return Mcentre;
}

double read_aei(
    int    *iid,
    double *imass,
    double *xc,
    double *yc,
    double *zc,
    double *vxc,
    double *vyc,
    double *vzc,
    int nbody)
{
  assert(nbody > 0);

  double mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);

  double Mcentre;
  std::cin >> Mcentre;
  Mcentre *= mass_scale;
  fprintf(stderr, " Mcentre= %g \n", Mcentre);

  for (int i = 0; i < nbody; i++)
  {
    int id;
    double mass;
    double a; // semi-major axis
    double e; // eccentricity
    double I; // inclination in degrees;
    double w; // argument of the pericentre in degrees
    double O; // longitude of the ascending node in degrees
    double M; // mean anomaly in degrees;
    std::cin >>
      id   >> 
      mass >> 
      a    >> 
      e    >> 
      I    >> 
      w    >> 
      O    >> 
      M;
    fprintf(stderr, " i =%d : index= %d mass= %g  a= %g  e= %g  I= %g  w= %g  O= %g  M= %g :: r_peri= %g  r_apo= %g\n",
        i, id,  mass,
        a, e, I, w, O, M,
        a*(1.0 -e), a*(1.0 +e));

    w *= M_PI/180.0;
    I *= M_PI/180.0;
    M *= M_PI/180.0;
    O *= M_PI/180.0;

    const double Mt   = Mcentre + mass;
    const double E    = Kepler::EccentricAnomaly(e, M);
    const double Edot = sqrt(Mt/(a*a*a))/(1 - e*cos(E));

    const double x    =  a*(cos(E) - e);
    const double y    =  a*sqrt(1.0 - e*e)*sin(E);
    const double vx   = -a*sin(E)*Edot;
    const double vy   =  a*sqrt(1.0 - e*e)*cos(E)*Edot;

    double r, f;

    r = sqrt(x*x+y*y);
    f = atan2(y,x);

    const vec3 pos(
        r*(cos(O)*cos(w+f) - sin(O)*sin(w+f)*cos(I)),
        r*(sin(O)*cos(w+f) + cos(O)*sin(w+f)*cos(I)),
        r*(sin(w+f)*sin(I)));

    r = sqrt(vx*vx+vy*vy);
    f = atan2(vy,vx);

    const vec3 vel(
        r*(cos(O)*cos(w+f) - sin(O)*sin(w+f)*cos(I)),
        r*(sin(O)*cos(w+f) + cos(O)*sin(w+f)*cos(I)),
        r*(sin(w+f)*sin(I)));

    fprintf(stderr, " m= %g  pos= %g %g %g  vel= %g %g %g \n",
        mass,
        pos.x, pos.y, pos.z,
        vel.x, vel.y, vel.z);

    iid[i] = id;
    imass[i] = mass;
    xc [i] = pos.x;
    yc[i] = pos.y;
    zc[i] = pos.z;
    vxc[i] = vel.x;
    vyc[i] = vel.y;
    vzc[i] = vel.z;
  }

  fprintf(stderr, "read= %d  particles \n", nbody);

  return Mcentre;
}

int main(int argc, char * argv[])
{
#ifndef _SSE_
  fpe_catch();
#endif

  char path[256] = "/tmp/out";
  if (argc > 1)
    sprintf(path, "%s", argv[1]);

  fprintf(stderr, " Writing output to %s \n", path);

  unsigned long long iteration = 0;

  int nbody, npl_out;
  std::cin >> nbody >> npl_out;
  assert(nbody  != 0);
  assert(npl_out >= 0);

  double time, Tscale;
  std::cin >> time >> Tscale;

  int dI;
  double Tend, dt_log, dt_out, dt_snap;
  std::cin >> Tend >> dt_out >> dt_log >> dt_snap >> dI;
  assert(dt_out > 0.0);
  assert(dt_log > 0.0);
  assert(dt_snap > 0.0);
  assert(Tscale > 0.0);
  const unsigned long long di_iter_max = 1LLU << 62;
  const unsigned long long di_iter = dI < 0 ? di_iter_max : dI;

  fprintf(stderr, " Time= %g   Tscale= %g \n", time, Tscale);
  fprintf(stderr, " Tend= %g \n", Tend);
  fprintf(stderr, " dTout= %g \n", dt_out);
  fprintf(stderr, " dTlog= %g \n", dt_log);
  fprintf(stderr, " dTsnap= %g \n", dt_snap);

  double dt;
  std::cin >> dt;
  fprintf(stderr, " dt= %g \n", dt);
  dt *= Tscale;

  double Mcentre = -1;

  int n = std::abs(nbody);
  int *id = new int[n];
  double *mass = new double[n];
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double *vx = new double[n];
  double *vy = new double[n];
  double *vz = new double[n];
  
  double *x1 = new double[n];
  double *y1 = new double[n];
  double *z1 = new double[n];
  double *vx1 = new double[n];
  double *vy1 = new double[n];
  double *vz1 = new double[n];

  
  if (nbody < 0)
  {
    fprintf(stderr, " Reading aei format ... \n");
    Mcentre = read_aei(id, mass, x,y,z, vx,vy,vz, -nbody);
  }
  else
  {
    fprintf(stderr, " Reading xyz format ... \n");
    Mcentre = read_xyz(id, mass, x,y,z, vx,vy,vz, nbody);
  }


  dh_open_(&Mcentre);

  const double E0 = dh_Ekin_(&n, mass, x,y,z, vx,vy,vz) + dh_Epot_(&n, mass, x,y,z, vx,vy,vz);
  fprintf(stderr, " time= %g   Etot= %g \n", time, E0);

  double Ep = E0;

  double t_log  = time/Tscale;
  double t_out  = time/Tscale;
  double t_snap = time/Tscale;

  double de_max = 0.0;

  dh_cvt2symp_(&n, &dt, mass, x,y,z, vx,vy,vz);

  double t0 = mytimer::get_wtime();
  const double t_start = t0;
  unsigned long long iteration1 = 0;
  while (time/Tscale < Tend)
  {
    int nstep = 32;
    int iter1 = 0;
    dh_iterate_(&n, &dt, &nstep, mass, x,y,z, vx,vy,vz, &time, &iter1);
    iteration1 += iter1;
    iteration++;

    const double t1 = mytimer::get_wtime() + 1.0/HUGE;
    if (time/Tscale > t_log || time/Tscale >= Tend || iteration % di_iter == 0)
    {
      for (int i = 0; i < n; i++)
      {
        x1[i] = x[i];
        y1[i] = y[i];
        z1[i] = z[i];
        vx1[i] = vx[i];
        vy1[i] = vy[i];
        vz1[i] = vz[i];
      }
      dh_cvt2phys_(&n, &dt, mass, x1,y1,z1, vx1,vy1,vz1);
      const double E1 = dh_Ekin_(&n, mass, x1,y1,z1, vx1,vy1,vz1) + dh_Epot_(&n, mass, x1,y1,z1, vx1,vy1,vz1);
      
      de_max = std::max(de_max, std::abs((E1 - E0)/E0));
      fprintf(stderr, "iter= %llu :: t= %g dt= %4.3g  dE= %4.3g ddE= %4.3g dEmax= %4.3g  Twall= %4.3g hr | <T>= %4.3g sec | iter1/iter= %g\n",
          iteration,
          time/Tscale, 
          dt/Tscale,
          (E1 - E0)/std::abs(E0),
          (E1 - Ep)/std::abs(Ep),
          de_max,
          (t1 - t_start) / 3600.0,
          t1 - t0,
          (double)iteration1/(double)iteration);

      t0 = t1;
      t_log += dt_log;
      Ep = E1;
    }

#if 0
    if (time/Tscale > t_out || time/Tscale >= Tend)
    {
      fprintf(stdout, "%g ", time/Tscale);
      for (int ipl = 0; ipl < npl_out; ipl++)
        fprintf(stdout,"%s ", s.print_orbit(ipl).c_str());
      fprintf(stdout, "\n");
      fflush(stdout);
      t_out += dt_out;
    }
#endif

#if 0
    if (s.time/Tscale >= t_snap || s.time/Tscale >= Tend)
    {
      t_snap += dt_snap;
      /*************/

      char fn[256];
      sprintf(fn, "%s/snap_%.8d.out", path, int(t_snap/dt_snap));
      std::ofstream fout(fn);


      char o1[256]; 
      sprintf(o1, "%lld \n", di_iter == di_iter_max ? -1 : di_iter);

      fout.precision(15);
      fout << s.ptcl.size()+1 << " " << npl << std::endl;
      fout << s.time << " " << Tscale << std::endl;
      fout << Tend << " " << dt_out << " " << dt_log  << " " << dt_snap << " " << o1;
      fout << dt << std::endl;
      fout << "1.0 1.0 1.0 \n";
      fout << s.print_output();

      fout.close();

      /*************/
    }
#endif
  }

  delete[] id;
  delete[] mass;
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] x1;
  delete[] y1;
  delete[] z1;
  delete[] vx1;
  delete[] vy1;
  delete[] vz1;
  dh_close_();

  return 0;

}
