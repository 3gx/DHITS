#include <cstdlib>
#include <fstream>

#if 1
#define _SSE_
#endif

#include "libdh.h"
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

real read_xyz(Particle::Vector &ptcl, const int nbody)
{
  assert(nbody > 0);
  real mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);
  
  real Mcentre;
  std::cin >> Mcentre;
  Mcentre *= mass_scale;
  fprintf(stderr, " Mcentre= %g \n", Mcentre);

  for (int i = 0; i < nbody; i++)
  {
    int id;
    real mass;
    vec3 pos, vel;
    std::cin >> id >> mass >> pos.x >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z;
    mass *= mass_scale;
    pos  *= pos_scale;
    vel  *= vel_scale;
    
    ptcl.push_back(Particle(mass, pos, vel, id));
  }

  fprintf(stderr, "read= %d  particles \n", nbody);

  return Mcentre;
}

real read_aei(Particle::Vector &ptcl, const int nbody)
{
  assert(nbody > 0);
  
  real mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);
  
  real Mcentre;
  std::cin >> Mcentre;
  Mcentre *= mass_scale;
  fprintf(stderr, " Mcentre= %g \n", Mcentre);

  for (int i = 0; i < nbody; i++)
  {
    int id;
    real mass;
    real a; // semi-major axis
    real e; // eccentricity
    real I; // inclination in degrees;
    real w; // argument of the pericentre in degrees
    real O; // longitude of the ascending node in degrees
    real M; // mean anomaly in degrees;
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

    const real Mt   = Mcentre + mass;
    const real E    = Kepler::EccentricAnomaly(e, M);
    const real Edot = sqrt(Mt/(a*a*a))/(1 - e*cos(E));

    const real x    =  a*(cos(E) - e);
    const real y    =  a*sqrt(1.0 - e*e)*sin(E);
    const real vx   = -a*sin(E)*Edot;
    const real vy   =  a*sqrt(1.0 - e*e)*cos(E)*Edot;

    real r, f;

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

    ptcl.push_back(Particle(mass, pos, vel, id));
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
  Particle::Vector ptcl;
  
  int nbody, npl_out;
  std::cin >> nbody >> npl_out;
  assert(nbody  != 0);
  assert(npl_out >= 0);

  real time, Tscale;
  std::cin >> time >> Tscale;

  int dI;
  real Tend, dt_log, dt_out, dt_snap;
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

  real dt;
  std::cin >> dt;
  fprintf(stderr, " dt= %g \n", dt);
  dt *= Tscale;

  real Mcentre = -1;
  if (nbody < 0)
  {
    fprintf(stderr, " Reading aei format ... \n");
    Mcentre = read_aei(ptcl, -nbody);
  }
  else
  {
    fprintf(stderr, " Reading xyz format ... \n");
    Mcentre = read_xyz(ptcl, nbody);
  }


  libDH dh(Mcentre);

  const real E0 = Etot(ptcl, Mcentre);
  fprintf(stderr, " time= %g   Etot= %g \n", time, E0);

  real Ep = E0;

  real t_log  = time/Tscale;
  real t_out  = time/Tscale;
  real t_snap = time/Tscale;

  real de_max = 0.0;

  ptcl = dh.cvt2symp(ptcl, dt);

  double t0 = mytimer::get_wtime();
  const double t_start = t0;
  unsigned long long iteration1 = 0;
  while (time/Tscale < Tend)
  {
    dh.iterate(ptcl, dt, time, iteration1);
    iteration++;
    const double t1 = mytimer::get_wtime() + 1.0/HUGE;
    if (time/Tscale > t_log || time/Tscale >= Tend || iteration % di_iter == 0)
    {
      const real E1 = Etot(dh.cvt2phys(ptcl, dt), Mcentre);
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
          (real)iteration1/(real)iteration);
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

  return 0;

}
