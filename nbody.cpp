#include <cstdlib>
#include <fstream>

#if 1
#define _SSE_
#endif

#include "DHITS.h"
#include "kepler.h"
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

typedef double real;

/* position are in AU, G = 1, time is in years */
/* Period @ R = 1 AU is 1 year */

Particle::Vector read_xyz(const int nbody, const real dt)
{
  Particle::Vector ptcl;
  assert(nbody > 1);

  real mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);

  for (int i = 0; i < nbody; i++)
  {
    int idummy;
    real mass, dt_scale;
    vec3 pos, vel;
    std::cin >> idummy >> dt_scale >> mass >> pos.x >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z;
    mass *= mass_scale;
    pos  *= pos_scale;
    vel  *= vel_scale;
    
    ptcl.push_back(Particle(mass, pos, vel, dt_scale*dt));
  }

  ptcl[0].dt = -1.0;
  
  fprintf(stderr, "read= %d  particles \n", nbody);
  return ptcl;
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
  
  int nbody, npl;
  std::cin >> nbody >> npl;
  assert(nbody != 0);
  assert(npl >= 0);

  real tepoch, Tscale;
  std::cin >> tepoch >> Tscale;

  int dI;
  real Tend, dt_log, dt_out, dt_snap;
  std::cin >> Tend >> dt_out >> dt_log >> dt_snap >> dI;
  assert(dt_out > 0.0);
  assert(dt_log > 0.0);
  assert(dt_snap > 0.0);
  assert(Tscale > 0.0);
  const unsigned long long di_iter_max = 1LLU << 62;
  const unsigned long long di_iter = dI < 0 ? di_iter_max : dI;

  fprintf(stderr, " Tepoch= %g   Tscale= %g \n", tepoch, Tscale);
  fprintf(stderr, " Tend= %g \n", Tend);
  fprintf(stderr, " dTout= %g \n", dt_out);
  fprintf(stderr, " dTlog= %g \n", dt_log);
  fprintf(stderr, " dTsnap= %g \n", dt_snap);

  real dt;
  std::cin >> dt;
  fprintf(stderr, " dt= %g \n", dt);

  assert(nbody > 0);
  fprintf(stderr, " Reading xyz format ... \n");
  ptcl = read_xyz(nbody, dt*Tscale);

  const real dt_max = 128.0;
  Nbody s(iteration, tepoch, dt_max, ptcl);

  const real E0 = s.Etot();
  fprintf(stderr, " tepoch= %g   Etot= %g \n", tepoch, E0);

  real Ep = E0;

  real t_log  = s.time/Tscale;
  real t_out  = s.time/Tscale;
  real t_snap = s.time/Tscale;


  real de_max = 0.0;

  double t0 = mytimer::get_wtime();
  const double t_start = t0;
  while (s.time/Tscale < Tend)
  {
    s.iterate();
    const double t1 = mytimer::get_wtime() + 1.0/HUGE;
    if (s.time/Tscale > t_log || s.time/Tscale >= Tend || s.iteration % di_iter == 0)
    {
      const real E1 = s.Epot() + s.Ekin();
      de_max = std::max(de_max, std::abs((E1 - E0)/E0));
      fprintf(stderr, "iter= %llu :: t= %g dt= %4.3g  dE= %4.3g ddE= %4.3g dEmax= %4.3g  Twall= %4.3g hr | <T>= %4.3g sec | iter1/iter= %g\n",
          s.iteration,
          s.time/Tscale, 
          s.dt/Tscale,
          (E1 - E0)/std::abs(E0),
          (E1 - Ep)/std::abs(Ep),
          de_max,
          (t1 - t_start) / 3600.0,
          t1 - t0,
          (real)s.iteration1/(real)s.iteration);
      t0 = t1;
      s.reset_counters();
      t_log += dt_log;
      Ep = E1;
    }

    if (s.time/Tscale > t_out || s.time/Tscale >= Tend)
    {
      fprintf(stdout, "%g ", s.time/Tscale);
      for (int ipl = 0; ipl < npl; ipl++)
        fprintf(stdout,"%s ", s.print_orbit(ipl).c_str());
      fprintf(stdout, "\n");
      fflush(stdout);
      t_out += dt_out;
    }

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
  }

  return 0;

}
