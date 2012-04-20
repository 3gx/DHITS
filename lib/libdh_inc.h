#ifndef _LIBDH_INC_H_
#define _LIBDH_INC_H_

struct libDH;

namespace libDH_cpp
{
  void open(const double);
  void close();
  void iterate(
      const int,
      const double,
      const int,
      const double[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[],
      double&,
      int&);
  
  void cvt2phys(
      const int,
      const double,
      const double[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[]);
  
  void cvt2symp(
      const int,
      const double,
      const double[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[]);
  
  double Ekin(
      const int,
      const double[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[]);
 
  double Epot(
      const int,
      const double[],
      double  x[], double  y[], double  z[],
      double vx[], double vy[], double vz[]);


}

extern "C"
{
  void dh_open_(double*);
  void dh_close_();
  void dh_iterate_(
      int *n,
      double *dt,
      int *nstep,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz,
      double *time,
      int *iteration);

  void dh_cvt2phys_(
      int *n,
      double *dt,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz);
  
  void dh_cvt2symp_(
      int *n,
      double *dt,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz);
  
  void dh_ekin_(
      double *Ekin,
      int *n,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz);
  
  
  void dh_epot_(
      double *Epot,
      int *n,
      double *m,
      double * x, double * y, double * z,
      double *vx, double *vy, double *vz);
  
  void get_wtime_(double *time);

}

#endif /* _LIBDH_INC_H_ */
