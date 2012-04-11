#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <iostream>
#include <vector>
#include <cassert>
#include "vector3.h"

struct Particle
{
  typedef std::vector<Particle> Vector;
  real  mass;
  vec3   pos;
  vec3   vel;
  real    id;

  Particle() {}
  Particle(const real v) : mass(v), pos(v), vel(v) {}
  Particle(const real m, const vec3 &_pos, const vec3 &_vel, const real _id = -1) :
    mass(m), pos(_pos), vel(_vel), id(_id) {}

  real Ekin() const {return mass*vel.norm2()*0.5;}
  vec3 momentum() const {return mass*vel;}
};

#ifdef _SSE_
typedef double v2df __attribute__((vector_size(16)));
typedef float  v4sf __attribute__((vector_size(16)));

inline void __divpd(double &a, double &b)
{
  const v2df one = {1.0, 1.0};
#if 1
  const v2df v = one/(v2df){a, b};
#else
  const v2df A = {a, b};
  const v4sf x = __builtin_ia32_rcpps((v4sf){a, b, a, b});
  v2df v = __builtin_ia32_cvtps2pd(x);
  const v2df h  = one - A*v;
  const v2df h2 = h*h;
  //  v *= (one + h) * (one + h2);
  //  v *= (one + (one + h2)*(h + h2));
  v *= (one + h)*(one + h2)*(one + h2*h2);
#endif

  a = __builtin_ia32_vec_ext_v2df(v, 0);
  b = __builtin_ia32_vec_ext_v2df(v, 1);
}

inline v2df __rsqrtpd(const v2df r2)
{
#if 0
  const v2df flag  = __builtin_ia32_cmpgtpd(r2, (v2df){0.0, 0.0});
  return __builtin_ia32_andpd(flag, (v2df){1.0, 1.0}/__builtin_ia32_sqrtpd(r2));
#endif

#if 1
  const v4sf x   = __builtin_ia32_cvtpd2ps(r2);
  const v4sf ysp = __builtin_ia32_rsqrtps(x);
  v2df y1 = __builtin_ia32_cvtps2pd(ysp);
  const v2df c1 = {-0.5, -0.5};
  const v2df c2 = {-3.0, -3.0};
  y1 = (c1 * y1) * (r2*y1*y1 + c2);
  y1 = (c1 * y1) * (r2*y1*y1 + c2);
  const v2df flag = __builtin_ia32_cmpgtpd(r2, (v2df){0.0, 0.0});
  return __builtin_ia32_andpd(flag, y1);
#endif

#if 0
  const v4sf x = __builtin_ia32_movlhps(__builtin_ia32_cvtpd2ps(r2_1),__builtin_ia32_cvtpd2ps(r2_2));
  const v4sf ysp = __builtin_ia32_rsqrtps(x);
  const v2df y1 = __builtin_ia32_cvtps2pd(ysp);
  const v2df y2 = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(ysp, ysp));
  const v2df c1 = {1.0, 1.0};
  const v2df c2 = {0.5, 0.5};
  const v2df c3 = {0.375, 0.375};
  const v2df c4 = {0.3125, 0.3125};
  const v2df z1 = c1 - r2_1*y1*y1;
  const v2df z2 = c1 - r2_2*y2*y2;
  const v2df flag_1  = __builtin_ia32_cmpgtpd(r2_1, (v2df){0.0, 0.0});
  const v2df flag_2  = __builtin_ia32_cmpgtpd(r2_2, (v2df){0.0, 0.0});
  r2_1  = __builtin_ia32_andpd(flag_1, y1*(c1 + z1*(c2 + z1*(c3 + z1*c4))));
  r2_2  = __builtin_ia32_andpd(flag_2, y2*(c1 + z2*(c2 + z2*(c3 + z2*c4))));
#endif
}

inline double reduce(const v2df val)
{
  return __builtin_ia32_vec_ext_v2df(__builtin_ia32_haddpd(val, val), 0);
}

#endif /* _SSE */

#endif //  __PARTICLE_H__

