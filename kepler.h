#ifndef __KEPLER_H__
#define __KEPLER_H__

#include <iostream>
#include <cassert>
#include "vector3.h"

#if 0
struct Kepler
{
  vec3 pos, vel;
  
  Kepler(const vec3 &_pos, const vec3 &_vel, const double Mo, const double dt = 0.0) : 
    pos(_pos), vel(_vel)
  {
    if (dt == 0.0) return;
    
    assert(step(Mo, dt));
  }

  bool step(const double Mo, const double dt)
  {
    const double    R = pos.abs();
    assert(R > 0.0);
    const double invR = 1.0/R;
    const double invM = 1.0/Mo;
    const vec3 Lvec = pos%vel;
    const vec3 Evec = (vel%Lvec)*invM - pos*invR; 
    const vec3 Bvec = Lvec%Evec;

    const double hs   = (vel*vel)*0.5 - Mo*invR;
    const double k    = -2.0*hs*invM;
    const double j    = std::sqrt(2.0*std::abs(hs));
    const double n    = j*j*j*invM;

    const double e = Evec.abs();
    if (e == 0.0)
    {
      assert(hs < 0.0);
      assert(false);
    }
    else if (e < 1.0)
    {
      const double inve = 1.0/e;
      assert(j > 0.0);
      const double invj = 1.0/j;
      assert(k > 0.0);
      const double invk = 1.0/k;

      const double x = inve * (1.0 - n*R*invj);
      const double y = n*invk*invM*inve * (pos*vel);
      const double E0 = std::atan2(y,x);
      const double nP = -(E0 - e*sin(E0));

      double M = n*dt - nP;
      while (M < 0.0)     M += 2*M_PI;
      while (M >= 2*M_PI) M -= 2*M_PI;

      const double E = EccentricAnomaly(e, M);
      const double cv = std::sin(E);
      const double cr = std::cos(E);

     
#if 0 
      assert((pos%vel).abs() > 0.0);
      
      fprintf(stderr, "pos= %g %g %g  vel= %g %g %g \n",
          pos.x, pos.y, pos.z,
          vel.x, vel.y, vel.z);
      
      pos = invk*(cr*inve - 1.0)*Evec + cv*inve*invj*Bvec;
      
      const double invr = 1.0/pos.abs();
      
      vel = inve*invr*(-Mo*cv*invj*Evec + cr*Bvec); 

      fprintf(stderr, "M= %g cr= %g  cv= %g  pos= %g %g %g  vel= %g %g %g \n",
          M,
          cr, cv,
          pos.x, pos.y, pos.z,
          vel.x, vel.y, vel.z);
      {
        const vec3 cr = Evec%Bvec;
        fprintf(stderr," Evec= %g %g %g \n", Evec.x, Evec.y, Evec.z);
        fprintf(stderr," Bvec= %g %g %g \n", Bvec.x, Bvec.y, Bvec.z);
        fprintf(stderr," cr= %g %g %g \n", cr.x, cr.y, cr.z);
      }
      assert((Evec%Bvec).abs() > 0.0);
      assert((pos%vel).abs() > 0.0);
#else
      const double invr = k/(1.0 - e*cr);
      pos = invk*(cr*inve - 1.0)*Evec + cv*inve*invj*Bvec;
      vel = inve*invr*(-Mo*cv*invj*Evec + cr*Bvec); 
#endif
    }
    else if (e == 1.0)
    {
      assert(0);
    }
    else  /* e > 1.0 */
    {
      assert(0);
    }

    return true;
  }

  static double sign(const double x) {return x < 0.0 ? -1.0 : (x > 0.0 ? 1.0 : 0.0);}

  static double EccentricAnomaly(const double e, const double M1)
  {

    const double M = M1 > M_PI ? 2*M_PI - M1 : M1;
    double E;

#if 1
    {
      int niter = 100;
      const double tol = 1.0e-10;
      double Ei = M + 0.85*e*sign(std::sin(M));
      E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));

#if 0

      while (std::abs(E - Ei) > tol*std::abs((E + Ei)*0.5))
      {
        if (niter < 10)
        {
          fprintf(stderr, "Ei= %g  E= %g  diff= %g \n",
              Ei, E, Ei-E);
        }
        assert(niter >= 0);
        Ei = E;
        E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));
        niter--;
      }
#else

      while (std::abs(E - Ei) > tol*E)
      {
        if (niter < 10)
        {
          fprintf(stderr, "Ei= %g  E= %g  diff= %g \n",
              Ei, E, (Ei-E)/E);
        }
        assert(niter >= 0);
        const double sinE = std::sin(E);
        const double cosE = std::cos(E);
        const double f3 = e*cosE;
        const double f2 = e*sinE;
        const double f1 = 1.0 - f3;
        const double f0 = E - e*sinE - M;
        const double d1 = -f0/ f1;
        const double d1h = 0.5*d1;
        const double d2 = -f0/(f1 + d1h* f2);
        const double d3 = -f0/(f1 + d1h*(f2 + (1.0/3.0)*d2*d2*f3));

        Ei = E;
        E  = Ei + d3;

        niter--;
      }
#endif
    }
#endif

#if 0
    {
      const int N = 23;
      static std::vector<real> coeff;
      static bool init = true;
      if (init)
      {
        init = false;

        unsigned int K = 2;
        real u = 0.5; 
        coeff.push_back(0.5);
        for (int n = 1; n <= N; n++)
        {
          K <<= 1;
          u  *= 0.5;
          for (unsigned int k = 1; k <= K; k += 2)
            coeff.push_back(k*u);
        }
      }

      size_t j = 1;
      for (int n = 2; n <= N; n++)
      {
        E = M_PI*coeff[j-1];
        const real g = E - e*std::sin(E) - M;
        j <<= 1;
        if (g < 0.0) 
          j++;
      }

    }
#endif

    if (M1 > M_PI) E = 2*M_PI - E;

    return E;
  }

};

#else


struct Kepler
{
  vec3 pos, vel;
#define TOLERANCE (1.0e-13)

  double G_func2(double q) 
  {
    double l = 3.0;
    double d = 15.0;
    double n = 0.0;
    double A, B, G;

    if(q==0.0) return 1.0;	/* this isn't necessary when first
                               Newt-Raph iteration is done by hand */

    A = B = G = 1.0;

    while (fabs(B) > TOLERANCE*fabs(G)) 
    {
      l += 2.0;
      d += 4.0*l;
      n += 10.0*l;

      A = d/(d-n*A*q);
      B *= A-1.0;
      G += B;

      l += 2.0;
      d += 4.0*l;
      n -= 8.0*l;

      A = d/(d-n*A*q);
      B *= A-1.0;
      G += B;
    }

    return G;
  };

  Kepler(const vec3 &_pos, const vec3 &_vel, const double Mo, const double dt = 0.0) : 
    pos(_pos), vel(_vel)
  {
    if (dt == 0.0) return;

    if (!step(Mo, dt))
    {
      assert(step(Mo, dt*0.5));
      assert(step(Mo, dt*0.5));
    }
  }

  bool step(const double Mo, const double dt)
  {
    double r0mag, v0mag2;
    double r0v0;	/* r dot v */
    double rcalc, dtcalc, terr;
    double u;	/* (?) universal variable */
    double beta;	/* (?) vis-a-vis integral */
    double P;	/* period (for elliptic orbits only) */
    double dU;
    int n;
    double q;
    double U0w2, U1w2;
    double U, U0, U1, U2, U3;
    double f, g, F, G;
    int no_iter;

    double du1, du2, du3, dqdu, d2qdu2, drdu, d2rdu2, fn, fnp, fnpp, fnppp;

    r0mag  = pos.abs();   // sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
    const double invr = 1.0/r0mag;
    v0mag2 = vel.norm2(); // v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
    r0v0   = pos*vel;     // r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2];
    beta   = 2*Mo*invr - v0mag2;

    if (beta > 0) 
    {	
      const double invb  = 1.0/beta;
      const double invbh = std::sqrt(invb); 
      /* elliptic orbit */
      P   = 2*M_PI*Mo*invbh*invb;
      n = floor((dt + P*0.5 + (-2.0)*r0v0*invb)/P);
      dU = 2*n*M_PI*invb*invb*invbh;
    } 
    else 
    {
      dU = 0.0;
    }

    u = 0;	/* a "better" guess is possible, see footnote at Battin p.219 */
    //	u = dt/(4*r0mag); 				  /* N-R step by hand */
    //	u = dt/(4*r0mag-2*dt*r0v0/r0mag);
    //	u = init_u(dt, mu, r0mag, r0v0, beta);

    no_iter = 0;
    do 
    {
      const real beta_uu = beta*u*u;
      const real invbuup = 1.0/(1.0 + beta_uu);
      q = beta_uu*invbuup; //beta*u*u/(1+beta*u*u);
      if (q > 0.5 || no_iter > 12) 
      {
        return false;
        assert(0);
      }
      //        return DRIFT_FAIL;
#if 0
      dqdu = 2*beta*u/(1+beta*u*u)/(1+beta*u*u);
      d2qdu2 = 2*beta/(1+beta*u*u) 
        - 8*beta*beta*u*u / (1+beta*u*u)/(1+beta*u*u)/(1+beta*u*u);
#else
      d2qdu2 = 2.0*beta*invbuup;
      dqdu   = d2qdu2*u*invbuup;
      d2qdu2 += (-4.0)*dqdu*beta*u * invbuup*invbuup;
#endif
      U0w2 = 1.0 + (-2.0)*q;
      U1w2 = 2.0*(1.0 - q)*u;
      const double U1w2_2 = U1w2*U1w2;
      const double U1w2_4 = U1w2_2*U1w2_2;
      U = (16.0/15.0) * U1w2_4*U1w2 * G_func2(q) + dU;
      U0 = 2.0*U0w2*U0w2 - 1.0;
      U1 = 2.0*U0w2*U1w2;
      U2 = 2.0*U1w2*U1w2;
      U3 = beta*U + U1*U2*(1.0/3.0);
      rcalc = r0mag*U0 + r0v0*U1 + Mo*U2;
#if 0
      drdu   = 4.0*(1.0-q)*(r0v0*U0 + (Mo-beta*r0mag)*U1);
      d2rdu2 = (-4.0)*dqdu*(r0v0*U0 + (Mo-beta*r0mag)*U1)
        + (4.0*(1.0-q)*4.0*(1.0-q))*(-beta*r0v0*U1 + (Mo-beta*r0mag)*U0);
      dtcalc = r0mag*U1 + r0v0*U2 + Mo*U3;

      fn    = dtcalc-dt;
      fnp   = 4.0*(1.0-q)*rcalc;
      fnpp  = 4.0*(drdu*(1.0-q) - rcalc*dqdu);
      fnppp = -8.0*drdu*dqdu - 4.0*rcalc*d2qdu2 + 4.0*(1.0-q)*d2rdu2;
#else
      const double C1 = Mo - beta*r0mag;
      const double C2 = r0v0*U0; 
      const double C3 = 4.0*(1.0-q);
      drdu   = C3*(C2 + C1*U1);
      d2rdu2 = (-4.0)*dqdu*(C2 + C1*U1)
        + C3*C3*(-beta*r0v0*U1 + C1*U0);
      dtcalc = r0mag*U1 + r0v0*U2 + Mo*U3;

      fn    = dtcalc-dt;
      fnp   = C3*rcalc;
      fnpp  = 4.0*(drdu*(1.0-q) - rcalc*dqdu);
      fnppp = -8.0*drdu*dqdu - 4.0*rcalc*d2qdu2 + C3*d2rdu2;
#endif

      du1  = -fn/ fnp;
      du2  = -fn/(fnp + du1*fnpp*0.5);
      du3  = -fn/(fnp + du2*fnpp*0.5 + du2*du2*fnppp*(1.0/6.0));

      u += du3;
      no_iter++;

      terr = fabs(dt-dtcalc);
    } 
    while (terr > TOLERANCE*fabs(dt));

    f = 1 - (Mo*invr)*U2;
    g = r0mag*U1 + r0v0*U2;
    const double invrc = 1.0/rcalc;
    F = -Mo*U1*invrc*invr;//(rcalc*r0mag);
    G = 1 - (Mo*invrc)*U2;

    const vec3 r1 = f*pos + g*vel;
    const vec3 v1 = F*pos + G*vel;

    pos = r1;
    vel = v1;
    return true;
  }


  static double EccentricAnomaly(const double e, const double M)
  {
    const int niter_max = 100;

    double Ei = M/2.0;
    double E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));

    const double tol = 1.0e-13;
    int iter = 0;
    while (std::abs(E - Ei) > tol*std::abs((E + Ei)*0.5))
    {
      Ei = E;
      E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));
      assert(iter < niter_max);
      iter++;
    }

    return E;
  }

#undef TOLERANCE 

};
#endif


#endif // __KEPLER_H__
