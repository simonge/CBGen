//
// TPhoton.cxx
//
// author:  Richard T. Jones  11/16/98
// version:  Dec. 12, 1998  v1.00
//
/*************************************************************************
 * Copyright(c) 1998, University of Connecticut, All rights reserved.    *
 * Author: Richard T. Jones, Asst. Prof. of Physics                      *
 *                                                                       *
 * Permission to use, copy, modify and distribute this software and its  *
 * documentation for non-commercial purposes is hereby granted without   *
 * fee, provided that the above copyright notice appears in all copies   *
 * and that both the copyright notice and this permission notice appear  *
 * in the supporting documentation. The author makes no claims about the *
 * suitability of this software for any purpose.                         *
 * It is provided "as is" without express or implied warranty.           *
 *************************************************************************/
//////////////////////////////////////////////////////////////////////////
//
// The TPhoton class represents a photon by its momentum and state of
// polarization.  The photon is not required to be on-shell; however
// only transverse polarization is represented.  Polarization is held
// as a Pauli matrix to allow for any degree of polarization, or none.
// The axis of quantization is the direction of the momentum vector.
// The TGhoston class is provided for testing gauge invariance in a
// code calculating observable quantities.
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "TPhoton.h"

ClassImp(TPhoton)

TThreeVectorReal TPhoton::Pol() const
{
   // The photon polarization is conventionally represented by the
   // trajectory of the tip of the electric field vector on the plane
   // perpendicular to the direction of propagation.  Two numbers are
   // required in that description: azimuthal angle phi of the major
   // axis of the ellipse, and eccentricity e of the ellipse, with a
   // sign indicating the direction of rotation around the ellipse.
   // For the more general case of mixed states an additional number
   // p between 0 and 1 is needed to describe the degree of polarization.
   // Here I pack all three numbers into a single 3-vector with:
   //   1) azimuthal angle is phi (gives plane of major axis);
   //   2) polar angle encodes eccentricity (r.h. circular is theta=0,
   //      l.h. circular is theta=pi, planar case is theta=pi/2);
   //   3) length is p.

   TThreeVectorReal pol;
   Double_t cons;
   fSpinDensity.Decompose(cons,pol);
   Double_t r(0),theta(0),phi(0);
   pol.GetPolar(r,theta,phi);
   r /= cons; phi /= 2;
   pol.SetPolar(r,theta,phi);
   return pol;
}

TPhoton &TPhoton::SetPol(const TThreeVectorReal &polar)
{
   // Photon polarization information must be packed into a single
   // 3-vector polar as follows:
   //   1) azimuthal angle is phi (gives plane of major axis);
   //   2) polar angle encodes eccentricity (r.h. circular is theta=0,
   //      l.h. circular is theta=pi, planar case is theta=pi/2);
   //   3) length of polar is degree of polarization between 0 and 1.

   TThreeVectorReal pol;
   Double_t r(0),theta(0),phi(0);
   polar.GetPolar(r,theta,phi);
   phi *= 2;
   pol.SetPolar(r,theta,phi);
   fSpinDensity.SetDensity(pol);
   return *this;
}

TFourVectorComplex TPhoton::Eps(const Int_t mode) const
{
   // Vector potential 4-vector for photon object is returned.  There
   // are four linearly-independent solutions (modes) which are:
   //    mode=1: right-circular polarization
   //    mode=2: left-circular polarization
   //    mode=3: longitudinal polarization
   //    mode=4: momentum vector (renormalized)
   // Modes 1-3 obey the Lorentz gauge condition k.Eps(mode)=0; mode 4 is
   // provided to complement the set.  These four modes are an orthonormal
   // set in the Minkowski sense (+--- metric), that is
   //    Eps(m1).Conj() . Eps(m2) = (m1==m2)*eta[m1]
   // where eta[1]=eta[2]=-1 and eta[3]=-eta[4].  The transverse modes 1-2
   // always have an invariant norm of -1, whereas the norm of modes 3 and
   // 4 are determined by whether the photon is time-like (eta[4]=+1) or
   // space-like (eta[4]=-1).  Both are zero if the photon is massless,
   // which means that modes 3 and 4 are only meaningful in the case of
   // virtual photons.  If the photon momentum has zero norm then these two
   // modes become generate and cannot be normalized.  The behaviour in this
   // case is to return the null vector for these modes.  Note that the
   // normalization defined above corresponds to the usual relativistic
   // convention with the particle density of a unit-amplitude plane wave
   // given by 2E.

   TThreeVectorComplex xhat,yhat,zhat(0,0,1);
   yhat.Cross(zhat,fMomentum);
   if (yhat.Length() > yhat.Resolution()) {
      xhat.Cross(yhat,fMomentum);
      xhat.Normalize(1);
      yhat.Normalize(1);
   }
   else {
      xhat[1]=1; xhat[2]=0; xhat[3]=0;
      yhat[1]=0; yhat[2]=1; yhat[3]=0;
   }

   const Double_t sqrt2=sqrt(2.);
   const Complex_t i_(0,1);
   TFourVectorComplex eps(fMomentum);
   switch (mode) {
   case 1:    // right-circular polarization (helicity frame)
      eps[0] = 0;
      (TThreeVectorComplex &)eps = xhat + i_*yhat;
      eps /= sqrt2;
      break;
   case 2:    // left-circular polarization (helicity frame)
      eps[0] = 0;
      (TThreeVectorComplex &)eps = xhat - i_*yhat;
      eps /= sqrt2;
      break;
   case 3:    // longitudinal mode (virtual photons only)
      eps[0] = fMomentum.Length();
      eps.Normalize(fMomentum[0]);
   case 4:    // ghost mode (absent in Lorentz gauge)
      Double_t norm=eps.InvariantSqr();
      if (norm > 0)
         eps /= sqrt(norm);
      else if (norm < 0)
         eps /= sqrt(-norm);
      else
         eps = TFourVectorComplex(0,0,0,0);
      break;
   }
   return eps;
}

TPhoton &TPhoton::SetPlanePolarization(const TThreeVectorReal &phat,
                                       Double_t pol)
{
   // This method is provided for setting a photon to a state of 
   // linear polarization along a particular physical plane, without
   // needing to remember what the Stokes parameters mean, or how
   // the somewhat confusing basis vectors in Eps() are defined for
   // the helicity basis states. Just set phat to the direction of
   // the polarization vector orthogonal to the photon momentum and
   // pol to the degree of linear polarization, and you are all set.
   // Remember that if you change the momentum, you need to come back
   // and invoke this method again to keep the polarization consistent.

   TThreeVectorComplex eps(Eps(1));
   TThreeVectorReal epsx(eps.RealPart());
   TThreeVectorReal epsy(eps.ImagPart());
   TThreeVectorReal s(epsx.Dot(phat), epsy.Dot(phat), 0);
   s *= pol / s.Length();
   return SetPol(s);
}

TPhoton &TPhoton::SetEllipticalPolarization(const TThreeVectorReal &phat,
                                            Double_t circ,
                                            Double_t pol)
{
   // This method is provided for setting a photon to a state of 
   // linear polarization along a particular physical plane, without
   // needing to remember what the Stokes parameters mean, or how
   // the somewhat confusing basis vectors in Eps() are defined for
   // the helicity basis states. Just set phat to the direction of
   // the major axis of the ellipse traced out by the tip of the
   // electric field vector orthogonal to the photon momentum, circ
   // to the circularity of the polarization ellipse, and pol to the
   // degree of polarization, and you are all set. The value of circ
   // should be +1 for right-circular, -1 for left-circular, 0 for
   // linear, and in between the cos(theta_Stokes) for the general case.
   // Remember that if you change the momentum, you need to come back
   // and invoke this method again to keep the polarization consistent.

   TThreeVectorComplex eps(Eps(1));
   TThreeVectorReal epsx = eps.RealPart();
   TThreeVectorReal epsy = eps.ImagPart();
   TThreeVectorReal s(epsx.Dot(phat), epsy.Dot(phat), 0);
   if (fabs(circ) < 1) {
      s *= sqrt(1 - circ*circ) / s.Length();
      s += TThreeVectorReal(0,0,circ);
      s *= pol;
   }
   else if (circ > 0) {
      s = TThreeVectorReal(0,0,1);
   }
   else {
      s = TThreeVectorReal(0,0,-1);
   }
   return SetPol(s);
}

TThreeVectorReal TPhoton::GetPolarizationPlane() const
{
   // This method is provided for checking the polarization plane
   // of a photon, without needing to remember what the Stokes
   // parameters mean, or how the somewhat confusing basis vectors
   // in Eps() are defined for the helicity basis states, so as to
   // be able to interpret what is returned by the Pol() method.

   TThreeVectorComplex eps(Eps(1));
   TThreeVectorReal epsx = eps.RealPart();
   TThreeVectorReal epsy = eps.ImagPart();
   TThreeVectorReal pol(Pol());
   return sqrt(2) * (pol[1] * epsx + pol[2] * epsy);
}

void TPhoton::Streamer(TBuffer &buf)
{
   // Put/get a TPhoton object to/from stream buffer buf.

   TPhoton *me=this;
   if (buf.IsReading()) {
      buf >> me;
   } else {
      buf << me;
   }
}

void TPhoton::Print(Option_t *option)
{
   // Output a Photon class member in ascii form.

   std::cout << "TPhoton is" << std::endl;
   std::cout << "Momentum: (" << fMomentum[0] << ", " << fMomentum[1] << ", "
        << fMomentum[2] << ", " << fMomentum[3] << " )" << std::endl;
   std::cout << "Helicity density matrix:" << std::endl;
   std::cout << fSpinDensity[0][0] << "  " << fSpinDensity[0][1] << std::endl;
   std::cout << fSpinDensity[1][0] << "  " << fSpinDensity[1][1] << std::endl;
}


ClassImp(TGhoston)

TFourVectorComplex TGhoston::Eps(const Int_t mode) const
{
   // This member function is what distinguishes a photon from a ghoston.
   // Instead of returning a vector-potential four-vector for the mode,
   // it returns the momentum four-vector.  This means that declaring
   // any external photon for a process as TGhoston instead of TPhoton 
   // should make the rate for the process vanish by gauge invariance.
   // The mode argument is ignored.

   return fMomentum;
}

void TGhoston::Print(Option_t *option)
{
   // Output a Ghoston class member in ascii form.

   std::cout << "TGhoston is" << std::endl;
   std::cout << "Momentum: (" << fMomentum[0] << ", " << fMomentum[1] << ", "
        << fMomentum[2] << ", " << fMomentum[3] << " )" << std::endl;
   std::cout << "Helicity density matrix:" << std::endl;
   std::cout << fSpinDensity[0][0] << "  " << fSpinDensity[0][1] << std::endl;
   std::cout << fSpinDensity[1][0] << "  " << fSpinDensity[1][1] << std::endl;
}
