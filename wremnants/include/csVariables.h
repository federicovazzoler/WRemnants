#ifndef WREMNANTS_CSVARIABLES_H
#define WREMNANTS_CSVARIABLES_H


#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>
#include <iostream>

namespace wrem {

typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;
typedef ROOT::Math::PxPyPzMVector PxPyPzMVector;
typedef ROOT::Math::PtEtaPhiMVector PtEtaPhiMVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Vector;

template <typename T>
double dot(const T& vec1, const T& vec2) {
    return vec1.x()*vec2.x()+vec1.y()*vec2.y()+vec1.z()*vec2.z();
}

template <typename T>
T cross(const T& vec1, const T& vec2) {
    auto cross = ROOT::Math::Cross<double>({vec1.x(), vec1.y(), vec1.z()}, {vec2.x(), vec2.y(), vec2.z()});
    T res(cross[0], cross[1], cross[2]);
    return res;
}

Vector unitBoostedVector(ROOT::Math::Boost& boostOp, PxPyPzEVector& vec) {
    PxPyPzEVector boostvec = boostOp(vec);
    return Vector(boostvec.x(), boostvec.y(), boostvec.z()).Unit();
}

struct CSVars {

  double sintheta;
  double costheta;
  double sinphi;
  double cosphi;

};

CSVars CalccsSineCosThetaPhi(const PtEtaPhiMVector& antilepton, const PtEtaPhiMVector& lepton) {
	// following https://arxiv.org/pdf/2006.11382.pdf pag. 18 eq. 2.56 and 2.57
  	PxPyPzEVector p1(lepton);
  	PxPyPzEVector p2(antilepton);
  	PxPyPzEVector Q = p1+p2;

	// arbitrarily chosen the positive orientation of the z axis by having hadron a move in the z direction in the lab frame. As a result, the negatively charged lepton moves into the same rest-frame hemisphere as hadron a for costheta > 0
  	const int zsign = std::copysign(1.0, Q.z());
  	const double energy = 6500.;
  	PxPyPzEVector pa(0., 0., zsign*energy, energy);
  	PxPyPzEVector pb(0., 0., -1.*zsign*energy, energy);

	double costheta = 1. / (Q.M() * Q.Mt()) * ((p1.E() + p1.z()) * (p2.E() - p2.z()) - (p1.E() - p1.z()) * (p2.E() + p2.z()));
	double sintheta = std::sqrt(1 - costheta * costheta);
	double cosphi = 1. / sintheta * (p1.Pt() * p1.Pt() - p2.Pt() * p2.Pt()) / (Q.Pt() * Q.Mt());
	double sinphi = 2. / sintheta * (p1.y() * p2.x() - p2.y() * p1.x()) / (Q.Pt() * Q.M());

#if 1
	costheta = costheta * zsign;
	double Phi = std::atan2(sinphi, cosphi) * zsign;
	cosphi = std::cos(Phi);
	sinphi = std::sin(Phi);
#endif

	CSVars angles = {sintheta, costheta, sinphi, cosphi};
    return angles;
}

}

#endif
