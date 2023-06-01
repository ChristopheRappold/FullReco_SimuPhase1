#ifndef CONSTANT_PARAMETER_HH
#define CONSTANT_PARAMETER_HH

#include <math.h>

namespace hitana {

const double mm2cm = 0.1;
const double cm2mm = 10.;

const double Deg2Rad = acos(-1)/180;
const double Rad2Deg = 1./Deg2Rad;

const double c = 299792458; // m/s
const double MassProton = 0.938272; // GeV
const double MassKaon   = 0.493677; // GeV
const double MassPion   = 0.139570; // GeV

}

#endif
