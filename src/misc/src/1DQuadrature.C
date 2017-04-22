//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <sstream>

#include <queso/1DQuadrature.h>

namespace QUESO {

//*****************************************************
// Base 1D quadrature class
//*****************************************************
Base1DQuadrature::Base1DQuadrature(
  double minDomainValue,
  double maxDomainValue,
  unsigned int order)
  : BaseQuadrature(),
  m_minDomainValue(minDomainValue),
  m_maxDomainValue(maxDomainValue),
  m_order         (order)
{
  queso_require_less_msg(m_minDomainValue, m_maxDomainValue, "min >= max");

  //                    UQ_UNAVAILABLE_RANK,
  //                    "Base1DQuadrature::constructor()",
  //                    "order = 0");
}

Base1DQuadrature::~Base1DQuadrature()
{
}

double
Base1DQuadrature::minDomainValue() const
{
  return m_minDomainValue;
}

double
Base1DQuadrature::maxDomainValue() const
{
  return m_maxDomainValue;
}

unsigned int
Base1DQuadrature::order() const
{
  return m_order;
}

//*****************************************************
// Generic 1D quadrature class
//*****************************************************
Generic1DQuadrature::Generic1DQuadrature(
  double minDomainValue,
  double maxDomainValue,
  const std::vector<double>& positions,
  const std::vector<double>& weights)
  :
  Base1DQuadrature(minDomainValue,maxDomainValue,positions.size()-1)
{
  m_positions = positions;
  m_weights   = weights;

  queso_require_not_equal_to_msg(m_positions.size(), 0, "invalid positions");

  queso_require_equal_to_msg(m_positions.size(), m_weights.size(), "inconsistent positions and weight");
}

Generic1DQuadrature::~Generic1DQuadrature()
{
}

//*****************************************************
// UniformLegendre 1D quadrature class
//*****************************************************
UniformLegendre1DQuadrature::UniformLegendre1DQuadrature(
  double       minDomainValue,
  double       maxDomainValue,
  unsigned int order,
  bool         densityIsNormalized)
  :
  Base1DQuadrature(minDomainValue,maxDomainValue,order)
{
  m_positions.resize(m_order+1,0.); // Yes, '+1'
  m_weights.resize  (m_order+1,0.); // Yes, '+1'

  // http://www.holoborodko.com/pavel/?page_id=679
  switch (m_order) { // eep2011
    case 1:
      m_weights  [0] =  1.;
      m_weights  [1] =  1.;

      m_positions[0] = -1./sqrt(3.);
      m_positions[1] =  1./sqrt(3.);
    break;

    case 2:
      m_weights  [0] =  5./9.;
      m_weights  [1] =  8./9.;
      m_weights  [2] =  5./9.;

      m_positions[0] = -sqrt(.6);
      m_positions[1] =  0.;
      m_positions[2] =  sqrt(.6);
    break;

    case 3:
      m_weights  [0] =  0.5 - sqrt(30.)/36.;
      m_weights  [1] =  0.5 + sqrt(30.)/36.;
      m_weights  [2] =  0.5 + sqrt(30.)/36.;
      m_weights  [3] =  0.5 - sqrt(30.)/36.;

      m_positions[0] = -sqrt(3.+2.*sqrt(1.2))/sqrt(7.);
      m_positions[1] = -sqrt(3.-2.*sqrt(1.2))/sqrt(7.);
      m_positions[2] =  sqrt(3.-2.*sqrt(1.2))/sqrt(7.);
      m_positions[3] =  sqrt(3.+2.*sqrt(1.2))/sqrt(7.);
    break;

    case 4:
      m_weights  [0] =  (322.-13.*sqrt(70.))/900.; // 0.236926885
      m_weights  [1] =  (322.+13.*sqrt(70.))/900.; // 0.478628670
      m_weights  [2] =  128./225.;                 // 0.568888889
      m_weights  [3] =  (322.+13.*sqrt(70.))/900.;
      m_weights  [4] =  (322.-13.*sqrt(70.))/900.;

      m_positions[0] = -sqrt(5.+2.*sqrt(10./7.))/3.;
      m_positions[1] = -sqrt(5.-2.*sqrt(10./7.))/3.;
      m_positions[2] =  0.;
      m_positions[3] =  sqrt(5.-2.*sqrt(10./7.))/3.;
      m_positions[4] =  sqrt(5.+2.*sqrt(10./7.))/3.;
    break;

    case 5:
      m_weights  [0] =  0.1713244923791703450402961;
      m_weights  [1] =  0.3607615730481386075698335;
      m_weights  [2] =  0.4679139345726910473898703;
      m_weights  [3] =  0.4679139345726910473898703;
      m_weights  [4] =  0.3607615730481386075698335;
      m_weights  [5] =  0.1713244923791703450402961;

      m_positions[0] = -0.9324695142031520278123016;
      m_positions[1] = -0.6612093864662645136613996;
      m_positions[2] = -0.2386191860831969086305017;
      m_positions[3] =  0.2386191860831969086305017;
      m_positions[4] =  0.6612093864662645136613996;
      m_positions[5] =  0.9324695142031520278123016;
    break;

    case 6:
      m_weights  [0] =  0.1294849661688696932706114;
      m_weights  [1] =  0.2797053914892766679014678;
      m_weights  [2] =  0.3818300505051189449503698;
      m_weights  [3] =  0.4179591836734693877551020;
      m_weights  [4] =  0.3818300505051189449503698;
      m_weights  [5] =  0.2797053914892766679014678;
      m_weights  [6] =  0.1294849661688696932706114;

      m_positions[0] = -0.9491079123427585245261897;
      m_positions[1] = -0.7415311855993944398638648;
      m_positions[2] = -0.4058451513773971669066064;
      m_positions[3] =  0.;
      m_positions[4] =  0.4058451513773971669066064;
      m_positions[5] =  0.7415311855993944398638648;
      m_positions[6] =  0.9491079123427585245261897;
    break;

    case 7:
      m_weights  [0] =  0.10122854;
      m_weights  [1] =  0.22238103;
      m_weights  [2] =  0.31370665;
      m_weights  [3] =  0.36268378;
      m_weights  [4] =  0.36268378;
      m_weights  [5] =  0.31370665;
      m_weights  [6] =  0.22238103;
      m_weights  [7] =  0.10122854;

      m_positions[0] = -0.96028986;
      m_positions[1] = -0.79666648;
      m_positions[2] = -0.52553241;
      m_positions[3] = -0.18343464;
      m_positions[4] =  0.18343464;
      m_positions[5] =  0.52553241;
      m_positions[6] =  0.79666648;
      m_positions[7] =  0.96028986;
    break;

    case 10:
      m_weights  [ 0] =  0.0556685671161736664827537;
      m_weights  [ 1] =  0.1255803694649046246346943;
      m_weights  [ 2] =  0.1862902109277342514260976;
      m_weights  [ 3] =  0.2331937645919904799185237;
      m_weights  [ 4] =  0.2628045445102466621806889;
      m_weights  [ 5] =  0.2729250867779006307144835;
      m_weights  [ 6] =  0.2628045445102466621806889;
      m_weights  [ 7] =  0.2331937645919904799185237;
      m_weights  [ 8] =  0.1862902109277342514260976;
      m_weights  [ 9] =  0.1255803694649046246346943;
      m_weights  [10] =  0.0556685671161736664827537;

      m_positions[ 0] = -0.9782286581460569928039380;
      m_positions[ 1] = -0.8870625997680952990751578;
      m_positions[ 2] = -0.7301520055740493240934163;
      m_positions[ 3] = -0.5190961292068118159257257;
      m_positions[ 4] = -0.2695431559523449723315320;
      m_positions[ 5] =  0.;
      m_positions[ 6] =  0.2695431559523449723315320;
      m_positions[ 7] =  0.5190961292068118159257257;
      m_positions[ 8] =  0.7301520055740493240934163;
      m_positions[ 9] =  0.8870625997680952990751578;
      m_positions[10] =  0.9782286581460569928039380;
    break;

    case 11:
      m_weights  [ 0] =  0.0471753363865118271946160;
      m_weights  [ 1] =  0.1069393259953184309602547;
      m_weights  [ 2] =  0.1600783285433462263346525;
      m_weights  [ 3] =  0.2031674267230659217490645;
      m_weights  [ 4] =  0.2334925365383548087608499;
      m_weights  [ 5] =  0.2491470458134027850005624;
      m_weights  [ 6] =  0.2491470458134027850005624;
      m_weights  [ 7] =  0.2334925365383548087608499;
      m_weights  [ 8] =  0.2031674267230659217490645;
      m_weights  [ 9] =  0.1600783285433462263346525;
      m_weights  [10] =  0.1069393259953184309602547;
      m_weights  [11] =  0.0471753363865118271946160;

      m_positions[ 0] = -0.9815606342467192506905491;
      m_positions[ 1] = -0.9041172563704748566784659;
      m_positions[ 2] = -0.7699026741943046870368938;
      m_positions[ 3] = -0.5873179542866174472967024;
      m_positions[ 4] = -0.3678314989981801937526915;
      m_positions[ 5] = -0.1252334085114689154724414;
      m_positions[ 6] =  0.1252334085114689154724414;
      m_positions[ 7] =  0.3678314989981801937526915;
      m_positions[ 8] =  0.5873179542866174472967024;
      m_positions[ 9] =  0.7699026741943046870368938;
      m_positions[10] =  0.9041172563704748566784659;
      m_positions[11] =  0.9815606342467192506905491;
    break;

    case 12:
      m_weights  [ 0] =  0.0404840047653158795200216;
      m_weights  [ 1] =  0.0921214998377284479144218;
      m_weights  [ 2] =  0.1388735102197872384636018;
      m_weights  [ 3] =  0.1781459807619457382800467;
      m_weights  [ 4] =  0.2078160475368885023125232;
      m_weights  [ 5] =  0.2262831802628972384120902;
      m_weights  [ 6] =  0.2325515532308739101945895;
      m_weights  [ 7] =  0.2262831802628972384120902;
      m_weights  [ 8] =  0.2078160475368885023125232;
      m_weights  [ 9] =  0.1781459807619457382800467;
      m_weights  [10] =  0.1388735102197872384636018;
      m_weights  [11] =  0.0921214998377284479144218;
      m_weights  [12] =  0.0404840047653158795200216;

      m_positions[ 0] = -0.9841830547185881494728294;
      m_positions[ 1] = -0.9175983992229779652065478;
      m_positions[ 2] = -0.8015780907333099127942065;
      m_positions[ 3] = -0.6423493394403402206439846;
      m_positions[ 4] = -0.4484927510364468528779129;
      m_positions[ 5] = -0.2304583159551347940655281;
      m_positions[ 6] =  0.;
      m_positions[ 7] =  0.2304583159551347940655281;
      m_positions[ 8] =  0.4484927510364468528779129;
      m_positions[ 9] =  0.6423493394403402206439846;
      m_positions[10] =  0.8015780907333099127942065;
      m_positions[11] =  0.9175983992229779652065478;
      m_positions[12] =  0.9841830547185881494728294;
    break;
#if 0
    case 13:
      m_weights  [ 0] =  ;
      m_weights  [ 1] =  ;
      m_weights  [ 2] =  ;
      m_weights  [ 3] =  ;
      m_weights  [ 4] =  ;
      m_weights  [ 5] =  ;
      m_weights  [ 6] =  ;
      m_weights  [ 7] =  ;
      m_weights  [ 8] =  ;
      m_weights  [ 9] =  ;
      m_weights  [10] =  ;
      m_weights  [11] =  ;
      m_weights  [12] =  ;
      m_weights  [13] =  ;

      m_positions[ 0] = -;
      m_positions[ 1] = -;
      m_positions[ 2] = -;
      m_positions[ 3] = -;
      m_positions[ 4] = -;
      m_positions[ 5] = -;
      m_positions[ 6] = -;
      m_positions[ 7] =  ;
      m_positions[ 8] =  ;
      m_positions[ 9] =  ;
      m_positions[10] =  ;
      m_positions[11] =  ;
      m_positions[12] =  ;
      m_positions[13] =  ;
    break;
#endif
#if 0
    14 ±0.1080549487073436620662447 0.2152638534631577901958764
    ±0.3191123689278897604356718 0.2051984637212956039659241
    ±0.5152486363581540919652907 0.1855383974779378137417166
    ±0.6872929048116854701480198 0.1572031671581935345696019
    ±0.8272013150697649931897947 0.1215185706879031846894148
    ±0.9284348836635735173363911 0.0801580871597602098056333
    ±0.9862838086968123388415973 0.0351194603317518630318329
#endif

    case 16:
      m_weights  [ 0] =  0.0241483028685479319601100;
      m_weights  [ 1] =  0.0554595293739872011294402;
      m_weights  [ 2] =  0.0850361483171791808835354;
      m_weights  [ 3] =  0.1118838471934039710947884;
      m_weights  [ 4] =  0.1351363684685254732863200;
      m_weights  [ 5] =  0.1540457610768102880814316;
      m_weights  [ 6] =  0.1680041021564500445099707;
      m_weights  [ 7] =  0.1765627053669926463252710;
      m_weights  [ 8] =  0.1794464703562065254582656;
      m_weights  [ 9] =  0.1765627053669926463252710;
      m_weights  [10] =  0.1680041021564500445099707;
      m_weights  [11] =  0.1540457610768102880814316;
      m_weights  [12] =  0.1351363684685254732863200;
      m_weights  [13] =  0.1118838471934039710947884;
      m_weights  [14] =  0.0850361483171791808835354;
      m_weights  [15] =  0.0554595293739872011294402;
      m_weights  [16] =  0.0241483028685479319601100;

      m_positions[ 0] = -0.9905754753144173356754340;
      m_positions[ 1] = -0.9506755217687677612227170;
      m_positions[ 2] = -0.8802391537269859021229557;
      m_positions[ 3] = -0.7815140038968014069252301;
      m_positions[ 4] = -0.6576711592166907658503022;
      m_positions[ 5] = -0.5126905370864769678862466;
      m_positions[ 6] = -0.3512317634538763152971855;
      m_positions[ 7] = -0.1784841814958478558506775;
      m_positions[ 8] =  0.;
      m_positions[ 9] =  0.1784841814958478558506775;
      m_positions[10] =  0.3512317634538763152971855;
      m_positions[11] =  0.5126905370864769678862466;
      m_positions[12] =  0.6576711592166907658503022;
      m_positions[13] =  0.7815140038968014069252301;
      m_positions[14] =  0.8802391537269859021229557;
      m_positions[15] =  0.9506755217687677612227170;
      m_positions[16] =  0.9905754753144173356754340;
    break;

    default:
      std::cerr << "In UniformLegendre1DQuadrature::constructor()"
                << ": m_order = " << m_order
                << std::endl;
      queso_error_msg("order not supported");
    break;
  }

  // Scale positions from the interval [-1, 1] to the interval [min,max]
  for (unsigned int j = 0; j < m_positions.size(); ++j) {
    m_positions[j] = .5*(m_maxDomainValue - m_minDomainValue)*m_positions[j] + .5*(m_maxDomainValue + m_minDomainValue);
    if (densityIsNormalized) {
      // Since \rho is "1/\Delta", we just multiply by ".5"
      m_weights[j] *= .5;
    }
    else {
      // Since \rho is "1", we multiply by ".5 * \Delta"
      m_weights[j] *= .5*(m_maxDomainValue - m_minDomainValue);
    }
  }
}

UniformLegendre1DQuadrature::~UniformLegendre1DQuadrature()
{
}

//*****************************************************
// GaussianHermite 1D quadrature class
//*****************************************************
GaussianHermite1DQuadrature::GaussianHermite1DQuadrature(
  double mean,
  double stddev,
  unsigned int order)
  :
  Base1DQuadrature(-INFINITY,INFINITY,order),
  m_mean  (mean),
  m_stddev(stddev)
{
  // FIX ME: prepare code for mean != 0 and stddev != 1
  m_positions.resize(m_order+1,0.); // Yes, '+1'
  m_weights.resize  (m_order+1,0.); // Yes, '+1'

  // http://www.efunda.com/math/num_integration/findgausshermite.cfm
  switch (m_order) { // eep2011
    case 1:
      m_weights  [0] =  sqrt(M_PI)/2.;
      m_weights  [1] =  sqrt(M_PI)/2.;

      m_positions[0] = -1./sqrt(2.);
      m_positions[1] =  1./sqrt(2.);
    break;

    case 2:
      m_weights  [0] =  sqrt(M_PI)/6.;
      m_weights  [1] =  2.*sqrt(M_PI)/3.;
      m_weights  [2] =  sqrt(M_PI)/6.;

      m_positions[0] = -sqrt(1.5);
      m_positions[1] =  0.;
      m_positions[2] =  sqrt(1.5);
    break;

    case 3:
      m_weights  [0] =  sqrt(M_PI)/4./(3.+sqrt(6.));
      m_weights  [1] =  sqrt(M_PI)/4./(3.-sqrt(6.));
      m_weights  [2] =  sqrt(M_PI)/4./(3.-sqrt(6.));
      m_weights  [3] =  sqrt(M_PI)/4./(3.+sqrt(6.));

      m_positions[0] = -sqrt(1.5+sqrt(1.5));
      m_positions[1] = -sqrt(1.5-sqrt(1.5));
      m_positions[2] =  sqrt(1.5-sqrt(1.5));
      m_positions[3] =  sqrt(1.5+sqrt(1.5));
    break;

    case 4:
      m_weights  [0] =  0.019953242049;
      m_weights  [1] =  0.393619323152;
      m_weights  [2] =  0.945308720483;
      m_weights  [3] =  0.393619323152;
      m_weights  [4] =  0.019953242059;

      m_positions[0] = -sqrt(2.5+sqrt(2.5));
      m_positions[1] = -sqrt(2.5-sqrt(2.5));
      m_positions[2] =  0.;
      m_positions[3] =  sqrt(2.5-sqrt(2.5));
      m_positions[4] =  sqrt(2.5+sqrt(2.5));
    break;

    case 5:
      m_weights   [0] = 0.00453000990551;
      m_weights   [1] = 0.157067320323;
      m_weights   [2] = 0.724629595224;
      m_weights   [3] = 0.724629595224;
      m_weights   [4] = 0.157067320323;
      m_weights   [5] = 0.00453000990551;

      m_positions [0] = -2.35060497367;
      m_positions [1] = -1.33584907401;
      m_positions [2] = -0.436077411928;
      m_positions [3] =  0.436077411928;
      m_positions [4] =  1.33584907401;
      m_positions [5] =  2.35060497367;
    break;

    case 6:
      m_weights   [0] = 0.0009717812451;
      m_weights   [1] = 0.0545155828191;
      m_weights   [2] = 0.42560725261;
      m_weights   [3] = 0.810264617557;
      m_weights   [4] = 0.42560725261;
      m_weights   [5] = 0.0545155828191;
      m_weights   [6] = 0.0009717812451;

      m_positions [0] = -2.65196135684;
      m_positions [1] = -1.67355162877;
      m_positions [2] = -0.816287882859;
      m_positions [3] =  0.;
      m_positions [4] =  0.816287882859;
      m_positions [5] =  1.67355162877;
      m_positions [6] =  2.65196135684;
    break;

    case 7:
      m_weights   [0] = 0.000199604072211;
      m_weights   [1] = 0.0170779830074;
      m_weights   [2] = 0.207802325815;
      m_weights   [3] = 0.661147012558;
      m_weights   [4] = 0.661147012558;
      m_weights   [5] = 0.207802325815;
      m_weights   [6] = 0.0170779830074;
      m_weights   [7] = 0.000199604072211;

      m_positions [0] = -2.93063742026;
      m_positions [1] = -1.9816567567;
      m_positions [2] = -1.15719371245;
      m_positions [3] = -0.381186990207;
      m_positions [4] =  0.381186990207;
      m_positions [5] =  1.15719371245;
      m_positions [6] =  1.9816567567;
      m_positions [7] =  2.93063742026;
    break;

    case 8:
      m_weights   [0] = 3.96069772633e-5;
      m_weights   [1] = 0.00494362427554;
      m_weights   [2] = 0.0884745273944;
      m_weights   [3] = 0.432651559003;
      m_weights   [4] = 0.720235215606;
      m_weights   [5] = 0.432651559003;
      m_weights   [6] = 0.0884745273944;
      m_weights   [7] = 0.00494362427554;
      m_weights   [8] = 3.96069772633e-5;

      m_positions [0] = -3.19099320178;
      m_positions [1] = -2.26658058453;
      m_positions [2] = -1.46855328922;
      m_positions [3] = -0.723551018753;
      m_positions [4] =  0.;
      m_positions [5] =  0.723551018753;
      m_positions [6] =  1.46855328922;
      m_positions [7] =  2.26658058453;
      m_positions [8] =  3.19099320178;
    break;

    case 9:
      m_weights   [0] = 7.64043285523e-6;
      m_weights   [1] = 0.00134364574678;
      m_weights   [2] = 0.0338743944555;
      m_weights   [3] = 0.240138611082;
      m_weights   [4] = 0.610862633735;
      m_weights   [5] = 0.610862633735;
      m_weights   [6] = 0.240138611082;
      m_weights   [7] = 0.0338743944555;
      m_weights   [8] = 0.00134364574678;
      m_weights   [9] = 7.64043285523e-6;

      m_positions [0] = -3.43615911884;
      m_positions [1] = -2.53273167423;
      m_positions [2] = -1.7566836493;
      m_positions [3] = -1.03661082979;
      m_positions [4] = -0.342901327224;
      m_positions [5] =  0.342901327224;
      m_positions [6] =  1.03661082979;
      m_positions [7] =  1.7566836493;
      m_positions [8] =  2.53273167423;
      m_positions [9] =  3.43615911884;
    break;

    case 19:
      m_weights   [0] =  2.22939364554e-13;
      m_weights   [1] =  4.39934099226e-10;
      m_weights   [2] =  1.08606937077e-7;
      m_weights   [3] =  7.8025564785e-6;
      m_weights   [4] =  0.000228338636017;
      m_weights   [5] =  0.00324377334224;
      m_weights   [6] =  0.0248105208875;
      m_weights   [7] =  0.10901720602;
      m_weights   [8] =  0.286675505363;
      m_weights   [9] =  0.462243669601;
      m_weights  [10] =  0.462243669601;
      m_weights  [11] =  0.286675505363;
      m_weights  [12] =  0.10901720602;
      m_weights  [13] =  0.0248105208875;
      m_weights  [14] =  0.00324377334224;
      m_weights  [15] =  0.000228338636017;
      m_weights  [16] =  7.8025564785e-6;
      m_weights  [17] =  1.08606937077e-7;
      m_weights  [18] =  4.39934099226e-10;
      m_weights  [19] =  2.22939364554e-13;

      m_positions [0] = -5.38748089001;
      m_positions [1] = -4.60368244955;
      m_positions [2] = -3.94476404012;
      m_positions [3] = -3.34785456738;
      m_positions [4] = -2.78880605843;
      m_positions [5] = -2.25497400209;
      m_positions [6] = -1.73853771212;
      m_positions [7] = -1.2340762154;
      m_positions [8] = -0.737473728545;
      m_positions [9] = -0.245340708301;
      m_positions[10] =  0.245340708301;
      m_positions[11] =  0.737473728545;
      m_positions[12] =  1.2340762154;
      m_positions[13] =  1.73853771212;
      m_positions[14] =  2.25497400209;
      m_positions[15] =  2.78880605843;
      m_positions[16] =  3.34785456738;
      m_positions[17] =  3.94476404012;
      m_positions[18] =  4.60368244955;
      m_positions[19] =  5.38748089001;
    break;

    default:
      std::stringstream ss;
      ss << "In GaussianHermite1DQuadrature::constructor()"
         << ": m_order = " << m_order
         << std::endl;
      queso_error_msg(ss.str());
    break;
  }
}

GaussianHermite1DQuadrature::~GaussianHermite1DQuadrature()
{
}

//*****************************************************
// WignerInverseChebyshev1st 1D quadrature class
//*****************************************************
WignerInverseChebyshev1st1DQuadrature::WignerInverseChebyshev1st1DQuadrature(
  double       minDomainValue,
  double       maxDomainValue,
  unsigned int order)
  :
  Base1DQuadrature(minDomainValue,maxDomainValue,order)
{
  m_positions.resize(m_order+1,0.); // Yes, '+1'
  m_weights.resize  (m_order+1,0.); // Yes, '+1'

  // http://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature
  switch (m_order) {
    default:
      queso_error_msg("order not supported");
    break;
  }

  // Scale positions from the interval [-1, 1] to the interval [min,max]
  for (unsigned int j = 0; j < m_positions.size(); ++j) {
    m_positions[j] = .5*(m_maxDomainValue - m_minDomainValue)*m_positions[j] + .5*(m_maxDomainValue + m_minDomainValue);
    m_weights[j] *= .5*(m_maxDomainValue - m_minDomainValue);
  }
}

WignerInverseChebyshev1st1DQuadrature::~WignerInverseChebyshev1st1DQuadrature()
{
}

//*****************************************************
// WignerChebyshev2nd 1D quadrature class
//*****************************************************
WignerChebyshev2nd1DQuadrature::WignerChebyshev2nd1DQuadrature(
  double       minDomainValue,
  double       maxDomainValue,
  unsigned int order)
  :
  Base1DQuadrature(minDomainValue,maxDomainValue,order)
{
  m_positions.resize(m_order+1,0.); // Yes, '+1'
  m_weights.resize  (m_order+1,0.); // Yes, '+1'

  // http://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature
  unsigned int n = m_order+1;
  for (unsigned int i = 0; i < n; ++i) {
    double angle = M_PI*((double)(i+1))/((double)(n+1));
    double cosValue = cos(angle);
    double sinValue = sin(angle);
    m_positions[i] = cosValue;
    m_weights[i] = ( M_PI/((double)(n+1)) )*sinValue*sinValue;
  }

  // Scale positions from the interval [-1, 1] to the interval [min,max]
  for (unsigned int j = 0; j < m_positions.size(); ++j) {
    m_positions[j] = .5*(m_maxDomainValue - m_minDomainValue)*m_positions[j] + .5*(m_maxDomainValue + m_minDomainValue);
    m_weights[j] *= .5*(m_maxDomainValue - m_minDomainValue);
  }
}

WignerChebyshev2nd1DQuadrature::~WignerChebyshev2nd1DQuadrature()
{
}

}  // End namespace QUESO
