/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file:
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <example_hyst.h>

struct cellStruct {
  double xu;
  double yu;
  double k;
  double r;
  double inter_y;
};

//------------------------------------------------------
// Calculates the restoring force, based on
// function [fs,rigidez,newparam] = restoringforce(u,uold,v,vold,fsold,param)
// Elastic plastic hysteresis model
//------------------------------------------------------
void restoringForce(
  double      u,
  double      uold,
  double      v,
  double      vold,
  double      fsold,
  cellStruct* inpParam,
  cellStruct* fsCell,   // output
  double*     fsDouble, // output
  double*     rigidez,  // output
  cellStruct* newParam) // output
{
  if (fsCell != NULL) {
    // Initialize hysteresis loop
    fsCell->xu      = 0;
    fsCell->yu      = 0;
    fsCell->k       = u;
    fsCell->r       = uold;
    fsCell->inter_y = v;
    return;
  }

  double xu      = inpParam->xu;
  double yu      = inpParam->yu;
  double k       = inpParam->k;
  double r       = inpParam->r;
  double inter_y = inpParam->inter_y;

  double k_r = k+r;

  *fsDouble = yu + k_r*(u-xu);

  *rigidez = k_r;

  if (v>0.) {
    double fyp = r*u + inter_y;
    if (*fsDouble > fyp) {
      *fsDouble = fyp;
      xu        = (inter_y + k_r*xu -yu)/k; // 1
      yu        = r*xu + inter_y;
      *rigidez  = r;
    }
    if (vold<0.) { //% && (v>0)
      if (fabs(fsold-r*uold+inter_y) < 1.e-3) {
        xu = uold;
        yu = fsold;
      }
      else {
        xu = (-inter_y + k_r*u - *fsDouble)/k; // 4
        yu = r*xu - inter_y;
      }
      *fsDouble = yu + k_r*(u-xu);
    }
  }
  else { //% (v<0)
    double fyn = r*u - inter_y;
    if (*fsDouble < fyn) {
      *fsDouble = fyn;
      xu        = (-inter_y + k_r*xu - yu)/k; // 3
      yu        = r*xu - inter_y;
      *rigidez  = r;
    }
    if (vold>0.) { //% && (v<0)
      if (fabs(fsold-r*uold-inter_y) < 1.e-3) {
        xu = uold;
        yu = fsold;
      }
      else {
        xu = (inter_y + k_r*u - *fsDouble)/k; // 2
        yu = r*xu + inter_y;
      }
      *fsDouble = yu + k_r*(u-xu);
    }
  }

  newParam->xu      = xu;
  newParam->yu      = yu;
  newParam->k       = k;
  newParam->r       = r;
  newParam->inter_y = inter_y;

  return;
}

//------------------------------------------------------
// Assembles: stiffness matrix K or dumping matrix C
//------------------------------------------------------
void ckmatrix(const QUESO::GslVector& vec, QUESO::GslMatrix& mat)
{
  UQ_FATAL_TEST_MACRO(vec.sizeLocal() != mat.numRowsLocal(),
                      vec.env().fullRank(),
                      "ckmatrix()",
                      "incompatible sizes");

  // Assemble C or K matrix
  unsigned int ndof = vec.sizeLocal();
  //C = spdiags([[-c(2:end) 0]; c+[c(2:end) 0]; -c]',-1:1,ndof,ndof);

  if (ndof == 4) {
    mat(0,0) =  vec[0]+vec[1]; mat(0,1) = -vec[1];        mat(0,2) = 0.;             mat(0,3) = 0.;
    mat(1,0) = -vec[1];        mat(1,1) =  vec[1]+vec[2]; mat(1,2) = -vec[2];        mat(1,3) = 0.;
    mat(2,0) = 0.;             mat(2,1) = -vec[2];        mat(2,2) =  vec[2]+vec[3]; mat(2,3) = -vec[3];
    mat(3,0) = 0.;             mat(3,1) = 0.;             mat(3,2) = -vec[3];        mat(3,3) =  vec[3];
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        vec.env().fullRank(),
                        "ckmatrix()",
                        "ndof != 4");
  }

  return;
}

//------------------------------------------------------
// The hysteretic model.
// function [u,ud,udd,ru,resfor,a]=hysteretic_model_example3(thi,a);
// a is the total acceleration at the base of the building
//------------------------------------------------------
void hystereticModel(
  const QUESO::BaseEnvironment& env,
  const QUESO::GslVector&       massInputVec,
  const QUESO::GslVector&       kInputVec,
  const QUESO::GslVector&       rInputVec,
  const QUESO::GslVector&       uInputVec,
  double                        rho,
  double                        gamma,
  const std::vector<double>&    a,
  std::vector<double>&          t,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& u,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& ud,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& udd,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& resfor,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& ru)
{
  bool         useLinear = false; // non linear model
  unsigned int ndof      = massInputVec.sizeLocal();
  unsigned int nsteps    = a.size();
  double       dt        = 0.01; // number of integration steps

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> dofSpace(env, "dof_", ndof, NULL);

  // Time interval, 0.01=10%*5thperiod
  QUESO::GslVector Mvec(massInputVec);
  QUESO::GslVector R   (kInputVec*rInputVec);
  QUESO::GslVector Kd  (kInputVec-R); // stiffness (N/m)
  QUESO::GslVector S_K (uInputVec); // s/k
  QUESO::GslVector S   (S_K*Kd);

  // Initial conditions
  QUESO::GslVector u0 (dofSpace.zeroVector());
  QUESO::GslVector ud0(dofSpace.zeroVector());

  // Restoring force
  std::vector<cellStruct> param(4);
  if (useLinear == false) {
    for (unsigned int w = 0; w < ndof; ++w) {
      restoringForce(Kd[w],R[w],S[w],0.,0.,NULL,
                     &param[w], NULL, NULL, NULL);
    }
  }

  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> p(dofSpace,nsteps,""); // momentum
  for (unsigned int i = 0; i < nsteps; ++i) {
    p.setPositionValues(i,-a[i]*Mvec);
  }

  QUESO::GslMatrix M(Mvec);

  //if (useLinear == false) {
  QUESO::GslVector kTdiag(dofSpace.zeroVector());
  //}
  //else {
  QUESO::GslMatrix kT(dofSpace.zeroVector());
  ckmatrix(Kd+R,kT);
  //}
  QUESO::GslMatrix tmp_mat(dofSpace.zeroVector());
  ckmatrix(kInputVec,tmp_mat);
  QUESO::GslMatrix C(rho*M + gamma*tmp_mat);

  double gammaHere = 0.5;
  double beta      = 0.25; // Average acceleration method

  double bdt      = beta*dt;
  double g_b      = gammaHere/beta;
  double g_2b     = g_b/2.;
  double g_bdt    = g_b/dt;
  double dt2      = dt*dt;
  double bdt2     = beta*dt2;
  double one_bdt2 = 1./bdt2;
  double one_2b   = 1./(2.*beta);

  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix> fs(dofSpace,nsteps,"");

  u.setPositionValues (0,u0 );
  ud.setPositionValues(0,ud0);

  QUESO::GslVector auxResforPosition(dofSpace.zeroVector());
  QUESO::GslVector auxFsPosition    (dofSpace.zeroVector());

  // Restoring Force (again)
  if (useLinear == false) {
    for (unsigned int w = 0; w < ndof; ++w) {
      restoringForce(u0[w],0.,ud0[w],0.,0.,&param[w],
                     NULL, &auxResforPosition[w], &kTdiag[w], &param[w]);
    }
    resfor.setPositionValues(0,auxResforPosition);
    auxResforPosition.matlabDiff(0,auxResforPosition[0],auxFsPosition);
    fs.setPositionValues(0,auxResforPosition);
    ckmatrix(kTdiag,kT);
  }
  else {
    //fs(:,1) = kT*u0;
    fs.setPositionValues(0,kT*u0);
  }

  //udd(:,1) = M\(p(:,1) - C*ud0 - fs(:,1));
  QUESO::GslVector auxPPosition  (dofSpace.zeroVector());
  QUESO::GslVector auxUddPosition(dofSpace.zeroVector());
  p.getPositionValues (0,auxPPosition);
  fs.getPositionValues(0,auxFsPosition);

  //auxUddPosition = M.invertMultiply(auxPPosition - C*ud0 - auxFsPosition);
  //using the fact that inv(M) is diagonal
  QUESO::GslVector invMvec(dofSpace.zeroVector());
  for (unsigned int ii=0; ii<4; ++ii){
      invMvec[ii]=1./Mvec[ii];
  }
  QUESO::GslMatrix invM (invMvec);
  auxUddPosition=invM*(auxPPosition-C*ud0-auxFsPosition);

  udd.setPositionValues(0,auxUddPosition);

  QUESO::GslMatrix A  ((1./bdt)*M + g_b*C);
  QUESO::GslMatrix B  ((1./(2.*beta))*M + dt*(g_2b-1.)*C);
  QUESO::GslMatrix kgT(kT + g_bdt*C + one_bdt2*M);

  QUESO::GslVector fsj    (dofSpace.zeroVector());

  QUESO::GslVector ui     (dofSpace.zeroVector());
  QUESO::GslVector udi    (dofSpace.zeroVector());
  QUESO::GslVector uddi   (dofSpace.zeroVector());
  QUESO::GslVector dpi    (dofSpace.zeroVector());
  QUESO::GslVector dpgi   (dofSpace.zeroVector());
  QUESO::GslVector auxVec1(dofSpace.zeroVector());
  QUESO::GslVector auxVec2(dofSpace.zeroVector());
  QUESO::GslVector resfor0(dofSpace.zeroVector());
  QUESO::GslVector fs0    (dofSpace.zeroVector());

  for (unsigned int i = 0; i < (nsteps-1); ++i) {
    u.getPositionValues  (i,ui);
    ud.getPositionValues (i,udi);
    udd.getPositionValues(i,uddi);
    p.getPositionValues  (i,auxVec1);
    p.getPositionValues  (i+1,auxVec2);
    dpi  = auxVec2 - auxVec1;
    dpgi = dpi + A*udi + B*uddi;

    // Solve for duj from kgi and dpgi using Newton-Raphson
    // BEGIN NEWTON RAPHSON
    resfor.getPositionValues(i,resfor0);
    fs.getPositionValues    (i,fs0);
    QUESO::GslVector dR(dpgi);
    unsigned int j = 0;
    QUESO::GslVector duj(dofSpace.zeroVector());

    QUESO::GslVector ru0(dofSpace.zeroVector());
    QUESO::GslVector rv0(dofSpace.zeroVector());
    ui.matlabDiff (1,ui[0], ru0); //ru0 = [ui(1);  diff(ui) ];
    udi.matlabDiff(1,udi[0],rv0); //rv0 = [udi(1); diff(udi)];

    unsigned int aOutput = 1;
    QUESO::GslVector ruj(dofSpace.zeroVector());
    QUESO::GslVector rvj(dofSpace.zeroVector());
    QUESO::GslVector uj (dofSpace.zeroVector());
    QUESO::GslVector vj (dofSpace.zeroVector());

    while (true) {
      j = j + 1;
      if (j > 10000) {
        // error('Newton Raphson is not converging');
        aOutput=0;
        break;
      }
      if (aOutput == 0) {
        break;
      }
      QUESO::GslVector dujold(duj);
      //duj = kgT.invertMultiply(dR);
      //replaced with exact solutions of linear eqs
      double c3=(dR[1]-kgT(1,0)*dR[0]/kgT(0,0))/kgT(1,2);
      double m3=(kgT(1,0)*kgT(0,1)/kgT(0,0)-kgT(1,1))/kgT(1,2);
      double c4=(dR[2]-kgT(2,2)*c3)/kgT(2,3);
      double m4=-(kgT(2,1)+kgT(2,2)*m3)/kgT(2,3);
      duj[1]=(dR[3]-kgT(3,2)*c3-kgT(3,3)*c4)/(kgT(3,2)*m3+kgT(3,3)*m4);
      duj[0]=dR[0]/kgT(0,0)-kgT(0,1)/kgT(0,0)*duj[1];
      duj[2]=c3+m3*duj[1];duj[3]=c4+m4*duj[1];

      QUESO::GslVector dudj(g_bdt*duj - g_b*udi + dt*(1.-g_2b)*uddi);
      uj = ui  + duj;
      vj = udi + dudj;

      if (useLinear == false) {
        uj.matlabDiff(1,uj[0],ruj); //ruj = [uj(1); diff(uj)];
        vj.matlabDiff(1,vj[0],rvj); //rvj = [vj(1); diff(vj)];

        for (unsigned int w = 0; w < ndof; ++w) {
          restoringForce(ruj[w],ru0[w],rvj[w],rv0[w],resfor0[w],&param[w],
                         NULL, &auxResforPosition[w], &kTdiag[w], &param[w]);
        }

        resfor.setPositionValues(i+1,auxResforPosition);
        auxResforPosition.matlabDiff(0,-auxResforPosition[ndof-1],fsj); //fsj = [-diff(resfor(:,i+1)); resfor(ndof,i+1)];
        fsj *= -1.;

        ckmatrix(kTdiag,kT);
        kgT = kT + g_bdt*C + one_bdt2*M;
      }
      else {
        fsj = kT*uj;
      }
      QUESO::GslVector df(fsj - fs0);
      df -= kT*duj;
      dR = dpgi - df;

      if (((duj-dujold).norm2() < 1.e-12) ||
          (df.norm2()           < 1.e-12)) {
        break;
      }
    } // end while

    fs.setPositionValues(i+1,fsj);
    if (useLinear == false) {
      ru.setPositionValues(i+1,ruj);
    }
    // END NEWTON RAPHSON

    t[i+1] = t[i] + dt;
    QUESO::GslVector duddi((1./bdt2)*duj - (1./bdt)*udi - one_2b*uddi);
    u.setPositionValues  (i+1,uj);
    ud.setPositionValues (i+1,vj);
    udd.setPositionValues(i+1,uddi+duddi);
  } // end for

  return;
}
