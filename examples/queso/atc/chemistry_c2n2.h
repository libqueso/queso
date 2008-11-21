#ifndef __UQ_ATC_CHEMISTRY_H__
#define __UQ_ATC_CHEMISTRY_H__
#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>


template <class V,class M>
void
atcLikelihoodRoutine<V,M>::reactionRoutine()
{
     double m_Ts       = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->temperature_s1;
     double m_NSPECI   = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->NSPECI;
     double m_NO_STEPS = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->NO_STEPS;
     const std::vector<double> m_lkc_s1 = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->lkc_s1;
     const std::vector<double> m_la_s1 = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->la_s1;
     const std::vector<double> m_p_s1 = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->p_s1;
     const std::vector<double> m_s_s1 = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->s_s1;
     const std::vector<double> m_reac = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->reac;
     const std::vector<double> m_prod = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->prod;
     const std::vector<double> m_conmo = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->conmo;
     const std::vector<double> m_conmo_old = ((atcLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->conmo_old;
     double m_WDOT[2,M_NSPECI]; 
     double m_WKIN[m_NO_STEPS];
     double m_REK[m_NO_STEPS],m_REB[m_NO_STEPS];
     double m_CO,m_CN,m_CC,m_CARG,m_CO2,m_CN2,m_CNO,m_CCO,m_CCN,m_CNCO,m_CNCN,m_CN2O,m_CC2N2;
     int i, NSP; NS, NR

     m_WKIN = 0.0;    

     for (i=0; i<m_NO_STEPS; i++){
         m_RKF[i]   = ( 10e0 ** m_la_s1[i]) * (m_Ts ** p_s[i]) * EXP(- m_s_s1[i] / m_Ts); 
         m_RKB[i]   = m_RKF[i] / (10e0 ** m_lkc_s1[i]);
     }

     for (NSP=0; NSP<2: NSP++) {
          m_CO    = m_conmo[0];
          m_CN    = m_conmo[1];
          m_CC    = m_conmo[2];
          m_CARG  = m_conmo[3];
          m_CO2   = m_conmo[4];
          m_CN2   = m_conmo[5];
          m_CNO   = m_conmo[6];
          m_CCO   = m_conmo[7];
          m_CCN   = m_conmo[8];
          m_CNCO  = m_conmo[9];;
          m_CNCN  = m_conmo[10];
          m_CN2O  = m_conmo[11];
          m_CC2N2 = m_conmo[12];

          m_WKIN[0]  = -m_RKF[0]  * m_CN2O  * m_CARG + m_RKB[0]  * m_CN2  * m_CO * m_CARG;
          m_WKIN[1]  = -m_RKF[1]  * m_CC2N2 * m_CO   + m_RKB[1]  * m_CCN  * m_CNCO;
          m_WKIN[2]  = -m_RKF[2]  * m_CNCO  * m_CO   + m_RKB[2]  * m_CCO  * m_CNO;
          m_WKIN[3]  = -m_RKF[3]  * m_CNCO  * m_CARG + m_RKB[3]  * m_CN   * m_CCO * m_CARG;
          m_WKIN[4]  = -m_RKF[4]  * m_CCN   * m_CO2  + m_RKB[4]  * m_CNCO * m_CO;
          m_WKIN[5]  = -m_RKF[5]  * m_CCN   * m_CO   + m_RKB[5]  * m_CCO  * m_CN;
          m_WKIN[6]  = -m_RKF[6]  * m_CN2O  * m_CO   + m_RKB[6]  * m_CNO  * m_CNO;
          m_WKIN[7]  = -m_RKF[7]  * m_CN2O  * m_CO   + m_RKB[7]  * m_CN2  * m_CO2;
          m_WKIN[8]  = -m_RKF[8]  * m_CN2   * m_CO   + m_RKB[8]  * m_CN   * m_CNO;
          m_WKIN[9]  = -m_RKF[9]  * m_CNO   * m_CO   + m_RKB[9]  * m_CN   * m_CO2;
          m_WKIN[10] = -m_RKF[10] * m_CNCO  * m_CN   + m_RKB[10] * m_CN2  * m_CCO;
          m_WKIN[11] = -m_RKF[11] * m_CNCO  * m_CN   + m_RKB[11] * m_CCN  * m_CNO;
          m_WKIN[12] = -m_RKF[12] * m_CNCO  * m_CC   + m_RKB[12] * m_CCN  * m_CCO;
          m_WKIN[13] = -m_RKF[13] * m_CCN   * m_CN   + m_RKB[13] * m_CC   * m_CN2;
          m_WKIN[14] = -m_RKF[14] * m_CN2O  * m_CCN  + m_RKB[14] * m_CNCN * m_CNO;
          m_WKIN[15] = -m_RKF[15] * m_CNCN  * m_CO   + m_RKB[15] * m_CCN  * m_CNO;
          m_WKIN[16] = -m_RKF[16] * m_CNCN  * m_CN   + m_RKB[16] * m_CCN  * m_CN2;
          m_WKIN[17] = -m_RKF[17] * m_CC2N2 * m_CARG + m_RKB[17] * m_CCN  * m_CCN * m_CARG;
          m_WKIN[18] = -m_RKF[18] * m_CC2N2 * m_CO   + m_RKB[18] * m_CCO  * m_CNCN;

          m_WDOT[NSP][:] = 0.0

          for(NS = 0; NS < m_NSPECI; NS++){
             for(NR = 0; NR < m_NO_STEPS; NR++1){
                m_WDOT[NSP][NS] = m_WDOT[NSP][NS] + (m_REAC[NR][NS]-m_PROD[NR][NS])*m_WKIN[NR];
             {
          }

          if(NSP == 1) {
             for(NS = 0; NS < m_NSPECI; NS++){
               m_conmo[NS] = m_conmo_old[NS] + DTIME * m_WDOT[1][NS];
             }
          }
          else if (NSP == 2) {
             for(NS = 0; NS < m_NSPECI; NS++){
               m_conmo[NS] = m_conmo_old[NS] + 0.5e0 *  DTIME *  ( WDOT[1][NS] + WDOT[2][NS] );
             }
          }
     }
  return ;
}




