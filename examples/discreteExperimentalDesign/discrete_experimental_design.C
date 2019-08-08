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

#include <queso/Environment.h>
#include <queso/VectorSpace.h>
#include <queso/VectorRV.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/ScalarFunction.h>
#include <queso/GridSearchExperimentalDesign.h>
#include <queso/ExperimentalLikelihoodInterface.h>
#include <queso/ExperimentalLikelihoodWrapper.h>
#include <queso/GaussianLikelihoodDiagonalCovariance.h>
#include <queso/ExperimentMetricMinVariance.h>
#include <queso/GaussianVectorRV.h>

#include <fstream>

/*
 *  This is a demonstration of the usage of the GridSearchExperimentalDesign infrastructure through a simple chemical kinetics example.
    The goal is to infer the 'activation energy' terms in Arrhenius reaction rate equations from synthetic noisy data.

 *  Requirements:
      1. Implementation of Likelihood class that is problem-specific, preferably one that derives from an existing class.
      2. Implementation of LikelihoodInterface to allow the scenario parameters to be passed from the experimental design to the likelihood.
      3. main() function that creates all necessary objects, provides (or calculates) a prior for the inverse problems,
            and runs the experimental design.

 *   This file provides a framework that can be readily modified to suit the user's needs for a given problem.
 */



/*
 * Likelihood evaluation function derived from GaussianLikelihoodDiagonalCovariance
 */
template<typename V = QUESO::GslVector, typename M = QUESO::GslMatrix>
class Likelihood :  public QUESO::GaussianLikelihoodDiagonalCovariance<V,M>
{
  public:

    Likelihood( const QUESO::VectorSet<V,M> & domain,
                const V & observations,
                const V & covariance,
                const V & T_vec,
                QUESO::GaussianVectorRV<V,M> & noise)
    : QUESO::GaussianLikelihoodDiagonalCovariance<V,M>("exp_like_",domain,observations,covariance),
      m_T_vec(T_vec),
      m_noise(noise)
    {}

    virtual void evaluateModel(const V & domainVector, V & modelOutput) const override
    {
      unsigned int num_data = this->m_observations.sizeGlobal() - 1;

      double delta_G1  = domainVector[0]; // 0.2
      double delta_Gm1 = domainVector[1]; // 0.3
      double delta_G2  = domainVector[2]; // 0.05
      double delta_Gm2 = domainVector[3]; // 0.1

      for (unsigned int i = 0; i < num_data; ++i)
      {
        double T = m_T_vec[i];
        modelOutput[i] = this->get_theta_c(delta_G1,delta_Gm1,delta_G2,delta_Gm2,T,this->m_env);
      }

      modelOutput[num_data] = this->get_theta_c(delta_G1,delta_Gm1,delta_G2,delta_Gm2,m_T,this->m_env);

    }

    // Need to implement this when using ML Sampling
    virtual double lnValue( const V & domainVector,
                            const V * domainDirection,
                            V * gradVector,
                            M * hessianMatrix,
                            V * hessianEffect) const override
    {
      return QUESO::GaussianLikelihoodDiagonalCovariance<V,M>::lnValue(domainVector);
    }

    // Set the scenario class variable and add a new synthetic noisy data point
    void setScenario(std::vector<double> & scenario_params)
    {
      m_T = scenario_params[0];

      unsigned int num_data = this->m_observations.sizeGlobal();

      // reference parameter values
      double delta_G1  = 0.2;
      double delta_Gm1 = 0.3;
      double delta_G2  = 0.05;
      double delta_Gm2 = 0.1;

      double value = this->get_theta_c(delta_G1,delta_Gm1,delta_G2,delta_Gm2,m_T,this->m_env);

      QUESO::GslVector sample( m_noise.imageSet().vectorSpace().zeroVector() );
      m_noise.realizer().realization(sample);
      double data_point = value + sample[0];

      const_cast<V &>(this->m_observations)[num_data-1] = data_point;
    }

    // Helper function to evaluate kinetic QoI
    static double get_theta_c(double delta_G1, double delta_Gm1, double delta_G2, double delta_Gm2, double T, const QUESO::BaseEnvironment & env)
    {
      double A1  = 1.0;
      double Am1 = 1.0;
      double A2  = 1.0;
      double Am2 = 1.0;

      double kB = 8.6173304e-5;

      double k1  = A1  * std::exp(-delta_G1/(kB*T));
      double km1 = Am1 * std::exp(-delta_Gm1/(kB*T));
      double k2  = A2  * std::exp(-delta_G2/(kB*T));
      double km2 = Am2 * std::exp(-delta_Gm2/(kB*T));

      QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> kinetic_space( env, "kinetics_", 3, NULL );

      M A(kinetic_space.zeroVector());
      A(0,0) = -k1; A(0,1) = km1;     A(0,2) = 0;
      A(1,0) =  k1; A(1,1) = -km1-k2; A(1,2) = km2;
      A(2,0) = 1.0; A(2,1) = 1.0;     A(2,2) = 1.0;

      V b(kinetic_space.zeroVector());
      b[0] = 0.0;
      b[1] = 0.0;
      b[2] = 1.0;

      V x = A.invertMultiply(b);

      return x[2]; // theta_c
    }

  private:
    // Temperature for current scenario
    double m_T;

    // Vector of previous scenarios run
    const V & m_T_vec;

    // Noise RV for generating synthetic data in setScenario()
    const QUESO::GaussianVectorRV<V,M> & m_noise;

};


/*
 * Interface class to pass the scenario temperature to Likelihood
 */
template<typename V = QUESO::GslVector, typename M = QUESO::GslMatrix>
class LikelihoodInterface : public QUESO::ExperimentalLikelihoodInterface<V,M>
{
public:

  LikelihoodInterface() {}

  virtual void reinit(std::vector<double> & scenario_params, QUESO::BaseScalarFunction<V,M> & likelihood) override
  {
    dynamic_cast<Likelihood<V,M> &>(likelihood).setScenario(scenario_params);
  }

};


/*
 * Separate likelihood class for the inverse problem to get the initial belief posterior
 */
template<typename V = QUESO::GslVector, typename M = QUESO::GslMatrix>
class SIPLikelihood :  public QUESO::GaussianLikelihoodDiagonalCovariance<V,M>
{
  public:

    SIPLikelihood( const QUESO::VectorSet<V,M> & domain,
                    const V & observations,
                    const V & covariance,
                    const V & T_vec)
    : QUESO::GaussianLikelihoodDiagonalCovariance<V,M>("exp_like_",domain,observations,covariance),
      m_T_vec(T_vec)
    {}

    virtual void evaluateModel(const V & domainVector, V & modelOutput) const override
    {
      unsigned int num_data = this->m_observations.sizeGlobal() - 1;

      double delta_G1  = domainVector[0]; // 0.2
      double delta_Gm1 = domainVector[1]; // 0.3
      double delta_G2  = domainVector[2]; // 0.05
      double delta_Gm2 = domainVector[3]; // 0.1

      for (unsigned int i = 0; i < num_data; ++i)
      {
        double T = m_T_vec[i];
        modelOutput[i] = Likelihood<V,M>::get_theta_c(delta_G1,delta_Gm1,delta_G2,delta_Gm2,T,this->m_env);
      }

    }

  private:
    // Vector of initial experiments
    const V & m_T_vec;

};


/*
 * Main function
 */
int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);

    QUESO::FullEnvironment env(MPI_COMM_WORLD,"./queso.in", "", NULL);

//////////////////////////////////////////////
/////////// Noise ////////////////////////////
//////////////////////////////////////////////

    // Set up the noise RV for generating synthetic data points

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Setting up the Noise parameters" <<std::endl
                <<"======================================" <<std::endl;    

    double sigma = 0.1;

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> noise_space( env, "noise_", 1, NULL );

    QUESO::GslVector noise_min( noise_space.zeroVector() );
    noise_min[0] = -1.0;

    QUESO::GslVector noise_max( noise_space.zeroVector() );
    noise_max[0] = 1.0;

    QUESO::GslVector noise_mean(noise_space.zeroVector());
    noise_mean.cwSet(0.0);

    QUESO::GslVector noise_var(noise_space.zeroVector());
    noise_var.cwSet(sigma*sigma);

    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> noise_domain( "noise_domain_", noise_space, noise_min, noise_max );
    QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> noise_rv( "noise_rv_", noise_domain, noise_mean, noise_var );


//////////////////////////////////////////////
/////////// Data /////////////////////////////
//////////////////////////////////////////////

    // These are the initial experiments that have already been run

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Setting up the Data" <<std::endl
                <<"======================================" <<std::endl;

    unsigned int n_obs = 4;
    unsigned int n_data = n_obs + 1;

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> data_space( env, "data_", n_data, NULL );

    QUESO::GslVector variance( data_space.zeroVector() );
    variance.cwSet(sigma*sigma);

    QUESO::GslVector T_vec( data_space.zeroVector() );
    T_vec[0] = 299.0;
    T_vec[1] = 312.0;
    T_vec[2] = 325.0;
    T_vec[3] = 351.0;

    QUESO::GslVector sample( noise_space.zeroVector() );
    QUESO::GslVector obs( data_space.zeroVector() );

    // ref param values
    double delta_G1  = 0.2;
    double delta_Gm1 = 0.3;
    double delta_G2  = 0.05;
    double delta_Gm2 = 0.1;

    for (unsigned int i=0; i<n_obs; ++i)
      {
        double thetaC = Likelihood<QUESO::GslVector,QUESO::GslMatrix>::get_theta_c(delta_G1,delta_Gm1,delta_G2,delta_Gm2,T_vec[i],env);
        noise_rv.realizer().realization(sample);
        obs[i] = thetaC + sample[0];
      }


//////////////////////////////////////////////
/////////// SIP Parameters ///////////////////
//////////////////////////////////////////////

    // Min/Max values for the inverse problem runs

    if (env.fullRank() == 0)
      std::cout <<"\n\n======================================" <<std::endl
                <<"Setting up the Inverse Problem parameters" <<std::endl
                <<"======================================" <<std::endl;

    unsigned int n_sip_params = 4;

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> param_space( env, "param_", n_sip_params, NULL );

    QUESO::GslVector min_vals( param_space.zeroVector() );
    min_vals[0] = 0.0;
    min_vals[1] = 0.0;
    min_vals[2] = 0.0;
    min_vals[3] = 0.0;

    QUESO::GslVector max_vals( param_space.zeroVector() );
    max_vals[0] = 1.0;
    max_vals[1] = 1.0;
    max_vals[2] = 1.0;
    max_vals[3] = 1.0;

    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> param_domain( "param_domain_", param_space, min_vals, max_vals );


//////////////////////////////////////////////
/////////// Run SIP for Current Belief ///////
//////////////////////////////////////////////

    // Based on the initial experiments (obs), run an inverse problem to get the current parameter belief distribution

    if (env.fullRank() == 0)
      std::cout <<"\n\n======================================" <<std::endl
                <<"Running a SIP to get current param beliefs" <<std::endl
                <<"======================================" <<std::endl;

    QUESO::GslVector prior_mean(param_space.zeroVector());
    prior_mean[0] = 0.4;
    prior_mean[1] = 0.4;
    prior_mean[2] = 0.2;
    prior_mean[3] = 0.2;

    QUESO::GslMatrix prior_cov(param_space.zeroVector());
    prior_cov(0,0) = 0.1;
    prior_cov(1,1) = 0.1;
    prior_cov(2,2) = 0.1;
    prior_cov(3,3) = 0.1;

    QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> sip_prior("sip_prior_", param_space,prior_mean,prior_cov);    

    QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> sip_post("sip_post_", param_space);

    QUESO::GslVector paramInitials( param_space.zeroVector() );
    paramInitials.cwSet(0.15);

    QUESO::GslMatrix propCovMatrix(param_space.zeroVector());
    propCovMatrix(0,0) = 0.01/20.0;
    propCovMatrix(1,1) = 0.01/20.0;
    propCovMatrix(2,2) = 0.01/20.0;
    propCovMatrix(3,3) = 0.01/20.0;

    SIPLikelihood<QUESO::GslVector,QUESO::GslMatrix> sip_likelihood( param_domain, obs, variance, T_vec );

    QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix> ip("", NULL,
                                                                          sip_prior,
                                                                          sip_likelihood,
                                                                          sip_post);

    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &propCovMatrix);

    // This will be the prior for all experimental scenarios below, set the type appropriately
    std::shared_ptr<QUESO::BaseVectorRV<QUESO::GslVector,QUESO::GslMatrix>> prior( new QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>(sip_post) );


//////////////////////////////////////////////
/////////// Scenario /////////////////////////
//////////////////////////////////////////////

    // Specifies the scenario space limits and discretization size
    // Our scenario space for this example is just the temperature T

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Setting up the Scenario parameters" <<std::endl
                <<"======================================" <<std::endl;

    unsigned int n_scenario_params = 1; // Ta

    QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> scenario_space( env, "scenario_", n_scenario_params, NULL );

    QUESO::GslVector scenario_min( scenario_space.zeroVector() );
    scenario_min[0] = 298.0;

    QUESO::GslVector scenario_max( scenario_space.zeroVector() );
    scenario_max[0] = 498.0;

    QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> scenario_domain( "scenario_domain_", scenario_space, scenario_min, scenario_max );

    // Number of discrete points for each scenario parameter
    std::vector<unsigned int> n_points(n_scenario_params);
    n_points[0] = 11;


//////////////////////////////////////////////
/////////// Likelihood ///////////////////////
//////////////////////////////////////////////

    // Create the likelihood and likelihood interface classes

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Setting up the Likelihood" <<std::endl
                <<"======================================" <<std::endl;    

    std::shared_ptr<QUESO::BaseScalarFunction<QUESO::GslVector,QUESO::GslMatrix>> likelihood( new Likelihood<QUESO::GslVector,QUESO::GslMatrix>(param_domain, obs, variance, T_vec, noise_rv) );

    std::shared_ptr<QUESO::ExperimentalLikelihoodInterface<QUESO::GslVector,QUESO::GslMatrix>> interface( new LikelihoodInterface<QUESO::GslVector,QUESO::GslMatrix>() );


//////////////////////////////////////////////
/////////// Experimental Design //////////////
//////////////////////////////////////////////

    // Create the ExperimentMetric and GridSearchExperimentalDesign objects

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Setting up the Experimental Design classes" <<std::endl
                <<"======================================" <<std::endl;  

    // MinVar metric uses Metropolis-Hastings to solve the SIP, EIG metric uses ML Sampling
    std::shared_ptr<QUESO::ExperimentMetricBase<QUESO::GslVector,QUESO::GslMatrix>> metric( new QUESO::ExperimentMetricMinVariance<QUESO::GslVector,QUESO::GslMatrix>(paramInitials,propCovMatrix) );

    std::shared_ptr<QUESO::ExperimentalLikelihoodWrapper<QUESO::GslVector,QUESO::GslMatrix>> wrapper( new QUESO::ExperimentalLikelihoodWrapper<QUESO::GslVector,QUESO::GslMatrix>(likelihood,interface) );

    std::shared_ptr<QUESO::ScenarioRunner<QUESO::GslVector,QUESO::GslMatrix>> runner( new QUESO::ScenarioRunner<QUESO::GslVector,QUESO::GslMatrix>(prior,wrapper,metric) );

    // If you want to output the rawchain for each scenario, set the output type    
    // runner->set_output_rawchain("h5");

    QUESO::GridSearchExperimentalDesign<QUESO::GslVector,QUESO::GslMatrix> exp_design(scenario_domain,n_points,runner);


//////////////////////////////////////////////
/////////// Get Results //////////////////////
//////////////////////////////////////////////

    // Iterate over the scenario space and run the SIP for each scenario

    if (env.fullRank() == 0)
      std::cout <<"======================================" <<std::endl
                <<"Running the Experimental Design" <<std::endl
                <<"======================================" <<std::endl;

    QUESO::GslVector experimental_params(param_space.zeroVector());

    // Give a non-empty output_prefix string to have the rawChain for each scenario SIP output
    std::string output_prefix = "";
    exp_design.run(experimental_params,output_prefix);


//////////////////////////////////////////////
/////////// Print Results ////////////////////
//////////////////////////////////////////////

    // Output the scenario with the highest metric to the screen

    if (env.fullRank() == 0)
      std::cout <<"\n==================================================" <<std::endl
                <<"Selected Experimental Design Scenario: " <<std::endl
                <<"T = " <<experimental_params[0] <<std::endl
                <<"==================================================" <<std::endl;

  MPI_Finalize();

  return 0;
}

