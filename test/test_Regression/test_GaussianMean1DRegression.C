#include <iostream>
#include <sstream>
#include <iterator>
#include <exception>
#include <vector>
#include <string>
#include <cmath>
#include <queso/Environment.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarSequence.h>

// Regression Test in One Dimension (Mean of a Gaussian Model)
// Model: Gaussian with known variance samplingVar, unknown mean.
// Objective: Use the Bayes multi-level sampler to estimate the posterior distribution of the unknown mean using an input data set  
// Output: The estimated mean of the posterior distribution and the standard deviation about this mean. Also the number of likelihood function calls. 

// User Input
// Everything in this test is controlled through a queso style input file processed by boost::program_options.  
// Options GaussianMean1DRegression_priorMean, GaussianMean1DRegression_priorVar, GaussianMean1DRegression_samplingVar, GaussianMean1DRegression_dataSet
// are added to the boost::program_options options descriptor and parsed from the input file.
// This allows the user to supply the mean and variance of a gaussian prior, the known model variance, and a string of samples from the "true" 
// model - i.e. the data. 

#define PI 3.14159265358979323846

// default input file name
const std::string inpFileDefault = "test_GaussianMean1DRegression.inp";

// default option values
const double priorMeanODV = 0.0;
const double priorVarODV  = 1.0;
const double samplingVarODV = 1.0;
const std::string dataSetODV = "0.0";

// information needed/provided by the likelihood function
struct likelihoodData {
  double samplingVar;                  // input
  std::vector<double> dataSet;         // input
};
unsigned long likelihoodCalls = 0;     // output

// evaluate log likelihood
template<class P_V,class P_M>
double LikelihoodFunc(
		      const P_V&  paramValues,
		      const P_V*  paramDirection,
		      const void* functionDataPtr,
		      P_V*        gradVector,
		      P_M*        hessianMatrix,
		      P_V*        hessianEffect)
{

  double samplingVar = ((likelihoodData*)functionDataPtr)->samplingVar;
  const std::vector<double> &dataSet = ((likelihoodData*)functionDataPtr)->dataSet;

  int n = dataSet.size();
  double loglh = -(n/2.)*std::log(2*PI*samplingVar);

  double currMu = paramValues[0];
  double diff,sumsq = 0.;
  for(int i=0; i<n; i++) {
    diff = dataSet[i] - currMu;
    sumsq += diff*diff; 
  }
  loglh -= sumsq/(2.*samplingVar);

  likelihoodCalls++;

  return loglh;
}

// QUESO work routine
template<class P_V,class P_M>
void GaussianMean1DRegressionCompute(const QUESO::BaseEnvironment& env, double priorMean, double priorVar, const likelihoodData& dat)
{
  // parameter space: 1-D on (-infinity, infinity)
  QUESO::VectorSpace<P_V,P_M> paramSpace(
					 env,            // queso environment
					 "param_",       // name prefix
					 1,              // dimensions
					 NULL);          // names
  
  P_V paramMin(paramSpace.zeroVector());
  P_V paramMax(paramSpace.zeroVector());
  paramMin[0] = -INFINITY;
  paramMax[0] = INFINITY;
  QUESO::BoxSubset<P_V,P_M> paramDomain(
					"paramBox_",        // name prefix
					paramSpace,         // vector space
					paramMin,           // min values
					paramMax);          // max values

  // gaussian prior with user supplied mean and variance
  P_V priorMeanVec(paramSpace.zeroVector());
  P_V priorVarVec(paramSpace.zeroVector());
  priorMeanVec[0] = priorMean;
  priorVarVec[0] = priorVar;
  QUESO::GaussianVectorRV<P_V,P_M> priorRv("prior_", paramDomain, priorMeanVec, priorVarVec); 

  // likelihood is important
  QUESO::GenericScalarFunction<P_V,P_M> likelihoodFunctionObj(
							      "like_",                        // name prefix 
							      paramDomain,                    // image set
							      LikelihoodFunc<P_V,P_M>,        // routine
							      (void *) &dat,                  // routine data ptr
							      true);                          // routineIsForLn

  QUESO::GenericVectorRV<P_V,P_M> postRv(
      "post_",          // name prefix
       paramSpace);     // image set
 

  // Initialize and solve the Inverse Problem with Bayes multi-level sampling
  QUESO::StatisticalInverseProblem<P_V,P_M> invProb(
      "",                     // name prefix
      NULL,                   // alt options
      priorRv,                // prior RV
      likelihoodFunctionObj,  // likelihood fcn
      postRv);                // posterior RV

  invProb.solveWithBayesMLSampling();

  // extract contents of posterior random variable (i.e. the MCMC samples)
  // store them in a QUESO::ScalarSequence and gather all of them to inter0 root proc
  int N = invProb.postRv().realizer().subPeriod();
  QUESO::ScalarSequence<double> subSamples(env, N, ""); 
  P_V sample(paramSpace.zeroVector());
  for(int i=0; i<N; i++) {   
    invProb.postRv().realizer().realization(sample);
    subSamples[i] = sample[0];
  }
  std::vector<double> unifiedSamples;
  subSamples.getUnifiedContentsAtProc0Only(true, unifiedSamples);

  // get the total number of likelihood calls at inter0 root proc
  unsigned long totalLikelihoodCalls = 0;
  MPI_Reduce(&likelihoodCalls, &totalLikelihoodCalls, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, env.inter0Comm().Comm());

  // compute posterior mean & std. output results on inter0 root proc
  if(env.inter0Rank() == 0) {
    double postMean = 0.0;
    double postVar = 0.0;
  
    // Knuth online mean/variance algorithm
    N = unifiedSamples.size();
    for(int n=1; n<=N; n++) {   
      double delta = unifiedSamples[n-1] - postMean;
      postMean += delta/n;
      postVar += delta*(unifiedSamples[n-1] - postMean);
    }
    postVar /= N;
  
    std::cout<<"Posterior Markov chain length: "<<N<<std::endl;
    std::cout<<"Posterior mean: "<<postMean<<" +/- "<<std::sqrt(postVar)<<std::endl;
    std::cout<<"Likelihood function calls: "<<totalLikelihoodCalls<<std::endl;
  }
}

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environments
  //************************************************
  MPI_Init(&argc,&argv);

  const char* configFile;
  if(argc > 1)
    configFile = argv[1];
  else
    configFile = inpFileDefault.c_str();

   // initilize environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD,       // MPI communicator
			     configFile,           // input file name
			     "",                   // name prefix
			     NULL );               // alt options

  // reset seed to value in input file + fullRank
  env.resetSeed(env.seed() + env.fullRank());

  double priorMean;
  double priorVar;
  std::string dataString;
  likelihoodData dat;

  // Define program options 
  po::options_description myOptDesc("GaussianMean1DRegression specific options");  
  myOptDesc.add_options()
    ("GaussianMean1DRegression_priorMean", po::value<double>(&priorMean)->default_value(priorMeanODV), "Prior Mean")
    ("GaussianMean1DRegression_priorVar", po::value<double>(&priorVar)->default_value(priorVarODV), "Prior Standard Deviation")
    ("GaussianMean1DRegression_samplingVar", po::value<double>(&dat.samplingVar)->default_value(samplingVarODV), "Data Sampling Standard Deviation")
    ("GaussianMean1DRegression_dataSet", po::value<std::string>(&dataString)->default_value(dataSetODV), "Calibration Data");

  // One more pass on the input file to get GaussianMean1DRegression option values
  po::variables_map myOptMap;
  std::ifstream config_stream(configFile);
  try {
    po::store(po::parse_config_file(config_stream, myOptDesc, true), myOptMap);
  } catch(std::exception& e) {
    if(env.fullRank() == 0)
      std::cout<<"Caught exception: "<<e.what()<<std::endl;
  }
  config_stream.close();
  po::notify(myOptMap);

  // parse data string into a vector of doubles stored in dat.dataSet  
  std::istringstream iss(dataString);
  std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(dat.dataSet));

  // output the option values that will be used by QUESO
  if(env.fullRank() == 0) {
    std::cout<<"Read option GaussianMean1DRegression_priorMean: "<<priorMean<<std::endl;
    std::cout<<"Read option GaussianMean1DRegression_priorVar: "<<priorVar<<std::endl;
    std::cout<<"Read option GaussianMean1DRegression_samplingVar: "<<dat.samplingVar<<std::endl;
    std::cout<<"Read option GaussianMean1DRegression_dataSet: ";
    for(unsigned int i=0; i<dat.dataSet.size(); i++) 
      std::cout<<dat.dataSet[i]<<" ";
    std::cout<<std::endl;
  }

  // Work routine, with GSL
  GaussianMean1DRegressionCompute<QUESO::GslVector,QUESO::GslMatrix>(env, priorMean, priorVar, dat);

  MPI_Finalize();

  return 0;
}
  

  
