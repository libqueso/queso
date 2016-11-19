#include <queso/Environment.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarSequence.h>

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

// Regression Test in One Dimension (Mean of a Gaussian Model)
//
// Model: Gaussian with known variance samplingVar, unknown mean.
//
// Objective: Use the Bayes multi-level sampler to estimate the posterior
// distribution of the unknown mean using an input data set
//
// Output: The estimated mean of the posterior distribution and the standard
// deviation about this mean. Also the number of likelihood function calls.

#define PI 3.14159265358979323846

// default input file name
std::string inputFileODV = "test_Regression/GaussianMean1DRegression_options";

// information needed/provided by the likelihood function
struct likelihoodData {
  double samplingVar;                  // input
  std::vector<double> dataSet;         // input
};
unsigned int likelihoodCalls = 0;     // output

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

  double samplingVar = ((likelihoodData *)functionDataPtr)->samplingVar;
  const std::vector<double> &dataSet = ((likelihoodData *)functionDataPtr)->dataSet;

  int n = dataSet.size();
  double loglh = -(n / 2.0) * std::log(2.0 * PI * samplingVar);

  double currMu = paramValues[0];
  double diff, sumsq = 0.0;
  for (int i = 0; i < n; i++) {
    diff = dataSet[i] - currMu;
    sumsq += diff * diff;
  }
  loglh -= sumsq / (2.0 * samplingVar);

  likelihoodCalls++;

  return loglh;
}

// QUESO work routine
template<class P_V, class P_M>
void GaussianMean1DRegressionCompute(const QUESO::BaseEnvironment& env,
    double priorMean, double priorVar, const likelihoodData& dat)
{
  // parameter space: 1-D on (-infinity, infinity)
  QUESO::VectorSpace<P_V, P_M> paramSpace(
					 env,       // queso environment
					 "param_",  // name prefix
					 1,         // dimensions
					 NULL);     // names

  P_V paramMin(paramSpace.zeroVector());
  P_V paramMax(paramSpace.zeroVector());
  paramMin[0] = -INFINITY;
  paramMax[0] = INFINITY;
  QUESO::BoxSubset<P_V, P_M> paramDomain(
					"paramBox_",  // name prefix
					paramSpace,   // vector space
					paramMin,     // min values
					paramMax);    // max values

  // gaussian prior with user supplied mean and variance
  P_V priorMeanVec(paramSpace.zeroVector());
  P_V priorVarVec(paramSpace.zeroVector());
  priorMeanVec[0] = priorMean;
  priorVarVec[0] = priorVar;
  QUESO::GaussianVectorRV<P_V, P_M> priorRv("prior_", paramDomain, priorMeanVec,
      priorVarVec);

  // likelihood is important
  QUESO::GenericScalarFunction<P_V, P_M> likelihoodFunctionObj(
							      "like_",                   // name prefix
							      paramDomain,               // image set
							      LikelihoodFunc<P_V, P_M>,  // routine
							      (void *) &dat,             // routine data ptr
							      true);                     // routineIsForLn

  QUESO::GenericVectorRV<P_V, P_M> postRv(
      "post_",       // name prefix
       paramSpace);  // image set


  // Initialize and solve the Inverse Problem with Bayes multi-level sampling
  QUESO::StatisticalInverseProblem<P_V, P_M> invProb(
      "",                     // name prefix
      NULL,                   // alt options
      priorRv,                // prior RV
      likelihoodFunctionObj,  // likelihood fcn
      postRv);                // posterior RV

  invProb.solveWithBayesMLSampling();

  // compute mean and second moment of samples on each proc via Knuth online mean/variance algorithm
  int N = invProb.postRv().realizer().subPeriod();
  double subMean = 0.0;
  double subM2 = 0.0;
  double delta;
  P_V sample(paramSpace.zeroVector());
  for (int n = 1; n <= N; n++) {
    invProb.postRv().realizer().realization(sample);
    delta = sample[0] - subMean;
    subMean += delta / n;
    subM2 += delta * (sample[0] - subMean);
  }

  // gather all Ns, means, and M2s to proc 0
  std::vector<int> unifiedNs(env.inter0Comm().NumProc());
  std::vector<double> unifiedMeans(env.inter0Comm().NumProc());
  std::vector<double> unifiedM2s(env.inter0Comm().NumProc());
  env.inter0Comm().Gather<int>(&N, 1, &(unifiedNs[0]), 1, 0, "", "");
  env.inter0Comm().Gather<double>(&subMean, 1, &(unifiedMeans[0]), 1, 0, "", "");
  env.inter0Comm().Gather<double>(&subM2, 1, &(unifiedM2s[0]), 1, 0, "", "");

  // get the total number of likelihood calls at proc 0
  unsigned int totalLikelihoodCalls = 0;
  env.inter0Comm().Allreduce<unsigned int>(&likelihoodCalls, &totalLikelihoodCalls, 1,
      RawValue_MPI_SUM, "", "");

  // compute global posterior mean and std via Chan algorithm, output results on proc 0
  if (env.inter0Rank() == 0) {
    int postN = unifiedNs[0];
    double postMean = unifiedMeans[0];
    double postVar = unifiedM2s[0];
    for (unsigned int i = 1; i < unifiedNs.size(); i++) {
      delta = unifiedMeans[i] - postMean;
      postMean = (postN * postMean + unifiedNs[i] * unifiedMeans[i]) /
        (postN + unifiedNs[i]);
      postVar += unifiedM2s[i] + delta * delta *
        (((double)postN * unifiedNs[i]) / (postN + unifiedNs[i]));
      postN += unifiedNs[i];
    }
    postVar /= postN;

    //compute exact answer - available in this case since the exact posterior is a gaussian
    N = dat.dataSet.size();
    double dataSum = 0.0;
    for (int i = 0; i < N; i++)
      dataSum += dat.dataSet[i];
    double datMean = dataSum / N;
    double postMeanExact = (N * priorVar / (N * priorVar + dat.samplingVar)) *
      datMean + (dat.samplingVar / (N * priorVar + dat.samplingVar)) * priorMean;
    double postVarExact = 1.0 / (N / dat.samplingVar + 1.0 / priorVar);

    std::cout << "Number of posterior samples: " << postN << std::endl;
    std::cout << "Estimated posterior mean: " << postMean << " +/- "
      << std::sqrt(postVar) << std::endl;
    std::cout << "Likelihood function calls: " << totalLikelihoodCalls
      << std::endl;
    std::cout << "\nExact posterior: Gaussian with mean " << postMeanExact
      << ", standard deviation " << std::sqrt(postVarExact) << std::endl;
  }
}

int main(int argc, char* argv[]) {
  //************************************************
  // Initialize environments
  //************************************************
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank = 0;
#endif

  // Find correct file for out-of-source builds
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileODV = test_srcdir + ('/' + inputFileODV);

  // default parameters
  double priorMean = 5.0;
  double priorVar = 1.0;
  std::string dataString("-0.1437   -0.0982   -1.3616   -0.0171   -1.0290    0.7133   -0.5496   -0.2415   -0.0933   -0.0669   -0.2311   -1.5371   -0.7418   -0.2220    0.0209    1.6905   -0.8379    0.4233   1.7065    1.5676    0.1265    1.0728    1.8126    2.7186   -0.6680    0.5251   -0.4768    0.6027   -1.9304   -2.1159   1.8698    0.3908   -1.2513    1.2597    0.1426   -0.6527   -0.5948   -0.6698    0.3874    0.3455   -0.1988   -0.4574   0.6930    1.3184   -0.2715   -0.9510    0.3593   -0.4184   0.5243    0.2846   -0.3556   -1.4877   -0.8829    1.0829   -0.7700   -0.1138    0.6430   -0.7929   -1.6493    0.1000");
  likelihoodData dat;
  dat.samplingVar = 1.0;

  // parse the data string into a vector of doubles and store them in dat.dataSet
  std::istringstream iss(dataString);
  std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(dat.dataSet));

  // output the default values that will be used by QUESO
  if (rank == 0) {
    std::cout << "Read option GaussianMean1DRegression_priorMean: "
      << priorMean << std::endl;
    std::cout << "Read option GaussianMean1DRegression_priorVar: "
      << priorVar << std::endl;
    std::cout << "Read option GaussianMean1DRegression_samplingVar: "
      << dat.samplingVar << std::endl;
    std::cout << "Read option GaussianMean1DRegression_dataSet: ";
    for (unsigned int i = 0; i < dat.dataSet.size(); i++)
      std::cout << dat.dataSet[i] << " ";
    std::cout << std::endl;
  }

   // initilize environment
#ifdef QUESO_HAS_MPI
  QUESO::FullEnvironment env(MPI_COMM_WORLD,  // MPI communicator
			     inputFileODV.c_str(),  // input file name
			     "",                 // name prefix
			     NULL );             // alt options
#else
  QUESO::FullEnvironment env(  // No MPI communicator
           inputFileODV.c_str(),  // input file name
           "",                 // name prefix
           NULL );             // alt options
#endif

  // reset seed to value in input file + fullRank
  env.resetSeed(env.seed() + env.fullRank());

  // Work routine, with GSL
  GaussianMean1DRegressionCompute<QUESO::GslVector, QUESO::GslMatrix>(env,
      priorMean, priorVar, dat);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
