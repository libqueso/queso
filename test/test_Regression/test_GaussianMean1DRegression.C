#include <iostream>
#include <sstream>
#include <iterator>
#include <exception>
#include <vector>
#include <string>
#include <cmath>

#include <boost/program_options.hpp>

#include <queso/Environment.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarSequence.h>

// Regression Test in One Dimension (Mean of a Gaussian Model)
//
// Model: Gaussian with known variance samplingVar, unknown mean.
//
// Objective: Use the Bayes multi-level sampler to estimate the posterior
// distribution of the unknown mean using an input data set
//
// Output: The estimated mean of the posterior distribution and the standard
// deviation about this mean. Also the number of likelihood function calls.

// User Input
//
// Everything in this test is controlled through a queso style input file
// processed by boost::program_options. The name of this file is passed to the
// program through the command line. If no name is passed in, the program tries
// the default file name specified by inputFileODV.  Options
// GaussianMean1DRegression_priorMean, GaussianMean1DRegression_priorVar,
// GaussianMean1DRegression_samplingVar, GaussianMean1DRegression_dataSet are
// added to the boost::program_options options descriptor and parsed from the
// input file. These options can also be specified on the command line, in
// which case the command line value will overwrite the input file value.  This
// allows the user to supply the mean and variance of a gaussian prior, the
// known model variance, and a string of samples from the "true" model - i.e.
// the data.  Running the program with the option '--help', or '-h', will show
// the proper usage of the program and list all the available options.

#define PI 3.14159265358979323846

// default input file name
const std::string inputFileODV = "test_Regression/GaussianMean1DRegression_options";

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
  MPI_Gather(&N, 1, MPI_INT, &(unifiedNs[0]), 1, MPI_INT, 0,
      env.inter0Comm().Comm());
  MPI_Gather(&subMean, 1, MPI_DOUBLE, &(unifiedMeans[0]), 1, MPI_DOUBLE, 0,
      env.inter0Comm().Comm());
  MPI_Gather(&subM2, 1, MPI_DOUBLE, &(unifiedM2s[0]), 1, MPI_DOUBLE, 0,
      env.inter0Comm().Comm());

  // get the total number of likelihood calls at proc 0
  unsigned long totalLikelihoodCalls = 0;
  MPI_Reduce(&likelihoodCalls, &totalLikelihoodCalls, 1, MPI_UNSIGNED_LONG,
      MPI_SUM, 0, env.inter0Comm().Comm());

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
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // variables to be filled by command line and/or input file values
  std::string inputFile;
  double priorMean;
  double priorVar;
  std::string dataString;
  likelihoodData dat;

  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "Produce help message");

  boost::program_options::options_description config("GaussianMean1DRegression specific options");
  config.add_options()
    ("GaussianMean1DRegression_priorMean", boost::program_options::value<double>(&priorMean)->default_value(priorMeanODV), "Prior Mean")
    ("GaussianMean1DRegression_priorVar", boost::program_options::value<double>(&priorVar)->default_value(priorVarODV), "Prior Standard Deviation")
    ("GaussianMean1DRegression_samplingVar", boost::program_options::value<double>(&dat.samplingVar)->default_value(samplingVarODV), "Data Sampling Standard Deviation")
    ("GaussianMean1DRegression_dataSet", boost::program_options::value<std::string>(&dataString)->default_value(dataSetODV), "Calibration Data");

  boost::program_options::options_description hidden("hidden options");
  hidden.add_options()
    ("input_file", boost::program_options::value<std::string>(&inputFile)->default_value(inputFileODV), "QUESO configuration file");
  boost::program_options::positional_options_description p;
  p.add("input_file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(hidden);
  boost::program_options::options_description visible("Allowed options");
  visible.add(generic).add(config);

  // parse command line for help, inputFile, and any config options present
  boost::program_options::variables_map vm;
  try {
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
	      options(cmdline_options).positional(p).run(), vm);
  } catch (std::exception& e) {
    if(rank == 0)
      std::cout<<"Caught exception: "<<e.what()<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  boost::program_options::notify(vm);

  // if help is requested, print out usgae and visible options, then exit cleanly
  if (vm.count("help")) {
    if (rank == 0) {
      std::cout << "Usage: " << argv[0] << " <input_file (=" << inputFileODV
        << ")>" << std::endl;
      std::cout << visible << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  // parse the input file to get GaussianMean1DRegression option values
  std::ifstream config_stream(inputFile.c_str());
  try {
    boost::program_options::store(boost::program_options::parse_config_file(config_stream, config, true), vm);
  } catch (std::exception& e) {
    if(rank == 0)
      std::cout << "Caught exception: " << e.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  boost::program_options::notify(vm);
  config_stream.close();

   // parse the data string into a vector of doubles and store them in dat.dataSet
  std::istringstream iss(dataString);
  std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(dat.dataSet));

  // output the option values that will be used by QUESO
  if (rank == 0) {
    std::cout << "QUESO options file: " << inputFile << std::endl;
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
  QUESO::FullEnvironment env(MPI_COMM_WORLD,  // MPI communicator
			     inputFile.c_str(),  // input file name
			     "",                 // name prefix
			     NULL );             // alt options

  // reset seed to value in input file + fullRank
  env.resetSeed(env.seed() + env.fullRank());

  // Work routine, with GSL
  GaussianMean1DRegressionCompute<QUESO::GslVector, QUESO::GslMatrix>(env,
      priorMean, priorVar, dat);

  MPI_Finalize();

  return 0;
}
