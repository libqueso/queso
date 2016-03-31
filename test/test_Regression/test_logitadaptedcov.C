#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/UniformVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/GaussianJointPdf.h>
#include <queso/InvLogitGaussianJointPdf.h>
#include <queso/StatisticalInverseProblem.h>

#include <cstdlib>
#include <iomanip>

#define TOL 1e-13

template<class V, class M>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain,
      const QUESO::GaussianJointPdf<V, M> & pdf)
  : QUESO::BaseScalarFunction<V, M>(prefix, domain),
    m_pdf(pdf)
  {
  }

  virtual ~Likelihood()
  {
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return this->m_pdf.lnValue(domainVector, NULL, NULL, NULL, NULL);
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }

  const QUESO::GaussianJointPdf<V, M> & m_pdf;
};

int main(int argc, char ** argv) {
  std::string inputFileName = "test_Regression/adaptedcov_input.txt";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, inputFileName, "", NULL);
#else
  QUESO::FullEnvironment env(inputFileName, "", NULL);
#endif

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 2, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(0.0);

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(1.0);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslVector, QUESO::GslMatrix> prior("prior_",
      paramDomain);

  QUESO::GslVector mean(paramSpace.zeroVector());
  mean.cwSet(0.5);

  QUESO::GslMatrix cov(paramSpace.zeroVector());
  cov(0, 0) = 0.05 * 0.05;
  cov(0, 1) = 0.001;
  cov(1, 0) = 0.001;
  cov(1, 1) = 0.05 * 0.05;

  QUESO::GaussianJointPdf<QUESO::GslVector, QUESO::GslMatrix> pdf("pdf_",
      paramDomain, mean, cov);

  Likelihood<QUESO::GslVector, QUESO::GslMatrix> likelihood("likelihood_",
      paramDomain, pdf);

  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix> posterior(
      "posterior_", paramDomain);

  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix> ip("",
      NULL, prior, likelihood, posterior);

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials.cwSet(0.5);

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0, 0) = 0.5;
  proposalCovMatrix(0, 1) = 0.0;
  proposalCovMatrix(1, 0) = 0.0;
  proposalCovMatrix(1, 1) = 0.5;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  QUESO::GslMatrix adaptedCovMatrix(
      dynamic_cast<const QUESO::InvLogitGaussianJointPdf<QUESO::GslVector,
          QUESO::GslMatrix> &>(
            ip.sequenceGenerator().transitionKernel().rv(0).pdf()).lawCovMatrix());

  std::cout << std::setprecision(15)
            << "Adapted covariance matrix is:" << std::endl
            << adaptedCovMatrix(0, 0) << " " << adaptedCovMatrix(0, 1) << std::endl
            << adaptedCovMatrix(1, 0) << " " << adaptedCovMatrix(1, 1) << std::endl;


  QUESO::GslMatrix regressionTestMatrix(proposalCovMatrix);
  regressionTestMatrix(0, 0) = 0.0162626079275191;
  regressionTestMatrix(0, 1) = 0.0065764502233059;
  regressionTestMatrix(1, 0) = 0.0065764502233059;
  regressionTestMatrix(1, 1) = 0.0154739306675151;

  QUESO::GslMatrix diff(regressionTestMatrix - adaptedCovMatrix);


  std::cout << "norm is " << diff.normFrob() / regressionTestMatrix.normFrob() << std::endl;

  unsigned int result;
  if (diff.normFrob() / regressionTestMatrix.normFrob() > TOL) {
    result = 1;
  }
  else {
    result = 0;
  }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return result;
}
