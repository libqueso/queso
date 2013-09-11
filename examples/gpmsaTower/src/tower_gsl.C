#include <tower_gsl.h>
#include <uqGpmsaComputerModel.h>
#include <uqJointPdf.h>
#include <uqVectorSpace.h>
#include <uqGslMatrix.h>

int main(int argc, char* argv[])
{
  //***********************************************************************
  // Initialize MPI
  //***********************************************************************
  MPI_Init(&argc,&argv);

  //***********************************************************************
  // Initialize QUESO environment
  //***********************************************************************
  UQ_FATAL_TEST_MACRO((argc < 2),
                      UQ_UNAVAILABLE_RANK,
                      "main()",
                      "run as <executable> 'inputFileName' 'useExperiments<default=0>' 'useML<default=0>'");
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);
  bool useExperiments = false;
  if (argc >=3) {
    useExperiments = !strcmp(argv[2],"1");
  }
  bool useML = false;
  if (argc >=4) {
    useML = !strcmp(argv[3],"1");
  }

  //***********************************************************************
  // Run program
  //***********************************************************************
  compute(*env,useExperiments,useML);

  //***********************************************************************
  // Finalize QUESO environment
  //***********************************************************************
  delete env;

  //***********************************************************************
  // Finalize MPI
  //***********************************************************************
  MPI_Finalize();

  return 0;
}

void compute(const uqFullEnvironmentClass& env, bool useExperiments, bool useML)
{
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Entering compute()..."
                          << std::endl;
  }

  //***********************************************************************
  // Step 01 of 09: Instantiate parameter space, parameter domain, and prior Rv
  //***********************************************************************
  unsigned int paper_p_t = 1; // parameter dimension; 'p_t' in paper; drag coefficient 'C' in tower example
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paper_p_t_space(env, "param_", paper_p_t, NULL);

  uqGslVectorClass paramMins(paper_p_t_space.zeroVector());
  uqGslVectorClass paramMaxs(paper_p_t_space.zeroVector());
  paramMins[0] = 0.;
  paramMaxs[0] = 1.;

  uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass>
    paramDomain("param_",paper_p_t_space,paramMins,paramMaxs);

  uqUniformVectorRVClass<uqGslVectorClass,uqGslMatrixClass>* paramPriorRvPtr = NULL;
//uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>* paramPriorRvPtr = NULL;

  //***********************************************************************
  // Step 02 of 09: Instantiate the 'scenario' and 'output' spaces
  //***********************************************************************
  unsigned int paper_p_x = 1; // scenario dimension; 'p_x' in paper; ball radius 'R' in tower example
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paper_p_x_space(env, "scenario_", paper_p_x, NULL);

  unsigned int paper_n_eta = 16; // simulation output dimension; 'n_eta' in paper; 16 in tower example
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paper_n_eta_space(env, "output_", paper_n_eta, NULL);

  //***********************************************************************
  // Step 03 of 09: Instantiate the simulation storage
  // Regarding simulation scenario input values, QUESO will standardize them so that
  //    they exist inside a hypercube: this will be done in the 'uqSimulationModelClass'
  //    constructor, step 04 of 09.
  // Regarding simulation output data, QUESO will transform it so that the mean is
  //    zero and the variance is one: this will be done in the 'uqSimulationModelClass'
  //    constructor, step 04 of 09.
  //***********************************************************************
  unsigned int paper_m = 25; // number of simulations; 'm' in paper; 25 in tower example
  uqSimulationStorageClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>
    simulationStorage(paper_p_x_space,paper_p_t_space,paper_n_eta_space,paper_m);

  // Add simulations: what if none?
  uqGslVectorClass extraSimulationGridVec(paper_n_eta_space.zeroVector());
  std::vector<uqGslVectorClass* > simulationScenarios(paper_m,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* > paramVecs          (paper_m,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* > outputVecs         (paper_m,(uqGslVectorClass*) NULL);

  for (unsigned int i = 0; i < paper_m; ++i) {
    simulationScenarios[i] = new uqGslVectorClass(paper_p_x_space.zeroVector());   // 'x_{i+1}^*' in paper
    paramVecs          [i] = new uqGslVectorClass(paper_p_t_space.zeroVector());   // 't_{i+1}^*' in paper
    outputVecs         [i] = new uqGslVectorClass(paper_n_eta_space.zeroVector()); // 'eta_{i+1}' in paper
  }

  // Populate all objects just instantiated
  std::set<unsigned int> tmpSetStep03;
  tmpSetStep03.insert(env.subId());

  extraSimulationGridVec.subReadContents("inputData/extraSimulHeights",
                                         "m",
                                         tmpSetStep03);

  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paper_m_space(env, "simulation_", paper_m, NULL);
  uqGslVectorClass inputSimulScenarioVector(paper_m_space.zeroVector());
  inputSimulScenarioVector.subReadContents("inputData/allSimulScenarios",
                                           "m",
                                           tmpSetStep03);
  uqGslVectorClass inputSimulParameterVector(paper_m_space.zeroVector());
  inputSimulParameterVector.subReadContents("inputData/allSimulParameters",
                                            "m",
                                            tmpSetStep03);
  uqGslMatrixClass inputSimulOutputsMatrix(env,paper_n_eta_space.map(),paper_m);
  inputSimulOutputsMatrix.subReadContents("inputData/allSimulOutputs",
                                          "m",
                                          tmpSetStep03);

  for (unsigned int i = 0; i < paper_m; ++i) {
    (*(simulationScenarios[i]))[0] = inputSimulScenarioVector [i];
    (*(paramVecs          [i]))[0] = inputSimulParameterVector[i];
    inputSimulOutputsMatrix.getColumn(i,*(outputVecs[i]));
  }

  // Finally, add information to the simulation storage
  for (unsigned int i = 0; i < paper_m; ++i) {
    simulationStorage.addSimulation(*(simulationScenarios[i]),*(paramVecs[i]),*(outputVecs[i]));
  }

  //***********************************************************************
  // Step 04 of 09: Instantiate the simulation model
  //***********************************************************************
  uqSimulationModelClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>
    simulationModel("",   // prefix
                    NULL, // options values
                    simulationStorage);
  unsigned int paper_p_eta = simulationModel.numBasis(); // number of simulation basis; 'p_eta' in paper; 2 in tower example

  //***********************************************************************
  // Step 05 of 09: Instantiate the experiment storage
  // Regarding experimental scenario input values, QUESO will standardize them so that
  //    they exist inside a hypercube: this will be done in the 'uqExperimentModelClass'
  //    constructor, step 06 of 09.
  // Regarding experimental data, the user has to provide it already in transformed
  //    format, that is, with mean zero and standard deviation one.
  //***********************************************************************
  uqExperimentStorageClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>* experimentStoragePtr = NULL;
  uqExperimentModelClass  <uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>* experimentModelPtr   = NULL;

  unsigned int paper_n = 0;

  std::vector<uqGslVectorClass* >                                      experimentScenarios_original(paper_n,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* >                                      experimentScenarios_standard(paper_n,(uqGslVectorClass*) NULL);
  std::vector<unsigned int>                                            experimentDims              (paper_n,0);
  std::vector<uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>* > experimentSpaces            (paper_n,(uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>*) NULL);
  std::vector<uqGslVectorClass* >                                      extraExperimentGridVecs     (paper_n,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* >                                      experimentVecs_original     (paper_n,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* >                                      experimentVecs_auxMean      (paper_n,(uqGslVectorClass*) NULL);
  std::vector<uqGslVectorClass* >                                      experimentVecs_transformed  (paper_n,(uqGslVectorClass*) NULL);
  std::vector<uqGslMatrixClass* >                                      experimentMats_original     (paper_n,(uqGslMatrixClass*) NULL);
  std::vector<uqGslMatrixClass* >                                      experimentMats_transformed  (paper_n,(uqGslMatrixClass*) NULL);

  std::vector<uqGslMatrixClass* >                                      DobsMats                    (paper_n, (uqGslMatrixClass*) NULL);

  std::vector<uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>* > Kmats_interp_spaces         (paper_n, (uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>*) NULL);
  std::vector<uqGslMatrixClass*                                      > Kmats_interp                (paper_n, (uqGslMatrixClass*) NULL); // Interpolations of 'Kmat_eta' = 'K_i's' in paper

  if (useExperiments) {
    paper_n = 3; // number of experiments; 'n' in paper; 3 in tower example

    experimentStoragePtr = new uqExperimentStorageClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>(paper_p_x_space, paper_n);

    // Add experiments
    experimentScenarios_original.resize(paper_n,(uqGslVectorClass*) NULL);
    experimentScenarios_standard.resize(paper_n,(uqGslVectorClass*) NULL);
    experimentDims.resize              (paper_n,0);
    experimentSpaces.resize            (paper_n,(uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>*) NULL);
    extraExperimentGridVecs.resize     (paper_n,(uqGslVectorClass*) NULL);
    experimentVecs_original.resize     (paper_n,(uqGslVectorClass*) NULL);
    experimentVecs_auxMean.resize      (paper_n,(uqGslVectorClass*) NULL);
    experimentVecs_transformed.resize  (paper_n,(uqGslVectorClass*) NULL);
    experimentMats_original.resize     (paper_n,(uqGslMatrixClass*) NULL);
    experimentMats_transformed.resize  (paper_n,(uqGslMatrixClass*) NULL);

    experimentDims[0] = 4; // 'n_y_1' in paper; 4 in tower example
    experimentDims[1] = 4; // 'n_y_2' in paper; 4 in tower example
    experimentDims[2] = 3; // 'n_y_3' in paper; 3 in tower example

    for (unsigned int i = 0; i < paper_n; ++i) {
      experimentScenarios_original[i] = new uqGslVectorClass                                     (paper_p_x_space.zeroVector());               // 'x_{i+1}' in paper
      experimentScenarios_standard[i] = new uqGslVectorClass                                     (paper_p_x_space.zeroVector());               // 'x_{i+1}' in paper
      experimentSpaces            [i] = new uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>(env, "expSpace", experimentDims[i], NULL);
      extraExperimentGridVecs     [i] = new uqGslVectorClass                                     (experimentSpaces[i]->zeroVector());          //
      experimentVecs_original     [i] = new uqGslVectorClass                                     (experimentSpaces[i]->zeroVector());          //
      experimentVecs_auxMean      [i] = new uqGslVectorClass                                     (experimentSpaces[i]->zeroVector());          //
      experimentVecs_transformed  [i] = new uqGslVectorClass                                     (experimentSpaces[i]->zeroVector());          // 'y_{i+1}' in paper
      experimentMats_original     [i] = new uqGslMatrixClass                                     (experimentSpaces[i]->zeroVector());          //
      experimentMats_transformed  [i] = new uqGslMatrixClass                                     (experimentSpaces[i]->zeroVector());          // 'W_{i+1}' in paper
    }

    //***********************************************************************
    // Populate information regarding experiment '0'
    //***********************************************************************
    (*(experimentScenarios_original[0]))[0] = 0.1; // 'x_1' in paper; Radius in tower example
    *(experimentScenarios_standard[0])  = *(experimentScenarios_original[0]);
    *(experimentScenarios_standard[0]) -= simulationModel.xSeq_original_mins();
    for (unsigned int j = 0; j < paper_p_x; ++j) {
      (*(experimentScenarios_standard[0]))[j] /= simulationModel.xSeq_original_ranges()[j];
    }

    (*(extraExperimentGridVecs[0]))[0] = 5.;
    (*(extraExperimentGridVecs[0]))[1] = 10.;
    (*(extraExperimentGridVecs[0]))[2] = 15.;
    (*(extraExperimentGridVecs[0]))[3] = 20.;

    (*(experimentVecs_original[0]))[0] = 1.2180;
    (*(experimentVecs_original[0]))[1] = 2.0126;
    (*(experimentVecs_original[0]))[2] = 2.7942;
    (*(experimentVecs_original[0]))[3] = 3.5747;

    experimentVecs_auxMean[0]->matlabLinearInterpExtrap(extraSimulationGridVec,simulationModel.etaSeq_original_mean(),*(extraExperimentGridVecs[0]));
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      *env.subDisplayFile() << "In compute(), step 05, experiment 0:"
                            << "\n  extraSimulationGridVec = "                 << extraSimulationGridVec
                            << "\n  simulationModel.etaSeq_original_mean() = " << simulationModel.etaSeq_original_mean()
                            << "\n  *(extraExperimentGridVecs[0]) = "          << *(extraExperimentGridVecs[0])
                            << "\n  *(experimentVecs_auxMean[0]) = "           << *(experimentVecs_auxMean[0])
                            << std::endl;
    }
    *(experimentVecs_transformed[0]) = (1./simulationModel.etaSeq_allStd()) * ( *(experimentVecs_original[0]) - *(experimentVecs_auxMean[0]) ); // 'y_1' in paper
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      *env.subDisplayFile() << "In compute(), step 05, experiment 0:"
                            << "\n  *(experimentVecs_original[0]) = "    << *(experimentVecs_original[0])
                            << "\n  *(experimentVecs_auxMean[0]) = "     << *(experimentVecs_auxMean[0])
                            << "\n  simulationModel.etaSeq_allStd() = "  << simulationModel.etaSeq_allStd()
                            << "\n  *(experimentVecs_transformed[0]) = " << *(experimentVecs_transformed[0])
                            << std::endl;
    }

    (*(experimentMats_original[0]))(0,0) = 1.;
    (*(experimentMats_original[0]))(1,1) = 1.;
    (*(experimentMats_original[0]))(2,2) = 1.;
    (*(experimentMats_original[0]))(3,3) = 1.;

    *(experimentMats_transformed[0]) = *(experimentMats_original[0]);

    //***********************************************************************
    // Populate information regarding experiment '1'
    //***********************************************************************
    (*(experimentScenarios_original[1]))[0] = 0.2; // 'x_2' in paper; Radius in tower example
    *(experimentScenarios_standard[1])  = *(experimentScenarios_original[1]);
    *(experimentScenarios_standard[1]) -= simulationModel.xSeq_original_mins();
    for (unsigned int j = 0; j < paper_p_x; ++j) {
      (*(experimentScenarios_standard[1]))[j] /= simulationModel.xSeq_original_ranges()[j];
    }

    (*(extraExperimentGridVecs[1]))[0] = 5.;
    (*(extraExperimentGridVecs[1]))[1] = 10.;
    (*(extraExperimentGridVecs[1]))[2] = 15.;
    (*(extraExperimentGridVecs[1]))[3] = 20.;

    (*(experimentVecs_original[1]))[0] = 1.1129;
    (*(experimentVecs_original[1]))[1] = 1.7225;
    (*(experimentVecs_original[1]))[2] = 2.2898;
    (*(experimentVecs_original[1]))[3] = 2.8462;

    experimentVecs_auxMean[1]->matlabLinearInterpExtrap(extraSimulationGridVec,simulationModel.etaSeq_original_mean(),*(extraExperimentGridVecs[1]));
    *(experimentVecs_transformed[1]) = (1./simulationModel.etaSeq_allStd()) * ( *(experimentVecs_original[1]) - *(experimentVecs_auxMean[1]) ); // 'y_2' in paper

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      *env.subDisplayFile() << "In compute(), step 05, experiment 1:"
                            << "\n  *(experimentVecs_original[1]) = "    << *(experimentVecs_original[1])
                            << "\n  *(experimentVecs_auxMean[1]) = "     << *(experimentVecs_auxMean[1])
                            << "\n  simulationModel.etaSeq_allStd() = "  << simulationModel.etaSeq_allStd()
                            << "\n  *(experimentVecs_transformed[1]) = " << *(experimentVecs_transformed[1])
                            << std::endl;
    }

    (*(experimentMats_original[1]))(0,0) = 1.;
    (*(experimentMats_original[1]))(1,1) = 1.;
    (*(experimentMats_original[1]))(2,2) = 1.;
    (*(experimentMats_original[1]))(3,3) = 1.;

    *(experimentMats_transformed[1]) = *(experimentMats_original[1]);

    //***********************************************************************
    // Populate information regarding experiment '2'
    //***********************************************************************
    (*(experimentScenarios_original[2]))[0] = 0.4; // 'x_3' in paper; Radius in tower example
    *(experimentScenarios_standard[2])  = *(experimentScenarios_original[2]);
    *(experimentScenarios_standard[2]) -= simulationModel.xSeq_original_mins();
    for (unsigned int j = 0; j < paper_p_x; ++j) {
      (*(experimentScenarios_standard[2]))[j] /= simulationModel.xSeq_original_ranges()[j];
    }

    (*(extraExperimentGridVecs[2]))[0] = 5.;
    (*(extraExperimentGridVecs[2]))[1] = 10.;
    (*(extraExperimentGridVecs[2]))[2] = 15.;

    (*(experimentVecs_original[2]))[0] = 1.0611;
    (*(experimentVecs_original[2]))[1] = 1.5740;
    (*(experimentVecs_original[2]))[2] = 2.0186;

    experimentVecs_auxMean[2]->matlabLinearInterpExtrap(extraSimulationGridVec,simulationModel.etaSeq_original_mean(),*(extraExperimentGridVecs[2]));
    *(experimentVecs_transformed[2]) = (1./simulationModel.etaSeq_allStd()) * ( *(experimentVecs_original[2]) - *(experimentVecs_auxMean[2]) ); // 'y_3' in paper

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      *env.subDisplayFile() << "In compute(), step 05, experiment 2:"
                            << "\n  *(experimentVecs_original[2]) = "    << *(experimentVecs_original[2])
                            << "\n  *(experimentVecs_auxMean[2]) = "     << *(experimentVecs_auxMean[2])
                            << "\n  simulationModel.etaSeq_allStd() = "  << simulationModel.etaSeq_allStd()
                            << "\n  *(experimentVecs_transformed[2]) = " << *(experimentVecs_transformed[2])
                            << std::endl;
    }

    (*(experimentMats_original[2]))(0,0) = 1.;
    (*(experimentMats_original[2]))(1,1) = 1.;
    (*(experimentMats_original[2]))(2,2) = 1.;

    *(experimentMats_transformed[2]) = *(experimentMats_original[2]);

    //***********************************************************************
    // Finally, add information to the experiment storage
    //***********************************************************************
    for (unsigned int i = 0; i < paper_n; ++i) {
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
        *env.subDisplayFile() << "In compute(), step 05"
                              << ": calling experimentStoragePtr->addExperiment() for experiment of id '" << i << "'..."
                              << std::endl;
      }
      experimentStoragePtr->addExperiment(*(experimentScenarios_standard[i]),*(experimentVecs_transformed[i]),*(experimentMats_transformed[i]));
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
        *env.subDisplayFile() << "In compute(), step 05"
                              << ": returned from experimentStoragePtr->addExperiment() for experiment of id '" << i << "'"
                              << std::endl;
      }
    }


    //***********************************************************************
    // Step 06 of 09: Instantiate the experiment model
    // User has to provide 'D' matrices
    // User has to interpolate 'K_eta' matrix in order to form 'K_i' matrices
    // 'K_eta' is 'Ksim' in the GPMSA tower example document (page 9)
    // 'K_i' is 'Kobs' in the GPMSA tower example document (page 9)
    //***********************************************************************
    unsigned int paper_p_delta = 13; // number of experiment basis; 'p_delta' in paper; 13 in tower example

    //***********************************************************************
    // Form and compute 'DsimMat'
    // Not mentioned in the paper
    // 'Dsim' in the GPMSA tower example document (page 11)
    //***********************************************************************
    double kernelSigma = 2.;
    uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paper_p_delta_space(env, "paper_p_delta_space", paper_p_delta, NULL);
    uqGslVectorClass kernelCenters(paper_p_delta_space.zeroVector());
    for (unsigned int i = 1; i < paper_p_delta; ++i) { // Yes, '1'
      kernelCenters[i] = kernelSigma*((double) i);
    }
    uqGslMatrixClass DsimMat(env,paper_n_eta_space.map(),paper_p_delta); // Important matrix (not mentioned on paper)
    uqGslVectorClass DsimCol(paper_n_eta_space.zeroVector());
    uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> oneDSpace(env, "oneDSpace", 1, NULL);
    uqGslVectorClass oneDVec(oneDSpace.zeroVector());
    uqGslMatrixClass oneDMat(oneDSpace.zeroVector());
    oneDMat(0,0) = kernelSigma*kernelSigma;
    for (unsigned int colId = 0; colId < DsimMat.numCols(); ++colId) {
      oneDVec[0] = kernelCenters[colId];
      uqGaussianJointPdfClass<uqGslVectorClass,uqGslMatrixClass> kernelPdf("",oneDSpace,oneDVec,oneDMat);
      for (unsigned int rowId = 0; rowId < paper_n_eta; ++rowId) {
        oneDVec[0] = extraSimulationGridVec[rowId];
        DsimCol[rowId] = kernelPdf.actualValue(oneDVec,NULL,NULL,NULL,NULL);
      }
      DsimMat.setColumn(colId,DsimCol);
    }

    //***********************************************************************
    // Populate information regarding experiment 'i'
    // 'D_{i+1}' in the paper
    // 'Dobs' in the GPMSA tower example document (page 11)
    //***********************************************************************
    DobsMats.resize(paper_n, (uqGslMatrixClass*) NULL); // Important matrices (D_i's on paper)
    for (unsigned int i = 0; i < paper_n; ++i) {
      DobsMats[i] = new uqGslMatrixClass(env,experimentSpaces[i]->map(),paper_p_delta); // 'D_{i+1}' in paper
      uqGslVectorClass DobsCol(experimentSpaces[i]->zeroVector());
      for (unsigned int colId = 0; colId < DobsMats[i]->numCols(); ++colId) {
        oneDVec[0] = kernelCenters[colId];
        uqGaussianJointPdfClass<uqGslVectorClass,uqGslMatrixClass> kernelPdf("",oneDSpace,oneDVec,oneDMat);
        for (unsigned int rowId = 0; rowId < DobsCol.sizeLocal(); ++rowId) {
          oneDVec[0] = (*(extraExperimentGridVecs[1]))[rowId];
          DobsCol[rowId] = kernelPdf.actualValue(oneDVec,NULL,NULL,NULL,NULL);
        }
        DobsMats[i]->setColumn(colId,DobsCol);
      }
    }
  
    //***********************************************************************
    // Normalize 'DsimMat' and all 'DobsMats'
    //***********************************************************************
    uqGslMatrixClass DsimMatTranspose(env,paper_p_delta_space.map(),paper_n_eta);
    DsimMatTranspose.fillWithTranspose(0,0,DsimMat,true,true);
    double Dmax = (DsimMat * DsimMatTranspose).max();
    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      *env.subDisplayFile() << "In compute()"
                            << ": Dmax = " << Dmax
                            << std::endl;
    }
    DsimMat /= std::sqrt(Dmax);

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
      DsimMat.setPrintHorizontally(false);
      *env.subDisplayFile() << "In compute()"
                            << ": 'DsimMat'"
                            << ", nRows = "      << DsimMat.numRowsLocal()
                            << ", nCols = "      << DsimMat.numCols()
                            << ", contents =\n " << DsimMat
                            << std::endl;
    }
    for (unsigned int i = 0; i < paper_n; ++i) {
      *(DobsMats[i]) /= std::sqrt(Dmax);
      DobsMats[i]->setPrintHorizontally(false);
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
        *env.subDisplayFile() << "In compute()"
                              << ": 'DobsMats', i = " << i
                              << ", nRows = "         << DobsMats[i]->numRowsLocal()
                              << ", nCols = "         << DobsMats[i]->numCols()
                              << ", contents =\n"     << *(DobsMats[i])
                              << std::endl;
      }
    }

    //***********************************************************************
    // Compute 'K_i' matrices
    // 'Kobs' in the GPMSA tower example document (page 9)
    //***********************************************************************
    Kmats_interp_spaces.resize(paper_n, (uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>*) NULL);
    Kmats_interp.resize       (paper_n, (uqGslMatrixClass*) NULL); // Important matrices (K_i's on paper)
    for (unsigned int i = 0; i < paper_n; ++i) {
      Kmats_interp_spaces[i] = new uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>(env,"Kmats_interp_spaces_",experimentStoragePtr->n_ys_transformed()[i],NULL);
      Kmats_interp       [i] = new uqGslMatrixClass(env,Kmats_interp_spaces[i]->map(),paper_p_eta);
      Kmats_interp       [i]->matlabLinearInterpExtrap(extraSimulationGridVec,simulationModel.Kmat_eta(),*(extraExperimentGridVecs[i])); // Important matrix (K_eta on paper) 
      Kmats_interp[i]->setPrintHorizontally(false);
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
        *env.subDisplayFile() << "In compute()"
                              << ": 'Kmats_interp', i = " << i
                              << ", nRows = "             << Kmats_interp[i]->numRowsLocal()
                              << ", nCols = "             << Kmats_interp[i]->numCols()
                              << ", contents =\n"         << *(Kmats_interp[i])
                              << std::endl;
      }
    }

    if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
      *env.subDisplayFile() << "In compute()"
                            << ": finished computing 'K_i' matrices"
                            << std::endl;
    }

    //***********************************************************************
    // The constructor below will read from the options file:
    // --> the 'G' values (a total of 'F' of them); 'Gs' in paper; G[0] = 13 in tower example
    // --> 'F' in paper = number of groups of experiment basis; 1 in tower example
    //***********************************************************************
    experimentModelPtr = new uqExperimentModelClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>("",   // prefix
                                                                                                                         NULL, // options values
                                                                                                                         *experimentStoragePtr,
                                                                                                                         DobsMats,
                                                                                                                         Kmats_interp);

    paramPriorRvPtr = new uqUniformVectorRVClass<uqGslVectorClass,uqGslMatrixClass>("prior_",paramDomain);
  //uqGslVectorClass meanVec(paper_p_t_space.zeroVector());
  //meanVec[0] = 0.5;
  //uqGslMatrixClass covaMat(paper_p_t_space.zeroVector());
  //covaMat(0,0) = 100.;
  //paramPriorRvPtr = new uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>("prior_",paramDomain/*paper_p_t_space*/,meanVec,covaMat);
  } // if (useExperiments)

  //***********************************************************************
  // Step 07 of 09: Instantiate the GPMSA computer model
  //***********************************************************************
  uqGpmsaComputerModelClass<uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass,uqGslVectorClass,uqGslMatrixClass>* gcm;
  gcm = new uqGpmsaComputerModelClass<uqGslVectorClass,uqGslMatrixClass,
                                      uqGslVectorClass,uqGslMatrixClass,
                                      uqGslVectorClass,uqGslMatrixClass,
                                      uqGslVectorClass,uqGslMatrixClass>
          ("",
           NULL,
           simulationStorage,
           simulationModel,
           experimentStoragePtr, // pass "NULL" if there is no experimental data available
           experimentModelPtr,   // pass "NULL" if there is no experimental data available
           paramPriorRvPtr);     // pass "NULL" if there is no experimental data available

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In compute()"
                          << ": finished instantiating 'gcm'"
                          << ", gcm->totalSpace().dimLocal() = " << gcm->totalSpace().dimLocal()
                          << std::endl;
  }

  //***********************************************************************
  // Step 08 of 09: Calibrate the computer model
  //***********************************************************************
  uqGslVectorClass totalInitialVec(gcm->totalSpace().zeroVector());
  gcm->totalPriorRv().realizer().realization(totalInitialVec); // todo_r

  uqGslVectorClass diagVec(gcm->totalSpace().zeroVector());
  diagVec.cwSet(0.25);

  if (useExperiments) {
#if 1
    // Use same initial values as matlab code, for debug
    totalInitialVec[ 0] = 2.9701e+04;          // lambda_eta = lamWOs
    totalInitialVec[ 1] = 1.;                  // lambda_w_1 = lamUz
    totalInitialVec[ 2] = 1.;                  // lambda_w_2 =
    totalInitialVec[ 3] = std::exp(-0.1*0.25); // rho_w_1_1  = exp(-model.betaU.*(0.5^2));
    totalInitialVec[ 4] = std::exp(-0.1*0.25); // rho_w_1_2  =
    totalInitialVec[ 5] = std::exp(-0.1*0.25); // rho_w_2_1  =
    totalInitialVec[ 6] = std::exp(-0.1*0.25); // rho_w_2_2  =
    totalInitialVec[ 7] = 1000.;               // lambda_s_1 = lamWs
    totalInitialVec[ 8] = 1000.;               // lambda_s_2 =
    totalInitialVec[ 9] = 999.999;             // lambda_y   = lamOs
    totalInitialVec[10] = 20.;                 // lambda_v_1 = lamVz
    totalInitialVec[11] = std::exp(-0.1*0.25); // rho_v_1_1  = betaV
    totalInitialVec[12] = 0.5;                 // theta_1    = theta
#endif

    diagVec[ 0] = 2500.; // lambda_eta = lamWOs
    diagVec[ 1] = 0.09;  // lambda_w_1 = lamUz
    diagVec[ 2] = 0.09;  // lambda_w_2 =
    diagVec[ 3] = 0.01;  // rho_w_1_1  = betaU
    diagVec[ 4] = 0.01;  // rho_w_1_2  =
    diagVec[ 5] = 0.01;  // rho_w_2_1  =
    diagVec[ 6] = 0.01;  // rho_w_2_2  =
    diagVec[ 7] = 2500.; // lambda_s_1 = lamWs
    diagVec[ 8] = 2500;  // lambda_s_2 = 
    diagVec[ 9] = 2500;  // lambda_y   = lamOs
    diagVec[10] = 2500;  // lambda_v_1 = lamVz
    diagVec[11] = 0.01;  // rho_v_1_1  = betaV
    diagVec[12] = 100.;  // theta_1    = theta
  }
  else {
#if 1
    // Use same initial values as matlab code, for debug
    totalInitialVec[ 0] = 2.9701e+04;          // lambda_eta = lamWOs
    totalInitialVec[ 1] = 1.;                  // lambda_w_1 = lamUz
    totalInitialVec[ 2] = 1.;                  // lambda_w_2 =
    totalInitialVec[ 3] = std::exp(-0.1*0.25); // rho_w_1_1  = exp(-model.betaU.*(0.5^2));
    totalInitialVec[ 4] = std::exp(-0.1*0.25); // rho_w_1_2  =
    totalInitialVec[ 5] = std::exp(-0.1*0.25); // rho_w_2_1  =
    totalInitialVec[ 6] = std::exp(-0.1*0.25); // rho_w_2_2  =
    totalInitialVec[ 7] = 1000.;               // lambda_s_1 = lamWs
    totalInitialVec[ 8] = 1000.;               // lambda_s_2 =
#endif

    diagVec[ 0] = 2500.; // lambda_eta = lamWOs
    diagVec[ 1] = 0.09;  // lambda_w_1 = lamUz
    diagVec[ 2] = 0.09;  // lambda_w_2 =
    diagVec[ 3] = 0.01;  // rho_w_1_1  = betaU
    diagVec[ 4] = 0.01;  // rho_w_1_2  =
    diagVec[ 5] = 0.01;  // rho_w_2_1  =
    diagVec[ 6] = 0.01;  // rho_w_2_2  =
    diagVec[ 7] = 2500.; // lambda_s_1 = lamWs
    diagVec[ 8] = 2500;  // lambda_s_2 = 
  }

  uqGslMatrixClass totalInitialProposalCovMatrix(diagVec); // todo_r

  if (useML) {
    gcm->calibrateWithBayesMLSampling();
  }
  else {
    gcm->calibrateWithBayesMetropolisHastings(NULL,totalInitialVec,&totalInitialProposalCovMatrix);
  }

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In compute()"
                          << ": finished calibrating 'gcm'"
                          << std::endl;
  }

  //***********************************************************************
  // Step 09 of 09: Make predictions with the calibrated computer model
  // Substep 09a: predict 'v' and 'u' weights
  // Substep 09b: predict 'w' weights
  // Substep 09c: predict new experiment results
  // Substep 09d: predict new simulation outputs
  //***********************************************************************
  unsigned int dim1 = 16;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> predictionGrid1Space(env, "predictionGrid1_", dim1, NULL);
  uqGslVectorClass predictionGrid1Vec(predictionGrid1Space.zeroVector());
  for (unsigned int i = 0; i < dim1; ++i) {
    predictionGrid1Vec[i] = ((double) i)/((double) (dim1-1));
  }
  
  unsigned int dim2 = 16;
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> predictionGrid2Space(env, "predictionGrid2_", dim2, NULL);
  uqGslVectorClass predictionGrid2Vec(predictionGrid2Space.zeroVector());
  for (unsigned int j = 0; j < dim2; ++j) {
    predictionGrid2Vec[j] = ((double) j)/((double) (dim2-1));
  }

  if (env.fullRank() == 0) { // Yes, only one processor in all 'full' communicator
    std::set<unsigned int> tmpSet;
    tmpSet.insert(0);

    predictionGrid1Vec.subWriteContents("towerPredictionGrid1", // varNamePrefix,
                                        "predictionGrid1",
                                        "m",
                                        tmpSet); // allowedSubEnvIds
    predictionGrid2Vec.subWriteContents("towerPredictionGrid2", // varNamePrefix,
                                        "predictionGrid2",
                                        "m",
                                        tmpSet); // allowedSubEnvIds
  }

  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paper_p_eta_space(env, "paper_p_eta_", paper_p_eta, NULL); // = gcm->unique_u_space() or gcm->unique_w_space()

  if (useExperiments) {
#if 1 // Might put '0' while code is not finished
    //***********************************************************************
    // Substep 09a of 09: Predict 'v' and 'u' weights
    //***********************************************************************
    uqGslVectorClass tmpExperimentScenarioVec (paper_p_x_space.zeroVector());
    uqGslVectorClass tmpSimulationParameterVec(paper_p_t_space.zeroVector());

    unsigned int paper_p_delta = 13; // todo_rr: repeated
    uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass> paper_p_delta_space(env, "paper_p_delta_space", paper_p_delta, NULL); // todo_rr: repeated // = gcm_unique_v_space()

    std::vector<uqGslMatrixClass* > predictionGridVMats(paper_p_delta,NULL);
    for (unsigned int k = 0; k < paper_p_delta; ++k) {
      predictionGridVMats[k] = new uqGslMatrixClass(predictionGrid1Space.zeroVector());
    }
    std::vector<uqGslMatrixClass* > predictionGridUMats(paper_p_eta,NULL);
    for (unsigned int k = 0; k < paper_p_eta; ++k) {
      predictionGridUMats[k] = new uqGslMatrixClass(predictionGrid1Space.zeroVector());
    }

    uqGslVectorClass vuMeanVec  (gcm->unique_vu_space().zeroVector());
    uqGslMatrixClass vuCovMatrix(gcm->unique_vu_space().zeroVector());
    uqGslVectorClass vMeanVec   (paper_p_delta_space.zeroVector());
    uqGslMatrixClass vCovMatrix (paper_p_delta_space.zeroVector());
    uqGslVectorClass uMeanVec   (paper_p_eta_space.zeroVector());
    uqGslMatrixClass uCovMatrix (paper_p_eta_space.zeroVector());
    for (unsigned int i = 0; i < dim1; ++i) {
      tmpExperimentScenarioVec[0] = predictionGrid1Vec[i];
      for (unsigned int j = 0; j < dim2; ++j) {
        tmpSimulationParameterVec[0] = predictionGrid2Vec[j];
#if 1
        gcm->predictVUsAtGridPoint(tmpExperimentScenarioVec,
                                   tmpSimulationParameterVec,
                                   vuMeanVec,
                                   vuCovMatrix,
                                   vMeanVec,
                                   vCovMatrix,
                                   uMeanVec,
                                   uCovMatrix);
#endif
        for (unsigned int k = 0; k < paper_p_delta; ++k) {
          (*predictionGridVMats[k])(i,j) = vMeanVec[k];
        }

        for (unsigned int k = 0; k < paper_p_eta; ++k) {
          (*predictionGridUMats[k])(i,j) = uMeanVec[k];
        }
      }
    }

    if (env.fullRank() == 0) { // Yes, only one processor in all 'full' communicator
      std::set<unsigned int> tmpSet;
      tmpSet.insert(0);

      char varNamePrefix[32+1];
      char fileName     [32+1];

      for (unsigned int i = 0; i < paper_p_delta; ++i) {
        sprintf(varNamePrefix,"towerV%dmat",i+1);
        sprintf(fileName,     "v%dmat",     i+1);
        predictionGridVMats[i]->subWriteContents(varNamePrefix,
                                                 fileName,
                                                 "m",
                                                 tmpSet); // allowedSubEnvIds
      }

      for (unsigned int i = 0; i < paper_p_eta; ++i) {
        sprintf(varNamePrefix,"towerU%dmat",i+1);
        sprintf(fileName,     "u%dmat",     i+1);
        predictionGridUMats[i]->subWriteContents(varNamePrefix,
                                                 fileName,
                                                 "m",
                                                 tmpSet); // allowedSubEnvIds
      }
    }

    for (unsigned int i = 0; i < paper_p_eta; ++i) {
      delete predictionGridUMats[i];
    }
    for (unsigned int i = 0; i < paper_p_delta; ++i) {
      delete predictionGridVMats[i];
    }
#endif // Might put '0'
  } // if (useExperiments)

#if 1 // Might put '0' while code is not finished
  //***********************************************************************
  // Substep 09b of 09: Predict 'w' weights
  //***********************************************************************
  uqGslVectorClass tmpSimulationScenarioVec (paper_p_x_space.zeroVector());
  uqGslVectorClass tmpSimulationParameterVec(paper_p_t_space.zeroVector());
  uqGslVectorClass wMeanVec  (paper_p_eta_space.zeroVector());
  uqGslMatrixClass wCovMatrix(paper_p_eta_space.zeroVector());

#if 0 // For debug only
  tmpSimulationScenarioVec [0] = 0.4;
  tmpSimulationParameterVec[0] = 0.4;

  uqGslVectorClass forcingSampleVecForDebug(gcm->totalSpace().zeroVector());
  forcingSampleVecForDebug[ 0] = 2.9700e+04;             // lambda_eta = lamWOs
  forcingSampleVecForDebug[ 1] = 1.0000;                 // lambda_w_1 = lamUz
  forcingSampleVecForDebug[ 2] = 1.0083;                 // lambda_w_2 =
  forcingSampleVecForDebug[ 3] = std::exp(-0.1004*0.25); // rho_w_1_1  = exp(-model.betaU.*(0.5^2));
  forcingSampleVecForDebug[ 4] = std::exp(-0.1022*0.25); // rho_w_1_2  =
  forcingSampleVecForDebug[ 5] = std::exp(-0.1013*0.25); // rho_w_2_1  =
  forcingSampleVecForDebug[ 6] = std::exp(-0.1011*0.25); // rho_w_2_2  =
  forcingSampleVecForDebug[ 7] = 999.5113;               // lambda_s_1 = lamWs
  forcingSampleVecForDebug[ 8] = 999.3138;               // lambda_s_2 =
  forcingSampleVecForDebug[ 9] = 996.6176;               // lambda_y   = lamOs
  forcingSampleVecForDebug[10] = 20.0194;                // lambda_v_1 = lamVz
  forcingSampleVecForDebug[11] = std::exp(-0.1020*0.25); // rho_v_1_1  = betaV
  forcingSampleVecForDebug[12] = 0.5000;                 // theta_1    = theta

  gcm->predictWsAtGridPoint(tmpSimulationScenarioVec,
                            tmpSimulationParameterVec,
                            &forcingSampleVecForDebug,
                            wMeanVec,
                            wCovMatrix);
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In compute()"
                          << ": for scenario = " << tmpSimulationScenarioVec [0]
                          << " and parameter = " << tmpSimulationParameterVec[0]
                          << ", wMeanVec = "     << wMeanVec
                          << ", wCovMatrix = "   << wCovMatrix
                          << std::endl;
  }
  env.fullComm().Barrier();
  sleep(1);
  queso_error();
  // exit(1);
#endif

  std::vector<uqGslMatrixClass* > predictionGridWMats(paper_p_eta,NULL);
  for (unsigned int k = 0; k < paper_p_eta; ++k) {
    predictionGridWMats[k] = new uqGslMatrixClass(env,predictionGrid1Space.map(),dim2);
  }

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In compute()"
                          << ": beginning loop on predictions"
                          << ", dim1 = " << dim1
                          << ", dim2 = " << dim2
                          << std::endl;
  }

  for (unsigned int i = 0; i < dim1; ++i) {
    tmpSimulationScenarioVec[0] = predictionGrid1Vec[i];
    for (unsigned int j = 0; j < dim2; ++j) {
      tmpSimulationParameterVec[0] = predictionGrid2Vec[j];

      struct timeval timevalBeforePredict;
      int iRC = gettimeofday(&timevalBeforePredict, NULL);
      if (iRC) {}; // just to remove compiler warning

      gcm->predictWsAtGridPoint(tmpSimulationScenarioVec,
                                tmpSimulationParameterVec,
                                NULL,
                                wMeanVec,
                                wCovMatrix);

      double predictTime = uqMiscGetEllapsedSeconds(&timevalBeforePredict);
      if ((env.subDisplayFile()) && (env.displayVerbosity() >= 3)) {
        *env.subDisplayFile() << "In compute()"
                              << ": looping on predictions"
                              << ", i = " << i << " < " << dim1
                              << ", j = " << j << " < " << dim2
                              << ", returned from 'gcm->predictWsAtGridPoint()' after " << predictTime << " seconds"
                              << std::endl;
      }

      for (unsigned int k = 0; k < paper_p_eta; ++k) {
        (*predictionGridWMats[k])(i,j) = wMeanVec[k];
      }
    }
  }

  if (env.fullRank() == 0) { // Yes, only one processor in all 'full' communicator
    std::set<unsigned int> tmpSet;
    tmpSet.insert(0);

    char varNamePrefix[32+1];
    char fileName     [32+1];
    for (unsigned int i = 0; i < paper_p_eta; ++i) {
      sprintf(varNamePrefix,"towerW%dmat",i+1);
      sprintf(fileName,     "w%dmat",     i+1);
      predictionGridWMats[i]->subWriteContents(varNamePrefix,
                                               fileName,
                                               "m",
                                               tmpSet); // allowedSubEnvIds
    }
  }

  for (unsigned int i = 0; i < paper_p_eta; ++i) {
    delete predictionGridWMats[i];
  }

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "In compute()"
                          << ": finished making predictions"
                          << std::endl;
  }
#endif // Might put '0'

#if 0 // Might put '0' while code is not finished
  //***********************************************************************
  // Substep 09c of 09: Predict new experiment results
  //***********************************************************************
  uqGslVectorClass newExperimentScenarioVec(paper_p_x_space.zeroVector());              // todo_rr
  uqGslMatrixClass newKmat_interp          (env,paper_n_eta_space.map(),paper_p_eta);   // todo_rr
  uqGslMatrixClass newDmat                 (env,paper_n_eta_space.map(),paper_p_delta); // todo_rr
  uqGslVectorClass simulationOutputMeanVec (paper_n_eta_space.zeroVector()); // Yes, size of simulation, since it is a prediction using the emulator
  uqGslVectorClass discrepancyMeanVec      (paper_n_eta_space.zeroVector());
  gcm->predictExperimentResults(newExperimentScenarioVec,newKmat_interp,newDmat,simulationOutputMeanVec,discrepancyMeanVec);

  //***********************************************************************
  // Substep 09d of 09: Predict new simulation outputs
  //***********************************************************************
  uqGslVectorClass newSimulationScenarioVec (paper_p_x_space.zeroVector  ()); // todo_rr
  uqGslVectorClass newSimulationParameterVec(paper_p_t_space.zeroVector  ()); // todo_rr
  uqGslVectorClass simulationOutputMeanVec2 (paper_n_eta_space.zeroVector());
  gcm->predictSimulationOutputs(newSimulationScenarioVec,newSimulationParameterVec,simulationOutputMeanVec2);
#endif // Might put '0'

  //***********************************************************************
  // Clean memory
  //***********************************************************************
  env.fullComm().Barrier();
  delete paramPriorRvPtr;
  delete experimentModelPtr;
  delete experimentStoragePtr;
  for (unsigned int i = 0; i < paper_n; ++i) {
    delete Kmats_interp       [i];
    delete Kmats_interp_spaces[i];
  }
  for (unsigned int i = 0; i < paper_m; ++i) {
    delete outputVecs         [i];
    delete paramVecs          [i];
    delete simulationScenarios[i];
  }
  for (unsigned int i = 0; i < paper_n; ++i) {
    delete experimentMats_transformed  [i];
    delete experimentMats_original     [i];
    delete experimentVecs_transformed  [i];
    delete experimentVecs_original     [i];
    delete experimentSpaces            [i];
    delete experimentScenarios_standard[i];
    delete experimentScenarios_original[i];
  }
  delete gcm;

  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Leaving compute()"
                          << std::endl;
  }

  return;
}
