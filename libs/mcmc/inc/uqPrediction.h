/* uq/libs/mcmc/inc/uqPrediction.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos

 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_PREDICTION_H__
#define __UQ_PREDICTION_H__

#include <uqPriorAndLikelihoodInterfaces.h>

template<class V, class M>
void
uqPrediction(
  const std::vector<const V*>&         chain,
  const std::vector<const V*>&         lrSigma2Chain,
  unsigned int                         numberOfSamples,
  tmp_QoIPredictionFunction_Class<V,M>& qoIPredictionFunction_Obj,
  const void*                          functionDataPtr,
  const std::string&                   outputFileName)
{
  std::vector<unsigned int> sampleIds(0,0);
  if (numberOfSamples < chain.size()) {
    sampleIds.resize(numberOfSamples,0);
  }
  else {
    sampleIds.resize(chain.size(),0);
  }
  gsl_rng* rng = chain[0]->env().rng();
  for (unsigned int i = 0; i < sampleIds.size(); ++i) {
    sampleIds[i] = (unsigned int) ( ceil(gsl_rng_uniform(rng) * chain.size()) - 1. );
  }

  std::vector<V*> predictions(0);//,NULL);
  const V* parameters;
  for (unsigned int i = 0; i < sampleIds.size(); ++i) {
    parameters = chain[sampleIds[i]];
    qoIPredictionFunction_Obj(*parameters,functionDataPtr,predictions);

    std::vector<V*> perturbedPredicitions(predictions.size(),NULL);
    unsigned int numTimeSteps  = predictions.size();
    unsigned int paramSpaceDim = predictions[0]->size();
    V            lrStdDev(*lrSigma2Chain[sampleIds[i]]);
    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      perturbedPredicitions[j] = new V(*predictions[j]);
      for (unsigned int k = 0; k < paramSpaceDim; ++k) {
        (*perturbedPredicitions[j])[k] = gsl_ran_gaussian(rng,lrStdDev[k]);
      }
    }

    for (unsigned int j = 0; j < numTimeSteps; ++j) {
      if (perturbedPredicitions[j]) delete perturbedPredicitions[j];
    }
    perturbedPredicitions.clear();
  }

#if 0
if ~iscell(data)
  d=data; data=[]; data{1}=d; clear d
end

  % random sample from the chain
  isample = ceil(rand(nsample,1)*nsimu);

  % modelfun   = @(d,th) algaefun(d(:,1),th   ,th(7:9),d    );
  % function y =         algaefun(time  ,theta,y0     ,xdata);
  % [t,y]      = ode15s(@algaesys,time,        y0,[],theta,xdata);
  % y = [num time steps, num state variables]
  % s2chain = [chain size, num output variables]

  for iisample = 1:nsample;
//theta(parind) = chain(isample(iisample),:)';
    th  = theta(local==0|local==i);
    y   = feval(modelfun,datai,th,varargin{:});
    ysave(iisample,:,:) = y;
    osave(iisample,:,:) = y + randn(size(y))*diag(sqrt(s2chain(isample(iisample),:)));
  end

  [my,ny] = size(y);
  for j=1:ny
    plim{i}{j} = plims(ysave(:,:,j),lims);
    olim{i}{j} = plims(osave(:,:,j),lims);
  end

out.predlims = plim;
out.obslims = olim;
out.data = data;
#endif

  return;
}
#endif // __UQ_PREDICTION_H__
