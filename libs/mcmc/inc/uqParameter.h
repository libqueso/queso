#ifndef __UQ_PARAMETER_H__
#define __UQ_PARAMETER_H__

#include <string>
#include <iostream>
#include <math.h>

class uqParameterClass
{
public:
  uqParameterClass(const std::string& name,
                   double             initialValue,
                   double             minValue = -INFINITY,
                   double             maxValue = INFINITY,
                   double             priorMu = 0.,
                   double             priorSigma = INFINITY);
 ~uqParameterClass();

  std::string name           () const;
  double      initialValue   () const;
  double      minValue       () const;
  double      maxValue       () const;
  double      priorMu        () const;
  double      priorSigma     () const;

  void        setName        (const std::string& name);
  void        setInitialValue(double initialValue);
  void        setMinValue    (double minValue);
  void        setMaxValue    (double maxValue);
  void        setPriorMu     (double priorMu);
  void        setPriorSigma  (double priorSigma);

  void        print          (std::ostream& os) const;

private:
  std::string m_name;
  double      m_initialValue;
  double      m_minValue;
  double      m_maxValue;
  double      m_priorMu;
  double      m_priorSigma;
};

std::ostream& operator<<(std::ostream& os, const uqParameterClass& param);

#endif // __UQ_PARAMETER_H__
