#ifndef __UQ_MISCELLANEOUS_H__
#define __UQ_MISCELLANEOUS_H__

#include <gsl/gsl_rng.h>
#include <vector>

void   uqMiscReadDoublesFromString (const std::string&   inputString,
                                    std::vector<double>& outputDoubles);
void   uqMiscReadWordsFromString   (const std::string&        inputString,
                                    std::vector<std::string>& outputWords);
double uqMiscGammar                (double   a,
                                    double   b,
                                    gsl_rng* rng);

#endif // __UQ_MISCELLANEOUS_H__
