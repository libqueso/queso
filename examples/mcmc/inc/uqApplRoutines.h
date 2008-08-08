#ifndef __UQ_APPL_ROUTINES_H__
#define __UQ_APPL_ROUTINES_H__

#include <uqEnvironment.h>

//*****************************************************
// User must provide an application routine, to be called by main().
//*****************************************************
template<class V, class M>
void
uqAppl(const uqEnvironmentClass& env);

#endif // __UQ_APPL_ROUTINES_H__
