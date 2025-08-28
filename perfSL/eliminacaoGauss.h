#ifndef __ELIM_GAUSS__
#define __ELIM_GAUSS__

#include "utils.h"
#include "sislin.h"

void triangulariza( SistLinear_t *C );
void retrosubst( SistLinear_t *C, real_t *X );

#endif