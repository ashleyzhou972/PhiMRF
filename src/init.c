#include "add.h"
#include "dm_call.h"
#include <R_ext/Rdynload.h>


void R_init_PhiMRF(DllInfo *info) {
  R_RegisterCCallable("PhiMRF", "add",  (DL_FUNC) &add_);
  R_RegisterCCallable("PhiMRF", "pmrf",  (DL_FUNC) &double_metropolis_cont);
  R_RegisterCCallable("PhiMRF", "preprocess_big",  (DL_FUNC) &preprocess);

}
