#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "EMPenModels.h"
#include "larsR.h"

static const R_CallMethodDef CallEntries[] = {
  {"cvEMfusedLasso1D",         (DL_FUNC) &cvEMfusedLasso1D,         12},
  {"cvEMfusedLasso2D",         (DL_FUNC) &cvEMfusedLasso2D,         11},
  {"cvEMlasso",                (DL_FUNC) &cvEMlasso,                10},
  {"cvEMlogisticFusedLasso1D", (DL_FUNC) &cvEMlogisticFusedLasso1D, 12},
  {"cvEMlogisticFusedLasso2D", (DL_FUNC) &cvEMlogisticFusedLasso2D, 11},
  {"cvEMlogisticLasso",        (DL_FUNC) &cvEMlogisticLasso,        10},
  {"EMfusedLasso",             (DL_FUNC) &EMfusedLasso,             10},
  {"EMlassoC",                 (DL_FUNC) &EMlassoC,                  9},
  {"EMlogisticFusedLasso",     (DL_FUNC) &EMlogisticFusedLasso,     10},
  {"EMlogisticLasso",          (DL_FUNC) &EMlogisticLasso,           9},
  {"fusion",                   (DL_FUNC) &fusion,                    7},
  {"lars",                     (DL_FUNC) &lars,                      7},
  {"cvlars",                   (DL_FUNC) &cvlars,                   11},
  {NULL, NULL, 0}
};

void R_init_HDPenReg(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
