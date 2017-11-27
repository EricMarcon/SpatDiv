#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _SpatDiv_parallelCountNbd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SpatDiv_parallelCountNbd",     (DL_FUNC) &_SpatDiv_parallelCountNbd,     6},
    {NULL, NULL, 0}
};

void R_init_SpatDiv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
