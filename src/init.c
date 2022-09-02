#include "survS.h"
#include "R_ext/Rdynload.h"
#include "Rversion.h"
#include "survproto.h"

static const R_CallMethodDef Callentries[] = {
  {"Ccoxscore2",  (DL_FUNC) &coxscore2,  6},
  {"Cmulticheck", (DL_FUNC) &multicheck, 6},
  {NULL, NULL, 0}
};

void R_init_survival(DllInfo *dll){
  R_registerRoutines(dll, NULL, Callentries, NULL, NULL);

  /* The following line makes only those routines defined above
   available to outside packages, i.e., internal things like
   dmatrix() are now invisible.
   */
  R_useDynamicSymbols(dll, FALSE);
  /*
   ** This line makes them only be available via the symbols above
   **  i.e., .Call("tmerge", ) won't work but .Call(Ctmerge, ) will.
   ** This feature was added in version 2.16
   */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
  R_forceSymbols(dll, TRUE);
#endif
}
