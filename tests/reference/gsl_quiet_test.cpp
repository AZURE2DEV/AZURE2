// GslQuiet (GSLException.h): quiet GSL calls from many OpenMP threads leave the
// process-wide GSL error handler as it was.
//
// A bare save / gsl_set_error_handler_off / restore from parallel threads can
// interleave so that a thread restores the "off" another one saved: the
// throwing handler AZURE2 installs is then gone for good.  With GslQuiet the
// first scope in turns the handler off and the last one out puts back the one
// that was installed, so after any number of concurrent scopes it is the
// original again -- and inside a scope it is off (a GSL error returns a status).

#include "GSLException.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_log.h>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif

static int handled = 0;
static void countingHandler(const char *, const char *, int, int) {
#pragma omp atomic
  handled++;
}

int main() {
  int fails = 0;
  gsl_set_error_handler(&countingHandler);

  // Inside a scope: off, so an error is a status, not a handler call.
  {
    GslQuiet quiet;
    gsl_sf_result r;
    int status = gsl_sf_log_e(-1.0, &r);  // domain error
    if (status == GSL_SUCCESS || handled != 0) {
      std::printf("FAIL  inside a scope the handler is off (status %d, handler calls %d)\n", status, handled);
      fails++;
    } else {
      std::printf("ok    inside a scope the handler is off\n");
    }
  }

  // Many threads, many short scopes, overlapping.
  long scopes = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : scopes)
  for (int i = 0; i < 20000; i++) {
    GslQuiet quiet;
    gsl_sf_result r;
    gsl_sf_hyperg_U_e(1.0 + (i % 7) * 0.1, 2.0, 0.5 + (i % 11), &r);
    scopes++;
  }
  gsl_error_handler_t *now = gsl_set_error_handler(&countingHandler);
  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif
  if (now != &countingHandler) {
    std::printf("FAIL  after %ld concurrent scopes (%d threads) the installed handler is lost\n", scopes, threads);
    fails++;
  } else {
    std::printf("ok    after %ld concurrent scopes (%d threads) the installed handler is back\n", scopes, threads);
  }

  // Outside every scope the installed handler is called again.
  gsl_sf_result r;
  gsl_sf_log_e(-1.0, &r);
  if (handled != 1) {
    std::printf("FAIL  outside a scope the installed handler is called (%d calls)\n", handled);
    fails++;
  } else {
    std::printf("ok    outside a scope the installed handler is called\n");
  }
  std::printf(fails ? "FAILED\n" : "PASSED\n");
  return fails ? 1 : 0;
}
