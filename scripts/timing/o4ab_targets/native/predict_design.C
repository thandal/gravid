// Form timing residuals and design derivatives at new TOAs without fitting.
#include "tempo2.h"
#include "t2fit.h"
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
int main(int argc, char **argv) {
  if (argc != 5)
    return 2;
  MAX_OBSN = 20000;
  pulsar *psr = (pulsar *)calloc(1, sizeof(pulsar));
  initialiseOne(psr, 1, 1);
  char par[1][MAX_FILELEN], tim[1][MAX_FILELEN];
  snprintf(par[0], MAX_FILELEN, "%s", argv[1]);
  snprintf(tim[0], MAX_FILELEN, "%s", argv[2]);
  snprintf(tempo2_clock_path, MAX_FILELEN, "%s", argv[3]);
  readParfile(psr, par, tim, 1);
  readTimfile(psr, tim, 1);
  char *options[] = {argv[0]};
  preProcess(psr, 1, 1, options);
  formBatsAll(psr, 1);
  formResiduals(psr, 1, 0);
  FILE *out = fopen(argv[4], "w");
  if (!out)
    return 3;
  fprintf(out, "# sat_mjd bat_mjd residual_s toaerr_s freq_MHz\n");
  for (int j = 0; j < psr->nobs; j++) {
    observation &o = psr->obsn[j];
    fprintf(out, "%.18Lg %.18Lg %.18Lg %.17g %.17g\n", o.sat, o.bat, o.residual,
            o.toaErr * 1e-6, o.freq);
  }
  fclose(out);
  std::vector<double> x(psr[0].nobs), y(psr[0].nobs), e(psr[0].nobs);
  std::vector<int> ip(psr[0].nobs);
  unsigned int nd =
      t2Fit_getFitData(psr, x.data(), y.data(), e.data(), ip.data());
  FitInfo globals;
  t2Fit_fillGlobalFitInfo(psr, 1, globals);
  t2Fit_fillFitInfo(psr, psr[0].fitinfo, globals, x.data(), ip.data(), nd);
  FitInfo &fi = psr[0].fitinfo;
  int p = fi.nParams;
  char path[4096];
  snprintf(path, sizeof path, "%s.meta", argv[4]);
  FILE *meta = fopen(path, "w");
  fprintf(meta, "%d %d\n", psr[0].nobs, p);
  for (int k = 0; k < p; k++)
    fprintf(meta, "%d %d\n", (int)fi.paramIndex[k], fi.paramCounters[k]);
  fclose(meta);
  snprintf(path, sizeof path, "%s.design", argv[4]);
  FILE *matrix = fopen(path, "wb");
  std::vector<double> row(p);
  for (int i = 0; i < psr[0].nobs; i++) {
    std::fill(row.begin(), row.end(), 0.);
    double xx =
        (double)(psr[0].obsn[i].bbat - psr[0].param[param_pepoch].val[0]);
    t2Fit_buildDesignMatrix(psr, 0, xx, i, row.data());
    fwrite(row.data(), sizeof(double), p, matrix);
  }
  fclose(matrix);
  destroyOne(psr);
  free(psr);
}
