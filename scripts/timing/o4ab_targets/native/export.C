// TEMPO2 output plugin for exporting the full design and the orbital fit.
#include "tempo2.h"
#include "t2fit.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
extern "C" int tempoOutput(int argc, char **argv, pulsar *psr, int npsr) {
  const param_label keys[] = {
      param_pb,    param_a1,    param_ecc,   param_om,      param_t0,
      param_pbdot, param_omdot, param_a1dot, param_tasc,    param_eps1,
      param_eps2,  param_m2,    param_sini,  param_eps1dot, param_eps2dot};
  const char *names[] = {"PB",    "A1",    "ECC",  "OM",      "T0",
                         "PBDOT", "OMDOT", "XDOT", "TASC",    "EPS1",
                         "EPS2",  "M2",    "SINI", "EPS1DOT", "EPS2DOT"};
  std::vector<int> indices;
  FILE *f = fopen("covariance.tsv", "w");
  if (!f)
    return 1;
  fprintf(f, "# binary covariance; native internal units; parameter rows "
             "include par-unit error\n");
  for (int j = 0; j < 15; j++)
    for (int k = 0; k < psr[0].fitinfo.output.totalNfit; k++)
      if (psr[0].fitinfo.output.indexParam[k] == keys[j] &&
          psr[0].fitinfo.output.indexCounter[k] == 0) {
        indices.push_back(k);
        fprintf(f, "parameter\t%s\t%d\t%.21Lg\t%.21Lg\t%.17g\n", names[j],
                (int)keys[j], psr[0].param[keys[j]].val[0],
                psr[0].param[keys[j]].err[0], sqrt(psr[0].covar[k][k]));
      }
  fprintf(f, "covariance\n");
  for (int i : indices) {
    for (int j : indices)
      fprintf(f, "%.17g\t", psr[0].covar[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);

  FitInfo &fi = psr[0].fitinfo;
  int n = psr[0].nobs, p = fi.nParams, nc = fi.nConstraints;
  FILE *meta = fopen("design.meta", "w");
  fprintf(meta, "%d %d %d\n", n, p, nc);
  for (int k = 0; k < p; k++)
    fprintf(meta, "%d %d\n", (int)fi.paramIndex[k], fi.paramCounters[k]);
  fclose(meta);
  FILE *design = fopen("design.bin", "wb"), *toas = fopen("toas.bin", "wb");
  std::vector<double> row(p);
  for (int i = 0; i < n; i++) {
    observation &o = psr[0].obsn[i];
    if (o.deleted)
      return 2;
    double x = (double)(o.bbat - psr[0].param[param_pepoch].val[0]);
    std::fill(row.begin(), row.end(), 0.);
    t2Fit_buildDesignMatrix(psr, 0, x, i, row.data());
    fwrite(row.data(), sizeof(double), p, design);
    double v[4] = {(double)o.sat, (double)o.residual, o.toaErr * 1e-6, o.freq};
    fwrite(v, sizeof(double), 4, toas);
  }
  fclose(design);
  fclose(toas);
  if (nc) {
    fprintf(stderr, "Constraints unsupported in robust audit\n");
    return 3;
  }
  textOutput(psr, npsr, 0, 0, 0, 1, (char *)"native.par");
  return 0;
}
const char *plugVersionCheck = TEMPO2_h_VER;
