// Evaluate the binary delay alone at equally spaced barycentric MJDs.
// Each input row offsets the fitted parameters in native TEMPO2 units.
#include "tempo2.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
int main(int argc, char **argv) {
  if (argc != 5)
    return 2;
  FILE *in = fopen(argv[2], "r"), *out = fopen(argv[3], "wb");
  if (!in || !out)
    return 3;
  int n, nd, k;
  long double start, end;
  if (fscanf(in, "%d %d %d %Lf %Lf", &n, &nd, &k, &start, &end) != 5)
    return 4;
  MAX_OBSN = n;
  pulsar *psr = (pulsar *)calloc(1, sizeof(pulsar));
  initialiseOne(psr, 1, 1);
  char par[1][MAX_FILELEN], tim[1][MAX_FILELEN];
  snprintf(par[0], MAX_FILELEN, "%s", argv[1]);
  tim[0][0] = 0;
  readParfile(psr, par, tim, 1);
  psr->nobs = n;
  std::vector<int> key(k);
  std::vector<long double> base(k);
  for (int j = 0; j < k; j++) {
    if (fscanf(in, "%d", &key[j]) != 1)
      return 5;
    base[j] = psr->param[key[j]].val[0];
  }
  for (int i = 0; i < n; i++) {
    memset(&psr->obsn[i], 0, sizeof(observation));
    psr->obsn[i].bbat = start + (end - start) * i / (n - 1);
  }
  for (int row = 0; row < nd; row++) {
    for (int j = 0; j < k; j++) {
      long double delta;
      if (fscanf(in, "%Lf", &delta) != 1)
        return 6;
      psr->param[key[j]].val[0] = base[j] + delta;
    }
    for (int i = 0; i < n; i++) {
      double d;
      if (!strcmp(psr->binaryModel, "DD"))
        d = DDmodel(psr, 0, i, -1);
      else if (!strcmp(psr->binaryModel, "ELL1"))
        d = ELL1model(psr, 0, i, -1, 0);
      else if (!strcmp(psr->binaryModel, "T2"))
        d = T2model(psr, 0, i, -1, 0);
      else
        return 7;
      fwrite(&d, sizeof(double), 1, out);
    }
  }
  fclose(in);
  fclose(out);
  destroyOne(psr);
  free(psr);
  return 0;
}
