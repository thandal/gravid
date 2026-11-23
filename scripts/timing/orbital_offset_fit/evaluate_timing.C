#include "tempo2.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
int main(int argc,char **argv){
 if(argc!=5)return 2;
 MAX_OBSN=20000;
 pulsar *psr=(pulsar*)calloc(1,sizeof(pulsar));initialiseOne(psr,1,1);
 char par[1][MAX_FILELEN],tim[1][MAX_FILELEN];snprintf(par[0],MAX_FILELEN,"%s",argv[1]);snprintf(tim[0],MAX_FILELEN,"%s",argv[2]);snprintf(tempo2_clock_path,MAX_FILELEN,"%s",argv[3]);
 readParfile(psr,par,tim,1);readTimfile(psr,tim,1);char *options[]={argv[0]};preProcess(psr,1,1,options);formBatsAll(psr,1);formResiduals(psr,1,0);
 FILE *out=fopen(argv[4],"w");if(!out)return 3;
 fprintf(out,"# sat_mjd bat_mjd residual_s toaerr_s freq_MHz\n");
 for(int j=0;j<psr->nobs;j++){observation &o=psr->obsn[j];fprintf(out,"%.18Lg %.18Lg %.18Lg %.17g %.17g\n",o.sat,o.bat,o.residual,o.toaErr*1e-6,o.freq);}
 fclose(out);destroyOne(psr);free(psr);
}
