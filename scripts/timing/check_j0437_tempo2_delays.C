// Read-only evaluation of the installed TEMPO2 T2 model. No fitting,
// covariance transformation, clock correction, or barycentering is performed.
// Inputs are shared TDB barycentric epochs and Earth/pulsar vectors.
#include "tempo2.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
int main(int argc,char **argv) {
    if(argc!=6) return 2;
    FILE *samples=fopen(argv[2],"r"), *draws=fopen(argv[3],"r"), *out=fopen(argv[4],"wb"), *meta=fopen(argv[5],"w");
    if(!samples || !draws || !out || !meta) return 3;
    int n,nd; if(fscanf(samples,"%d",&n)!=1 || fscanf(draws,"%d",&nd)!=1) return 4;
    MAX_OBSN=n;
    pulsar *psr=(pulsar*)calloc(1,sizeof(pulsar)); initialiseOne(psr,1,1);
    char par[1][MAX_FILELEN],tim[1][MAX_FILELEN];
    snprintf(par[0],MAX_FILELEN,"%s",argv[1]);tim[0][0]=0;
    readParfile(psr,par,tim,1);psr->nobs=n;
    if(strcmp(psr->binaryModel,"T2") || psr->nCompanion!=1) return 5;
    const param_label keys[10]={param_pb,param_a1,param_ecc,param_om,param_t0,param_pbdot,param_omdot,param_kin,param_kom,param_m2};
    longdouble base[10];
    const longdouble rad2deg=180.L/acosl(-1.L);
    const longdouble scale[10]={1.L/86400,1,1,rad2deg,1.L/86400,1e-12L,1,rad2deg,rad2deg,1};
    for(int j=0;j<10;j++) {base[j]=psr->param[keys[j]].val[0];fprintf(meta,"%d %.21Lg\n",j,base[j]);}
    for(int i=0;i<n;i++) {
        observation &o=psr->obsn[i];memset(&o,0,sizeof(observation));
        if(fscanf(samples,"%Lf",&o.bbat)!=1) return 6;
        for(int j=0;j<3;j++) if(fscanf(samples,"%lf",&o.earth_ssb[j])!=1)return 7;
        for(int j=0;j<3;j++) if(fscanf(samples,"%lf",&o.psrPos[j])!=1)return 8;
    }
    for(int k=0;k<nd;k++) {
        for(int j=0;j<10;j++) {
            longdouble delta;if(fscanf(draws,"%Lf",&delta)!=1)return 9;
            psr->param[keys[j]].val[0]=base[j]+delta*scale[j];
        }
        psr->param[param_sini].val[0]=sinl(psr->param[param_kin].val[0]/rad2deg);
        for(int i=0;i<n;i++) {double delay=T2model(psr,0,i,-1,0);fwrite(&delay,sizeof(double),1,out);}
    }
    fclose(samples);fclose(draws);fclose(out);fclose(meta);destroyOne(psr);free(psr);
    return 0;
}
