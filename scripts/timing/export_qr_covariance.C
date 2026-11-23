#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "tempo2.h"
#include "t2fit.h"
#include <algorithm>
#include <vector>

extern "C" int tempoOutput(int argc, char **argv, pulsar *psr, int npsr) {
    const param_label labels[5] = {param_pb, param_a1, param_ecc, param_om, param_t0};
    const char *names[5] = {"PB", "A1", "ECC", "OM", "T0"};
    int indices[5];
    for (int j=0;j<5;j++) {
        indices[j]=-1;
        for (int k=0;k<psr[0].fitinfo.output.totalNfit;k++) {
            if (psr[0].fitinfo.output.indexParam[k]==labels[j] && psr[0].fitinfo.output.indexCounter[k]==0) indices[j]=k;
        }
        if(indices[j]<0) { fprintf(stderr,"Missing orbital fit parameter %s\n",names[j]); return 1; }
    }
    const char *file=getenv("AUDIT_COV_OUTPUT");
    FILE *out=fopen(file ? file : "native-covariance.tsv","w");
    if(!out) return 1;
    fprintf(out,"NATIVECOV1\n");
    fprintf(out,"total_fit\t%d\n",psr[0].fitinfo.output.totalNfit);
    fprintf(out,"chi2\t%.17g\ndof\t%d\n",psr[0].fitChisq,psr[0].fitNfree);
    for(int j=0;j<5;j++) {
        fprintf(out,"parameter\t%s\t%d\t%.21Lg\t%.21Lg\t%.17g\n",names[j],indices[j],psr[0].param[labels[j]].val[0],psr[0].param[labels[j]].err[0],sqrt(psr[0].covar[indices[j]][indices[j]]));
    }
    fprintf(out,"covariance_internal\n");
    for(int j=0;j<5;j++) {
        for(int k=0;k<5;k++) fprintf(out,"%s%.17g", k ? "\t" : "",psr[0].covar[indices[j]][indices[k]]);
        fprintf(out,"\n");
    }
    fclose(out);
    // Optional read-only export for an independent QR check. Derivatives are
    // rebuilt at the post-fit parameters, so this is a relinearization check.
    const char *prefix=getenv("QR_DESIGN_PREFIX");
    if(prefix) {
        FitInfo &fi=psr[0].fitinfo;
        int n=psr[0].nobs, p=fi.nParams, nc=fi.nConstraints;
        char path[4096];
        snprintf(path,sizeof path,"%s.meta",prefix);
        FILE *meta=fopen(path,"w"); if(!meta) return 1;
        fprintf(meta,"%d %d %d\n",n,p,nc);
        for(int k=0;k<p;k++)
            fprintf(meta,"%d %d\n",(int)fi.paramIndex[k],fi.paramCounters[k]);
        fclose(meta);
        std::vector<double> row(p);
        snprintf(path,sizeof path,"%s.design",prefix);
        FILE *design=fopen(path,"wb"); if(!design) return 1;
        snprintf(path,sizeof path,"%s.toas",prefix);
        FILE *toas=fopen(path,"wb"); if(!toas) { fclose(design); return 1; }
        for(int i=0;i<n;i++) {
            observation &o=psr[0].obsn[i];
            if(o.deleted) { fprintf(stderr,"Deleted TOAs unsupported by QR check\n"); return 1; }
            double x=(double)(o.bbat-psr[0].param[param_pepoch].val[0]);
            std::fill(row.begin(),row.end(),0.0);
            t2Fit_buildDesignMatrix(psr,0,x,i,row.data());
            fwrite(row.data(),sizeof(double),p,design);
            double values[2]={(double)o.prefitResidual,o.toaErr*1e-6};
            fwrite(values,sizeof(double),2,toas);
        }
        fclose(design); fclose(toas);
        snprintf(path,sizeof path,"%s.constraints",prefix);
        FILE *constraints=fopen(path,"wb"); if(!constraints) return 1;
        for(int c=0;c<nc;c++) {
            std::fill(row.begin(),row.end(),0.0);
            t2Fit_buildConstraintsMatrix(psr,0,c,row.data());
            fwrite(row.data(),sizeof(double),p,constraints);
        }
        fclose(constraints);
        snprintf(path,sizeof path,"%s.orbit",prefix);
        FILE *orbit=fopen(path,"w"); if(!orbit) return 1;
        for(int j=0;j<5;j++) {
            int index=-1;
            for(int k=0;k<p;k++)
                if(fi.paramIndex[k]==labels[j] && fi.paramCounters[k]==0) index=k;
            if(index<0) return 1;
            fprintf(orbit,"%d\n",index);
        }
        fclose(orbit);
    }
    // Standard output is unchanged apart from the additional covariance file.
    textOutput(psr,npsr,0,0,0,1,(char*)"native.par");
    return 0;
}
const char *plugVersionCheck = TEMPO2_h_VER;
