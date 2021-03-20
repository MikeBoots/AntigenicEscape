#include <stdio.h>
#include <math.h>

#define XMIN    0.0
#define XMAX    80.0
#define NX      1600     // number of antigenicity strains
#define REC_NX_STEP   10
#define DEL_X   ((XMAX-XMIN)/NX)
#define DIFFUSION_COEFFICIENT_X 0.01
//#define MUT     0.01    // mutation rate
#define MUT     (2*DIFFUSION_COEFFICIENT_X / (DEL_X*DEL_X))

#define ALFMIN  0.0
#define ALFMAX  20.0
#define NALF    100     // number of virulence genotypes
#define DEL_ALF ((ALFMAX-ALFMIN)/NALF)
#define DIFFUSION_COEFFICIENT_ALF 0.01 //(2.0 e-4)
//#define MUT_ALF 0.01    // mutation rate between adjacent virulence genoetypes
#define MUT_ALF (2*DIFFUSION_COEFFICIENT_ALF / (DEL_ALF*DEL_ALF))

#define U       0.001   // birth rate/death rate
#define GAMMA   0.5     // recovery rate

#define C       5.0     // coefficient of beta-alpha tradeoff: beta = C sqrt(alpha)
#define DT      0.001

//#define SIG     5.0     // width of partial cross immunity (porality immunity)
#define S0      1.0     // total host population density

#define EPS     0.01    // initial infected density at x = 0
#define ABS(x)  ((x) >= 0 ? (x) : (-(x)))


#define GEND        150000
#define GSTEP_MEAN     100

#define NREC_MARGINAL  500
#define GSTEP_MARGINAL ((int) GEND/NREC_MARGINAL)

#define GS_JOINT         0
#define GE_JOINT    100000
#define GSTEP_JOINT   5000


double beta_func(double alpha) {
    return C*sqrt(alpha);
}

int main()
{
    double S[NX];         // S(x): the number of people susceptible to strain x
    double I[NX][NALF];   // I(x, alpha): the number of people infectious with strain (x,alpha)
    double SS[NX], II[NX][NALF];
    //double disease_induced_death[NX];
    //double sigma[NX];   // sigma(x-y): the chance that the host acquires the immunity to strain x by the infection of strain y
    double xval[NX];
    double beta[NALF];  // beta[i]: transmission rate of i-th virulence genotype
    double alpha[NALF]; // alpha[i]: viruelnce of i-th virulence genotype
    int i, j, k, g;
    FILE *fpx, *fpa, *fpt, *fq, *fp_joint;
    
    fpx = fopen("dat/marginal-x.d", "w");
    fpa = fopen("dat/marginal-a.d", "w");
    fpt = fopen("dat/total-I.d", "w");
    fq = fopen("dat/mean.d", "w");
    fp_joint = fopen("dat/joint.d","w");
    
    
    // Partial cross immunity
    /*
    for (i=0; i<NX; i++) {
        sigma[i] = exp(-0.5*pow(i/SIG, 2.0));
    }
     */
    // antigenicity
    for (i=0; i<NX; i++) {
        xval[i] = XMIN + DEL_X * i;
    }
    
    // virulence
    for (k=0; k<NALF; k++) {
        alpha[k] = ALFMIN + (ALFMAX-ALFMIN) * k / (NALF - 1);
    }
    // transmission rate
    for (k=0; k<NALF; k++) {
        beta[k] = beta_func(alpha[k]);
    }
    
    // initial distribution
    for (i=0; i<NX; i++) {
        S[i] = S0;
        for (k=0; k<NALF; k++)
            I[i][k] = 0;
    }
    I[0][0] = EPS;
    
    // iteration
    for (g=0; g<GEND; g++) {
        for (i=0; i<NX; i++) {
            // susceptile
            // force of infection against S(i)
            /*
            double lambda = 0;
            for (j=0; j<NX; j++) {
                for (k=0; k<NALF; k++) {
                    lambda += sigma[ABS(i-j)] * beta[k] * I[j][k];
                }
            }
            */
            double lambda = 0;
            for (k=0; k<NALF; k++) {
                lambda += beta[k] * I[i][k];
            }
            
            /*
            disease_induced_death[i] = 0;
            for (k=0; k<NALF; k++) {
                disease_induced_death[i] += alpha[k] * I[i][k];
            }
             */
            //printf("lambda = %f\n", lambda);
            
            SS[i] = S[i] +
                    DT*(
                        //U*(1-S[i]) + disease_induced_death[i]
                        - lambda * S[i]
                    );
            // infectious
            double diff, diff_alpha;

            for (k=0; k<NALF; k++) {
                switch (i) {
                    case 0:
                        diff = 0.5 * (I[1][k] - I[0][k]);
                        break;
                    case NX-1:
                        diff = 0.5 * (I[NX-2][k] - I[NX-1][k]);
                        break;
                    default:
                        diff = 0.5 * (I[i+1][k] + I[i-1][k]) - I[i][k];
                }
                
                switch (k) {
                    case 0:
                        diff_alpha = 0.5 * (I[i][1] - I[i][0]);
                        break;
                    case NALF-1:
                        diff_alpha = 0.5 * (I[i][NALF-2] - I[i][NALF-1]);
                        break;
                    default:
                        diff_alpha = 0.5 * (I[i][k+1] + I[i][k-1]) - I[i][k];
                }

                II[i][k] = I[i][k] +
                            DT * (
                                  (beta[k] * S[i] - (U + GAMMA + alpha[k])) * I[i][k]
                                  + MUT     * diff
                                  + MUT_ALF * diff_alpha
                            );
            }
        }
        for (i=0; i<NX; i++) {
            S[i] = SS[i];
            for (k=0; k<NALF; k++)
                I[i][k] = II[i][k];
        }
    
        if (g % GSTEP_MEAN == 0) {
            double itotal=0, ex=0, exx = 0, ealf=0, eaa=0;
            double vx, valf;
            for (i=0; i<NX; i++) {
                for (k=0; k<NALF; k++) {
                    itotal += I[i][k];
                    ex     += xval[i] * I[i][k];
                    exx    += xval[i] * xval[i] * I[i][k];
                    ealf   += alpha[k] * I[i][k];
                    eaa    += alpha[k] * alpha[k] * I[i][k];
                }
            }
            ex /= itotal;
            ealf /= itotal;
            exx /= itotal;
            eaa /= itotal;
            vx = exx - ex*ex;
            valf = eaa - ealf*ealf;
            printf("t = %15.7f IT = %15.7e ex = %15.7f vx = %15.7f ealf = %15.7f valf = %15.7f\n",
                   g*DT, itotal, ex, vx, ealf, valf);
            fprintf(fq, "%f %e %f %f %f %f\n", g*DT, itotal, ex, vx, ealf, valf);
        }
        if (g % GSTEP_MARGINAL== 0) {
            double Ix, Ia, Itotal;

            for (i=0; i<NX; i++) {
                if (i % REC_NX_STEP==0) {
                    Ix = 0;
                    for (k=0; k<NALF; k++) {
                        Ix += I[i][k];
                    }
                    fprintf(fpx, "%f %f %f\n", g*DT, xval[i], Ix);
                }
            }
            
            Itotal = 0;
            for (k=0; k<NALF; k++) {
                Ia = 0;
                for (i=0; i<NX; i++) {
                    Ia += I[i][k];
                }
                fprintf(fpa, "%f %f %f\n", g*DT, alpha[k], Ia);
                Itotal += Ia;
            }
            fprintf(fpt, "%f %f\n", g*DT, Itotal);
            
        }
        if (g % GSTEP_JOINT== 0 && g>=GS_JOINT && g<=GE_JOINT) {
            double Itotal;
            Itotal = 0;
            for (i=0; i<NX; i++) {
                //if (i % REC_NX_STEP==0) {
                    for (k=0; k<NALF; k++) {
                        Itotal += I[i][k];
                    }
                //}
            }
            for (i=0; i<NX; i++) {
                //if (i % REC_NX_STEP==0) {
                    for (k=0; k<NALF; k++) {
                        fprintf(fp_joint, "%f %f %f %f\n", g*DT, xval[i], alpha[k], I[i][k]/Itotal);
                    }
                //}
            }
        }
    }
    fclose(fq);
    fclose(fpx);
    fclose(fpa);
    fclose(fpt);
    fclose(fp_joint);
    return 0;
}
