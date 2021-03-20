#include <stdio.h>
#include <math.h>

//#define BETA    5.0     // transmission rate
#define U       0.0 //0.001
    // birth rate/death rate
    // --> set to zero to prevent the reemergence of antigenicity
//#define ALF     0.1     // virulence
#define GAMMA   0.5     // recovery rate
#define MUT     0.01    // mutation rate
#define NX      300    // number of antigenicity strains
#define SIG     5.0     // width of partial cross immunity (porality immunity)
#define S0      1.0     // total host population density
#define GEND    200000
#define GSTEP_MEAN   10
#define NREC_MARGINAL 500
#define GSTEP_MARGINAL ((int) GEND/NREC_MARGINAL)
#define GS_JOINT    100000
#define GE_JOINT    120000
#define GSTEP_JOINT    50// 500
#define XMIN_JOINT  100
#define XMAX_JOINT  150
#define AMIN_JOINT  2
#define AMAX_JOINT  6
#define DT      0.001
#define EPS     0.01    // initial infected density at x = 0
#define ABS(x)  ((x) >= 0 ? (x) : (-(x)))

#define ALFMIN  0.025 //0.01 -> TO MAKE Initial R0 = 1.5
#define ALFMAX  10.0 // 20.0
#define C       5.0    // coefficient of beta-alpha tradeoff: beta = C sqrt(alpha)
#define NALF    100     // number of virulence genotypes
#define MUT_ALF 0.04 //0.01
    // mutation rate between adjacent virulence genoetypes
    // --> To maintain the same value for D_ALF = 0.0002
double beta_func(double alpha) {
    return C*sqrt(alpha);
}

int main()
{
    double S[NX];         // S(x): the number of people susceptible to strain x
    double I[NX][NALF];   // I(x, alpha): the number of people infectious with strain (x,alpha)
    double SS[NX], II[NX][NALF];
    double disease_induced_death[NX];
    double sigma[NX];   // sigma(x-y): the chance that the host acquires the immunity to strain x by the infection of strain y
    double beta[NALF];  // beta[i]: transmission rate of i-th virulence genotype
    double alpha[NALF]; // alpha[i]: viruelnce of i-th virulence genotype
    int i, j, k, g;
    FILE *fpx, *fpa, *fpt, *fq, *fp_joint, *fps;
    
    fpx = fopen("sim_data/marginal-x.d", "w");
    fpa = fopen("sim_data/marginal-a.d", "w");
    fps = fopen("sim_data/susceptibility.d","w");
    fpt = fopen("sim_data/total-I.d", "w");
    fq = fopen("sim_data/mean.d", "w");
    fp_joint = fopen("sim_data/joint.d","w");
    
    
    // Partial cross immunity
    for (i=0; i<NX; i++) {
        sigma[i] = exp(-0.5*pow(i/SIG, 2.0));
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
            double lambda = 0;
            for (j=0; j<NX; j++) {
                for (k=0; k<NALF; k++) {
                    lambda += sigma[ABS(i-j)] * beta[k] * I[j][k];
                }
            }
            disease_induced_death[i] = 0;
            for (k=0; k<NALF; k++) {
                disease_induced_death[i] += alpha[k] * I[i][k];
            }
            //printf("lambda = %f\n", lambda);
            
            SS[i] = S[i] +
                    DT*(
                        //U*(1-S[i]) + disease_induced_death[i]
                        - lambda * S[i]
                    );
            // infectious
            double diff, diff_alpha;

            for (k=0; k<NALF; k++) {
                if (i==0) {
                    diff = MUT * (I[1][k] - I[0][k]);
                } else if (i==NX-1) {
                    diff = MUT * (I[NX-2][k] - I[NX-1][k]);
                } else {
                    diff = MUT * (I[i+1][k] + I[i-1][k] - 2 * I[i][k]);
                }
                if (k==0) {
                    diff_alpha = MUT_ALF * (I[i][1] - I[i][0]);
                } else if (k==NALF-1) {
                    diff_alpha = MUT_ALF * (I[i][NALF-2] - I[i][NALF-1]);
                } else {
                    diff_alpha = MUT_ALF * (I[i][k+1] + I[i][k-1] - 2 * I[i][k]);
                }
                II[i][k] = I[i][k] +
                            DT * (
                                  (beta[k] * S[i] - (U + GAMMA + alpha[k])) * I[i][k]
                                  + diff + diff_alpha
                            );
            }
        }
        for (i=0; i<NX; i++) {
            S[i] = SS[i];
            for (k=0; k<NALF; k++)
                I[i][k] = II[i][k];
        }
    
        if (g % GSTEP_MEAN == 0) {
            double itotal=0, ex=0, exx = 0, ealf=0, eaa=0, exa=0;
            double vx, valf, cov_xalf;
            for (i=0; i<NX; i++) {
                for (k=0; k<NALF; k++) {
                    itotal += I[i][k];
                    ex     += i * I[i][k];
                    exx    += i * i * I[i][k];
                    ealf   += alpha[k] * I[i][k];
                    eaa    += alpha[k] * alpha[k] * I[i][k];
                    exa    += i * alpha[k] * I[i][k];
                }
            }
            ex /= itotal;
            ealf /= itotal;
            exx /= itotal;
            eaa /= itotal;
            exa /= itotal;
            vx = exx - ex*ex;
            valf = eaa - ealf*ealf;
            cov_xalf = exa - ex * ealf;
            printf("t = %7.2f IT = %12.5e ex = %12.5f ealf = %12.5f vx = %12.5e valf = %12.5e cxalf = %12.5e\n",
                   g*DT, itotal, ex, ealf, vx, valf, cov_xalf);
            fprintf(fq, "%f %e %f %f %f %f %f\n", g*DT, itotal, ex, ealf, vx, valf, cov_xalf);
        }
        if (g % GSTEP_MARGINAL== 0) {
            double Ix, Ia, Itotal;
            Itotal = 0;
            for (i=0; i<NX; i++) {
                Ix = 0;
                for (k=0; k<NALF; k++) {
                    Ix += I[i][k];
                }
                fprintf(fpx, "%f %d %f\n", g*DT, i, Ix);
                Itotal += Ix;
            }
            fprintf(fpt, "%f %f\n", g*DT, Itotal);
            for (k=0; k<NALF; k++) {
                Ia = 0;
                for (i=0; i<NX; i++) {
                    Ia += I[i][k];
                }
                fprintf(fpa, "%f %f %f\n", g*DT, alpha[k], Ia);
            }
            for (i=0; i<NX; i++) {
                fprintf(fps, "%f %d %f\n", g*DT, i, S[i]);
            }
        }
        if (g % GSTEP_JOINT== 0 && g>=GS_JOINT && g<=GE_JOINT) {
            double Itotal;
            Itotal = 0;
            for (i=0; i<NX; i++) {
                for (k=0; k<NALF; k++) {
                    Itotal += I[i][k];
                }
            }
            for (i=0; i<NX; i++) {
                if (i>=XMIN_JOINT && i<=XMAX_JOINT) {
                    for (k=0; k<NALF; k++) {
                        if (alpha[k]>=AMIN_JOINT && alpha[k]<=AMAX_JOINT) {
                            fprintf(fp_joint, "%f %d %f %f\n", g*DT, i, alpha[k], I[i][k]/Itotal);
                        }
                    }
                }
            }
        }
    }
    fclose(fq);
    fclose(fpx);
    fclose(fpa);
    fclose(fps);
    fclose(fpt);
    fclose(fp_joint);
    return 0;
}
