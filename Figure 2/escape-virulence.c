#include <stdio.h>
#include <math.h>

#define BETA    2.0     // transmission rate
#define U       0.001    // birth rate/death rate
#define ALF     0.1     // virulence
#define GAMMA   0.5     // recovery rate
#define DX      0.2     // GLID WIDTH
#define MUT     (0.001/(DX*DX))    // mutation rate
#define NX      (40*5) // number of antigenicity strains
#define SIG     2.0     // width of partial cross immunity (porality immunity)
#define S0      1.0     // total host population density
#define GEND    30000
#define GSTEP   100
#define GSTEP_MEAN 10
#define DT      0.01
#define EPS     0.01    // initial infected density at x = 0
#define ABS(x)  ((x) >= 0 ? (x) : (-(x)))

#define GREC_START  100
#define XC      2.79    // For x > xc, BETA S(x) > U + GAMMA + ALF after the first epidemijc
#define GREC_MOMENT 9000    // The central moments for within-morph distributions (for x<xc and for x>xc) are recorded up to this time step

#define XMIN    (-2.0)  // to start from the middle

int main()
{
    double S[NX];   // S(x): the number of people susceptible to strain x
    double I[NX];   // I(x): the number of people infectious with strain x
    double xval[NX], ex;
    double SS[NX], II[NX];
    double sigma[NX];   // sigma(x-y): the chance that the host acquires the immunity to strain x by the infection of strain y
    int i, j, g;
    FILE *fp, *fq, *fr;
    char key;
    double xmax;
    
    fp = fopen("dat/traj.d", "w");
    fq = fopen("dat/mean.d", "w");
    fr = fopen("dat/within-morph-moments", "w");
    
    //
    for (i=0; i<NX; i++) {
        xval[i] = XMIN + DX*i;
    }
    xmax = XMIN + DX*(NX-1);
    
    // Partial cross immunity
    for (i=0; i<NX; i++) {
        sigma[i] = exp(-0.5*pow(DX*i/SIG, 2.0));
    }
    
    // initial distribution
    for (i=0; i<NX; i++) {
        S[i] = S0;
        I[i] = 0;
    }
    //I[0] = EPS;
    I[10] = EPS;
    //printf("%f\n", xval[10]);
    //scanf("%c", &key);
    
    // iteration
    for (g=0; g<GEND; g++) {
        for (i=0; i<NX; i++) {
            // susceptile
            // force of infection
            double lambda = 0;
            for (j=0; j<NX; j++) {
                lambda += sigma[ABS(i-j)] * I[j];
            }
            lambda *= BETA;
            //printf("lambda = %f\n", lambda);
            
            //SS[i] = S[i] + DT*(U*(1-S[i]) - lambda * S[i] + ALF * I[i]);
            SS[i] = S[i] + DT*(
                               //U*(1-S[i])
                               - lambda * S[i]
                               //+ ALF * I[i]
                               );
            // infectious
            double diff;
            if (i==0) {
                diff = MUT * (I[1] - I[0]);
            } else if (i==NX-1) {
                diff = MUT * (I[NX-2] - I[NX-1]);
            } else {
                diff = MUT * (I[i+1] + I[i-1] - 2 * I[i]);
            }
            II[i] = I[i] + DT * ( (BETA * S[i] - (U + GAMMA + ALF)) * I[i] + diff);
        }
        for (i=0; i<NX; i++) {
            S[i] = SS[i];
            I[i] = II[i];
        }

        if (g > GREC_START && g % GSTEP == 0) {
            double eps = 0.0001;
            //printf("t = %f\n", g*DT);
            for (i=0; i<NX; i++) {
                fprintf(fp, "%g %f %e %e\n", DT*g, xval[i], S[i], I[i]);
                if (I[i] > eps) {
                    printf("x = %15.7f S = %15.7e I = %15.7e\n", xval[i], S[i], I[i]);
                }
            }
        }
        if (g > GREC_START && g % GSTEP_MEAN == 0) {
            double itotal=0;
            ex = 0;
            for (i=0; i<NX; i++) {
                itotal += I[i];
                ex     += xval[i] * I[i];
            }
            ex /= itotal;
            
            double vx = 0;
            for (i=0; i<NX; i++) {
                vx += pow(xval[i] - ex, 2.0) * I[i];
            }
            vx /= itotal;
            printf("t = %15.7f IT = %15.7e ex = %15.7e vx = %15.7e\n",
               g*DT, itotal, ex, vx);
            fprintf(fq, "%f %e %e %e\n", g*DT, itotal, ex, vx);
        }
        if (ex > 0.95*xmax)
            break;
        
        if (g > GREC_START && g < GREC_MOMENT) {
            if (g % GSTEP_MEAN == 0) {
                double i0, i1, ex0, ex1, vx0, vx1, tx0, tx1, qx0, qx1;
                i0 = i1 = ex0 = ex1 = 0;
                for (i=0; i<NX; i++) {
                    if (xval[i] < XC) {
                        i0 += I[i];
                        ex0 += xval[i] * I[i];
                    } else {
                        i1 += I[i];
                        ex1 += xval[i] * I[i];
                    }
                }
                ex0 /= i0;
                ex1 /= i1;
                vx0 = vx1 = tx0 = tx1 = qx0 = qx1 = 0;
                for (i=0; i<NX; i++) {
                    if (xval[i] < XC) {
                        vx0 += pow(xval[i] - ex, 2.0) * I[i];
                        tx0 += pow(xval[i] - ex, 3.0) * I[i];
                        qx0 += pow(xval[i] - ex, 4.0) * I[i];
                    } else {
                        vx1 += pow(xval[i] - ex, 2.0) * I[i];
                        tx1 += pow(xval[i] - ex, 3.0) * I[i];
                        qx1 += pow(xval[i] - ex, 4.0) * I[i];
                    }
                }
                vx0 /= i0; tx0 /= i0; qx0 /= i0;
                vx1 /= i1; tx1 /= i1; qx1 /= i1;
                printf("t = %15.7f:\n",g*DT);
                printf("i0 = %15.7e ex0 = %15.7e vx0 = %15.7e tx0 = %15.7e qx0 = %15.7e\n", i0, ex0, vx0, tx0, qx0);
                printf("i1 = %15.7e ex1 = %15.7e vx1 = %15.7e tx1 = %15.7e qx1 = %15.7e\n", i1, ex1, vx1, tx1, qx1);
                fprintf(fr, "%f %e %e %e %e %e %e %e %e %e %e\n",
                        g*DT, i0, ex0, vx0, tx0, qx0, i1, ex1, vx1, tx1, qx1);
            }
        }
    }
    fclose(fp);
    fclose(fq);
    fclose(fr);
    return 0;
}
