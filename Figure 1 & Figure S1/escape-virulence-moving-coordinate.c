#include <stdio.h>
#include <math.h>

#define EPS_QE  (1.0e-48) //THRESHOLD DENSITY I FOR "EXTINCTION"

#define BETA    2.0     // transmission rate
#define U       0.001    // birth rate/death rate
#define ALF     0.1     // virulence
#define GAMMA   0.5     // recovery rate
//#define SIG     2.0     // width of partial cross immunity (porality immunity)
#define SIGMIN  0.1
#define SIGMAX  2.0
#define NPARA   20

#define S0      1.0     // total host population density

#define DX      0.025//0.2     // GLID WIDTH
// Diffusion coefficient = 0.001
#define MUT     (0.001/(DX*DX))    // mutation rate
// MUT * DT must be smaller than 1
#define NX      (60*40)//(30*20)//(40*5) // number of antigenicity strains
                // XMAX = XMIN + NX*DX = 0 + 60*40*0.025 = 60

#define DT                  0.1
#define GEND                (600*10)
#define G_TRAJ_REC_END      (150*10)
#define TRAJ_REC_NPARA_A    1   // SIGMA = 0.2
#define TRAJ_REC_NPARA_B    5
#define GSTEP               (1*10)
#define GSTEP_MEAN          (1)

#define EPS     0.01    // initial infected density at x = 0
#define ABS(x)  ((x) >= 0 ? (x) : (-(x)))

//#define GREC_START  100
//#define XC      2.79    // For x > xc, BETA S(x) > U + GAMMA + ALF after the first epidemijc
//#define GREC_MOMENT 9000    // The central moments for within-morph distributions (for x<xc and for x>xc) are recorded up to this time step

#define XMIN    0.0//(-2.0)  // to start from the middle

int drift(double sig, int npara)
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
    char fname[80];
    
    double Itotal;
    
    double beta, u, alf, gamma, mut;
    beta = BETA*DT;
    u = U*DT;
    alf = ALF*DT;
    gamma = GAMMA*DT;
    mut = MUT*DT;
    
    // The range of i to update densities I[i] and S[i]
    // imin++ if I[imin] < EPS_QE
    // imax++ if I[imax] > EPS_QE
    // imax does not increase if imax = NX-1
    int imin=0, imax=10;
    
    
    //fp = fopen("traj.d", "w");
    if (npara == TRAJ_REC_NPARA_A || npara == TRAJ_REC_NPARA_B) {
        sprintf(fname, "dat/traj-%d.d", npara);
        fp = fopen(fname, "w");
    }
    
    //fq = fopen("mean.d", "w");
    sprintf(fname, "dat/mean-%d.d", npara);
    fq = fopen(fname, "w");
    
    //fr = fopen("within-morph-moments", "w");
    
    //
    for (i=0; i<NX; i++) {
        xval[i] = XMIN + DX*i;
    }
    xmax = XMIN + DX*(NX-1);
    
    // Partial cross immunity
    for (i=0; i<NX; i++) {
        //sigma[i] = exp(-0.5*pow(DX*i/SIG, 2.0));
        sigma[i] = exp(-0.5*pow(DX*i/sig, 2.0));
    }
    
    // initial distribution
    for (i=0; i<NX; i++) {
        S[i] = S0;
        I[i] = 0;
    }
    //I[0] = EPS;
    // XNIN + 40*DX = -2 + 40*0.05 = 0
    //I[40] = EPS;
    I[0] = EPS;
    
    Itotal = 0;
    for (i=0; i<NX; i++) {
        Itotal += I[i];
    }
    
    

    
    // iteration
    for (g=0; g<GEND; g++) {
        //*****************
        if (I[imin] < EPS_QE) {
            imin++;
            printf("xmin-xmax: %f-%f\n",xval[imin],xval[imax]);
            //scanf("%c", &key);
        }
        if (I[imax] > EPS_QE && imax <NX-1) { imax++;
            printf("xmin-xmax: %f-%f\n",xval[imin],xval[imax]);
            //scanf("%c", &key);
        }
        
        // TRUNCTION ONLY IN CROSS IMMUNITY
        //for (i=imin; i<=imax; i++) {
        //*****************
        for (i=0; i<NX; i++) {
            // susceptile
            // force of infection
            double lambda = 0;
            //for (j=0; j<NX; j++) {
            for (j=imin; j<=imax; j++) {
                lambda += sigma[ABS(i-j)] * I[j];
            }
            lambda *= beta;
            //printf("lambda = %f\n", lambda);
            
            //SS[i] = S[i] + DT*(U*(1-S[i]) - lambda * S[i] + ALF * I[i]);
            SS[i] = S[i] - lambda * S[i];
            // infectious
            double diff;
            if (i==0) {
                diff = mut * (I[1] - I[0]);
            } else if (i==NX-1) {
                diff = mut * (I[NX-2] - I[NX-1]);
            } else {
                diff = mut * (I[i+1] + I[i-1] - 2 * I[i]);
            }
            //II[i] = I[i] + DT * ( (BETA * S[i] - (U + GAMMA + ALF)) * I[i] + diff);
            II[i] = I[i] + (beta * S[i] - (u + gamma + alf)) * I[i] + diff;
        }
        Itotal = 0;
        //for (i=imin; i<=imax; i++) {
        for (i=0; i<NX; i++) {
            S[i] = SS[i];
            I[i] = II[i];
            Itotal += I[i];
        }

        if (g % GSTEP == 0) {
            double eps = 0.0001;
            if (npara == TRAJ_REC_NPARA_A || npara == TRAJ_REC_NPARA_B) {
                //printf("t = %f\n", g*DT);
                if (g < G_TRAJ_REC_END) {
                    for (i=0; i<NX; i++) {
                        fprintf(fp, "%g %f %f %f\n", DT*g, xval[i], S[i], I[i]);
                        //if (I[i] > eps) {
                        if (I[i] > EPS_QE) {
                            printf("x = %15.7f S = %15.7e I = %15.7e\n", xval[i], S[i], I[i]);
                        }
                    }
                }
            }
        }
        if (g % GSTEP_MEAN == 0) {
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
        

    }
    if (npara == TRAJ_REC_NPARA_A || npara == TRAJ_REC_NPARA_B) {
        fclose(fp);
    }
    fclose(fq);
    //fclose(fr);
    return 0;
}


int main()
{
    int npara;
    double sig;
    char key;
    
    printf("D*DT/(DX*DX) = %f\n", MUT*DT);
    printf("Hit return to continue\n");
    scanf("%c", &key);
    
    for (npara=0; npara<NPARA;npara++) {
        sig = SIGMIN + (SIGMAX-SIGMIN)*npara/(NPARA-1);
        drift(sig, npara);
    }
}
