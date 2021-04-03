#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int I4B, I2B, I1B;
double SP, DP, SPC, DPC;
bool LGT;
int i,j,ncells,Nx,Ny,Np,pretime,simtime,framecount,idum,incells,equil,Ncity;
double L,rm,rminwca,rcutoff,dx,dy,D_T,dt,time,epst,epst_attr,mindestabstand,rinf,inf_prob,maxdisease_time,availvacc,impftime,epst_attr_eff,epst_eff,A_city,A_citytmp,rmincity,rn,mincitydist,distx,disty,dist,dummy;
FILE *infile;
double *dxpart, *dypart, *SIR, *disease_time, *k_city, *R_city, *rnrw, *Pnrw;
double **rp, **inc, **eff_inc, **rcity;
int *ncounter, *ncounter_new;
int **cellcounter, **nindex, **nindex_new, **ccounter_normal, **ccounter_infected;
int ***cellindex, ***cindex_normal;

double rando(int*);

int main(int argc, char *argv[]) {
    if (argc != 2){
        printf("You must supply the input file\n");
        exit(1);
    }
    infile = fopen(argv[1], "r");
    fscanf(infile, "%d %*[^\n]\n", &Np);
    fscanf(infile, "%lf %*[^\n]\n", &L);
    fscanf(infile, "%lf %*[^\n]\n", &dt);
    fscanf(infile, "%lf %*[^\n]\n", &epst);
    fscanf(infile, "%lf %*[^\n]\n", &epst_attr);
    fscanf(infile, "%d %*[^\n]\n", &idum);
    fscanf(infile, "%d %*[^\n]\n", &pretime);
    fscanf(infile, "%d %*[^\n]\n", &simtime);
    fscanf(infile, "%lf %*[^\n]\n", &D_T);
    fscanf(infile, "%lf %*[^\n]\n", &mindestabstand);
    fscanf(infile, "%lf %*[^\n]\n", &inf_prob);
    fscanf(infile, "%lf %*[^\n]\n", &maxdisease_time);
    fscanf(infile, "%lf %*[^\n]\n", &rinf);
    fscanf(infile, "%lf %*[^\n]\n", &impftime);
    fscanf(infile, "%lf %*[^\n]\n", &availvacc);
    fscanf(infile, "%d %*[^\n]\n", &incells);
    fscanf(infile, "%lf %*[^\n]\n", &A_city);
    fscanf(infile, "%d %*[^\n]\n", &Ncity);
    fscanf(infile, "%lf %*[^\n]\n", &mincitydist);

    rm = mindestabstand + 0.1;
    rcutoff = rm;
    ncells = 100;
    rminwca = 0.8*rm;

    epst_eff = 0;
    epst_attr_eff = 0;
    rmincity = 0;

    dx = 1.0;
    dy = dx;
    Nx = (int)round(L/dx);
    Ny = Nx;

    dxpart = malloc(sizeof(double)*Np);
    dypart = malloc(sizeof(double)*Np);
    cellcounter = malloc(sizeof(int*)*ncells);
    for (int idx = 0; idx < ncells; idx++){
        cellcounter[idx] = malloc(sizeof(int*)*ncells);
    }
    rp = malloc(sizeof(double*)*2);
    rp[0] = malloc(sizeof(double)*Np);
    rp[1] = malloc(sizeof(double)*Np);
    cellindex = malloc(sizeof(int**)*ncells);
    for (int idx = 0; idx < ncells; idx++){
        cellindex[idx] = malloc(sizeof(int*)*ncells);
        for (int jdx = 0; jdx < ncells; jdx++){
            cellindex[idx][jdx] = malloc(sizeof(int)*Np);
        }
    }
    ncounter = malloc(sizeof(int)*Np);
    nindex = malloc(sizeof(int*)*Np);
    for (int idx = 0; idx < Np; idx++){
        nindex[idx] = malloc(sizeof(int)*Np);
    }
    SIR = malloc(sizeof(double)*Np);
    disease_time = malloc(sizeof(double)*Np);
    ncounter_new = malloc(sizeof(int)*Np);
    nindex_new = malloc(sizeof(int*)*Np);
    for (int idx = 0; idx < Np; idx++){
        nindex_new[idx] = malloc(sizeof(int)*Np);
    }
    cindex_normal = malloc(sizeof(int**)*incells);
    for (int idx = 0; idx < incells; idx++){
        cindex_normal[idx] = malloc(sizeof(int*)*incells);
        for (int jdx = 0; jdx < incells; jdx++){
            cindex_normal[idx][jdx] = malloc(sizeof(int)*incells);
        }
    }
    ccounter_infected = malloc(sizeof(int*)*incells);
    ccounter_normal = malloc(sizeof(int*)*incells);
    inc = malloc(sizeof(double*)*incells);
    eff_inc = malloc(sizeof(double*)*incells);
    for (int idx = 0; idx < incells; idx++){
        ccounter_infected[idx] = malloc(sizeof(int)*incells);
        ccounter_normal[idx] = malloc(sizeof(int)*incells);
        inc[idx] = malloc(sizeof(int)*incells);
        eff_inc[idx] = malloc(sizeof(int)*incells);
    }

    rcity = malloc(sizeof(double*)*2);
    rcity[0] = malloc(sizeof(double)*Ncity);
    rcity[1] = malloc(sizeof(double)*Ncity);
    R_city = malloc(sizeof(double)*Ncity);
    k_city = malloc(sizeof(double)*Ncity);
    rnrw = malloc(sizeof(double)*26);
    Pnrw = malloc(sizeof(double)*26);

    for (int idx = 0; idx < Np; idx++) {
        disease_time[idx] = 0;
        dxpart[idx] = 0;
        dypart[idx] = 0;
        ncounter[idx] = 0;
        for (int jdx = 0; jdx < Np; jdx++){
            nindex[idx][jdx] = 0;
            nindex_new[idx][jdx] = 0;
        }
        ncounter_new[idx] = 0;
    }

    dist = 0;
    while (dist < mincitydist){
        for (i = 0; i < Ncity; i++){
            rcity[0][i] = (L-2*dx)*rando(idum)+dx;
            rcity[1][i] = (L-2*dx)*rando(idum)+dx;
            if (distx > L/2) {
                distx = L - distx;
            }
            if (disty > L/2) {
                disty = L - disty;
            }
            if (sqrt(distx*distx + disty*disty) <= dist) {
                dist = sqrt(pow(rcity[0][i] - rcity[0][j],2) + pow(rcity[1][i] - rcity[1][j],2));
            }
        }
    }

    for (i = 0; i < 26; i++){
        rnrw[i] = 20+i*2;
    }
    for (i = 0; i < 26; i++){
        Pnrw[i] = 1/rnrw[i]*1/log(rnrw[25]/rnrw[0]);
    }

    for (int idx = 0; idx < Ncity; idx++){
        R_city[idx] = 0;
    }
    for (j = 0; j < Ncity; j++){
        while (R_city[j] == 0){
            rn = rando(&idum);
            for (i = 0; i < 26; i++){
                if (rn < Pnrw[i] && rn > Pnrw[25]) {
                    R_city[j] = rnrw[i];
                }
            }
        }
        k_city[j] = 1/(2*pow(R_city[j],2));
    }

    for (i = 1; i < Np; i++){
        rp[0][i] = (L-2*dx)*rando(&idum)+dx;
        rp[1][i] = (L-2*dx)*rando(&idum)+dx;
    }

    time = 0;
    equil = 1;
    A_citytmp = A_city;
    while (time < pretime){
        if (time < pretime*1.0/4){
            A_city = 10*A_citytmp;
        } else {
            A_city = A_citytmp;
        }
        assigncells();
        calcforces();
        diffusion();
        move();
        boundarycondition();
        time += dt;
    }
    equil = 0;

    for (int idx = 0; idx < Np; idx++){
        SIR[idx] = 1;
    }
    for (int idx = 0; idx < Np/2000; idx++){
        SIR[idx] = 2;
    }
    framecount = 0;
    time = 0.0;

    while (time < simtime){
        if (time > 0){
            epst_eff = epst;
            epst_attr_eff = epst_attr;
        }
        assigncells();
        calcforces();
        diffusion();
        disease_progression();
        if (time > impftime){
            vaccination();
        }
        move();
        boundarycondition();
        if (fmod(time,2) < dt){
            writetrajectory(framecount);
            framecount++;
        }
        time += dt;
    }

    return 0;
}

double rando(int *idum){
    const int IA=16807;
    const int IM=2147483647;
    const int IQ=127773;
    const int IR=2836;
    static double am;
    static int ix=-1;
    static int iy = -1;
    static int k;
    if (idum <= 0 || iy < 0){
        am = 1.0/IM;
        iy = (888889999 ^ abs(*idum)) | 1;
        ix = 777755555 ^ abs(*idum);
        *idum = abs(*idum) + 1;
    }
    ix ^= ix << 13;
    ix ^= ix >> 17;
    ix ^= ix << 5;
    k = iy/IQ;
    iy = IA*(iy-k*IQ)-IR*k;
    if (iy < 0)
        iy += IM;
    return am * ((IM & (ix ^ iy)) | 1);
}
