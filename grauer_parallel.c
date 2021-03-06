#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#ifdef __linux__ 
#include <sys/time.h>
#endif
#include <time.h>
#include <omp.h>

int THREADS = 1;
int I4B, I2B, I1B;
double SP, DP, SPC, DPC;
bool LGT;
int i,j,ncells,Nx,Ny,Np,pretime,simtime,framecount,idum,incells,equil,Ncity;
double L,rm,rminwca,rcutoff,dx,dy,D_T,dt,main_time,epst,epst_attr,mindestabstand,rinf,inf_prob,maxdisease_time,availvacc,impftime,epst_attr_eff,epst_eff,A_city,A_citytmp,rmincity,rn,mincitydist,distx,disty,dist,dummy;
FILE *infile;
double *dxpart, *dypart, *SIR, *disease_time, *k_city, *R_city, *rnrw, *Pnrw;
double **rp, **inc, **eff_inc, **rcity;
int *ncounter, *ncounter_new;
int **cellcounter, **nindex, **nindex_new, **ccounter_normal, **ccounter_infected;
int ***cellindex, ***cindex_normal;
float g;
bool gaus_stored = false;

double rando(int*);
void assigncells();
void calcforces();
void diffusion();
void move();
void boundarycondition();
void disease_progression();
void writetrajectory(int framecount);
void vaccination();

int main(int argc, char *argv[]) {
    if (argc != 3){
        printf("You must supply the input file and # of OMP threads to use!\n");
        exit(1);
    }

    char* tc = argv[2];
    THREADS = atoi(tc);

    if (THREADS <= 0) {
        printf("You must supply a # of OMP threads that is greater than 0!\n");
        exit(1);
    }

    // Disable dynamic teams and use THREADS # of threads
    omp_set_dynamic(0);
    omp_set_num_threads(THREADS);
    printf("New total max threads: %d\n", omp_get_max_threads());

    // Get the starting time of our program's calculations
    struct timeval t0, t1;
    gettimeofday(&t0, 0);

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
            cindex_normal[idx][jdx] = malloc(sizeof(int)*Np);
        }
    }
    ccounter_infected = malloc(sizeof(int*)*incells);
    ccounter_normal = malloc(sizeof(int*)*incells);
    inc = malloc(sizeof(double*)*incells);
    eff_inc = malloc(sizeof(double*)*incells);
    for (int idx = 0; idx < incells; idx++){
        ccounter_infected[idx] = malloc(sizeof(int)*incells);
        ccounter_normal[idx] = malloc(sizeof(int)*incells);
        inc[idx] = malloc(sizeof(double)*incells);
        eff_inc[idx] = malloc(sizeof(double)*incells);
    }

    rcity = malloc(sizeof(double*)*2);
    rcity[0] = malloc(sizeof(double)*Ncity);
    rcity[1] = malloc(sizeof(double)*Ncity);
    R_city = malloc(sizeof(double)*Ncity);
    k_city = malloc(sizeof(double)*Ncity);
    rnrw = malloc(sizeof(double)*26);
    Pnrw = malloc(sizeof(double)*26);

		memset(disease_time, 0, sizeof(double) * Np);
		memset(dxpart, 0, sizeof(double) * Np);
		memset(dypart, 0, sizeof(double) * Np);
		memset(ncounter, 0, sizeof(int) * Np);
		memset(ncounter_new, 0, sizeof(int) * Np);
    for (int idx = 0; idx < Np; idx++) {
			memset(nindex[idx], 0, sizeof(int) * Np);
			memset(nindex_new[idx], 0, sizeof(int) * Np);
    }

    dist = 0;
    while (dist < mincitydist){
        for (i = 0; i < Ncity; i++) {
            rcity[0][i] = (L-2*dx)*rando(&idum)+dx;
            rcity[1][i] = (L-2*dx)*rando(&idum)+dx;
        }
        dist = L;
        for (i = 0; i < Ncity; i++){
            for (j = 0; j < Ncity; j++){
                if (i != j){
                    distx = fabs(rcity[0][i] - rcity[0][j]);
                    disty = fabs(rcity[1][i] - rcity[1][j]);
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
        }
    }

    for (i = 0; i < 26; i++){
        rnrw[i] = 20+i*2;
    }
    for (i = 0; i < 26; i++){
        Pnrw[i] = 1/rnrw[i]*1/log(rnrw[25]/rnrw[0]);
    }
		memset(R_city, 0, sizeof(double) * Ncity);

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

    for (i = 0; i < Np; i++){
        rp[0][i] = (L-2*dx)*rando(&idum)+dx;
        rp[1][i] = (L-2*dx)*rando(&idum)+dx;
    }

    main_time = 0;
    equil = 1;
    A_citytmp = A_city;
    while (main_time < pretime){
        if (main_time < pretime/4){
            A_city = 10*A_citytmp;
        } else {
            A_city = A_citytmp;
        }

//        printf("Pretime: time=%d\n", time);
//        printf("Before assigncells: %f %f\n", rp[0][0], rp[1][0]);
        assigncells();
//        printf("Before calcforces: %f %f\n", rp[0][0], rp[1][0]);
        calcforces();
//        printf("Before diffusion: %f %f\n", rp[0][0], rp[1][0]);
//        printf("dxpart before diffusion: %.50f\n", dxpart[0]);
        diffusion();
//        printf("Before move: %f %f\n", rp[0][0], rp[1][0]);
//        printf("dxpart before move: %.50f\n", dxpart[0]);
        move();
//        printf("Before boundarycondition: %f %f\n", rp[0][0], rp[1][0]);
        boundarycondition();
//        printf("Before time: %f %f\n", rp[0][0], rp[1][0]);
        main_time += dt;
    }
    equil = 0;

		memset(SIR, 1, sizeof(double) * Np);
		memset(SIR, 2, sizeof(double) * Np/2000);
    framecount = 0;
    main_time = 0.0;

    while (main_time < simtime){
        if (main_time > 0){
            epst_eff = epst;
            epst_attr_eff = epst_attr;
        }
//        printf("Timestep: time=%f\n", time);
//        printf("Before assigncells: %f %f\n", rp[0][0], rp[1][0]);
        assigncells();
//        printf("Before calcforces: %f %f\n", rp[0][0], rp[1][0]);
        calcforces();
//        printf("Before diffusion: %f %f\n", rp[0][0], rp[1][0]);
        diffusion();
//        printf("Before disease_progression: %f %f\n", rp[0][0], rp[1][0]);
        disease_progression();
        if (main_time > impftime){
//            printf("Before vaccination: %f %f\n", rp[0][0], rp[1][0]);
            vaccination();
        }
//        printf("Before move: %f %f\n", rp[0][0], rp[1][0]);
        move();
//        printf("Before boundarycondition: %f %f\n", rp[0][0], rp[1][0]);
        boundarycondition();
//        printf("Before writetrajectory: %f %f\n", rp[0][0], rp[1][0]);
        if (fmod(main_time,2) < dt){
            writetrajectory(framecount);
            framecount++;
        }
        main_time += dt;
//        exit(0);
    }

    // Get the time elapsed since starting the program timer, print that value
    gettimeofday(&t1, 0);
    double elapsed = (t1.tv_sec-t0.tv_sec) * 1.0f + (t1.tv_usec - t0.tv_usec) / 1000000.0f;
    printf("Total time elapsed: %f\n\n", elapsed);
    return 0;
}

double rando(int *idum){
    const int IA=16807;
    const int IM=2147483647;
    const int IQ=127773;
    const int IR=2836;

    const float NEAREST_BELOW_ONE = 0.999999940395355224609375f;
    static float am;
    static int ix=-1;
    static int iy = -1;
    static int k;
    if (*idum <= 0 || iy < 0){
        am = NEAREST_BELOW_ONE/IM;
        iy = (888889999 ^ abs(*idum)) | 1;
        ix = 777755555 ^ abs(*idum);
        *idum = abs(*idum) + 1;
    }
    ix ^= (unsigned int)ix << 13;
    ix ^= (unsigned int)ix >> 17;
    ix ^= (unsigned int)ix << 5;
    k = iy/IQ;
    iy = IA*(iy-k*IQ)-IR*k;
    if (iy < 0)
        iy += IM;
    return am * ((IM & (ix ^ iy)) | 1);
}

/**
 * Assign particles to cells based on their coordinates;
 * particle indices in each cell are sorted from low to high.
 */
void assigncells() {

    for (int i = 0; i < ncells; i++){
        memset(cellcounter[i], 0, sizeof(int) * ncells);
        for (int j = 0; j < ncells; j++){
            memset(cellindex[i][j], 0, sizeof(int) * Np);
        }
    }

    for (int i = 0; i < incells; i++){
        memset(ccounter_infected[i], 0, sizeof(int) * incells);
        memset(ccounter_normal[i], 0, sizeof(int) * incells);
        for (int j = 0; j < incells; j++){
            memset(cindex_normal[i][j], 0, sizeof(int) * incells);
        }
    }    

    #pragma omp parallel
    {

    #pragma omp for
    for (int ip = 0; ip < Np; ip++){
        int i, j, ii, jj;

        ii = (int) ceil(rp[0][ip] * ncells/(L)) - 1;
        jj = (int) ceil(rp[1][ip] * ncells/(L)) - 1;
        cellcounter[ii][jj]++;
        cellindex[ii][jj][cellcounter[ii][jj]-1] = ip;
        i = (int) ceil(rp[0][ip] * incells/(L)) - 1;
        j = (int) ceil(rp[1][ip] * incells/(L)) - 1;

        if ((SIR[ip] == 2) || (SIR[ip] == 3)){
            ccounter_infected[i][j] = ccounter_infected[i][j] + 1;
        }
        else if (SIR[ip] == 1){
            ccounter_normal[i][j]++;
            cindex_normal[i][j][ccounter_normal[i][j]-1] = ip;
        }
    }
    }
}

// Convention: xij = xi-xj, yij = yi-yj, zij = zi-zj; Fij = -Fji.
void calcforces() {
    // Get the starting time of our program's calculations
    //struct timeval t0, t1;
    //gettimeofday(&t0, 0);
    memset(dxpart, 0, sizeof(double) * Np);
    memset(dypart, 0, sizeof(double) * Np);
    #pragma omp parallel
    {

    int icell, jcell, ip, jp, jcount, icity;
    double r, r2, xij, yij, xforce, yforce, rforce, rmin;

    #pragma omp for
    for (int i = 0; i < ncells; i++){
        for (int j = 0; j < ncells; j++){
            for (int icount = 0; icount < cellcounter[i][j]; icount++){
                rmin = L;
                ip = cellindex[i][j][icount];
                for (int ii = i-1; ii <= i+1; ii++){
                    for (int jj = j-1; jj <= j+1; jj++){
                        if (ii == -1)
                            icell = ncells - 1;
                        else if (ii == ncells)
                            icell = 0;
                        else
                            icell = ii;

                        if (jj == -1)
                            jcell = ncells - 1;
                        else if (jj == ncells)
                            jcell = 0;
                        else
                            jcell = jj;

                        for (jcount = 0; jcount < cellcounter[icell][jcell]; jcount++){
                            jp = cellindex[icell][jcell][jcount];
                            if (ip != jp){
                                xij = rp[0][ip] - rp[0][jp];
                                yij = rp[1][ip] - rp[1][jp];

                                if (i == 0 && icell == ncells-1)
                                    xij = xij + L;
                                else if (i == ncells-1 && icell == 0)
                                    xij = xij - L;
                                if (j == 0 && jcell == ncells-1)
                                    yij = yij + L;
                                else if (j == ncells-1 && jcell == 0)
                                    yij = yij - L;

                                r2 = pow(xij, 2) + pow(yij, 2);
                                r = sqrt(r2);

                                // count neighbors - social distancing
                                if (r <= 1.5*mindestabstand){
                                    ncounter_new[ip] = ncounter_new[ip] + 1;
                                    nindex_new[ip][ncounter_new[ip]] = jp;
                                    if (r <= rmin) {
                                        nindex_new[ip][ncounter_new[ip]] = nindex_new[ip][0];
                                        nindex_new[ip][0] = jp;
                                        rmin = r;
                                    }
                                }

                                // Infection
                                if ((equil == 0) && (r <= rinf) && ((SIR[ip] == 2) || (SIR[ip] == 3)) && (SIR[jp] == 1)){
                                    if (rando(&idum) <= inf_prob){
                                        if (rando(&idum) <= 0.75) SIR[jp] = 2;
                                        else SIR[jp] = 3;
                                    }
                                }

                                // Potential
                                if ((r <= rcutoff) && ((jp != nindex[ip][0]) || (SIR[jp] == 3))){
                                    // WCA-pot
                                    if (r >= rminwca)
                                        rforce=-12*(pow(rm, 6)/pow(r, 7) - pow(rm, 12)/pow(r, 13));
                                    else
                                        rforce=-12*(pow(rm, 6)/pow(rminwca, 7) - pow(rm, 12)/pow(rminwca, 13));
                                    xforce = xij * rforce / r;
                                    yforce = yij * rforce / r;
                                    dxpart[ip] += dt*epst_eff*xforce;
                                    dypart[ip] += dt*epst_eff*yforce;
                                }

                                if ((r <= mindestabstand) && (jp == nindex[ip][0]) && (SIR[jp] != 6)){
                                    // attractive part of WCA-pot for pair attraction
                                    if (r >= rminwca) rforce = -12*(pow(rm, 6)/pow(r, 7));
                                    else rforce = -12*(pow(rm, 6)/pow(rminwca, 7));
                                    xforce = xij*rforce/r;
                                    yforce = yij*rforce/r;
                                    dxpart[ip] = dxpart[ip] + dt*epst_attr_eff*xforce;
                                    dypart[ip] = dypart[ip] + dt*epst_attr_eff*yforce;
                                }
                            }
                        }
                    }
                }
                // attraction from City
                for (icity = 0; icity < Ncity; icity++){
                    xij = rp[0][ip] - rcity[0][icity];
                    yij = rp[1][ip] - rcity[1][icity];
                    if (xij > L/2)
                        xij = L - xij;
                    else if (xij < -L/2)
                        xij = L + xij;
                    if (yij > L/2)
                        yij = L - yij;
                    else if (yij < -L/2)
                        yij = L + yij;
                    r2 = xij*xij + yij*yij;
                    r = sqrt(r2);
                    if (r >= rmincity){
                        rforce=-2*A_city*k_city[icity]*r*exp(-k_city[icity]*(r*r));
                    } else {
                        rforce=-2*A_city*k_city[icity]*rmincity*exp(-k_city[icity]*(rmincity*rmincity));
                    }

//                    printf("rforce: %.50f\n", rforce);
//                    exit(0);
                    xforce = xij*rforce/r;
                    yforce = yij*rforce/r;
                    dxpart[ip] += dt*xforce;
                    dypart[ip] += dt*yforce;
                }
            }
        }
    }
    }

    // Get the time elapsed since starting the program timer, print that value
    //gettimeofday(&t1, 0);
    //double elapsed = (t1.tv_sec-t0.tv_sec) * 1.0f + (t1.tv_usec - t0.tv_usec) / 1000000.0f;
    //printf("Total time elapsed in calcforces(): %f\n\n", elapsed);
//    printf("dxpart: %.50f\n", dxpart[0]);
//    exit(0);
}

void gasdev(float *harvest){
    float rsq, v1, v2;
    if (gaus_stored) {
        *harvest = g;
        gaus_stored = false;
    }
    else {
        while (true) {
            v1 = rando(&idum);
            v2 = rando(&idum);
            v1 = (float)2.0 * v1 - (float)1.0;
            v2 = (float)2.0 * v2 - (float)1.0;
            rsq = v1*v1 + v2*v2;
            if (rsq > 0.0 && rsq < 1.0) {break;}
        }
        rsq = sqrtf((float)-2.0 * logf(rsq) / rsq);
        *harvest = v1 * rsq;
        g = v2 * rsq;
        gaus_stored = true;
    }
    if (harvest == 0) *harvest = pow(10, -10);
}

void diffusion() {
    int ip;
    float harvest;

    #pragma omp parallel for
    for (ip = 0; ip < Np; ip++){
        gasdev(&harvest);
        dxpart[ip] = dxpart[ip] + sqrt(2 * D_T) * sqrt(dt) * harvest;
        gasdev(&harvest);
        dypart[ip] = dypart[ip] + sqrt(2 * D_T) * sqrt(dt) * harvest;
    }
}

void move() {
    ncounter = ncounter_new;
    nindex = nindex_new;

    #pragma omp parallel for
    for (int idx = 0; idx < Np; idx++){
        rp[0][idx] += dxpart[idx];
        rp[1][idx] += dypart[idx];

        dxpart[idx] = 0;
        dypart[idx] = 0;

        ncounter_new[idx] = 0;
        for (int jdx = 0; jdx < Np; jdx++){
            nindex_new[idx][jdx] = 0;
        }
    }
}

void boundarycondition() {
    #pragma omp parallel for
    for (int i = 0; i < Np; i++){
        if (rp[0][i] > L) rp[0][i] = rp[0][i] - L;
        else if (rp[0][i] < 0) rp[0][i] = rp[0][i] + L;
        if (rp[1][i] > L) rp[1][i] = rp[1][i] - L;
        else if (rp[1][i] < 0) rp[1][i] = rp[1][i] + L;
    }
}

void disease_progression() {
    #pragma omp parallel for
    for (int i = 0; i < Np; i++){
        if (SIR[i] == 2){
            disease_time[i] += dt;
            if (disease_time[i] >= maxdisease_time){
                if (rando(&idum) <= 0.01) SIR[i] = 6;
                else SIR[i]=4;
            }
        } else if (SIR[i] == 3){
            disease_time[i] += dt;
            if (disease_time[i] >= maxdisease_time){
                if (rando(&idum) <= 0.035) SIR[i] = 6;
                else SIR[i] = 4;
            }

        } else if (SIR[i] == 4 || SIR[i] == 6){
            disease_time[i] += dt;
        }

        if ((SIR[i] == 3 && disease_time[i] >= maxdisease_time/3) || SIR[i] == 6){
            dxpart[i] = 0;
            dypart[i] = 0;
        }
    }
}

void vaccination() {
    double sum_inc = 0;
    double rest = 0;

    #pragma omp parallel
    {

    //#pragma omp single
    //printf("Before first vacc block\n");

    #pragma omp for
    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            if (ccounter_infected[i][j] > 0)
                inc[i][j] = ccounter_normal[i][j] * ccounter_infected[i][j];
            else
                inc[i][j] = 0;
            sum_inc += inc[i][j];
        }
    }

    //# pragma omp single
    //printf("Before second vacc block\n");

    # pragma omp for
    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            eff_inc[i][j] = inc[i][j] / sum_inc * availvacc;
        }
    }

    //# pragma omp single
    //printf("Before third vacc block\n");

    # pragma omp for
    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            if (eff_inc[i][j] > ccounter_normal[i][j])
                rest += eff_inc[i][j] - ccounter_normal[i][j];
            for (int icount = 0; icount < ccounter_normal[i][j]; icount++){
                int idx = cindex_normal[i][j][icount];
                if (rando(&idum) <= eff_inc[i][j] / ccounter_normal[i][j])
                    SIR[idx] = 5;
            }
        }
    }

    //#pragma omp single
    //printf("Before fourth vacc block\n");

    if (rest > 0){
        int sum_ccounter_normal = 0;
        #pragma omp for
        for (int i = 0; i < incells; i++){
            for (int j = 0; j < incells; j++){
                sum_ccounter_normal += ccounter_normal[i][j];
            }
        }
        double prob = rest / sum_ccounter_normal;
        #pragma omp for
        for (int i = 0; i < Np; i++){
            if (SIR[i] == 1){
                if (rando(&idum) <= prob)
                    SIR[i] = 5;
            }
        }
    }
    }

    //printf("End of vacc\n");

}

void writetrajectory(int framecount) {

    printf("writetrajectory(%d)\n", framecount);

    char trajFileName[14];
    snprintf(trajFileName, 14, "traj_%04d.dat", framecount);
    FILE *trajFile = fopen(trajFileName, "w+");
    for (int idx = 0; idx < Np; idx++){
        char line[40];
        snprintf(line, 50, "%.18g %.18g\n", rp[0][idx], rp[1][idx]);
        fputs(line, trajFile);
    }
    fclose(trajFile);

    char dataFileName[14];
    snprintf(dataFileName, 14, "data_%04d.dat", framecount);
    FILE *dataFile = fopen(dataFileName, "w+");
    for (int idx = 0; idx < Np; idx++){
        char line[40];
        snprintf(line, 50, "%.18g %.18g\n", SIR[idx], disease_time[idx]);
        fputs(line, dataFile);
    }
    fclose(dataFile);
}
