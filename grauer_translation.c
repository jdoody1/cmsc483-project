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

    for (i = 0; i < Np; i++){
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
        exit(0);
    }

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
        for (int j = 0; j < ncells; j++){
            cellcounter[i][j] = 0;
            for (int k = 0; k < Np; k++){
                cellindex[i][j][k] = 0;
            }
        }
    }

    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            ccounter_infected[i][j] = 0;
            ccounter_normal[i][j] = 0;
            for (int k = 0; k < incells; k++){
                cindex_normal[i][j][k] = 0;
            }
        }
    }

    int i, j, ii, jj, ip;

    for (ip = 0; ip++; ip < Np){
        ii = (int) ceil(rp[0][ip] * ncells/(L));
        jj = (int) ceil(rp[1][ip] * ncells/(L));
        cellcounter[ii][jj] = cellcounter[ii][jj] + 1;
        cellindex[ii][jj][cellcounter[ii][jj]] = ip;
        i = (int) ceil(rp[0][ip] * incells/(L));
        j = (int) ceil(rp[1][ip] * incells/(L));

        if ((SIR[ip] == 2) || (SIR[ip] == 3)){
            ccounter_infected[i][j] = ccounter_infected[i][j] + 1;
        }
        else if (SIR[ip] == 1){
            ccounter_normal[i][j] = ccounter_normal[i][j] + 1;
            cindex_normal[i][j][ccounter_normal[i][j]] = ip;
        }
    }
}

// Convention: xij = xi-xj, yij = yi-yj, zij = zi-zj; Fij = -Fji.
void calcforces() {
    for (int i = 0; i < Np; i++){
        dxpart[i] = 0;
        dypart[i] = 0;
    }

    int i, j, ii, jj, icell, jcell, ip, jp, icount, jcount, icity;
    double r, r2, xij, yij, xforce, yforce, rforce, rmin;

    for (i = 0; i < ncells; i++){
        for (j = 0; j < ncells; j++){
            for (icount = 0; icount < cellcounter[i][j]; icount++){
                rmin = L;
                ip = cellindex[i][j][icount];

                for (ii = i-1; ii < i+1; ii++){
                    for (jj = j-1; jj < j+1; jj++){
                        if (ii == -1) icell = ncells - 1;
                        else if (ii == ncells) icell = 0;
                        else icell = ii;

                        if (jj == -1) jcell = ncells - 1;
                        else if (jj == ncells) jcell = 0;
                        else jcell = jj;

                        for (jcount = 0; jcount < cellcounter[icell][jcell]; jcount++){
                            jp = cellindex[icell][jcell][jcount];
                            if (ip != jp){
                                xij = rp[0][ip] - rp[0][jp];
                                yij = rp[1][ip] - rp[1][jp];

                                if (i == 0 && icell == ncells-1) xij = xij + L;
                                else if (i == ncells-1 && icell == 0) xij = xij - L;
                                if (j == 0 && jcell == ncells-1) yij = yij + L;
                                else if (j == ncells-1 && jcell == 0) yij = yij - L;

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
                                    if (r >= rminwca) rforce=-12*(pow(rm, 6)/pow(r, 7) - pow(rm, 12)/pow(r, 13));
                                    else rforce=-12*(pow(rm, 6)/pow(rminwca, 7) - pow(rm, 12)/pow(rminwca, 13));
                                    xforce = xij * rforce / r;
                                    yforce = yij * rforce / r;
                                    dxpart[ip] = dxpart[ip] + dt*epst_eff*xforce;
                                    dypart[ip] = dypart[ip] + dt*epst_eff*yforce;
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
                    if (xij > L/2) xij = L - xij;
                    else if (xij < -L/2) xij = L + xij;
                    if (yij > L/2) yij = L - yij;
                    else if (yij < -L/2) yij = L + yij;

                    xforce = xij*rforce/r;
                    yforce = yij*rforce/r;
                    dxpart[ip] = dxpart[ip] + dt*xforce;
                    dypart[ip] = dypart[ip] + dt*yforce;
                }
            }
        }
    }
}

float gasdev(float harvest){
    float rsq, v1, v2;
    if (gaus_stored) {
        harvest = g;
        gaus_stored = false;
    }
    else {
        while (true) {
            v1 = rando(&idum);
            v2 = rando(&idum);
            v1 = 2.0 * v1 - 1.0;
            v2 = 2.0 * v2 - 1.0;
            rsq = pow(v1, 2) + pow(v2, 2);
            if (rsq > 0.0 && rsq < 1.0) {break;}
        }
        rsq = sqrt(-2.0 * log(rsq) / rsq);
        harvest = v1 * rsq;
        g = v2 * rsq;
        gaus_stored = true;
    }
    if (harvest == 0) harvest = pow(10, -10);
    return harvest;
}

void diffusion() {
    int ip;
    float harvest;

    for (ip = 0; ip < Np; ip++){
        gasdev(harvest);
        dxpart[ip] = dxpart[ip] + sqrt(2 * D_T) * sqrt(dt) * harvest;
        gasdev(harvest);
        dypart[ip] = dypart[ip] + sqrt(2 * D_T) * sqrt(dt) * harvest;
    }
}

void move() {
    ncounter = ncounter_new;
    nindex = nindex_new;

    for (int i = 0; i < Np; i++){
        rp[0][i] = rp[0][i] + dxpart[i];
        rp[1][i] = rp[1][i] + dypart[i];

        dxpart[i] = 0;
        dypart[i] = 0;

        ncounter_new[i] = 0;
        for (int j = 0; j < Np; j++){
            nindex_new[i][j] = 0;
        }
    }
}

void boundarycondition() {
    for (int i = 0; i < Np; i++){
        if (rp[0][i] > L) rp[0][i] = rp[0][i] - L;
        else if (rp[0][i] < 0) rp[0][i] = rp[0][i] + L;
        if (rp[1][i] > L) rp[1][i] = rp[1][i] - L;
        else if (rp[1][i] < 0) rp[1][i] = rp[1][i] + L;
    }
}

void disease_progression() {
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

    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            if (ccounter_infected[i][j] > 0) inc[i][j] = ccounter_normal[i][j] * ccounter_infected[i][j];
            else inc[i][j] = 0;
            sum_inc += inc[i][j];
        }
    }

    double rest = 0;

    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            eff_inc[i][j] = inc[i][j] / sum_inc * availvacc;
        }
    }

    for (int i = 0; i < incells; i++){
        for (int j = 0; j < incells; j++){
            if (eff_inc[i][j] > ccounter_normal[i][j]) rest += eff_inc[i][j] - ccounter_normal[i][j];
            for (int k = 0; k < ccounter_normal[i][j]; k++){
                int idx = cindex_normal[i][j][k];
                if (rando(&idum) <= eff_inc[i][j] / ccounter_normal[i][j]) SIR[idx] = 5;
            }
        }
    }

    if (rest > 0){
        int sum_ccounter_normal = 0;
        for (int i = 0; i < incells; i++){
            for (int j = 0; j < incells; j++){
                sum_ccounter_normal += ccounter_normal[i][j];
            }
        }
        double prob = rest / sum_ccounter_normal;
        for (int i = 0; i < Np; i++){
            if (SIR[i] == 1){
                if (rando(&idum) <= prob) SIR[i] = 5;
            }
        }
    }

}

void writetrajectory(int framecount) {

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
