#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***
        CONSTANTS OF THE SYSTEM
                                    ***/

#define NA 54
#define NB 54
#define epsBB 1
#define epsAB 1
#define sigBB 1
#define sigAB 1
#define dt 0.001
#define tmax 5000
#define teq 3000
#define rho 0.8442
#define nhis 108
#define T0 0.728
#define PI 3.14159

typedef struct{
    char type;
    double pos[3];
    double pos_msd[3];
    double pos_ref[3];
    double vel[3];
    double vel_ref[1];
    double force[3];
}Particle;


/***
        GLOBAL VARIABLES
                            ***/

int N = NA + NB;
double L, T, P;
double U = 0, K = 0, E = 0;

void init(Particle *particles);
void force(Particle *particles);
void integrate(int stage, Particle *particles);
void savePositions(int time, Particle *particles);

/***
        MAIN PROGRAM
                        ***/

int main(){

    printf("\nThis is a Molecular Dynamic Simulation that uses Lennard-Jones Potential and Velocity-Verlet Algorithm.\n");
    printf("It positions a determined number of particles in a box, calculates the forces between them and integrate these forces.\n");
    printf("Energies, positions, temperature and pressure are saved in files.\n");
    printf("This version does not make measures as radial distribution function and mean squared displacement.\n");
    printf("Author: Vitoria Tesser Henkes. Coauthor: Camila Raupp da Luz. Both from Federal University of Rio Grande do Sul.\n\n");
    printf("\t|Initial Conditions|\n");
    printf("\tTemperature: %.3f\n\tNumber of particles: %i\n\tDensity: %.4f\n\tMaximum simulation step: %i\n\tTime step: %.3f\n", T0, N, rho, tmax, dt);

    int t, n;
    Particle *particles;
    particles = (Particle*) malloc(N * sizeof(Particle));

    FILE *file_energies = fopen("energies.txt", "w");
    FILE *file_temp = fopen("temp.txt", "w");
    FILE *file_press = fopen("press.txt", "w");

    fprintf(file_energies, "# t U K E\n");
    fprintf(file_temp, "# t T\n");
    fprintf(file_press, "# t P\n");

    init(particles);
    force(particles);

    E = U + K;
    fprintf(file_energies, "%i %f %f %f\n", 0, U, K, E);
    fprintf(file_temp, "%i %f\n", 0, T);
    fprintf(file_press, "%i %f\n", 0, P);

    savePositions(0, particles);

    for(t = 1; t <= tmax; t++){

        integrate(1, particles);
        force(particles);
        integrate(2, particles);

        E = U + K;
        fprintf(file_energies, "%i %f %f %f\n", t, U, K, E);
        fprintf(file_temp, "%i %f\n", t, T);
        fprintf(file_press, "%i %f\n", t, P);

        if(t % 100 == 0){
            savePositions(t, particles);
            printf("\nArrived on time %d", t);
        }
    }

    printf("\n\n\t|Final Conditions|");
    printf("\n\tPotential Energy per particle: %.4f", U);
    printf("\n\tTotal Energy per particle: %.4f", E);
    printf("\n\tKinetic Energy per particle: %.4f", K);
    printf("\n\tTemperature: %.4f", T);
    printf("\n\tPressure: %.4f\n", P);

    fclose(file_temp);
    fclose(file_energies);
    free(particles);

    return 0;
}

/***
        INITIALIZATION
                            ***/

void init(Particle *particles){

    int n, n3;
    int i = 0, j = 0, k = 0;
    double space, scaleFactor;
    double sumvx = 0, sumvy = 0, sumvz = 0, sumv2 = 0;

    L = cbrt((double) N/rho); //cubic box length
    n3 = ceil(cbrt(N)); //maximum number of particles per dimension
    space = (double) L/n3; //space between particles

    printf("\n\t|Box Characteristics|\n");
    printf("\tCubic\n\tSide: %.4f\n", L);

    for(n = 0; n < N; n++){

        //differentiation of particles
        if(n < NA)
            particles[n].type = 'A';
        else
            particles[n].type = 'B';

        //initial positions
        particles[n].pos[0] = i * space;
        particles[n].pos[1] = j * space;
        particles[n].pos[2] = k * space;

        i++;
        if(i == n3){
            i = 0;
            j++;
            if(j == n3){
                j = 0;
                k++;
            }
        }

        //random velocities betweeen -0.5 and 0.5
        particles[n].vel[0] = (double) rand()/RAND_MAX - 0.5;
        particles[n].vel[1] = (double) rand()/RAND_MAX - 0.5;
        particles[n].vel[2] = (double) rand()/RAND_MAX - 0.5;

        //sum of velocities
        sumvx += particles[n].vel[0];
        sumvy += particles[n].vel[1];
        sumvz += particles[n].vel[2];

        //sum of squared velocities
        sumv2 += (pow(particles[n].vel[0], 2) + pow(particles[n].vel[1], 2) + pow(particles[n].vel[2], 2));
    }

    //mean of velocities
    sumvx /= N;
    sumvy /= N;
    sumvz /= N;

    //kinetic energy and temperature
    K = sumv2/(2*N);
    T = (2/3.)*K;

    //correction of velocities
    scaleFactor = sqrt(T0/T);
    for(n = 0; n < N; n++){
        particles[n].vel[0] = (particles[n].vel[0] - sumvx) * scaleFactor;
        particles[n].vel[1] = (particles[n].vel[1] - sumvy) * scaleFactor;
        particles[n].vel[2] = (particles[n].vel[2] - sumvz) * scaleFactor;
    }
}

/***
        FORCE CALCULATION
                                ***/

void force(Particle *particles){

    int n, ni, nj;
    double r, rc2, ecut;
    double xr, yr, zr;
    double r2, r2i, r6i, ff, s6;

    //potential energy and forces null
    U = 0;
    for(n = 0; n < N; n++)
        particles[n].force[0] = particles[n].force[1] = particles[n].force[2] = 0;
    P = rho*T;

    rc2 = pow(L/2, 2);

    for(ni = 0; ni < N - 1; ni++){

        for(nj = ni + 1; nj < N; nj++){

            //distance between particles
            xr = particles[ni].pos[0] - particles[nj].pos[0];
            xr = xr - (L * round(xr/L)); //periodic boundary condition

            yr = particles[ni].pos[1] - particles[nj].pos[1];
            yr = yr - (L * round(yr/L));

            zr = particles[ni].pos[2] - particles[nj].pos[2];
            zr = zr - (L * round(zr/L));

            r2 = pow(xr, 2) + pow(yr, 2) + pow(zr, 2);

            //force calculation only if particles are close enough
            if(r2 < rc2){

                r2i = 1/r2;
                r6i = pow(r2i, 3);

                //differentiation of interactions according to the types
                if(particles[ni].type == 'A' && particles[nj].type == 'A'){
                    ecut = 4 * (1/pow(rc2, 6) - 1/pow(rc2, 3));
                    ff = 48 * r2i * r6i * (r6i - 0.5);
                    U += (4 * r6i * (r6i - 1)) - ecut;
                }else if(particles[ni].type == 'B' && particles[nj].type == 'B'){
                    s6 = pow(sigBB, 6);
                    ecut = 4 * epsBB * s6 * (s6/pow(rc2, 6) - 1/pow(rc2, 3));
                    ff = 48 * epsBB * r2i * r6i * s6 * (r6i * s6 - 0.5);
                    U += (4 * epsBB* r6i * s6 * (r6i * s6 - 1)) - ecut;
                }else{
                    s6 = pow(sigAB, 6);
                    ecut = 4 * epsAB * s6 * (s6/pow(rc2, 6) - 1/pow(rc2, 3));
                    ff = 48 * epsAB * r2i * r6i * s6 * (r6i * s6 - 0.5);
                    U += (4 * epsAB* r6i * s6 * (r6i * s6 - 1)) - ecut;
                }

                //force at each direction
                particles[ni].force[0] += ff * xr;
                particles[nj].force[0] -= ff * xr;

                particles[ni].force[1] += ff * yr;
                particles[nj].force[1] -= ff * yr;

                particles[ni].force[2] += ff * zr;
                particles[nj].force[2] -= ff * zr;

                //pressure
                P += (sqrt(r2) * ff) * rho/(3*N*T);
            }
        }
    }

    U /= N;
}

/***
        INTEGRATION OF MOTION
                                 ***/

void integrate(int stage, Particle *particles){

   int n;
   double xx = 0, yy = 0, zz = 0;
   double sumvx = 0, sumvy = 0, sumvz = 0, sumv2 = 0;

    switch(stage){
        //first integration
        case 1:
            for(n = 0; n < N; n++){
                //first middle time interval displacement
                xx = particles[n].vel[0]*dt/2;
                yy = particles[n].vel[1]*dt/2;
                zz = particles[n].vel[2]*dt/2;

                //middle time interval position
                particles[n].pos[0] += xx;
                particles[n].pos[1] += yy;
                particles[n].pos[2] += zz;

                //periodic boundary condition
                particles[n].pos[0] -= floor(particles[n].pos[0]/L)*L;
                particles[n].pos[1] -= floor(particles[n].pos[1]/L)*L;
                particles[n].pos[2] -= floor(particles[n].pos[2]/L)*L;
            }
            break;
        //second integration
        case 2:
            for(n = 0; n < N; n++){
                //next velocities
                particles[n].vel[0] += particles[n].force[0]*dt;
                particles[n].vel[1] += particles[n].force[1]*dt;
                particles[n].vel[2] += particles[n].force[2]*dt;

                //second middle time interval displacement
                xx = particles[n].vel[0]*dt/2;
                yy = particles[n].vel[1]*dt/2;
                zz = particles[n].vel[2]*dt/2;

                //next positions
                particles[n].pos[0] += xx;
                particles[n].pos[1] += yy;
                particles[n].pos[2] += zz;

                //sum of velocities
                sumvx += particles[n].vel[0];
                sumvy += particles[n].vel[1];
                sumvz += particles[n].vel[2];

                //sum of squared velocities
                sumv2 += (pow(particles[n].vel[0], 2) + pow(particles[n].vel[1], 2) + pow(particles[n].vel[2], 2));
            }

            //check on center of mass movements
            if((-1 > sumvx || sumvx > 1) || (-1 > sumvy || sumvy > 1) || (-1 > sumvz || sumvz > 1))
                printf("\nCAUTION: the center of mass velocity is not zero.\nvx = %f vy = %f vz = %f\n", sumvx, sumvy, sumvz);

            //kinetic energy and temperature
            K = sumv2/(2*N);
            T = (2/3.)*K;

            break;
    }
}

void savePositions(int time, Particle *particles){

    int n;
    FILE *file_pos;
    char filename[30] = {0};
    sprintf(filename, "%s//pos%i.txt", "./positions", time);
    file_pos = fopen(filename, "w");
    fprintf(file_pos, "%i\n", N);
    fprintf(file_pos, "Box at time %i\n", time);
    for(n = 0; n < N; n++)
        fprintf(file_pos, "%c %f %f %f\n", particles[n].type, particles[n].pos[0], particles[n].pos[1], particles[n].pos[2]);
    fclose(file_pos);
}
