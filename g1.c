#include <SDL2/SDL.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*******************************************
  predefining stuff (Konstanten, typedef)
********************************************/


#define G 6.6743e-11
#define M_sol 20		// Sonnenmasse
#define M_bh1 6e8		// Masse Schwarzes Loch (1)
#define M_bh2 8e7		// Masse Schwarzes Loch (2)
#define N 500			// num_particles
#define HEIGHT 600		// Fensterhöhe
#define WIDTH 800		// Fensterbreite

typedef struct {
   double x;
   double y;
} Vector2D;

typedef struct {
   double m;
   double x;
   double y;
   double vx;
   double vy;
   double fx;
   double fy;
} Particle;



/**************
  Subroutines
***************/


// Kräfte zurücksetzen in Integration

void clear_forces(Particle* p) {
   for (int i=0; i < N; i++) {
      p[i].fx = 0;
      p[i].fy = 0;
   }
}


// Kraft berechnen

double calc_force(Particle* p, Particle* b) {
   for (int i=0; i < N; i++) {
      double dist_x = p[i].x - b->x;
      double dist_y = p[i].y - b->y;
      double r = sqrt(dist_x * dist_x + dist_y * dist_y);
      double fg = G * p[i].m * b->m / (r*r);
      p[i].fx += -fg * dist_x / r;
      p[i].fy += -fg * dist_y / r;
   }
   return 0;
}


// Runge-Kutta4 Schritt

void rk4_step(Particle* p) {
   int i;
   double dt = 0.05;
   double t;
   Particle k1[N], k2[N], k3[N], k4[N], k_g[N];

   // k1
   for (i = 0; i < N; i++) {
      k1[i].x = p[i].vx;
      k1[i].y = p[i].vy;
      k1[i].vx = p[i].fx / p[i].m;
      k1[i].vy = p[i].fy / p[i].m;
   }

   // k2
   for (i = 0; i < N; i++) {
      k2[i].x = p[i].vx + 0.5 * dt * k1[i].vx;
      k2[i].y = p[i].vy + 0.5 * dt * k1[i].vy;
      k2[i].vx = (p[i].fx + 0.5 * dt * k1[i].vx) / p[i].m;
      k2[i].vy = (p[i].fy + 0.5 * dt * k1[i].vy) / p[i].m;
   }

   // k3
   for (i = 0; i < N; i++) {
      k3[i].x = p[i].vx + 0.5 * dt * k2[i].vx;
      k3[i].y = p[i].vy + 0.5 * dt * k2[i].vy;
      k3[i].vx = (p[i].fx + 0.5 * dt * k2[i].vx) / p[i].m;
      k3[i].vy = (p[i].fy + 0.5 * dt * k2[i].vy) / p[i].m;
   }

   // k4
   for (i = 0; i < N; i++) {
      k4[i].x = p[i].vx + dt * k3[i].vx;
      k4[i].y = p[i].vy + dt * k3[i].vy;
      k4[i].vx = (p[i].fx + dt * k3[i].vx) / p[i].m;
      k4[i].vy = (p[i].fy + dt * k3[i].vy) / p[i].m;
   }

   // Endgültiger Schritt
   for (i = 0; i < N; i++) {
      p[i].x += (dt / 6.0) * (k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x);
      p[i].y += (dt / 6.0) * (k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y);
      p[i].vx += (dt / 6.0) * (k1[i].vx + 2 * k2[i].vx + 2 * k3[i].vx + k4[i].vx);
      p[i].vy += (dt / 6.0) * (k1[i].vy + 2 * k2[i].vy + 2 * k3[i].vy + k4[i].vy);
   }
}


// Anfangspositionen für Partikel

void generate_particles(Particle** p){
	*p =malloc(sizeof(Particle) * N);
	Vector2D bh_pos = {WIDTH / 2, HEIGHT / 2}; // position of the black hole
	double r0 = 120.0; // maximum distance of particles from black hole
	double a = 20.0; // scale length of the Burbidge equation distribution
	double rho0 = 1.1; // central density of the Burbidge equation distribution

for (int i = 0; i < N; i++) {
    double r = 5 + r0 * drand48(); // distance from black hole
    double theta = 2.0 * M_PI * drand48();

    double rho = rho0 * exp(-pow(r/a, 0.25)); // density at distance r

    double m = rho * (4.0/3.0)*M_PI*pow(r+1e-10,3.0); // mass of the particle

    (*p)[i].x = bh_pos.x + r * cos(theta);
    (*p)[i].y = bh_pos.y + r * sin(theta);

    double dx = (*p)[i].x - bh_pos.x;
    double dy = (*p)[i].y - bh_pos.y;
    double v = sqrt(G * M_bh1 / r );
    double phi = atan2(dy, dx);
    (*p)[i].vx = -v*sin(phi);
    (*p)[i].vy = v*cos(phi);

    (*p)[i].m = m;
}
}


// Anfangsbedingung zweites Schwarzes Loch

void generate_bh2(Particle** p){
    *p =malloc(sizeof(Particle) * 16);

    Vector2D bh_pos = {WIDTH / 2, HEIGHT / 2};

    double r = 150; // distance from black hole
    (*p)->x = bh_pos.x + r;
    (*p)->y = bh_pos.y + r;

    double dx = (*p)->x - bh_pos.x;
    double dy = (*p)->y - bh_pos.y;
    double v = sqrt(G * M_bh1 / r );
    double phi = atan2(dy, dx);
    (*p)->vx = -v*sin(phi);
    (*p)->vy = v*cos(phi);

    (*p)->m = M_bh2;
}


// Schwarzes Loch1

void generate_bh1(Particle** bh_1){
    *bh_1 =malloc(sizeof(Particle) * 16);

    (*bh_1)->x = WIDTH / 2;
    (*bh_1)->y = HEIGHT / 2;
    (*bh_1)->vx = 0.0;
    (*bh_1)->vy = 0.0;
    (*bh_1)->m = M_bh1;
}


// Partikel dartellen

void draw_particle(SDL_Renderer* renderer, Particle* p) {
   SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
   SDL_Rect rect = {p->x, p->y, 1, 1};
   SDL_RenderFillRect(renderer, &rect);
}




/******************************************
---------------MAINPART--------------------
******************************************/

int main() {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Particles", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    Particle* p;
    Particle* bh_2;
    Particle* bh_1;

    generate_particles(&p);
    generate_bh2(&bh_2);
    generate_bh1(&bh_1);
    double dt = 2e-5;

    while (1) {
        SDL_Event e;
        if (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                break;
            }
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        Vector2D bh_pos = {WIDTH / 2, HEIGHT / 2}; // position of the black hole

        for (int i = 0; i < N; i++) {
            clear_forces(p);
            clear_forces(bh_2);
            clear_forces(bh_1);

            // Update position
            calc_force(p, bh_2);
            calc_force(p, bh_1);
            rk4_step(p);

	    // Update BH2
	    calc_force(bh_2, bh_1);
            calc_force(bh_2, p);
            rk4_step(bh_2);

            // Update BH1
            calc_force(bh_1, bh_2);
            calc_force(bh_1, p);
            rk4_step(bh_1);

            // draw particles
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_Rect rect = {(int)p[i].x, (int)p[i].y, 1, 1};
            SDL_RenderFillRect(renderer, &rect);
          }


	// Draw blach_hole2
	SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        SDL_Rect bh_rect2 = {(int)bh_2->x, (int)bh_2->y, 2, 2};
        SDL_RenderFillRect(renderer, &bh_rect2);


        // Draw black hole
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        SDL_Rect bh_rect = {(int)bh_1->x, (int)bh_1->y, 2, 2};
        SDL_RenderFillRect(renderer, &bh_rect);

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
