#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

// Parameters
const double L = 300;                       // System length
const double rho = 0.05;                     // Particle density
const int N = static_cast<int>(L * rho);    // Number of particles
const double vlim = 0.1;                    // Absolute value of the limits of initial velocities distribution
const double v_0 = 0.1;                     // Desired average velocity
const double eta = 2.0;                     // Noise amplitude
const double delta = 10.0;                    // Vision distance (neighborhood)
const int t_f = 3000;                       // Number of simulation steps
const double d_min = 0.1;                   // Minimum distance between particles to avoid overlap

// G function based on the article
double G(double u) {
    if (u >= 0) {
        return (u + 1) / 2;
    } else {
        return (u - 1) / 2;
    }
}

// Apply periodic boundary conditions
double apply_periodic_boundary(double x, double L) {
    if (x >= L) return x - L;
    if (x < 0) return x + L;
    return x;
}

// Compute the average velocity in the neighborhood of particle i
double average_velocity(int i, const vector<double>& positions, const vector<double>& velocities, double delta, double L) {
    vector<double> neighbors;
    for (int j = 0; j < N; ++j) {
        double dist = fabs(positions[i] - positions[j]);
        dist = min(dist, L - dist);  // Account for periodic boundary
        if (dist < delta) {
            neighbors.push_back(velocities[j]);
        }
    }
    double sum = accumulate(neighbors.begin(), neighbors.end(), 0.0);
    return sum / neighbors.size();
}

// Ensure no overlap of particles
void avoid_overlap(vector<double>& X, int i, double d_min, double L) {
    if (i > 0) {
        double dist_prev = X[i] - X[i - 1];
        if (dist_prev < 0) dist_prev += L;  // Correct with periodic boundary
        if (dist_prev < d_min) {
            X[i] = X[i - 1] + d_min;
            X[i] = apply_periodic_boundary(X[i], L);
        }
    }
}

int main() {
    // Random number generation setup
    mt19937 rng(2906);  // Seed
    uniform_real_distribution<double> pos_dist(0, L);
    uniform_real_distribution<double> vel_dist(-2 * vlim, 2 * vlim);
    uniform_real_distribution<double> noise_dist(-eta / 2, eta / 2);

    // Initialize position and velocity arrays
    vector<vector<double>> X(t_f, vector<double>(N));
    vector<vector<double>> U(t_f, vector<double>(N));

    // Initial positions are random and sorted
    vector<double> initial_positions(N);
    for (int i = 0; i < N; ++i) {
        initial_positions[i] = pos_dist(rng);
    }
    sort(initial_positions.begin(), initial_positions.end());
    X[0] = initial_positions;

    // Initial velocities are random
    for (int i = 0; i < N; ++i) {
        U[0][i] = vel_dist(rng);
    }

    // Simulation loop
    for (int k = 1; k < t_f; ++k) {
        for (int j = 0; j < N; ++j) {
            // Compute the average velocity of neighbors
            double u_prom = average_velocity(j, X[k-1], U[k-1], delta, L);

            // Update velocity using G(u) and adding uniform noise
            double noise = noise_dist(rng);
            U[k][j] = G(u_prom) + noise;

            // Update the position based on the velocity
            X[k][j] = X[k-1][j] + v_0 * U[k][j];

            // Apply periodic boundary conditions to the positions
            X[k][j] = apply_periodic_boundary(X[k][j], L);

            // Ensure no particle overlap
            avoid_overlap(X[k], j, d_min, L);
        }
    }

    // Save the positions and times to a file for Python plotting
    ofstream outFile("data/delta/rho0p5/10.txt");    
    for (int k = 0; k < t_f; ++k) {
        for (int j = 0; j < N; ++j) {
            outFile << X[k][j] << " ";
        }
        outFile << "\n";
    }
    outFile.close();

    return 0;
}