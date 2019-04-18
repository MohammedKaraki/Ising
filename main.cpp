#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::array;
using std::exp;
using std::tanh;
using std::log;
using std::abs;


// for random number generator
#define FIXED_SEED 139933

struct coordinates
{
    int x, y, t;

    coordinates(int x_p, int y_p, int t_p) 
        : x(x_p), y(y_p), t(t_p)
    {
    }

    friend std::ostream& operator<< (std::ostream& out, const coordinates& c)
    {
        return out << "(" << c.x << ", " << c.y << ", " << c.t << ")";
    }
};


class square_lattice
{
private:
    const int side_length;  /* number of lattice sites on each side of
                               the square */
    const int replica_count;    /* Divisions of 'time' dimension */
    double beta; // inverse temperature
    double mag_field;

    double Gamma,   // = \eta / \delta
           delta; // = \beta / N



    vector<int> config; 


private:
    std::mt19937 rand_gen;
    std::uniform_int_distribution<> bin_dist; /* 0 or 1 (to initialize 
                                                 'config' randomly) */
    std::uniform_int_distribution<> site_dist; /* 0, 1, ..., site_count-1
                                                  (to sample one site) */
    std::uniform_real_distribution<> prob_dist; /* Sampling from [0., 1.] 
                                       to make probabilistic decisions */


    // an attempt to increase efficiency by eliminating the need to calculate 
    // the indices of the 6 neighbors every time needed.
    // It stores a list of nearest neighbors for each site.
    vector< array<int, 6> > nbr_table;

public:
    const int site_count; // = side_length^2 * replica_count

    square_lattice(int side_length_p, int replica_count_p)
        : side_length(side_length_p), 
        replica_count(replica_count_p),
        site_count(side_length_p * side_length_p * replica_count_p)
    {
        rand_gen = std::mt19937(FIXED_SEED);
        bin_dist = std::uniform_int_distribution<> (0, 1);
        site_dist = std::uniform_int_distribution<> (0, site_count - 1);
        prob_dist = std::uniform_real_distribution<> (0., 1.);




        // Filling the table of neighbors. It uses PBC in all three directions.
        nbr_table.resize(site_count);
        for (int i_x = 0; i_x < side_length; ++i_x) {
            for (int i_y = 0; i_y < side_length; ++i_y) { 
                for (int i_t = 0; i_t < replica_count; ++i_t) {

                    int index = coord_to_index( coordinates(i_x, i_y, i_t) );

                    auto relative_index = [&, this] (int dx, int dy, int dt) {
                        int j_x = i_x + dx;
                        int j_y = i_y + dy;
                        int j_t = i_t + dt;

                        if (j_x < 0) { j_x += side_length; }
                        if (j_y < 0) { j_y += side_length; }
                        if (j_t < 0) { j_t += replica_count; }

                        j_x %= side_length;
                        j_y %= side_length;
                        j_t %= replica_count;

                        return coord_to_index( coordinates(j_x, j_y, j_t) );
                    };

                    nbr_table[index][0] = relative_index(+1, 0, 0);
                    nbr_table[index][1] = relative_index(-1, 0, 0);
                    nbr_table[index][2] = relative_index(0, +1, 0);
                    nbr_table[index][3] = relative_index(0, -1, 0);
                    nbr_table[index][4] = relative_index(0, 0, +1);
                    nbr_table[index][5] = relative_index(0, 0, -1);

                }
            }
        }


    }

    void fill_random()
    {
        config.resize(site_count);

        for (auto it = config.begin(); it != config.end(); ++it) {
            *it = bin_dist(rand_gen) * 2 - 1;
        }
    }

    void set_params(double beta_p, double mag_field_p) 
    {
        beta = beta_p;
        mag_field = mag_field_p;

        delta = beta / replica_count;

        double eta = -.5 * log(tanh(delta * mag_field));
        Gamma = eta / delta;
    }

    // Wolff's algorithm
    int cluster_update()
    {
        vector<int> cluster, old_set, new_set;

        int initial_site = site_dist(rand_gen);
        cluster.push_back(initial_site);
        old_set.push_back(initial_site);


        vector<bool> on_cluster; /* using which, it takes O(1) to check
                                    if a given site belongs to the cluster; no
                                    need to do a search every time. */
        on_cluster.resize(site_count);
        for (auto it = on_cluster.begin(); it != on_cluster.end(); ++it) {
            *it = false;
        }
        on_cluster[initial_site] = true;


        // Return true with probability p
        auto prob_true = [this](double p) {
            if (prob_dist(rand_gen) < p) { 
                return true;
            }
            return false;
        };


        int seed_spin = config[initial_site];
        const double p_xy = 1. - exp(-2. * delta);
        const double p_t = 1. - exp(-2. * delta * Gamma);


        while (!old_set.empty()) {
            new_set.clear();
            for (auto it = old_set.begin(); it != old_set.end(); ++it) {
                for (auto j = 0; j < 6; ++j){

                    int nbr_index = nbr_table[*it][j];
                    if ( config[nbr_index] == seed_spin &&
                            on_cluster[nbr_index] == false) {
                        double p = j <= 3 ? p_xy : p_t; 
                        if (prob_true(p)) { 
                            new_set.push_back(nbr_index);
                            cluster.push_back(nbr_index);
                            on_cluster[nbr_index] = true;
                        }
                    }
                }
            }

            old_set = new_set;
        }

        if (cluster.size() > 1) {
            for (auto it = cluster.begin(); it != cluster.end(); ++it) {
                config[*it] *= -1;
            }
        } else {
            array<int, 6>& nbrs = nbr_table[initial_site];
            int nbr_spin_xy = config[nbrs[0]] + config[nbrs[1]]
                + config[nbrs[2]] + config[nbrs[3]];
            int nbr_spin_t = config[nbrs[4]] + config[nbrs[5]];

            int site_spin = config[initial_site];

            double delta_E = 2. * nbr_spin_xy * site_spin 
                + 2. * Gamma * nbr_spin_t * site_spin;

            double p = exp(-delta * delta_E);
            if (prob_true(p)) {
                config[initial_site] *= -1;
            }
        }
        return cluster.size();
    }


    double energy_density()
    {
        double total_E = 0.;
        for (int index = 0; index < site_count; ++index) {
            const array<int, 6>& nbrs = nbr_table[index];
            total_E += -config[index] * (config[nbrs[0]] + config[nbrs[2]]
                    + Gamma * config[nbrs[4]]);
        }
        return total_E / site_count;
    }

    double magnetization_density()
    {
        return abs(std::accumulate(config.begin(), config.end(), 0.))
            / site_count;
    }


    int coord_to_index(const coordinates& c)
    {
        return 
            c.x + 
            c.y * side_length + 
            c.t * side_length * side_length;
    }

    coordinates index_to_coord(const int i)
    {
        int x, y, t;
        x = i % side_length;
        y = ((i - x) / side_length) % side_length;
        t = i / (side_length * side_length);

        return coordinates(x, y, t);
    }

};


double mean(const vector<double>& v) 
{
    return std::accumulate(v.begin(), v.end(), 0.) / v.size();
}

double std_dev(const vector<double>& v)
{
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.);
    double avg = mean(v);
    return std::sqrt(sq_sum / v.size() - avg*avg);
}

#define FORMATTED(a,b,c,d,e,f) a << std::left << std::setprecision(5) << \
    std::setw(10) << b << std::setw(10) << c << std::setw(15) << d << \
    std::setw(10) << e << std::setw(10) << f << endl

int main()
{
    square_lattice test(10, 10);

    int sample_count = 50;
    int measure_count = 100;

    double beta = 2.;


    auto t1 = std::chrono::high_resolution_clock::now();

    FORMATTED(cout, "field", "E/V","\\sigma_{E/V}", "m", "\\sigma_{m}");

    for (double mag_field = .1; mag_field < sample_count * .1; 
            mag_field += .1) {
        test.set_params(beta, mag_field);
        test.fill_random();
        for (int j = 0; j < 2000; j++) {
            test.cluster_update();
        }

        vector<double> energies;
        vector<double> magnetizations;
        for (int j = 0; j < measure_count; ++j) {
            energies.push_back(test.energy_density());
            magnetizations.push_back(test.magnetization_density());

            int cluster_size = test.cluster_update();
            int interv = 2 * int(test.site_count / cluster_size);
            for (int k = 0; k < interv; ++k) {
                test.cluster_update();
            }
        }


        FORMATTED(cout, mag_field, mean(energies), std_dev(energies),
                mean(magnetizations), std_dev(magnetizations));
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    cout << std::chrono::duration<double>(t2-t1).count() << "s" << endl;

    cout << "... done" << endl;
    return 0;
}

