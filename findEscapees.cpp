/*============================================================================
 *
 *  findEscapees.cpp:  a tool to determine which stars have sufficient energy
 *                    to escape the cluster.
 *                    The program accepts a stream of N-body snapshots, and
 *                    outputs a list of escaping stars for each snapshot.
 *
 *  Note: in this first version, all functions are included in one file,
 *        without any use of a special library or header files.
 *_____________________________________________________________________________
 *
 *  Usage: findEscapees [-h (for help)]
 *
 *         Input/output are read/written from/to the standard i/o streams.
 *         Since there are no options, the code can simply be run by
 *         specifying the input file for the N-body snapshots:
 *
 *            findEscapees < data.in
 *
 *         This will produce data on the screen.  In order to capture the data
 *         in an output file "escapinglist.out", use:
 *
 *            findEscapees < data.in > escapinglist.out
 *_____________________________________________________________________________

 *
 *  External data format:
 *
 *     The program expects input of a single snapshot of an N-body system,
 *     in the following format: the number of particles in the snapshot n;
 *     the time t; mass mi, position ri and velocity vi for each particle i,
 *     with position and velocity given through their three Cartesian
 *     coordinates, divided over separate lines as follows:
 *
 *                      n
 *                      t
 *                      m1 r1_x r1_y r1_z v1_x v1_y v1_z
 *                      m2 r2_x r2_y r2_z v2_x v2_y v2_z
 *                      ...
 *                      mn rn_x rn_y rn_z vn_x vn_y vn_z
 *
 *     Output for each snapshot consists of a list of escaping stars, each one
 *     defined by the identity {id1,id2} of its members, as follows:
 *
 *                      time = t
 *                      star1 = id1  star2 = id2
 *
 *     Each line "time = t" is followed by zero, one, or more escaping star
 *     listings, depending on how many are found in the corresponding
 *     snapshot.
 *
 *  Internal data format:
 *
 *     The data for an N-body system is stored internally as 1-dimensional
 *     arrays for the masses and densities, and 2-dimensional arrays for the
 *     positions and velocities of all particles.
 *_____________________________________________________________________________
 *
 *    Version 1:  Mar 2019 R Norman
 *_____________________________________________________________________________
 */

#include  <iostream>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <unistd.h>                       // for getopt()
using namespace std;

typedef double  real;                      // "real" as a general name for the
                                           // standard floating-point data type

const int NDIM = 3;                        // number of spatial dimensions

bool read_options(int argc, char *argv[]);
int get_snapshot(real * * mass, real (* * pos)[NDIM], real (* * vel)[NDIM],
		 int & n, real & t);
void delete_snapshot(const real mass[], const real pos[][NDIM],
		     const real vel[][NDIM]);
bool find_escapees(real mass[], real pos[][NDIM], real vel[][NDIM], int n,
		   real t);

/*-----------------------------------------------------------------------------
 *  main  --  reads option values, and starts a loop; in each round of the loop
 *            a new shapshot is read, the escaping stars are found and reported.
 *-----------------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
    if (! read_options(argc, argv))
        return 1;                // halt criterion detected by read_options()

    real * mass;                 // masses for all particles
    real (* pos)[NDIM];          // positions for all particles
    real (* vel)[NDIM];          // velocities for all particles

    int n;                       // N, number of particles in the N-body system
    real t;                      // time

    while(get_snapshot(&mass, &pos, &vel, n, t)){
        if (! find_escapees(mass, pos, vel, n, t))
	    return 1;
	delete_snapshot(mass, pos, vel);
    }
    return 0;
}

/*-----------------------------------------------------------------------------
 *  read_options  --  reads the command line options, and implements them.
 *
 *  note: when the help option -h is invoked, the return value is set to false,
 *        to prevent further execution of the main program; similarly, if an
 *        unknown option is used, the return value is set to false.
 *-----------------------------------------------------------------------------
 */

bool read_options(int argc, char *argv[])
{
    int c;
    while ((c = getopt(argc, argv, "h")) != -1)
        switch(c){
            case 'h':
            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << endl;
                      return false;      // execution stops after help or error
            }

    return true;                         // ready to continue program execution
}

/*-----------------------------------------------------------------------------
 *  get_snapshot  --  reads a single snapshot from the input stream cin
 *
 *  note: memory allocation for masses, positions and velocities is done here
 *        after reading in the number of particles (n).
 *        If the end of file is reached, get_snapshot() returns 0;
 *        after successful completion, get_snapshot() returns 1.
 *-----------------------------------------------------------------------------
 */

int get_snapshot(real * * mass, real (* * pos)[NDIM], real (* * vel)[NDIM],
		 int & n, real & t)
{
    cin >> n;
    if (cin.fail())
        return 0;
    cin >> t;

    *mass = new real[n];                  // masses for all particles
    *pos = new real[n][NDIM];             // positions for all particles
    *vel = new real[n][NDIM];             // velocities for all particles

    for (int i = 0; i < n ; i++){
        cin >> (*mass)[i];                       // mass of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> (*pos)[i][k];                 // position of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> (*vel)[i][k];                 // velocity of particle i
    }
    return 1;
}

/*-----------------------------------------------------------------------------
 *  delete_snapshot  --  frees up the memory that was allocated to the masses,
 *                       positions, and velocities for the particles in a
 *                       snapshot.
 *-----------------------------------------------------------------------------
 */

void delete_snapshot(const real mass[], const real pos[][NDIM],
		     const real vel[][NDIM])
{
    delete[] mass;
    delete[] pos;
    delete[] vel;
}

/*-----------------------------------------------------------------------------
 *  find_escapees  --   For each star, calculate the total gravitational potential
 *                      at the position of the star due to the rest of the cluster.
 *                      If the KE of the star is greater than the GPE at that
 *                       point then the star will escape.
 *-----------------------------------------------------------------------------
 */

bool find_escapees(real mass[], real pos[][NDIM], real vel[][NDIM], int n,
		   real t)
{
    cout << "time = " << t << endl;
    for (int i = 0; i < n; i++){    //Check each star
        //Calculate the grav potential at star due to rest of cluster
        real grav_pot = 0.0;
        for (int j = 0; j < n ; j++){
            if(j!=i){
                real delr[NDIM];
                real delr_sq = 0.0;
                for (int k = 0; k < NDIM ; k++){
                    delr[k] = pos[j][k]-pos[i][k];
                    delr_sq += delr[k] * delr[k];
                }
                real r = sqrt(delr_sq);
                grav_pot += mass[j]/r;
            }
        }
        //Calculate KE/m of particle
        real v_sq = 0.0;
        for(int k=0; k<NDIM; k++){
            v_sq += vel[i][k]*vel[i][k];
        }
        real ke = 0.5*v_sq;
        //Determine if star is escaping
        if(ke>grav_pot){
            cout << "Star" << i << endl;
        }
    }
    return true;
}

/*-----------------------------------------------------------------------------
 *                                                                    \\   o
 *  end of file:  find_escapees2.C                                    /\\'  O
 *                                                                   /\     |
 *=============================================================================
 */
