/*=============================================================================
 *
 *  createCluster.cpp:  generates a file with the initial conditions for a
 *                      sphere of stars.
 *
 *_____________________________________________________________________________
 *
 *  usage: createCluster
 *                  -n number_of_stars
 *                  -m Type of mass distribution
 *                          1   -[Default] All stars have the equal mass = 1/n
 *                          2   -Stars have a random, but equal probability, mass
 *                               between mass_min and mass_max
 *                          3   -Salpeter (one stage) distribution
 *                          4   -Salpeter (two stage) distribution
 *                          5   -Kroupa, Tout & Gilmore distribution
 *                          6   -Unused at this stage
 *                  -p Type of position distribution
 *                          1   -[Default] All positions are scattered
 *                               homogeneously across the cluster.
 *                  -v  Type of velocity distribution
 *                          1   -[Default] All velocities are zero
 *                    [-s random_number_generator_seed]
 *                    [-r turn off mass rescaling]
 *                    [-h (for help)]
 *
 *         The number of particles has to be specified, since there is no
 *         natural default value.  If no seed is specified for the random
 *         number generator, a random value for the seed is chosen which
 *         depends on the unix clock and will be different every second.
 *
 *         Example:  "./createCluster -n 3 > data.out"  will produce a file in the
 *         following format, for particle number, time, and masses mi,
 *         positions ri, and velocities vi for particles i:
 *
 *                      3
 *                      0
 *                      m1 r1_x r1_y r1_z v1_x v1_y v1_z
 *                      m2 r2_x r2_y r2_z v2_x v2_y v2_z
 *                      m3 r3_x r3_y r3_z v3_x v3_y v3_z
 *_____________________________________________________________________________
 *
 *  Our sphere has unit radius and density determined by the chosen distribution type.
 *  The total mass is unity unless the 'rescaled (-r)' option is chosen in which case
 *  the stars have masses in multiples of the Sun's mass - this is for diagnostic
 *  purposes as the nbody routine expects the total mass of the cluster to be unity.
 *
 *  Each star will be sprinkled somewhere randomly within the unit sphere,
 *  with initial velocity determined by the velocity distribution option (-v).
 *_____________________________________________________________________________
 *
 *    version 1:    Jan 2002   Piet Hut, Jun Makino
 *    version 2:    Oct 2019   Robert Norman
 *____________________________________________________________________________*/

#include  <iostream>
#include  <cmath>                   // to include sqrt(), etc.
#include  <cstdlib>                 // formerly <cstdlib.h>; for atoi(), atof()
                                    // and rand(), srand()
#include  <unistd.h>                // for getopt()
#include  <time.h>                  // for time()
using namespace std;

typedef double  real;                      // "real" as a general name for the
                                           // standard floating-point data type

const int NDIM = 3;                        // number of spatial dimensions

bool read_options(int argc, char *argv[], int &n, int &massDistribType, int &seed, bool &rescale,
                    int &posDistribType, int &velDistribType);
real randunit(int seed);
real randinter(real a, real b);
void masses(real mass[], int distType, bool rescale, real mass_min, real mass_max, real *mass_cluster, int n);
void massDistrib_equal(int n, real mass[], real mass_min, real mass_max, real *mass_cluster);
void massDistrib_uniform(int n, real mass[], real mass_min, real mass_max, real *mass_cluster);
void massDistrib_Salpeter1(int n, real mass[], real mass_min, real mass_max, real *mass_cluster);
void massDistrib_Salpeter2(int n, real mass[], real *mass_cluster);
void massDistrib_Kroupa_etal(int n, real mass[], real mass_min, real mass_max, real *mass_cluster);
void massDistrib_Kroupa2(int n, real mass[], real mass_min, real mass_max, real *mass_cluster);
void positions(int n, real pos[][NDIM], int posDistribType);
void posDistrib_Uniform(int n, real pos[][NDIM]);
void velocities(int n, real vel[][NDIM], int velDistribType);
void velDistrib_Zero(int n, real vel[][NDIM]);
void velDistrib_Random(int n, real vel[][NDIM]);
void put_snapshot(real mass[], real pos[][NDIM],
                  real vel[][NDIM], int n, real t, real *mass_cluster);

/*-----------------------------------------------------------------------------
 *  main  --  read in option values, invoke the model builder
 *-----------------------------------------------------------------------------
 */
int main(int argc, char *argv[])
{
    int n = 0;            // N, number of particles; start with an
    				      // unphysical value, to allow us to force
    				      // the user to make a definite choice.

    int seed = 0;         // seed for the random number generator;
    				      // a default of zero will be replaced by
    				      // the time taken from the unix clock.

    int massDistribType = 0;    //Type of mass distribution required
    bool rescale = true;        //Rescale the star masses to cluster mass of unity
    int posDistribType = 1;     //Type of position distribution required.
    int velDistribType = 1;     //Type of velocity distribution required.

    real mass_min = 0.0;   //Default values for the mass limits of the stars.
    real mass_max = 1.0;
    real *mass_cluster = new real;      //Total mass of the cluster in solar masses.

    if (! read_options(argc, argv, n, massDistribType, seed, rescale, posDistribType, velDistribType))
        return 1;

    if (n <= 0){
	cerr << "a value of N = " << n << " is not allowed." << endl;
	return 1;
    }

    if (seed == 0)       /* no particular positive seed provided?            */
        seed = time(0);  /* then give a random value, different every second */

    cerr << "seed = " << seed << endl;

    randunit(seed);

    real * mass = new real[n];                  // masses for all particles
    real (* pos)[NDIM] = new real[n][NDIM];     // positions for all particles
    real (* vel)[NDIM] = new real[n][NDIM];     // velocities for all particles

    masses(mass, massDistribType, rescale, mass_min, mass_max, mass_cluster, n);
    positions(n, pos, posDistribType);
    velocities(n, vel, velDistribType);

    put_snapshot(mass, pos, vel, n, 0, mass_cluster);

    delete mass;
    delete pos;
    delete vel;

}

/*==============================================================================
 *  read_options  --  reads the command line options, and implements them.
 *
 *  note: when the help option -h is invoked, the return value is set to false,
 *        to prevent further execution of the main program; similarly, if an
 *        unknown option is used, the return value is set to false.
 *------------------------------------------------------------------------------
 */

bool read_options(int argc, char *argv[], int &n, int &massDistribType, int &seed, bool &rescale, int &posDistribType, int &velDistribType)
{
    int c;
    while ((c = getopt(argc, argv, "hn:m:s:rp:v:")) != -1)       //colon indicates the option has an argument following
        switch(c){
            case 'n': n = atoi(optarg);
                      break;
            case 'm': massDistribType = atoi(optarg);
                      break;
            case 's': seed = atoi(optarg);
                      break;
            case 'r': rescale = false;
                      break;
            case 'p': posDistribType = atoi(optarg);
                      break;
            case 'v': velDistribType = atoi(optarg);
                      break;
            case 'h':
            case '?': cerr << "usage: " << argv[0] << endl
                           << "         -n number_of_particles\n"
                           << "         -m mass distribution_type\n"
                           << "         -p position distribution type\n"
                           << "         [-s random_number_generator_seed]\n"
                           << "         [-r turn off mass rescaling]\n"
                           << "         [-h (for help)]\n"
                           << endl;
                      return false;      // execution stops after help or error
            }
    if(n==0 || massDistribType==0 || massDistribType>6){
        cerr << endl
             << "You must enter the number of stars in the system,\n"
             << "      -n num_of_stars\n"
             << "and, a valid distribution type\n"
             << "      -d dist_type ( 1 .. 5 )"
             << endl << endl;
        return false;
    }

    return true;                         // ready to continue program execution
}

/*==============================================================================
 *  randunit  --  returns a random real number within the unit interval
 *                note: based on      @(#)rand.c   4.1 (Berkeley) 12/21/80,
 *                      but returning a positive number smaller than unity.
 *
 *  note: to initialize the random number generator, invoke it with an nonzero
 *        argument, which will then become the seed;
 *        to run the random number generator, invoke it with argument 0.
 *-----------------------------------------------------------------------------
 */
real randunit(int seed)
{
    const real MAXN = 2147483647;  // the maximum value which rand() can return

    static int randx;

    if (seed)
        {
        randx = seed;
        return(0.0);        // to make the compiler happy, we return a value,
        }                   // even though it will not be used in this case
    else
        return((real)((randx = randx * 1103515245 + 12345) & 0x7fffffff)/MAXN);
}

/*=================================================================================
 *  randinter  --  returns a random real number within an interval [a,b]
 *                 by invoking  randunit() .
 *-----------------------------------------------------------------------------
 */
real  randinter(real a, real b)
    {
    return(a + (b-a)*randunit(0));
    }

/*=================================================================================
*   masses  -- calculates the masses of the stars.
*
*   The calculation of the masses is a two step process:
*      1.  The N masses are calculated as a random number between 'mass_min' and
*          'mass_max'.  These units of these masses are multiples of the Sun's mass.
*      2.  All masses are rescaled to: m[i]/mass_cluster.  This makes the
*          total mass of the cluster equal to unity.
*
*   The total mass of the cluster (mass_cluster) can be recovered and used to
*   rescale the masses back to real world units.
-------------------------------------------------------------------------------
*/
void masses(real mass[], int distType, bool rescale, real mass_min, real mass_max, real *mass_cluster, int n)
{
    //All masses are initially calculated in multiples of the sun's mass and then
    //rescaled at the end.
    switch(distType) {
        case 1: //Equal masses
                massDistrib_equal(n, mass, mass_min, mass_max, mass_cluster);
                break;

        case 2: //Uniformly random between mass_min and mass_max
                //Calculate the masses of the particles in multiples of the sun's mass
                massDistrib_uniform(n, mass, mass_min, mass_max, mass_cluster);
                break;
        case 3: //Salpeter distribution
                massDistrib_Salpeter1(n, mass, mass_min, mass_max, mass_cluster);
                break;
        case 4: //Salpeter two stage distribution
                massDistrib_Salpeter2(n, mass, mass_cluster);
                break;
        case 5: //Kroupa, et al distribution
                massDistrib_Kroupa_etal(n, mass, mass_min, mass_max, mass_cluster);
                break;
        case 6:
                break;
    }
        if(rescale) {
            //Rescale the masses so the total cluster mass is unity
            for(int i = 0; i < n; i++)
                mass[i] = mass[i]/(*mass_cluster);
        }

    cerr << "Total mass of cluster = " << *mass_cluster << " solar masses" << endl;
}

// Mass distribution functions --------------------------------------------------------

//Tyoe 1 - Equal masses
void massDistrib_equal(int n, real mass[], real mass_min, real mass_max, real *mass_cluster)
{
    real mass_tot = 0.0;

    for (int i = 0; i < n; i++) {
        mass[i] = 1.0;
        mass_tot += mass[i];
    }
    *mass_cluster = mass_tot;
}

//Type 2 - Uniform masses
void massDistrib_uniform(int n, real mass[], real mass_min, real mass_max, real *mass_cluster)
{
    real mass_tot = 0.0;

    for (int i = 0; i < n; i++) {
        mass[i] = randinter(mass_min, mass_max);
        mass_tot += mass[i];
    }
    *mass_cluster = mass_tot;
}

//Type 3 - Salpeter distribution
//Uses Salpeter's classical power law IMF p(m) = 0.53*m^2.3
void massDistrib_Salpeter1(int n, real mass[], real mass_min, real mass_max, real *mass_cluster)
{
    real mass_tot = 0.0;
    const double alpha = 2.3;
    const double norm_const = 0.53;

    for (int i = 0; i < n; i++) {
        real rand_X = randunit(0);
        mass[i] = pow(rand_X*(1.0-alpha)/norm_const + pow(0.5, 1.0-alpha), 1.0/(1.0-alpha));
//        cerr << rand_X << "\t" << mass[i] << endl;
        mass_tot += mass[i];
    }
    *mass_cluster = mass_tot;
}

//Type 4 - Salpeter Two Stage distribution
void massDistrib_Salpeter2(int n, real mass[], real *mass_cluster)
{
    real mass_tot = 0.0;
    const double m0 = 0.08;
    const double m1 = 0.5;
    const double m2 = 150.0;
    const double alpha1 = 1.3;
    const double alpha2 = 2.3;

    //Calculate the constants
    double A = (pow(m1, 1.0-alpha1) - pow(m0, 1.0-alpha1))/(1.0-alpha1);
    double B = (pow(m2, 1.0-alpha2) - pow(m1, 1.0-alpha2))/(1.0-alpha2);
    double k1 = pow(m1, alpha1 - alpha2);
    double k2 = 1.0/(A*k1 + B);
    real p = k1*k2*A + k2*B;  //Check normality of probability function p(m). Should be unity
    cerr << "k1 =\t" << k1 << endl
         << "k2 =\t" << k2 << endl
         << "A =\t" << A << endl
         << "B =\t" << B << endl
         << "p =\t" << p << endl << endl;

    real bound = k1*k2*A;       //Boundary between the two stages
    cerr << "boundary = " << bound << endl << endl;

    //Generate the masses of the stars
    for (int i = 0; i < n; i++) {
        real X = randunit(0);
        if(X<bound)
            mass[i] = pow(X*(1.0-alpha1)/(k1*k2) + pow(m0, 1.0-alpha1), 1.0/(1.0-alpha1));
        else
            mass[i] = pow((1.0-alpha2)*(X-k1*k2*A)/k2 + pow(m1, 1.0-alpha2), 1.0/(1.0-alpha2));

        mass_tot += mass[i];
    }
    *mass_cluster = mass_tot;
}

//Type 5 - Kroupa, Tout, Gilmore
//See "Gravitational N-Body Simulations"  Sverre Aarseth 2003 p121
//  Uses a mass generating function
void massDistrib_Kroupa_etal(int n, real mass[], real mass_min, real mass_max, real *mass_cluster)
{
    real mass_tot = 0.0;
    real rand_X;
    real g1 = 0.19;     //gamma1
    real g2 = 1.55;     //gamma2
    real g3 = 0.05;     //gamma3
    real g4 = 0.6;      //gamma4

    for (int i = 0; i < n; i++) {
        rand_X = randunit(0);
        mass[i] = (mass_max - mass_min)*(0.08 + (g1*pow(rand_X,g2) + g3*pow(rand_X,g4))/pow((1.0 - rand_X),0.58));
        mass_tot += mass[i];
    }
    *mass_cluster = mass_tot;

}

//Type 6 - Kroupa2 - Five part power law
void massDistrib_Kroupa2(int n, real mass[], real mass_min, real mass_max, real *mass_cluster)
{

}

// End of mass distribution functions ------------------------------------------

/*==============================================================================
*   positions
*       Constructs the positions of each star in the cluster.
*
*       The distribution of the positions within the cluster depends upon the
*       type and shape chosen by the user.  See the top of the file.
*
*-------------------------------------------------------------------------------
*/
void positions(int n, real pos[][NDIM], int posDistribType)
{
    switch(posDistribType) {
        case 1: //Uniform distribution (Homogenous)
                posDistrib_Uniform(n, pos);
                break;
        case 2: //Central cluster
                break;
    }
}

/*------------------------------------------------------------------------------
*   posDistrib_Uniform
*       Constructs an homogeneous sphere, with unit radius
*
*       Note: we use a rejection technique, in which we choose particle positions
*       at random within a cube that encloses our sphere.  After each choice
*       of a new position, we check whether that position lies within our
*       sphere.  If it does, we accept that particle; if not, we reject it.
--------------------------------------------------------------------------------
*/
void posDistrib_Uniform(int n, real pos[][NDIM])
{
    for (int i = 0; i < n; i++){
        real rsq;
        do{
            rsq = 0;
            for (int k = 0; k < NDIM; k++){
                pos[i][k] = randinter(-1.0, 1.0);
                rsq += pos[i][k] * pos[i][k];
            }
        }while (rsq > 1);
    }
}

/*------------ End of Position distribution functions --------------------------*/

/*==============================================================================
*  velocities
*       Calculation of the velocities of the stars.
*
*-----------------------------------------------------------------------------
*/
void velocities(int n, real vel[][NDIM], int velDistribType)
{
    switch(velDistribType) {
        case 1: //All velocities are zero
                velDistrib_Zero(n, vel);
                break;
        case 2: //Random velocities
                velDistrib_Random(n, vel);
                break;
    }
}

/*------------------------------------------------------------------------------
*   velDistrib_Zero
*       All the stars have zero velocity [Default]
*
*-------------------------------------------------------------------------------
*/
void velDistrib_Zero(int n, real vel[][NDIM])
{
    for (int i = 0; i < n; i++)
        for (int k = 0; k < NDIM; k++)
            vel[i][k] = 0;
}

/*-------------------------------------------------------------------------------
*   velDistrib_Random
*--------------------------------------------------------------------------------
*/
void velDistrib_Random(int n, real vel[][NDIM])
{
    for (int i = 0; i < n; i++)
        for (int k = 0; k < NDIM; k++)
            vel[i][k] = randinter(-1.0,1.0);
}


/*===============================================================================
 *  put_snapshot  --  write a single snapshot on the output stream cout, in
 *                    the same format as described above for get_snapshot.
 *
 *  note: we use "const" here for the arguments, since they are not intended
 *        to be altered by a call to put_snapshot.
 *-----------------------------------------------------------------------------
 */

void put_snapshot(real mass[], real pos[][NDIM],
                  real vel[][NDIM], int n, real t, real *mass_cluster)
{
    cout.precision(16);

    cout << n << endl;
    cout << t << endl;
    for (int i = 0; i < n ; i++){
        cout << mass[i];
        for (int k = 0; k < NDIM; k++)
            cout << ' ' << pos[i][k];
        for (int k = 0; k < NDIM; k++)
            cout << ' ' << vel[i][k];
        cout << endl;
    }
    cout << endl;
    cout << *mass_cluster << endl;
}


/*-----------------------------------------------------------------------------
 *                                                                    \\   o
 *  end of file:  createCluster.cpp                                         /\\'  O
 *                                                                   /\     |
 *=============================================================================
 */

