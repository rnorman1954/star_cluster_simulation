/*==========================================================================
*    convert.cpp
*        Takes a file in the 'createCluster' or 'nbody' format, strips away the
*        the mass, time and velocities to leave a file (stars.dat) of positions only.
*        This can be easily opened with gnuplot to display the positions
*        of the stars.
*        This is used for diagnostic purposes.
*
*---------------------------------------------------------------------------
*    Usage:
*        convert < file.in
*
*---------------------------------------------------------------------------
*   Version: 1  Mar 2019    Robert Norman
*
*-------------------------------------------------------------------------*/

#include <iostream>
#include <stdio.h>

using namespace std;

int main()
{
    int n; //Number of stars
    int startTime;

    cin >> n;
    cin >> startTime;

    double *mass = new double[n];
    double (*pos)[3] = new double[n][3];
    double vel[3];

    for (int i = 0; i < n ; i++){
        cin >> mass[i];                    // mass of particle i
        for (int k = 0; k < 3; k++)
            cin >> pos[i][k];              // position of particle i
        for (int k = 0; k < 3; k++)
            cin >> vel[k];                 // velocity of particle i. Note these are overwritten.
    }


    // Save data in file, ready for reading by an external plotting program.
    FILE *output;
    output= fopen("stars.dat","w");

    for (int i=0; i<n; i++)  {
        for (int k = 0; k < 3; k++)
            fprintf(output, "%f\t", pos[i][k]);        // position
        fprintf(output,"\n");
    }

    fclose(output);

    cout << "Data stored in stars.dat" << endl;

    delete[] mass;
    delete[] pos;

    return 0;
}
