/*  Takes a file that is output by nbody_sh1.exe and creates separate
    gnuplot files (png format) that can be combined in GIMP to create a gif animation.
    The png files are stored in the subfolder 'images'.
    Invoke by:
        createAnimationFiles [-h -d -i] < files.out

        Parameters:
            h   Help.  This lists the valid parameters and their meaning.

            t   Determines the time of the last snapshot for which an image is to be
                created.  Default is 10.0 which is the same as in nbody_sh1.

            d   Determines how many iterations are passed before an image is created.
                eg -d 5  This would produce an image for every 5th data snapshot.

            i   Flag to indicate that only the first snapshot is to be created.

            f   Basis for the filename of of all images.  eg. "-f cluster" will produce
                cluster0.png, cluster1.png, cluster2.png, etc
*/
#include <iostream>
#include <stdio.h>          //for fprintf, FILE
#include <unistd.h>         //for getopt()
#include <cstdlib>          //for atoi()
#include <time.h>           //for clock(), clock_t, CLOCKS_PER_SEC
#include <cstring>

using namespace std;

const int NDIM = 3; //Number of spatial dimensions

bool read_options(int argc, char *argv[], int &display_freq, double &t_max, bool &initial_flag, char * filename);
void wait(int seconds);

int main(int argc, char *argv[])
{
    int n;                      //Number of stars
    int dummy;                  //Used to store subsequent values of n
    int display_freq = 1;       //Number of iterations between each image created.
    double t;                   //Iteration time
    double t_max = 10.0;        //Time of the last image created
    bool initial_flag = false;  //Flag that detemines if only the initial snapshot is required
    char filename[200];

    if (! read_options(argc, argv, display_freq, t_max, initial_flag, filename))
        return 1;                // halt criterion detected by read_options()

    cin >> n;
    cin >> t;

    double mass;    //Not used by gnuplot so continually overwritten
    double (*pos)[NDIM] = new double[n][NDIM];    //Store the positions for one iteration
    double vel[NDIM];  //Not used by gnuplot so continually overwritten

    int iter = 0;       //Loop through the sequence of results in the file
    do{
        if(iter!=0){
            cin >> dummy;
            cin >> t;
        }
        for (int i = 0; i < n ; i++){
            cin >> mass;                    // mass of particle i
            for (int k = 0; k < NDIM; k++)
                cin >> pos[i][k];              // position of particle i
            for (int k = 0; k < NDIM; k++)
                cin >> vel[k];                 // velocity of particle i. Note these are overwritten.
        }

        //Set standard configuration for the plot
        double xMax = 3.0;
        double yMax = xMax;
        double zMax = xMax;

//        char imagename[] = "clusterCold26";

        // Plot the position of all stars for each time period.
        if(iter%display_freq==0){  //Only plot every 'display_freq' set of data

            //#define GNUPLOT_NAME "pgnuplot.exe -persist"  //This version for window
//            #define GNUPLOT_NAME "pgnuplot.exe"         //This version for file creation in Win
            #define GNUPLOT_NAME "gnuplot"              //This version for file creation in Linux
            FILE *pipe = popen(GNUPLOT_NAME, "w");

            if (pipe != NULL)
            {
                fprintf(pipe, "set term png size 800,600\n");     //set the terminal to png image
                if (iter<10) {
                    fprintf(pipe, "set out './images/%s_0%d.png'\n", filename, iter);  //Save the image file
                    } else {
                    fprintf(pipe, "set out './images/%s_%d.png'\n", filename, iter);
                    }
                fprintf(pipe, "set key off\n");
                fprintf(pipe, "set xrange [%f:%f]\n", -xMax, xMax);
                fprintf(pipe, "set yrange [%f:%f]\n", -yMax, yMax);
                fprintf(pipe, "set zrange [%f:%f]\n", -zMax, zMax);
                fprintf(pipe, "set view 360,360\n");              //Set to 360,360 for 2D, 57,36 for normal view
                for(int i=0; i<n; i++){                         //Create variables and assign values
                    fprintf(pipe, "ax%d=%f\n", i, pos[i][0]);  //for input into the plot routine.
                    fprintf(pipe, "ay%d=%f\n", i, pos[i][1]);
                    fprintf(pipe, "az%d=%f\n", i, pos[i][2]);
                }
                fprintf(pipe, "splot ");
                for(int i=0; i<n; i++) {
                    fprintf(pipe, "'+' using (ax%d):(ay%d):(az%d) ps 2 pt 3 title ''", i,i,i); // plot data
                    if(i!=(n-1)){
                        fprintf(pipe, ",");
                    }
                }
                fprintf(pipe, "\n");                    //End the plot command



//                fprintf(pipe, "set output \n");          //Set output to nothing - required for windows

                fflush(pipe);                           // flush and
                pclose(pipe);                          // close the pipe
            }//end if pipe
            else
                printf("Could not open the pipe\n");
        }//end if iter
        cout << ".";    //Show progress
        wait(2);        //Wait while gnuplot catches up
        iter++;
    } while((t<t_max)&&(!initial_flag));   //Stop just before the end of the file which is t = t_max.
                                           //Stop if only the first snapshot is requested.
                                           //End of 'iter' loop
    cout << endl;
    delete[] pos;

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

bool read_options(int argc, char *argv[], int &display_freq, double &t_max, bool &initial_flag, char *filename)
{
    int c;
    while ((c = getopt(argc, argv, "ht:f:d:i")) != -1)
        switch(c){
            case 'd': display_freq = atoi(optarg);  //Number of iterations between each image
                      break;
            case 't': t_max = atof(optarg);
                      break;
            case 'i': initial_flag = true;
                      break;
            case 'f': strcpy(filename, optarg);
                      cout << filename << endl;
                      break;  //Base file name for images (future usage??)
            case 'h':
            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-t Time of last snapshot to make into an image (Default=10.0)]\n"
                           << " [-d No. of iterations between each image (Default=1)]\n"
                           << " [-f Base filename for images (Images named 'baseFileName###')]\n"
                           << endl;
                      return false;      // execution stops after help or error
            }

    return true;                         // ready to continue program execution
}

//Routine to halt the execution for a specified time.
//Used to give other threads/processes time to catch up.
void wait(int seconds){
      clock_t endwait;
      endwait = clock () + seconds*CLOCKS_PER_SEC;
      while (clock() < endwait) {}
}

