#include <iostream>
#include <math.h>

/* a code for solving the radiative transfer equation in a plane-parallel
   atmosphere using variable Eddington factors (see article "Difference
   Equations and Linearization Methods" by Lawrence Auer
   in "Methods in Radiative Transfer" ed. Wolfgang Kalkofen (1984) */

#define N_DEPTH_PTS 100
#define N_MU_PTS 20

int main()
{
  // thermalization parameter for source function in isotropic, monochromatic
  // scattering (Milne-Eddington problem)
  const double therm_parm = 1.0e-4;
  // optical depth grid
  double *tau = new double [ N_DEPTH_PTS ];
  // direction cosine grid
  double *mu = new double [ N_MU_PTS ];
  // Eddington factor f_K = K / J
  double **f_K = new double *[ N_DEPTH_PTS ];
  // Eddington factor f_H = \int_0^1 j(\mu) \mu d\mu / J
  double **f_H = new double *[ N_DEPTH_PTS ];
  // source function
  double **S = new double *[ N_DEPTH_PTS ];
  // mean intensity
  double **J = new double *[ N_DEPTH_PTS ];

  for ( int i = 0; i < N_DEPTH_PTS; i++ )
  {
    f_K [ i ] = new double [ N_MU_PTS ];
    f_H [ i ] = new double [ N_MU_PTS ];
    S   [ i ] = new double [ N_MU_PTS ];
    J   [ i ] = new double [ N_MU_PTS ];
  }

  // optical depth points spaced logarithmically
  std::cout << "OPTICAL DEPTH GRID:" << std::endl;
  tau [ 0 ] = 1.0e-8;
  std::cout << "tau [ 0 ] = " << tau [ 0 ] << std::endl;
  for ( int i = 1; i < N_DEPTH_PTS; i++ )
  {
    tau [ i ] = tau [ i-1 ] * 1.27;
    std::cout << "tau [ " << i << " ] = " << tau [ i ] << std::endl;
  }

  // Eddington approximation becomes initial guess for Eddington factors
  for ( int i = 0; i < N_DEPTH_PTS; i++ )
  {
    for ( int j = 0; j < N_MU_PTS; j++ )
    {
      f_K [ i ][ j ] = 1.0 / 3.0;
      f_H [ i ][ j ] = 1.0 / sqrt( 3.0 );
    }
  }

  delete [] J;
  delete [] S;
  delete [] f_H;
  delete [] f_K;
  delete [] mu;
  delete [] tau;

  return 0;
}
