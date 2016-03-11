#include<iostream>  // Writing
#include<vector>    // Gives us vectors
#include<random>    // So we can make random data for the knapsack problem
#include"tpm.h"     // The two phase solver
#include"CplexModel.h" // The CplexModel



typedef std::mt19937_64 G; // Random number generator based on the 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000
int main(int argc, char** argv)
{
    try
    {
        int     n   = 50,   // Number of items in the knapsack problem
                cap = 0.0;  // Capacity of the knapsack. Later set to sum_i w[i] / 2
        double weightSum = 0.0; // The sum of the weights
        std::vector<int> w, // Weights of the items
                         p1,// Profit of items in first objective
                         p2;// Profit of items in second objective


        /*================================================================*/
        /*          Generating random data for KP-problem                 */
        /*================================================================*/
        std::mt19937_64::result_type TheSeed = 0; // Seed for random number generator
        std::mt19937_64 generator ( TheSeed );    // Random number generator
        // Distributions
        std::uniform_int_distribution<>  Weight  = std::uniform_int_distribution<> ( 10 , 50 );
        std::uniform_int_distribution<>  Profit  = std::uniform_int_distribution<> ( 1 , 100 );

        for ( int i=0; i<n; ++i )
        {
            w.push_back ( Weight ( generator ) );
            p1.push_back ( Profit ( generator ) );
            p2.push_back ( Profit ( generator ) );
            weightSum += w[i];
        }

        // Set the capacity
        cap = int ( weightSum / 2.0 );

        /*================================================================*/
        /*          Construct the CplexModel for the KP                   */
        /*================================================================*/
        CplexModel theModel = CplexModel ( );
        theModel.buildBOKP ( n , cap , w , p1, p2 );

        /*================================================================*/
        /*          Run the two phase method                              */
        /*================================================================*/
        tpm twoPhaseMethod = tpm ( ); // Construct the object
        twoPhaseMethod.printProgress(); // Let the algorithm print the progress of phase two to the screen
        const std::string FileName = "TheOutputFile.txt"; // Choose a file name for printing results
        twoPhaseMethod.printToFile( FileName ); // Set the file name in tpm
        twoPhaseMethod.RUN ( theModel );

        /*================================================================*/
        /*          Print statistics                                      */
        /*================================================================*/
        testStatistics* TS = twoPhaseMethod.getTestStatistics ( );
        std::cout   << "Number of branching nodes        : " << TS->NumberOfBranchingNodes << "\n"
                    << "Number of phase one solutions    : " << TS->NumberOfPhaseOneSolutions << "\n"
                    << "Number of Phase two solutions    : " << TS->NumberOfPhaseTwoSolutions << "\n"
                    << "Total number of solutions        : " << TS->TotalNumberOfSolutions << "\n"
                    << "Time used in phase one (seconds) : " << TS->PhaseOneTime << "\n"
                    << "Time used in phase two (seconds) : " << TS->PhaseTwoTime << "\n"
                    << "Total time consumption (sedonds) : " << TS->TotalTime << std::endl;

        // Return 0, as all seems to be in good order
        return 0;
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception : " << e.what ( ) << std::endl;
        return 1;
    }
}
