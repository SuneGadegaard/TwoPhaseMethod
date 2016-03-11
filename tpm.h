/** \mainpage An introduction to two phase method
* \author Sune Lauth Gadegaard
* \version 1.0.0
*
* \section License
*
* Copyright 2015, Sune Lauth Gadegaard.
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*
\section Description
*
* This program provides a two--phase method for combinatorial optimization problems. A long with the implementations, a main file giving an example of how to use the program
* is provided. The two phase method relies solely on CPLEX as a solver for the subproblems, so if you do not have CLPEX installed, this program will for sure not run!
*
* The class tpm implements the two phase method, while the class CplexModel implements a class used to hold the cplex objects. The classes NDS and solutions implements a non--dominated
* set and a solution, respectively.
*
* \section Compiling
* The codes were compiled using the GNU GCC compiler on a Linux Ubuntu 14.04 machine.
* The following flags were used: -Wall -O3 -std=c++11 -DIL_STD.
* The Code::blocks IDE was used as well. At \href{http://www-01.ibm.com/support/docview.wss?uid=swg21449771}{this page} you can find a guide to how one can configure Code::Block with
* CPLEX on a linux machine. The programs has not been tested on any other operating systems or with any other IDEs.
*
* \latexonly
* \section{Change log for tpm.h, NDS.h, CplexModel.h, and solution.h}
* \begin{center}
*     \begin{tabularx}{\textwidth}{llr X}\toprule
*        FILE:          &   \multicolumn{3}{l}{tpm.h, NDS.h, CplexModel.h, and solution.h}\\
*        Version:       &   \multicolumn{3}{l}{1.0.0}\\
*        \multicolumn{4}{l}{- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -}\\
*        CHANGE LOG:    &   DATE    &   VER.-NO.    &   CHANGES MADE\\ \midrule
*                       &   2016--03--10    &   1.0.0       & First implementation\\ \bottomrule
*     \end{tabularx}
* \end{center}
* \endlatexonly
*/

#ifndef TPM_H_INCLUDED
#define TPM_H_INCLUDED

//! C++ includes
#include<iostream>
#include<set>
#include<vector>
#include<ilcplex/ilocplex.h>
#include<stdexcept>
#include<chrono>

//! My own C++ includes
#include"NDS.h" //! Implementation of a non domminated set
#include"solution.h" //! Implementation of a solution class to hold a solution
#include"CplexModel.h" //! Implememntation of the class holding the cplex model

typedef IloArray<IloNumVarArray>    IloVarMatrix;
using namespace std::chrono;
typedef std::chrono::high_resolution_clock CPUclock;

struct testStatistics{
    unsigned long NumberOfBranchingNodes;       //!< Total number of branching nodes
    unsigned long NumberOfPhaseOneSolutions;    //!< Number of solutions generated in phase one
    unsigned long NumberOfPhaseTwoSolutions;    //!< Number of solutions generated in phase two
    unsigned long TotalNumberOfSolutions;       //!< Total number of solutions (in objective space) on the non--dominated frontier
    double PhaseOneTime;                        //!< Time in seconds used in phase one
    double PhaseTwoTime;                        //!< Time in seconds used in phase two
    double TotalTime;                           //!< Total time used on the entire algorithm
}; //!< Struct used to gather test statistics


class tpm{
    private:

        struct BOUND{
            double UB;
            double LB;
        }; //!< Struct used to specify bounds on an objective function

        struct BOUNDS{
            BOUND f1;
            BOUND f2;
        }; //!< Struct used to specify bounds on both objective functions


        CPUclock::time_point StartTime;

        /**
         * @name Parameters and flags
         * This section contains a list of parameters and flags used internally in the SSCFLPsolver class.
         */
            double myZero;      //!< Below this value, and you are effectually equal to zero
            double myOne;       //!< Above this value and you are effectually equal to one
            double myTol;       //!< Tolerance for equality.
            //double TimeLeft;    //!< Time left to use by cplex!
            int problemtype;    //!< Indicating the problem type
            bool PrintProgress; //!< Prints the progress in the Ranking based two phase method
            double totalTime;   //!< Variable holding the time limit for the whole solve
            bool PrintToFile;     //!< True if solutions should be printed to file. Default is false
            std::string FileName; //!< Name of the file, which should printed to
            bool DoRanking;     //!< If true, the rannking based two phase method is used. Default is false, meaning the perpendicular search method is used in phase two as default.
            bool TakeFront;     //!< If true, the stack of problems generated in the PSM method is taken on a FIFO principle. Otherwise, FILO principle
        ///@}

        /**
         * @name Test statistics
         * This section contains all test statistics reported by the algorithm. All variables should have self explanatory names.
         */
        ///@{
            testStatistics* theStatistics;
        ///@}

        /**
         * @name Data structures
         * This section contains all the data structures used in the two phase method
         */
         NDS NonDomSet;
        ///@{

        /*! \brief Runs the Non-Inferior Set Estimation algorithm as a first phase.
         *
         * This function runs the NISE algorithm as a first phase used to generate all extreme supported non dominated outcomes
         * It uses a NISE algorithmic framework and calls the function
         * \param theModel reference to a CplexModel object. The CplexModel object contains a bi-objective combinatorial optimization problem
         */
        void RunPhaseOne ( CplexModel &theModel );

        /*! \brief Runs a perpendicular search method as a phase two.
         * This function runs a perpendicular search method algorithm for each triangle created by the first phase. This algorithm is default.
         * \param theModel reference to a CplexModel object. The CplexModel object contains a bi-objective combinatorial optimization problem
         */
        void RunPhaseTwo ( CplexModel &theModel );

        /*! \brief Runs a ranking algorithm as phase two
         * Generates all the solutions which are not found in phase one by ranking the solutions in the
         */
        void RunPhaseTwoRanking ( CplexModel &theModel );
    public:
        /*!
         * Default constructor setting default values for parameters and flags.
         */
        tpm ( );

        /*!
         * Destructor releasing allocated memory and cleaning up
         */
        ~tpm ( );

        /*! \brief Works as the public API for the two--phase method.
         *
         * Works as the public API for the tpm class. It runs a two phase method in order to generate all non--dominated outcomes of
         * bi--objective combinatorial optimization problem specified in the CplexModel object theModel parsed as an argument.
         * \param theModel reference to a CplexModel object. Contains the BOCO problem which should be solved by the two phase algorithm
         */
        int RUN( CplexModel &theModel );

        /*! \brief Sets a time limit for the whole algorithm.
         *
         * Sets a time limit for the entire algorithm. It does, however, not terminate the cplex search when time is up, but
         * terminates the algorithm if, after cplex finishes an optimization, the time limit is reached. Therefore, the procedure
         * will most likely use more time than the time limit.
         * \param timeLimit double. The time limit in seconds.
         */
        void setTimeLimit ( double timeLimit ){ totalTime = timeLimit; }

        /*!
         * Tells the algorithm to print progress of phase two to the screen
         */
        void printProgress ( ) { PrintProgress = true; }

        /*! \brief Sets a file for printing solution info
         * Sets the internal falg PrintToFile = true and stores the fileName so that results are printet to this file.
         * \param fileName const reference to a string. Contains the file name to where results are printed
         */
        void printToFile ( const std::string & fileName );

        /*! \brief Sets the phase two algorithm to a ranking based algorithm
         *
         * Sets the phase two algorithm to a ranking based algorithm which ranks all solutions between two supported non--dominated solutions until
         * a solution value exceeds the value of the worst local Nadir point in the corresponding triangle. Can potentially work better than default
         * PSM method on problems with a totally unimodular constraint matrix.
         */
        void doRanking ( ) { DoRanking = true; }

        /*! \brief Specifies the order in which subproblems are processed in the PSM method
         *
         * Specifies to do a depth first search in the PSM method. If doRanking () has been called prior to calling setFILO () it has no effect.
         * Default is a best bound search.
         */
        void setDepthFirst ( ) { TakeFront = false; }

        /*! \brief Returns the test statistics
         *
         * Returns a point to a testStatistics struct. The struct contains test statistics obtained throughout the algorithm.
         */
        testStatistics* getTestStatistics ( ) { return theStatistics; }

};

#endif // TPM include guard ends here
