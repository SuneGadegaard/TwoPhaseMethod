#include"tpm.h"

using namespace std::chrono;

/********************************************************************************************/
inline
double dist( const solution &sol1, const solution &sol2 )
{
    double x = sol1.getFirst() - sol2.getFirst();
    double y = sol1.getSecond() - sol2.getSecond();
    double dist = pow ( x , 2 ) + pow ( y , 2 );
    dist = sqrt ( dist );
    return dist;
}

/********************************************************************************************/
tpm::tpm ( ):
    myZero ( 0.0001 ),
    myOne ( 0.9999 ),
    myTol ( 0.001 ),
    PrintProgress( false ),
    totalTime ( 1e+75 ),
    PrintToFile ( false ),
    DoRanking ( false ),
    TakeFront ( true )
{
    theStatistics = new testStatistics;
}

/********************************************************************************************/
tpm::~tpm ( ){
    delete theStatistics;
}

/********************************************************************************************/
int tpm::RUN( CplexModel &theModel )
{
    try{
        /*================================================*/
        /*      Initialize test statistics                */
        /*================================================*/
        theStatistics->TotalTime =
        theStatistics->PhaseOneTime =
        theStatistics->PhaseTwoTime =
        theStatistics->NumberOfBranchingNodes =
        theStatistics->TotalNumberOfSolutions =
        theStatistics->NumberOfPhaseOneSolutions =
        theStatistics->NumberOfPhaseTwoSolutions = 0;

        // Tell cplex not to reduce the problem!
        theModel.cplex.setParam( IloCplex::Reduce , 0);

        /*================================================*/
        /*      Phase one starts here                     */
        /*================================================*/
        auto Start_1 = CPUclock::now ( );
        StartTime = Start_1;
        RunPhaseOne ( theModel );
        auto End_1 = CPUclock::now ( );
        // Gather statistics
        theStatistics->PhaseOneTime = duration_cast< duration < double > > (End_1 -Start_1 ).count ( );
        // Print out time consumption to the screen
        std::cout << "Phase one time : " << theStatistics->PhaseOneTime << " seconds\n";


        /*================================================*/
        /*      Phase two starts here                     */
        /*================================================*/
        auto Start_2 = CPUclock::now ( );
        if ( DoRanking )
        {
            // Ranking based two phase method
            RunPhaseTwoRanking( theModel );
        }
        else
        {
            // Perpendicular search method based TPM
            RunPhaseTwo ( theModel );
        }
        auto End_2 = CPUclock::now ( );

        // Print out time consumption to the screen
        theStatistics->PhaseTwoTime    = duration_cast< duration < double > > (End_2 -Start_2 ).count ( );
        theStatistics->TotalTime       = theStatistics->PhaseOneTime + theStatistics->PhaseTwoTime;
        std::cout << "Phase two time : " << theStatistics->PhaseTwoTime << " seconds\n";
        std::cout << "Total time     : " << theStatistics->TotalTime    << " seconds\n";



        /*================================================*/
        /*      If PrintToFile = True, we print           */
        /*================================================*/
        if ( PrintToFile ){
            std::vector<double> VarVals;
            std::ofstream OutputFile;
            OutputFile.open ( FileName , std::ofstream::out | std::ofstream::app);

            if ( !OutputFile )
            {
                throw 101;
            }
            // First we print the test statistics
            OutputFile  << "NumVar \t TotalTime \t PhaseOne time \t PhaseTwo time \t BranchingNodes \t Number of Solutions \t PhaseOneSol \t PhaseTwoSol \n";
            OutputFile  << theModel.AllVars.getSize ( ) << "\t"
                        << theStatistics->TotalTime << "\t"
                        << theStatistics->PhaseOneTime << "\t"
                        << theStatistics->PhaseTwoTime << "\t"
                        << theStatistics->NumberOfBranchingNodes << "\t"
                        << theStatistics->TotalNumberOfSolutions << "\t"
                        << theStatistics->NumberOfPhaseOneSolutions << "\t"
                        << theStatistics->NumberOfPhaseTwoSolutions << "\n";

            // Then we print the solution
            OutputFile  << "Fontier and corresponding solution\n";
            for (auto it = NonDomSet.NDs.begin ( ); it != NonDomSet.NDs.end ( ); ++it )
            {
                OutputFile  << it->getFirst ( ) << "\t"
                            << it->getSecond ( );

                it->getVarValues ( VarVals );
                for ( auto varIt = VarVals.begin ( ); varIt != VarVals.end ( ); ++varIt )
                {
                    OutputFile << "\t" << *varIt;
                }
                OutputFile << "\n";
            }
            OutputFile.close ( );

        }
        return 0;
    }
    catch ( int i )
    {
        if ( 101 == 1 )
        {
            std::cerr << "Could not open the specified file for printing!\n";
        }
        return 101;
    }
    catch(std::exception &e)
    {
        std::cerr << "Exception in the RUN in tpm class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
    catch(IloException &ie){
        std::cerr << "IloException in RUN in tpm class : " << ie.getMessage ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/********************************************************************************************/
void tpm::RunPhaseOne( CplexModel &theModel){
    try
    {
        bool OnlyOneNonDomSol = false;  // If phase one found only one solution, we are done, as that is the only non--dominated solution. Assume false.
        int maxiterations = 1000;       // Iterations performed by the NISE algorithm. Method is still exact, as phase two finds any missing solutions
        double  lambda1,                // Weight for the first objective
                lambda2;                // Weight for the second objective
        std::pair<double,double> p;     // Pair used to store the outcome vector of a solution
        solution solUL, solLR;          // Solutions used to hold the solutions marking the current search dirrection
        CPUclock::time_point endTime;   // Time point used to measure time from start until now
        std::chrono::duration<double> TimeUntilNow; // Difference between time points
        std::vector<double> VarVals;    // Vector used to hold the current solution

        auto PlusIt = NonDomSet.SupNDs.begin ( ); // Iterator used to point to solUL
        auto MinusIt= NonDomSet.SupNDs.begin ( ); // Iterator used to point to solLR

        // Set cplex parameters
        theModel.cplex.setParam ( IloCplex::EpAGap , 0.0 ); // Absolute gap is zero as we can assume nothing about the integrality of coeficients
        theModel.cplex.setParam ( IloCplex::EpGap  , 0.0 ); // Relative gap is zero as we can assume nothing about the magnitude of the solutions
        theModel.cplex.setOut( theModel.env.getNullStream ( ) );    // Tell cplex not to print the log
        theModel.cplex.setWarning( theModel.env.getNullStream ( ) );// Tell cplex not to print warnings

        /*==========================================================*/
        /*      Start by finding the upper left point               */
        /*==========================================================*/
        theModel.OBJ.setLinearCoef( theModel.f1 , 1.0 ); // Full weight on objective one
        theModel.OBJ.setLinearCoef( theModel.f2 , 0.0 ); // No weight on objective two

        if( theModel.cplex.solve ( ) ) // Solve the problem. If we cant, we need to terminate, and an error is thrown
        {
            theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( ); // Collect the number of branching nodes

            /*=============================================================*/
            /*      Test if time limit is reached                          */
            /*          If so, go to the end!                              */
            /*=============================================================*/
            endTime = CPUclock::now ( );
            TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
            if ( TimeUntilNow.count ( ) > 3600.0 ) //
            {
                std::cout << "=========== Timeout ===========\n";
                goto END;
            }
            // We now have the objective function value of the first objective
            p.first = theModel.cplex.getValue( theModel.f1 );
            // Now change objective function coefficients and bounds on f1
            theModel.f1.setUB( p.first + myTol ); // Add a little to ensure nothing goes wrong
            theModel.OBJ.setLinearCoef( theModel.f1 , 0 );
            theModel.OBJ.setLinearCoef( theModel.f2 , 1 );

            // Resolve to get the lexicographic minimizer
            theModel.cplex.setParam(IloCplex::Param::TimeLimit , 36000 );
            if ( theModel.cplex.solve ( ) )
            {
                theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( );
                endTime = CPUclock::now ( );
                TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
                if ( TimeUntilNow.count ( ) > 3600.0 )
                {
                    std::cout << "=========== Timeout ===========\n";
                    goto END;
                }
                p.second = theModel.cplex.getValue( theModel.f2 );
                // Retrieve the solution corresponding to this outcome
                for ( int s = 0; s < theModel.AllVars.getSize (); ++s )
                {
                    VarVals.push_back( theModel.cplex.getValue( theModel.AllVars[s] ) );
                }
                // Create a new solution, and add it to the frontier
                solution sol = solution ( true , p , VarVals );
                NonDomSet.SupNDs.push_back( sol );
                // Clear VarVals for reuse
                VarVals.clear ( );
                solUL = sol;
                // Remember to set the upper bound of f1 back to IloInfinity
                theModel.f1.setUB ( IloInfinity );
            }else
            {
                throw std::runtime_error ( "Could not solve the model second time in order to find z^ul" );
            }
        }else
        {
            throw std::runtime_error ( "Could not solve the model first time in order to find z^ul" );
        }

        /*==========================================================*/
        /*      Continue by finding the lower right point           */
        /*==========================================================*/
        theModel.OBJ.setLinearCoef( theModel.f1 , 0 );  // No weight on first objective
        theModel.OBJ.setLinearCoef( theModel.f2 , 1 );  // Full weight on second objective

        if ( theModel.cplex.solve ( ) ) // Solve the problem. If we cant, we need to terminate, and an error is thrown
        {
            theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( ); // Collect branching nodes

            /*=============================================================*/
            /*      Test if time limit is reached                          */
            /*          If so, go to the end!                              */
            /*=============================================================*/
            endTime = CPUclock::now ( );
            TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
            if ( TimeUntilNow.count ( ) > 3600.0 )
            {
                std::cout << "=========== Timeout ===========\n";
                goto END;
            }
            // We now have the objective function value of the second objective
            p.second = theModel.cplex.getValue( theModel.f2 );
            // Now change objective function coefficients and bounds on f2
            theModel.f2.setUB( p.second+0.1 ); // Add a little to ensure nothing goes wrong
            theModel.OBJ.setLinearCoef( theModel.f1 , 1 );
            theModel.OBJ.setLinearCoef( theModel.f2 , 0 );

            // Resolve to get the lexicographic minimizer
            if ( theModel.cplex.solve ( ) )
            {
                theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( ); // Collext branching nodes
                /*=============================================================*/
                /*      Test if time limit is reached                          */
                /*          If so, go to the end!                              */
                /*=============================================================*/
                endTime = CPUclock::now ( );
                TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
                if ( TimeUntilNow.count ( ) > 3600.0 )
                {
                    std::cout << "=========== Timeout ===========\n";
                    goto END;
                }
                for ( int s = 0; s < theModel.AllVars.getSize (); ++s )
                {
                    VarVals.push_back( theModel.cplex.getValue( theModel.AllVars[s] ) );
                }
                p.first = theModel.cplex.getValue( theModel.f1 );
                solution sol = solution ( true , p , VarVals );
                VarVals.clear ( ); // Clear VarVals for reuse
                // Check if we have produced the same solution again!
                if ( sol.getFirst() != NonDomSet.SupNDs.begin()->getFirst() ){
                    NonDomSet.SupNDs.push_back( sol );
                }
                else OnlyOneNonDomSol = true;
                solLR = sol;
                // Remember to set the ubber buund of f2 back to iloinfinity
                theModel.f2.setUB ( IloInfinity );
            }else
            {
                throw std::runtime_error ( "Could not solve the model in order to find z^ul" );
            }
        }else
        {
            throw std::runtime_error ( "Could not solve the model in order to find z^lr" );
        }

        /*==========================================================*/
        /*      Continue with the main algorithm,                   */
        /*          Iterate as in a NISE algorithm                  */
        /*==========================================================*/
        PlusIt = NonDomSet.SupNDs.begin ( );
        MinusIt= std::next ( NonDomSet.SupNDs.begin ( ) );

        if (!OnlyOneNonDomSol){ // If we have only one solution, we must terminate her, as there is only one non--dominated outcome
            for ( int iteration = 1; (iteration<maxiterations) && dist( solLR , *PlusIt )>=0.1; ++iteration )
            { // Loop as long as PlusIt!=lower right corner and as long as we are below maxiterations
                //Update the scalars
                lambda1 = ( PlusIt->getSecond ( ) - MinusIt->getSecond ( ) );   // Calculate the weight of the first objective
                lambda2 = ( MinusIt->getFirst ( ) - PlusIt->getFirst ( ) );     // Calculate the weight of the second objective

                // Update the objective coefficients
                theModel.OBJ.setLinearCoef( theModel.f1 , lambda1 );
                theModel.OBJ.setLinearCoef( theModel.f2 , lambda2 );

                // If we could not solve, throw a runtime error
                if ( !( theModel.cplex.solve ( ) ) ) throw std::runtime_error ( "Could not solve the subproblem in phase on of the two phase method " );
                theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( ); // Collect branching nodes

                /*=============================================================*/
                /*      Test if time limit is reached                          */
                /*          If so, go to the end!                              */
                /*=============================================================*/
                endTime = CPUclock::now ( );
                TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
                if ( TimeUntilNow.count ( ) > 3600.0 )
                {
                    std::cout << "=========== Timeout ===========\n";
                    goto END;
                }

                // Check if we have found a new solution
                if ( theModel.cplex.getObjValue ( ) <= ( lambda1*PlusIt->getFirst() + lambda2*PlusIt->getSecond() -myTol ) )
                {
                    p.first = theModel.cplex.getValue ( theModel.f1 );
                    p.second= theModel.cplex.getValue ( theModel.f2 );
                    for ( int s = 0; s < theModel.AllVars.getSize ( ); ++s )
                    {
                        VarVals.push_back( theModel.cplex.getValue( theModel.AllVars[s] ) );
                    }
                    solution sol = solution ( true , p , VarVals );
                    VarVals.clear ( );
                    NonDomSet.SupNDs.insert ( MinusIt , sol );
                }else{
                    //The solution is not new and we go to the next one
                    PlusIt = MinusIt;
                }
                MinusIt = std::next ( PlusIt );
            }
        }
        END:
        // Copy the solutions found in phase one into the set of non-dominated solutions
        NonDomSet.copySupToNonDom ( );

        theStatistics->NumberOfPhaseOneSolutions = NonDomSet.NDs.size ( );
    }catch(std::exception &e){
        std::cerr << "Exception in the RunPhaseOne in tpm class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }catch(IloException &ie){
        std::cerr << "IloException in RunPhaseOne in tpm class : " << ie.getMessage ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/********************************************************************************************/
void tpm::RunPhaseTwo ( CplexModel &theModel )
{
    try
    {
        CPUclock::time_point endTime;
        std::chrono::duration<double> TimeUntilNow;
        IloBoolVarArray branchVar(theModel.env);
        bool OnlyOneNonDomSol = (NonDomSet.NDs.size ( ) == 1); // Check if only one solution was found in phase one
        int triangle        = 0,   // Variable used to count the triangles
            NumOfTriangles  = NonDomSet.NDs.size ( ) - 1,  // Variable holding the number of triangles we should process
            NumOfVars       = theModel.AllVars.getSize ( );
        double  lambda1=0.0, // Weight of first objective
                lambda2=0.0; // Weight of second objective
        std::pair<double,double> p; // Pair used to store outcome vector of a solution
        BOUNDS CurrentBounds;   // BOUNDS variable to hold the current bounds on the objective functions
        std::list<BOUNDS> bounds;   // List of BOUNDS
        std::list<solution>::iterator nextSol; // Iterator to the "next" solution on the list of non--dominated solutions
        std::vector<double> VarValues;
        // Tell cplex not to print to the console
        theModel.cplex.setOut( theModel.env.getNullStream ( ) );

        // loop over all supported non-dominated points
        for ( auto it = NonDomSet.SupNDs.begin ( ); !OnlyOneNonDomSol && std::next ( it )!=NonDomSet.SupNDs.end ( ); ++it )
        {
            if ( PrintProgress )
            {
                std::cout << "Processing triangle " << ++triangle << " of " << NumOfTriangles << "\n";
            }
            // Retreive the solution following it
            nextSol = std::next ( it );

            // Set the bounds on the objectives based on the current triangle
            theModel.f2.setUB ( it->getSecond ( ) );
            theModel.f1.setUB ( nextSol->getFirst ( ) );

            // Calculate the slope of the search direction
            lambda1 = it->getSecond ( ) - nextSol->getSecond ( );
            lambda2 = nextSol->getFirst ( ) - it->getFirst ( );

            // Set the objective function coefficients according to it and nextSol
            theModel.OBJ.setLinearCoef( theModel.f1 , lambda1 );
            theModel.OBJ.setLinearCoef( theModel.f2 , lambda2 );

            // Initialize the stack of subproblems
            BOUNDS FirstBounds;
            FirstBounds.f1.UB = nextSol->getFirst ( ) ;
            FirstBounds.f1.LB = it->getFirst ( ) ;
            FirstBounds.f2.UB = it->getSecond ( ) ;
            FirstBounds.f2.LB = nextSol->getSecond ( ) ;
            bounds.push_back ( FirstBounds );

            while ( !bounds.empty ( ) )
            {
                if ( TakeFront )
                {   // Get the first element on the
                    CurrentBounds = bounds.front ( ); // Get the back of the vector
                    bounds.erase ( bounds.begin( ) ); // Pop the back, so we do not need to inspect it again
                }
                else
                {   // Get the last element on the
                    CurrentBounds = bounds.back ( );
                    bounds.pop_back ( );
                }


                theModel.f1.setBounds ( CurrentBounds.f1.LB , CurrentBounds.f1.UB );
                theModel.f2.setBounds ( CurrentBounds.f2.LB , CurrentBounds.f2.UB );

                theModel.cplex.setParam(IloCplex::Param::TimeLimit , 3600 );
                endTime = CPUclock::now ( );
                TimeUntilNow = std::chrono::duration_cast<duration<double>>(endTime - StartTime);
                if ( TimeUntilNow.count ( ) > 3600.0 )
                {
                    std::cout << "=========== Timeout ===========\n";
                    goto END;
                }
                if ( theModel.cplex.solve ( ) )
                {
                    // Update the Time left
                    theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes();

                    // If the current model has a solution, get it!
                    p.first = theModel.cplex.getValue( theModel.f1 );
                    p.second = theModel.cplex.getValue( theModel.f2 );
                    // Create a new solution, and insert it into the non-dominated set
                    for ( int  i = 0; i < NumOfVars; ++i )
                    {
                        VarValues.push_back( theModel.cplex.getValue ( theModel.AllVars[i] ) );
                    }
                    solution sol = solution ( false , p , VarValues );
                    NonDomSet.updateNDS( sol );
                    VarValues.clear ( );
                    // Create two new subproblems:
                    {  // First subproblem, to the left of the current outcome vector
                        if ( p.first -1 < CurrentBounds.f1.LB || p.second +1 > CurrentBounds.f2.UB )
                        {}// The left subproblem is infeasible, and should not be added!
                        else
                        { // The left subproblem might be feasible, we create it, and se what happens
                           BOUNDS leftBounds;
                           leftBounds.f1.LB = CurrentBounds.f1.LB;
                           leftBounds.f1.UB = p.first - 1.0;
                           leftBounds.f2.LB = p.second + 1.0;
                           leftBounds.f2.UB = CurrentBounds.f2.UB;
                           bounds.push_back( leftBounds );
                        }

                        // Now create the subproblem to the right
                        if ( ( p.second -1.0 < CurrentBounds.f2.LB ) || ( p.first + 1.0 >= CurrentBounds.f1.UB ) )
                        {} // The right subproblem is infeasible, and we can descard it
                        else
                        { // The right subproblem might be feasible, we create it, and se what happens.
                            BOUNDS rightBounds;
                            rightBounds.f1.LB = p.first + 1.0;
                            rightBounds.f1.UB = CurrentBounds.f1.UB;
                            rightBounds.f2.LB = CurrentBounds.f2.LB;
                            rightBounds.f2.UB = p.second - 1.0;
                            bounds.push_back( rightBounds );
                        }
                    }
                }
            }
        }
        END:
        std::cout << "Number of supported efficient solutions     : " << NonDomSet.SupNDs.size ( ) << std::endl;
        std::cout << "Number of non supported efficient solutions : " << (NonDomSet.NDs.size ( ) - NonDomSet.SupNDs.size ( ) ) << std::endl;

        theStatistics->TotalNumberOfSolutions = NonDomSet.NDs.size ( );
        theStatistics->NumberOfPhaseTwoSolutions = theStatistics->TotalNumberOfSolutions - theStatistics->NumberOfPhaseOneSolutions;


    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in RunTwoPhase in the tpm class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in RunTwoPhase in the tpm class : " << ie.getMessage ( )  << std::endl;
        exit ( EXIT_FAILURE );
    }
}

/********************************************************************************************/
void tpm::RunPhaseTwoRanking( CplexModel &theModel )
{
    try
    {
        IloInt NumOfVars = theModel.AllVars.getSize ( ); // Variable used to store the number of variables in the current model
        IloBoolVarArray branchVar(theModel.env);
        bool OnlyOneNonDomSol = (NonDomSet.NDs.size ( ) == 1);
        int     triangle = 0,
                NumOfTriangles = NonDomSet.NDs.size ( ) - 1;
        double  f1_bound=0.0,// Upper bound on objective 1
                f2_bound=0.0,// Upper bound on objective 2
                lambda1=0.0, // Weight of first objective
                lambda2=0.0, // Weight of second objective
                LNP    =0.0, // Local Nadir point value wrt the current weight vector
                WLNP   =0.0, // Worst local Nadir point wrt the current weight vector
                ObjV   =0.0; // Objective function value of cplex.
        unsigned long iterations = 0;
        IloExpr NoGood = IloExpr( theModel.env ); // IloExpression used to build the no good inequalities
        std::pair<double,double> p; // Pair used to store outcome vector of a solution
        std::vector<double> Sol(NumOfVars), oldSol(NumOfVars);



        if ( !OnlyOneNonDomSol )
        {
            auto startTime = CPUclock::now ( );
            for ( auto SupIt = NonDomSet.SupNDs.begin (); std::next( SupIt ) != NonDomSet.SupNDs.end ( ); ++SupIt )
            {
                std::cout << "Triangle " << ++triangle << " of " << NumOfTriangles << std::endl;
                // Retrieve the bound of the current triangle
                f1_bound = std::next( SupIt )->getFirst ( ) - 1.0;
                f2_bound = SupIt->getSecond ( ) - 1.0;
                // Set the bound in the cplex model
                theModel.f1.setBounds( SupIt->getFirst ( ) , f1_bound  );
                theModel.f2.setBounds( std::next ( SupIt )->getSecond ( ) , f2_bound );

                // Calculate the slope of the search direction
                lambda1 = SupIt->getSecond ( ) - std::next ( SupIt )->getSecond ( );
                lambda2 = std::next ( SupIt )->getFirst ( ) - SupIt->getFirst ( );
                // Set the objective function coefficients according to it and nextSol
                theModel.OBJ.setLinearCoef( theModel.f1 , lambda1 );
                theModel.OBJ.setLinearCoef( theModel.f2 , lambda2 );

                // Set time limit
                theModel.cplex.setParam(IloCplex::ClockType , 2);
                theModel.cplex.setParam( IloCplex::Param::TimeLimit , totalTime);

                // As long as cplex solves the problem, we continue to rank
                while ( theModel.cplex.solve ( ) )
                {
                    ++iterations; // Iterations counter is incremented

                    /*=====================================================*/
                    /*      Calculate the total time consumption and       */
                    /*      to the time limit                              */
                    /*=====================================================*/
                    auto nowTime = CPUclock::now ( );
                    double time = duration_cast<std::chrono::duration<double>>( nowTime - startTime ).count ( );
                    if ( totalTime < time ) break;

                    theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( );
                    // Retrieve the info of the current solution!
                    p.first = theModel.cplex.getValue( theModel.f1 );
                    p.second = theModel.cplex.getValue( theModel.f2 );
                    ObjV = theModel.cplex.getObjValue ( );

                    // Calculate the worst local NAdir point in the current triangle
                    WLNP = 0.0;
                    for ( auto it = NonDomSet.NDs.begin (); std::next(it)!= NonDomSet.NDs.end ( ); ++it )
                    {
                        if ( std::next(it)->getFirst ( ) > f1_bound ) break;
                        else if ( it->getFirst ( ) >= SupIt->getFirst ( ) )
                        {
                            LNP = lambda1 * std::next(it)->getFirst ( ) + lambda2 * (it->getSecond ( ) );
                            if ( LNP > WLNP ) WLNP = LNP;
                        }
                    }

                    // Retreive the current solution, and build the no good inequality
                    for ( int var = 0; var<NumOfVars; ++ var )
                    {
                        oldSol[var] = Sol[var];
                        if ( theModel.cplex.getValue ( theModel.AllVars[var] ) >= 0.5 )
                        {
                            Sol[var] = 1;
                            NoGood += (1 - theModel.AllVars[var] );
                        }
                        else
                        {
                            Sol[var] = 0;
                            NoGood += theModel.AllVars[var];
                        }
                    }
                    // Create a new solution, and update the non-dominated set
                    solution sol = solution ( false , p , Sol);
                    NonDomSet.updateNDS( sol );

                    // Calculate the hamming distance between the current and the previous solutions
                    int Diff = 0;
                    for ( int i=0; i<NumOfVars; ++i )
                    {
                        Diff += std::max ( Sol[i] - oldSol[i] , oldSol[i] - Sol[i] );
                    }

                    // Add the no-good inequality and  clear the iloexpr
                    theModel.model.add ( NoGood >= 1 );
                    NoGood.clear ( );

                    // If PrintProgress is true, print the progress:
                    if ( PrintProgress )
                    {
                        std::cout   << "it : " << iterations
                                    << "\t UB : " << WLNP
                                    << "\t LB : " << ObjV
                                    << "\t Gap : " << ( WLNP - ObjV ) / ObjV
                                    << "\t Diff : " << Diff
                                    << "\t Time : " << totalTime << std::endl;
                    }
                    // If the value of the worst local Nadir point exceeds the current objective function value, we can stop the search in the current trianle
                    if ( ObjV >= WLNP ) break;
                }

            }
        }
        theStatistics->NumberOfBranchingNodes += theModel.cplex.getNnodes ( );
        theStatistics->TotalNumberOfSolutions = NonDomSet.NDs.size ( );
        theStatistics->NumberOfPhaseTwoSolutions = theStatistics->TotalNumberOfSolutions - theStatistics->NumberOfPhaseOneSolutions;

        NoGood.end ( );
    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in RunPhaseTwoRanking in the tpm class : " << ie.getMessage ( ) << std::endl;
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in RunPhaseTwoRanking in the tpm class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Some crazy error I do not know happened and now I am terminating!\n";
        exit ( EXIT_FAILURE );
    }

}

/********************************************************************************************/
void tpm::printToFile( const std::string& fileName )
{
    try
    {
        // Copy file name using assignment
        FileName = fileName ;
        // Set the PrintToFile flag to true
        PrintToFile = true;
    }
    catch ( std::exception &e )
    {
        std::cerr << "Could not set the file name. Error : " << e.what ( ) << std::endl;
        PrintToFile = false;
    }
    catch ( ... )
    {
        std::cerr << "Could not set the file name." << std::endl;
        PrintToFile = false;
    }
}












