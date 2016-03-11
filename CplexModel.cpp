#include"CplexModel.h"

/*****************************************************************************************/
CplexModel::CplexModel()
{
    model = IloModel ( env );
    cplex = IloCplex ( model );
    cplex.setParam ( IloCplex::Param::TimeLimit , 3600 );
    cplex.setParam(IloCplex::EpAGap , 0.0 );
    cplex.setParam(IloCplex::EpGap , 0.0 );
    f1 = IloNumVar ( env, 0 , IloInfinity, ILOFLOAT , "f1" );
    f2 = IloNumVar ( env, 0 , IloInfinity, ILOFLOAT , "f2" );
    OBJ = IloMinimize ( env , f1 + f2 );
    AllVars = IloNumVarArray( env );
}

/*****************************************************************************************/
CplexModel::~CplexModel()
{
    if ( AllVars.getImpl() != nullptr )
    {   // End all variables and end the array used to store them
        AllVars.endElements ( );
        AllVars.end ( );
    }
    if ( f1.getImpl ( ) != nullptr )
    {   // End f1
        f1.end ( );
    }
    if ( f2.getImpl ( ) != nullptr )
    {   // end f2
        f2.end ( );
    }
    if ( OBJ.getImpl ( ) )
    {   // End the objective function
        OBJ.end ( );
    }
    if ( cplex.getImpl ( ) )
    {   // End the cplex instance
        cplex.end ( );
    }
    if ( model.getImpl ( ) )
    {   // End the model instance
        model.end ( );
    }
    if ( env.getImpl ( ) )
    {   // end the environment
        env.end ( );
    }

}

/*****************************************************************************************/
void CplexModel::buildBOKP( int n, int cap, const std::vector<int>& w, const std::vector<int>& p1, const std::vector<int> p2 )
{
    try{
        // Initialize IloExprs to hold objective and the constraint.
        IloExpr Obj1 ( env ), Obj2 ( env ), cst ( env );
        // Initialize binary variables
        IloNumVarArray x = IloNumVarArray ( env, n , 0 , 1 , ILOBOOL );

        // Build  Objectives and constraints. Note that we use -p1[i] and -p2[i] as the tpm algorithm is implemented as a minimization problem
        for ( int i=0; i<n; ++i )
        {
            AllVars.add ( x[i] );
            Obj1 += p1[i]*x[i];
            Obj2 += p2[i]*x[i];
            cst += w[i]*x[i];
        }

        // Add the objective function to the model
        model.add ( OBJ );

        // Add the individual objective functions as constraints
        model.add ( Obj1 == f1 );
        model.add ( Obj2 == f2 );

        // Add the knapsack constraint
        model.add ( cst >= cap );

        Obj1.end ( );
        Obj2.end ( );
        cst.end ( );
    }catch(IloException &ie){
        std::cerr << "IloException in the buildSSCFLP of the CplexModelClass : " << ie.getMessage ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }catch ( std::exception &e){
        std::cerr << "Exception in the buildSSCFLP of the CplexModelClass : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}

