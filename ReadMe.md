# TwoPhaseMethod

This repository contains an implementation of a two phase method for bi-objective combinatorial optimization problems. The main routine is implemented in tpm.cpp and tpm.h ((T)wo (P)hase (M)ethod). It implements two two phase algorithms differing in the algorithms used in phase two. The first, default algorithm uses a perpendicular search method while the second uses a ranking algorithm. If the ranking algorithm is used, the problem solved MUST be a pure binary optimization problem, as no-good constraints are used in the ranking. Note that the tpm methods are implemented as minimization problems, implying if one objective should be maximized instead, one needs to use the "negative" of the objective function.

The following subclases are needed for the program to run

CplexModel -> Implements the boco problem which should be solved by the algorithm implemented in tmp. It is up to the user ti implement this class. The only requirements are the following public members:

    IloEnv env, holds the environment on which the model is build
    IloModel model holds the actual IloModel that needs to be solved, that is model is initialize as model ( env )
    IloCplex cplex is initialized as cplex ( model )
    IloNumVar f1 holds the first objective.
    IloNumVar f2 holds the second objective.
    IloObjective OBJ holds the objective function minimize( f1 + f2 ) on start up.
    Note that all of the above is taken care of in the constructor of the CplexModel class
    IloNumVarArray AllVars contains all DECISION variables. That is, it does not contain f1 and f2
    If AllVars contains variables which are not binary, the tpm algorithm with doRanking turned on, will not behave as it should!
    CplexModel is implemented in CplexModel.h and CplexModel.cpp

NDS -> Implements a non dominated set class consisting of solutions. It should be fully functional, but you are more than welcome to report bugs. The NDS class is implemented in NDS.h and NDS.cpp solution -> Implements the solution class. It is implemented in solution.h and solution.cpp

