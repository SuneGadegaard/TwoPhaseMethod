# A general two phase method for bi-objective combinatorial optimization problems.

This repository contains an implementation of a two phase method for bi-objective combinatorial optimization problems. The main routine is implemented in tpm.cpp and tpm.h ((T)wo (P)hase (M)ethod). It implements two two phase algorithms differing in the algorithms used in phase two. The first, default algorithm uses a perpendicular search method while the second uses a ranking algorithm. If the ranking algorithm is used, the problem solved MUST be a pure binary optimization problem, as no-good constraints are used in the ranking. Note that the tpm methods are implemented as minimization problems, implying if one objective should be maximized instead, one needs to use the "negative" of the objective function. Furthermore, it is assumed that the objective functions of the problems takes only integer values. If this assumtion is not met, the program does not guarantee an optimal solution!

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

# How to use the program
The programs have been tested on a Linux Ubuntu 14.04 LTS machine. All codes have been compiled using the GNU gcc compilers with optimization options O3 and the C++11 flag enabled. You need to link CPLEX to the programs, and a guide to this is using the Code::blocks IDE is found here: https://www-304.ibm.com/support/docview.wss?uid=swg21449771 .

In order to use the program, you need a method for reading in your data, and you need to implement your linear integer optmization problem in the CplexModel class. After this, you simply hand the instance of the CplexModel class to the tpm class' run function and the problem is solved.

# An example
In the main.cpp file an example solving the bi-objective knapsack problem is given. First data for the problem is generated. Then an instance of the CplexModel class is created and the self-implemented buildBOKP function is called to build the bi-objective knapsak problem. Then an instance of the tpm class is created. Af the instance is created we set the "printProgress" falg to true by calling printProgress() and we tell the tpm instance that we want the solution printet to the file "TheOutputFile.txt". The we run the two phase algorithm by calling the RUN () function. Finally, the test statistics are printet to screen.

