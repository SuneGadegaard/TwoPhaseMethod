#ifndef CPLEXMODEL_H_INCLUDED
#define CPLEXMODEL_H_INCLUDED

#include<vector>
#include<iostream>
#include<ilcplex/ilocplex.h>
#include<string>
typedef IloArray<IloNumVarArray>    IloVarMatrix;

class CplexModel{
    public:
        IloEnv env;             //!< IloEnv used to build the model on
        IloModel model;         //!< IloModel used to hold the model that should be solved
        IloCplex cplex;         //!< IloCplex used to solve all the subproblems arising in the tpm algorihm
        IloNumVar f1;           //!< Variable used to hold the first objective
        IloNumVar f2;           //!< Variable used to hold the second objective
        IloObjective OBJ;       //!< IloObjective used to hold the objective f1 + f2. It is needed so that objective function coefficients can later be changed
        IloNumVarArray AllVars; //!< Array of variables. Used to hold all the variables so that no-good inequalities can be generated in the tpm algorithm

        /*! \brief Constructor of the CplexModel class
         *
         * Constructor of the CplexModel class. Initializes all of the above public members, that is all the cplex stuff.
         */
        CplexModel();

        /*! \brief Destructor of the CplexModel class.
         *
         * Destructor deallocating all memory allocated in the CplexModel object's life time
         */
        ~CplexModel();

        /**
         * Implementing an instance of the bi-objective binary knapsack problem
         * \param n integer. Number of items in the KP problem
         * \param cap integer. The capacity of the knapsack
         * \param w constant reference to a vector of integers. The weight of each item
         * \param p1 onstant reference to a vector of integers. The profit of each item in the first objective
         * \param p2 onstant reference to a vector of integers. The profit of each item in the second objective
         * \note The size of w, p1, and p2 must be atleast n
         */
        void buildBOKP( int n, int cap, const std::vector<int>& w, const std::vector<int>& p1, const std::vector<int> p2 );

};

#endif // CPLEXMODEL_H_INCLUDED
