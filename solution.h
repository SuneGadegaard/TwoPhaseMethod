#ifndef SOLUTION_H_INCLUDED
#define SOLUTION_H_INCLUDED

#include<vector>
#include<iostream>
#include<stdexcept>

class solution{
    private:
        bool supported;     //! Flag indicating if a point is supported
        std::pair<double,double> p; //! pair representing the outcome vector
        std::vector<double> var;    //! vector holding the variable value of the solution
    public:
        solution ( );   //! Empty constructor
        solution ( bool supported );    //! Constructor only setting the support-flag
        solution ( bool supported, const std::pair<double,double>& p ); //! Constructor initializing the sopport flag and the outcome vector
        solution ( bool supported, const std::pair<double,double>& p, const std::vector<double>& var );//! Constructor initializing the sopport flag, the outcome vector and the var vector

        /*!
         * Function querying if the solution is supported or not
         * \return true iff the solution is known to be supported.
         */
        inline
        bool isSupported ( ) const { return supported; }

        /*!
         * Function returning the firrst objective function value of the solution
         * \return double The value of the first objective function value
         */
        inline
        double getFirst ( ) const { return p.first; }

        /*!
         * Function returning the second objective function value of the solution
         * \return double The value of the first objective function value
         */
        inline
        double getSecond ( ) const { return p.second; }

        /*!
         * Function returning the vector of variable values
         * \param VarVector vector of doubles. Equals var on output
         */
        void getVarValues( std::vector<double> &VarVector ) const;
};

#endif // SOLUTION_H_INCLUDED
