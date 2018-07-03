#ifndef NDS_H_INCLUDED
#define NDS_H_INCLUDED


/**
 * Class declaring a non dominated frontier for bi-objective integer programming problems
 */


//! C++ includes
#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<stdexcept>
#include<limits>

//! My own C++ includes
#include"solution.h"

class NDS{

       double TheWorstLocalNadirPoint;

    public:
        std::list< solution > NDs; //! Set of non-dominated solutions
        std::list< solution > SupNDs; //! List of supported non--dominated solutions

        /*!
         * Default constructor
         */
        NDS ( );
        /*!
         * Explicit copy constructor
         */
        NDS ( const NDS& other);

        /*!
         * Function returning an iterator to the first element on the list of supported non dominated points
         */
        inline
        std::list< solution >::iterator getFirstSupportedPoints ( ) { return SupNDs.begin(); }

        /*!
         * Function returning the end iterator for the list of supported non dominated solutions
         */
        inline
        std::list< solution >::iterator getSupportedPointsEnd ( ) { return SupNDs.end(); }

        /*!
         * Function returning true if the list of supported non-dominated solutions is empty
         */
        inline
        bool isEmpty ( ) { return SupNDs.empty(); }

        /*!
         * Function adding a solution to the list of supported non dominated solutions
         */
        void addToSupportedNDs ( const solution& p );

        /*!
         * Function overwriting the SupNDs list with list other of pairs of doubles
         * \param other List of pairs of doubles.
         * \note The function should be used if one knows the list of supported non dominated points from some other method and just wants to store these in the NDS class
         */
        inline
        void createSupportedNDs( const std::list< solution > & other ){ SupNDs.clear(); SupNDs = other; }

        /*!
         * Function copying the set of supported non dominated solutions into the set of non dominated solutions
         */
        inline
        void copySupToNonDom ( ){ for ( auto it = SupNDs.begin(); it!=SupNDs.end(); ++it ) NDs.push_back ( *it ); }

        /*!
         * Adds a point p to NDs if it is non--dominated by all points in UBset.
         * The set UBset is updated if it happens that the new point p dominated solutions in UBset.
         * \param sol solution. Contains the solution we want to test for non dominancy
         */
        void updateNDS(const solution &sol);

};

#endif // NDS include guard ends here
