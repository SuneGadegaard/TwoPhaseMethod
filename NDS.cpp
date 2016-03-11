#include"NDS.h"

/*
 * C++ implementation of NDs.h
 * author Sune Lauth Gadegaard
 * version 1.0
 */

/********************************************************************************************/
NDS::NDS(){}

/********************************************************************************************/
NDS::NDS( const NDS& other ):NDs ( other.NDs ), SupNDs ( other.SupNDs ) {}

/********************************************************************************************/
void NDS::addToSupportedNDs( const solution &s ){
    try
    {
        bool insertP    = true;
        auto insertIt   = SupNDs.begin();

        for( auto it = SupNDs.begin( ); it!=SupNDs.end() && insertP; ++it ){
            insertP = true;//!( (it->getFirst()-0.1<=s.getFirst()) && (it->getSecond()-0.1<=s.getSecond()) );
            s.getFirst();
            if (    (it->getFirst() < s.getFirst() && std::next(it)==SupNDs.end() ) ||
                    (it->getFirst() < s.getFirst() && std::next(it)->getFirst() >= s.getFirst() ) ){
                        insertIt = it;
                    }
        }
        if ( insertP ){
            if ( std::next(insertIt)==SupNDs.end() ) SupNDs.push_back( s ); // If insertIt incremented once is the end, push back
            else SupNDs.insert (std::next(insertIt),s); // Else, insert after insertIT
        }
    }catch( std::exception &e ){
        std::cerr << "Exception in addToSupportedNDs in the NDS class : " << e.what ( ) << std::endl;
        exit ( 1 );
    }
}

/********************************************************************************************/
void NDS::updateNDS( const solution &sol)
{
    try
    {
        bool pShouldInsert = true;

        if ( NDs.size( ) == 0 ){
            NDs.push_back( sol );
            return;
        }
        auto it = NDs.begin();
        while( it != NDs.end() ){
            // Check each element if p is dominated
            if ( ( it->getFirst( )-0.1<=sol.getFirst( ) ) && ( it->getSecond( )-0.1<=sol.getSecond( ) ) ){
                // sol is dominated by *it and nothing more should be done!
                //std::cout << "Should not be added\n";
                pShouldInsert = false;
                break;
            // Now check if sol dominates *it
            }else if ( (( it->getFirst( )>sol.getFirst( ) ) && (it->getSecond( )>=sol.getSecond() )) || ( (it->getFirst()>=sol.getFirst()) && (it->getSecond()>sol.getSecond()) ) ){
                // sol dominates *it. We erase it, insert sol instead and check the following points for domination
                auto insertIt = NDs.erase ( it );
                NDs.insert ( insertIt , sol );
                // We insert p instead of the one it dominates. Therefore, we should not insert it again
                pShouldInsert = false;
                // We loop through the rest of the vector of solutions to see if more solutions are dominated.
                // If so, we delete the solution. The erase function returns an iterator to the next element on the list
                // Therefore we automatically increment the iterator. Furthermore, we loop only as long as the solution is dominated by the new solution
                for ( auto itt = std::next(it); itt!=NDs.end(); ){
                    std::cout << "Starting to erase\n";
                    if ( (sol.getFirst()<=itt->getFirst()+0.1) && (sol.getSecond() <= itt->getSecond()+0.1) )
                    {
                        itt = NDs.erase(itt);
                    }
                    else ++itt;
                }
                if ( NDs.size() >= 2 ){
                    TheWorstLocalNadirPoint =  std::numeric_limits<double>::min( ) ;
                    for ( auto it = std::next(NDs.begin()); it!=NDs.end(); ++it ){
                        double NadirPoint = it->getFirst( ) + std::prev(it)->getSecond( );
                        if ( NadirPoint > TheWorstLocalNadirPoint ) TheWorstLocalNadirPoint = NadirPoint;
                    }
                }
            }else ++it;
        }
        // If p survives the check, it is added at the end and the pareto front is sorted
        // Definitely not the most efficient way of doing this. Consider using a list instead
        if ( pShouldInsert ){
            auto it = NDs.begin ( );
            while ( it != NDs.end() )
            {
                if ( it->getFirst( ) > sol.getFirst( ) )
                { // As sol parsed the tests above, we know that sol should be inserted before the first solution which has a first coordinate larger than sol itself
                    NDs.insert( it , sol );
                    break;
                }
                ++it;
            }

            TheWorstLocalNadirPoint = std::numeric_limits<double>::min( ) ;
            for ( auto it = std::next(NDs.begin()); it!=NDs.end(); ++it ){
                double NadirPoint = it->getFirst( ) + std::prev(it)->getSecond( );
                if ( NadirPoint > TheWorstLocalNadirPoint ) TheWorstLocalNadirPoint = NadirPoint;
            }
        }

    }catch(std::exception &e)
    {
        std::cerr << "Exception in updateNDS in the NDS class : " << e.what ( ) << std::endl;
        exit ( EXIT_FAILURE );
    }
}
