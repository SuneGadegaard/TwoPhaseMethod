#include"solution.h"
/************************************************************************************/
solution::solution():supported(false),p(std::pair<double,double>(0.0,0.0)){}

/************************************************************************************/
solution::solution ( bool supported ):supported(supported){}

/************************************************************************************/
solution::solution ( bool supported, const std::pair<double,double>& p): supported(supported), p(p) {}

/************************************************************************************/
solution::solution ( bool supported, const std::pair<double,double>& p, const std::vector<double>& var ): supported(supported), p(p), var(var) {}

/************************************************************************************/
void solution::getVarValues ( std::vector<double> &VarVector ) const {
    try{
        if ( !VarVector.empty() ) VarVector.clear(); // Clear the incomning vector to be sure

        if ( var.empty ( ) ) return; // If var-vector is empty, just return
        else for( auto it= var.begin(); it!=var.end(); ++it ) VarVector.push_back( *it ); // Else, copy the content of var to VarVector

    }catch( std::exception &e ){
        std::cerr << "Exception in getVarVector in the solution class : " << e.what() << std::endl;
    }
}
