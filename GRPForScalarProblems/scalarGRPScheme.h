/*
The interface of GRP scheme in the cases of scalar conservation laws,
where we assume that the flux function is strictly convex. The GRP
algorithm in scalar cases has a general form for which the numerical
procedures needn't to be significantly changed with different equations.
*/

#ifndef SCALARGRPSCHEME_H
#define SCALARGRPSCHEME_H

#include <vector>
#include <functional>
#include <algorithm>

using std::vector;
using fluxFun = std::function<double(double)>;

class scalarGRPScheme
{
public:
    void setUp(const fluxFun &flux, const fluxFun &derivativeOfFlux, const vector<double> &U_0, double delta_t, double endTime, double delta_x, double xRange_min, double xRange_max);
    vector<vector<double>> getNumericalResult();

protected:
    // Prequisitions for the launch of the GRP scheme
    fluxFun flux_;
    fluxFun derivativeOfFlux_;
    vector<double> U_0_;
    double delta_t_ = 0.01, endTime_ = 1.0, delta_x_ = 0.01, xRange_min_ = -1.0, xRange_max_ = 1.0;

    /*
    For different equations, one should define the boundary conditions accordingly.
    It's necessary to override this function as one's need.
    */  
    virtual vector<double> setBoundaryCondition(const vector<double> incompleteU, int timeLevel);

    // The main result is stored here
    vector<vector<double>> s_Limited_, U_;

    // The first step of GRP scheme is to linearize the discrete initial data
    void linearizeU_0();

    // minmod function
    double minmod(double a, double b);

    /*
    The line search method returns the optimal point of a function in a given range,
    which is required in the sonic cases.
    */
    double directLineSearch(const fluxFun &f, double start, double end);

    // The common slope limiter algorithm that is used in the GRP scheme.
    vector<double> bySlopeLimiter(const vector<double> &s, const vector<double> &U);
};

#endif