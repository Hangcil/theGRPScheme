#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include "./GRPForScalarProblems/scalarGRPScheme.h"
using std::ofstream;

/*
Here we consider the Burgers equation u_t+u*u_x=0 with periodic initial data u(x,t=0)=sin(2*\pi*x).
The calculation range of x is [0.0,1.0] and the evolution time interval is [0.0,1.0]. The grid is setup
as \Delta x=1/22 and \Delta t=0.5*\Delta x.
*/

class example3_11 : public scalarGRPScheme
{
protected:
    /*
    The bisection method used for calculating the root of a function.
    It's simple yet satisfyingly robust.
    Complementing the boundary values that are not calculated by the GRP scheme finally reduces
    to a problem of finding the root of a function.
    */
    double bisection(const fluxFun &f, double start, double end)
    {
        double a = start, b = end, c = 0.0;
        int i = 0;
        while (i <= 10)
        {
            i++;
            c = (a + b) / 2.0;
            if (f(c) >= 0.0001)
            {
                b = c;
            }
            else if (f(c) <= -0.0001)
            {
                a = c;
            }
            else
            {
                return c;
            }
        }
        return c;
    }

    /*
    The overrided function setBoundaryCondition() of the original in the base class scalarGRPSheme.
    As for specific problems, one should set different boundary conditions corespondingly.
    5 points of values (the first 2 and last 3 points) should be complemented in the scalar GRP problem.
    */
    vector<double> setBoundaryCondition(const vector<double> incompleteU, int timeLevel) override
    {
        double t = (timeLevel + 1) * delta_t_;
        double dx = delta_x_;
        vector<double> ret;
        ret.push_back(0.0);
        auto auxFun1 = [&](double x) -> double
        {
            return t * std::sin(2 * M_PI * x) + x - dx;
        };
        auto auxFun2 = [&](double x) -> double
        {
            return t * std::sin(2 * M_PI * x) + x + 2 * dx - 1.0;
        };
        auto auxFun3 = [&](double x) -> double
        {
            return t * std::sin(2 * M_PI * x) + x + dx - 1.0;
        };
        ret.push_back(std::sin(2 * M_PI * bisection(auxFun1, 0.0, 0.1)));
        ret.insert(ret.end(), incompleteU.begin(), incompleteU.end());
        ret.push_back(std::sin(2 * M_PI * bisection(auxFun2, 0.9, 1.0)));
        ret.push_back(std::sin(2 * M_PI * bisection(auxFun3, 0.9, 1.0)));
        ret.push_back(0.0);
        return ret;
    }
};

int main()
{
    fluxFun f = [](double u) -> double
    {
        return 0.5 * u * u;
    };
    fluxFun derivativeOfF = [](double u) -> double
    {
        return u;
    };
    example3_11 solver3_11;
    double dx = 1.0 / 22.0;
    vector<double> U_0;
    for (double i = 0.0; i <= 1.0; i += dx)
    {
        U_0.push_back(std::sin(2 * M_PI * i));
    }
    solver3_11.setUp(f, derivativeOfF, U_0, dx / 2.0, 1.0, dx, 0.0, 1.0);
    auto result = solver3_11.getNumericalResult();

    // Output the numeric results, which can be read and plotted by other scientific programs.
    ofstream oFile("./example3_11.txt");
    if(oFile)
    {
        for(auto i=0;i<45;i++)
        {
            for(auto j=0;j<23;j++)
            {
                oFile<<result[i][j]<<' ';
            }
            oFile<<'\n';
        }
        oFile.close();
    }

    return 0;
}
