#define _USE_MATH_DEFINES
#include <fstream>
#include "./GRPForScalarProblems/scalarGRPScheme.h"
using std::ofstream;

/*
Here we consider the Burgers equation u_t+u*u_x=0 with the initial data of the step function.
The calculation range of x is [0.0,1.0] and the evolution time interval is [0.0,1.0]. The grid is setup
as \Delta x=0.04 and \Delta t=0.5*\Delta x.
*/

class example3_12 : public scalarGRPScheme
{
protected:
    /*
    The overrided function setBoundaryCondition() of the original in the base class scalarGRPSheme.
    As for specific problems, one should set different boundary conditions corespondingly.
    5 points of values (the first 2 and last 3 points) should be complemented in the scalar GRP problem.
    */
    vector<double> setBoundaryCondition(const vector<double> incompleteU, int timeLevel) override
    {
        vector<double> ret;
        ret.push_back(1.0);
        ret.push_back(1.0);
        ret.insert(ret.end(), incompleteU.begin(), incompleteU.end());
        ret.push_back(0.0);
        ret.push_back(0.0);
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
    example3_12 solver3_12;
    double dx = 0.04;
    vector<double> U_0;
    for (double i = 0.0; i <= 1.0; i += dx)
    {
        if (i <= 0.22)
        {
            U_0.push_back(1.0);
        }
        else
        {
            U_0.push_back(0.0);
        }
    }
    solver3_12.setUp(f, derivativeOfF, U_0, dx / 2.0, 0.8, dx, 0.0, 1.0);
    auto result = solver3_12.getNumericalResult();

    // Output the numeric results, which can be read and plotted by other scientific programs.
    ofstream oFile("./example3_12.txt");
    if (oFile)
    {
        for (auto i = 0; i < 40; i++)
        {
            for (auto j = 0; j < 25; j++)
            {
                oFile << result[i][j] << ' ';
            }
            oFile << '\n';
        }
        oFile.close();
    }

    return 0;
}
