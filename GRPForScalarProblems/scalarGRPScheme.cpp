/*
The implementation of GRP scheme in the cases of scalar conservation laws,
where we assume that the flux function is strictly convex. The GRP
algorithm in scalar cases has a general form for which the numerical
procedures needn't to be significantly changed with different equations.
*/

#include "scalarGRPScheme.h"

void scalarGRPScheme::setUp(const fluxFun &flux, const fluxFun &derivativeOfFlux, const vector<double> &U_0, double delta_t, double endTime, double delta_x, double xRange_min, double xRange_max)
{
    flux_ = flux;
    derivativeOfFlux_ = derivativeOfFlux;
    U_0_ = U_0;
    U_.push_back(U_0_);
    if (delta_t > 0)
    {
        delta_t_ = delta_t;
    }
    if (endTime >= delta_t)
    {
        endTime_ = endTime;
    }
    if (delta_x > 0)
    {
        delta_x_ = delta_x;
    }
    if (xRange_min < xRange_max)
    {
        xRange_min_ = xRange_min;
        xRange_max_ = xRange_max;
    }
}

// Always override this function in one's specific problems
vector<double> scalarGRPScheme::setBoundaryCondition(const vector<double> incompleteU, int timeLevel)
{
    if (incompleteU.empty())
    {
        return {};
    }
    vector<double> ret;
    ret.push_back(incompleteU[0]);
    ret.push_back(incompleteU[0]);
    ret.insert(ret.end(), incompleteU.begin(), incompleteU.end());
    ret.push_back(*(incompleteU.end() - 1));
    ret.push_back(*(incompleteU.end() - 1));
    ret.push_back(*(incompleteU.end() - 1));
    return ret;
}

void scalarGRPScheme::linearizeU_0()
{
    vector<double> s_0_Limited;
    for (auto i = 1; i < U_0_.size() - 1; i++)
    {
        s_0_Limited.push_back(minmod(U_0_[i + 1] - U_0_[i], U_0_[i] - U_0_[i - 1]) / delta_x_);
    }
    s_Limited_.push_back(s_0_Limited);
}

double scalarGRPScheme::minmod(double a, double b)
{
    if (a * b <= 0.0)
    {
        return 0.0;
    }
    if (a < 0)
    {
        return std::max(a, b);
    }
    return std::min(a, b);
}

double scalarGRPScheme::directLineSearch(const fluxFun &f, double start, double end)
{
    auto g = [&](double x) -> double
    {
        double temp = (1 - x) * start + x * end;
        return flux_(temp);
    };
    double a1 = 0.0, b1 = 1.0, l = 0.001,
           lambda = a1 + 0.382 * (b1 - a1), mu = a1 + 0.618 * (b1 - a1);
    int iter = 0;
    while (true)
    {
        iter++;
        if (b1 - a1 <= l || iter >= 50)
        {
            double ret = (1 - (b1 + a1) / 2) * start + (b1 + a1) / 2 * end;
            return ret;
        }
        else
        {
            if (g(lambda) > g(mu))
            {
                a1 = lambda;
                lambda = mu;
                mu = a1 + 0.618 * (b1 - a1);
                continue;
            }
            else
            {
                b1 = mu;
                mu = lambda;
                lambda = a1 + 0.382 * (b1 - a1);
                continue;
            }
        }
    }
}

vector<double> scalarGRPScheme::bySlopeLimiter(const vector<double> &s, const vector<double> &U)
{
    if (delta_x_ <= 0 || U.empty())
    {
        return {};
    }
    auto minmod3 = [&](double a, double b, double c) -> double
    {
        double min = std::min(std::min(a, b), c), max = std::max(std::max(a, b), c);
        if (min * max <= 0.0)
        {
            return 0.0;
        }
        if (min < 0)
        {
            return max;
        }
        return min;
    };
    vector<double> s_ret;
    for (auto i = 1; i < U.size() - 1; i++)
    {
        s_ret.push_back(minmod3(2 * (U[i + 1] - U[i]), 2 * (U[i] - U[i - 1]), delta_x_ * s[i - 1]) / delta_x_);
    }
    return s_ret;
}

vector<vector<double>> scalarGRPScheme::getNumericalResult()
{
    double curTime = delta_t_;
    int curTimeLevel = 0;
    auto gridNum = U_0_.size();
    linearizeU_0();
    while (curTime <= endTime_)
    {
        vector<double> vec_U_RiemannSol;
        vector<double> vec_U_derivatives;
        vector<double> vec_f_midTime_atStagered;
        vector<double> vec_SlopesNextTimeLevel;
        vector<double> vec_UNextTimeLevel;
        vector<double> vec_U_atStageredNext;
        for (auto i = 1; i < gridNum - 2; i++)
        {
            double U_atStagered_l = U_[curTimeLevel][i] + delta_x_ / 2 * s_Limited_[curTimeLevel][i - 1];
            double U_atStagered_r = U_[curTimeLevel][i + 1] - delta_x_ / 2 * s_Limited_[curTimeLevel][i];
            double RiemannSol = 0.0;
            double U_derivative = 0.0;
            double U_midTime_atStagered = 0.0;
            double f_midTime_atStagered = 0.0;
            if (U_atStagered_l > U_atStagered_r + 0.00001)
            {
                double shockDir = (U_atStagered_l - U_atStagered_r) / (flux_(U_atStagered_l) - flux_(U_atStagered_r));
                if (abs(shockDir) <= 0.00001)
                {
                    double temp1 = derivativeOfFlux_(U_atStagered_r);
                    double temp2 = derivativeOfFlux_(U_atStagered_l);
                    shockDir = (-temp1 * temp1 * s_Limited_[curTimeLevel][i - 1] + temp2 * temp2 * s_Limited_[curTimeLevel][i - 2]) / (U_atStagered_r - U_atStagered_l);
                }
                if (shockDir <= -0.00001)
                {
                    RiemannSol = U_atStagered_r;
                }
                else if (shockDir >= 0.00001)
                {
                    RiemannSol = U_atStagered_l;
                }
            }
            else if (U_atStagered_r > U_atStagered_l + 0.00001)
            {
                double minSpeed = derivativeOfFlux_(U_atStagered_l);
                double maxSpeed = derivativeOfFlux_(U_atStagered_r);
                if (minSpeed >= 0.00001)
                {
                    RiemannSol = U_atStagered_l;
                }
                else if (maxSpeed <= -0.00001)
                {
                    RiemannSol = U_atStagered_r;
                }
                else
                {
                    RiemannSol = directLineSearch(flux_, U_atStagered_l, U_atStagered_r);
                }
            }
            else
            {
                RiemannSol = U_atStagered_l;
            }
            vec_U_RiemannSol.push_back(RiemannSol);
            if (derivativeOfFlux_(RiemannSol) > 0)
            {
                U_derivative = -1 * derivativeOfFlux_(RiemannSol) * s_Limited_[curTimeLevel][i - 1];
            }
            else if (derivativeOfFlux_(RiemannSol) < 0)
            {
                U_derivative = -1 * derivativeOfFlux_(RiemannSol) * s_Limited_[curTimeLevel][i];
            }
            else
            {
                U_derivative = 0.0;
            }
            vec_U_derivatives.push_back(U_derivative);
            U_midTime_atStagered = RiemannSol + delta_t_ / 2 * U_derivative;
            f_midTime_atStagered = flux_(RiemannSol) + delta_t_ / 2 * derivativeOfFlux_(RiemannSol) * U_derivative;
            vec_f_midTime_atStagered.push_back(f_midTime_atStagered);
            vec_U_atStageredNext.push_back(RiemannSol + delta_t_ / 2 * U_derivative);
        }
        for (auto i = 2; i < gridNum - 3; i++)
        {
            vec_UNextTimeLevel.push_back(U_[curTimeLevel][i] - delta_t_ / delta_x_ * (vec_f_midTime_atStagered[i - 1] - vec_f_midTime_atStagered[i - 2]));
            vec_SlopesNextTimeLevel.push_back((vec_U_atStageredNext[i - 1] - vec_U_atStageredNext[i - 2]) / delta_x_);
        }
        vec_UNextTimeLevel = setBoundaryCondition(vec_UNextTimeLevel, curTimeLevel);
        vec_SlopesNextTimeLevel.insert(vec_SlopesNextTimeLevel.begin(), minmod(vec_UNextTimeLevel[1] - vec_UNextTimeLevel[0], vec_UNextTimeLevel[2] - vec_UNextTimeLevel[1]) / delta_x_);
        vec_SlopesNextTimeLevel.insert(vec_SlopesNextTimeLevel.begin(), minmod(vec_UNextTimeLevel[gridNum - 1] - vec_UNextTimeLevel[gridNum - 2], vec_UNextTimeLevel[gridNum - 2] - vec_UNextTimeLevel[gridNum - 3]) / delta_x_);
        vec_SlopesNextTimeLevel.insert(vec_SlopesNextTimeLevel.begin(), minmod(vec_UNextTimeLevel[gridNum - 2] - vec_UNextTimeLevel[gridNum - 3], vec_UNextTimeLevel[gridNum - 3] - vec_UNextTimeLevel[gridNum - 4]) / delta_x_);
        s_Limited_.push_back(bySlopeLimiter(vec_SlopesNextTimeLevel, vec_UNextTimeLevel));
        U_.push_back(vec_UNextTimeLevel);
        curTimeLevel++;
        curTime += delta_t_;
    }
    return U_;
}
