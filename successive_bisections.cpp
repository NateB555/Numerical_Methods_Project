#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

float simple_power(float x, int exponent)
{
    float answer = x;
    for (int i = 2; i <= exponent; i++)
    {answer *= x;}
    return answer;
}


// float secondOrderEquation(float a1, float a2, float a3, float x)
// {    
//     return a1*pow(x,2) + a2*x + a3;
// }

//! need to integrate a way to try different guesses (p14) 

double solve_2nd_Order_equation(std::function<double(double)> f, double precision)
{
    double f1; double f2;
    // first guesses:
    double x1 = 10;
    double x2 = 15;
    int trial = 1;
    do
    {
        // std::cout << "x1: " << x1 << ", f1: " << f(x1) << std::endl;
        f1 = f(x1);
        // std::cout << "x2: " << x2 << ", f2: " << f(x2) << std::endl;
        f2 = f(x2);
        // f1 = secondOrderEquation(a1, a2, a3, x1);
        // f2 = secondOrderEquation(a1, a2, a3, x2);
        if ((f1 == 0) || (f2 == 0)) {break;} // condition in case x1 or x2 happens to be the right answer
        if (std::abs(f2) > std::abs(f1)) {x2 -= 0.25;} // 
        else {x2 += 0.25;}
        trial++; 
    } while (f1*f2 > 0);    // makes sure one of the trial functions is negative; i.e: there is one on each side of 0;
    
    if (f1 == 0) {return x1;}
    if (f2 == 0) {return x2;}

    double x3; double f3;
    int iteration = 1;

    do
    {
        x3 = (x1+x2)/2;
        f3 = f(x3);
        // f3 = secondOrderEquation(a1, a2, a3, x3);
        std::cout << iteration << ".  " << "x = " << x3 << " : f(x) = " << f3 << endl;
        if ((abs(f3)>precision) && (f1*f3>0)) {x1=x3;}
        else if ((abs(f3)>precision) && (f1*f3<0)) {x2=x3;}
        iteration++;
        // cout << std::abs(f3) << endl;

    } while(std::abs(f3)>=precision);

    // && (abs(x3-x2)>=0.00001))
    return x3;
}

double derivative(std::function<double(double)> f)
{
    
}




int main()
{
    // 
    std::cout << "f(x) = 4x^2 + 2x = 5" << endl;
    auto f = [](double x) 
    {
        // return 4*pow(x,2) + 2*x - 5;
        return x - 50*exp(-0.1*x);
    };
    double result = solve_2nd_Order_equation(f, 0.00001);
    std::cout << "\n x = " << result << endl;

    return 0;
}