#ifndef DISTRIBUTION_MATH_H
#define DISTRIBUTION_MATH_H

long double Beta_Continued_Fraction( long double x, long double a,
                                                               long double b);
long double xBeta_Distribution(double x, double a, double b);

double Beta_Function(double a, double b);
long double xBeta_Function(long double a, long double b);
double Entire_Incomplete_Gamma_Function(double x, double nu);
long double xEntire_Incomplete_Gamma_Function(long double x, long double nu);

long double xSmall_x(long double x, long double nu);
long double xMedium_x(long double x, long double nu);
long double xLarge_x(long double x, long double nu);

double Factorial(int n);
long double xFactorial(int n);
int Factorial_Max_Arg( void );


double Gamma_Function(double x);
long double xGamma_Function(long double x);
double Gamma_Function_Max_Arg( void );
long double xGamma_Function_Max_Arg( void );

double Gamma_Distribution(double x, double nu);

long double xGamma(long double x);
long double Duplication_Formula( long double two_x );

double Ln_Beta_Function(double a, double b);
long double xLn_Beta_Function(long double a, long double b);
double Ln_Gamma_Function(double x);
long double xLn_Gamma_Function(long double x);

long double xLnGamma_Asymptotic_Expansion( long double x );


double Chi_Square_Distribution(double x, int n);

double F_Distribution(double f, int v1, int v2);

double F_Density(double x, int v1, int v2 );


#endif // DISTRIBUTION_MATH_H
