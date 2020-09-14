using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;

namespace Aproksymacja_funkcjami_sklejanymi
{
    public class SplineApproximation
    {
        public double[] Approximate(double a, double b, int n, Dictionary<double,double> points)
        {
            //przygotowanie tablic i przedziałów
            double h = (b - a) / n;
            double[] x = new double[n + 3];
            double sum = 0;
            x[0] = a - h; x[n + 2] = b + h;
            for (int i = 1; i < n + 2; i++)
                x[i] = x[i - 1] + h;
            double[,] factors = new double[n + 3, n + 3];
            //współczynniki
            for(int i = 0; i < n + 3; i++)
            {
                sum = 0;
                for (int j = 0; j < n + 3; j++)
                {
                    sum = 0;
                    foreach (var ptk in points)
                        sum += BaseSpline(x[i], h, ptk.Key) * BaseSpline(x[j], h, ptk.Key);
                    factors[i, j] = sum;
                }                   
            }
            //wyrazy wolne
            //sum = 0;
            double[] constantTerms = new double[n + 3];
            for(int k = 0; k<n+3; k++)
            {
                sum = 0;
                foreach (var ptk in points)
                    sum += ptk.Value * BaseSpline(x[k], h, ptk.Key);
                constantTerms[k] = sum;
            }
            //obliczanie ukladu rownan;
            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            var wynik = A.SolveIterative(B, new MlkBiCgStab()); ;

            return wynik.AsArray();
        }
        
        public double BaseSpline(double baseptk, double h, double x)
        {
            double constant = (1 / Math.Pow(h, 3));
            if ((x >= baseptk - 2 * h) && (x <= baseptk - h))
                return constant * Math.Pow(x - (baseptk - 2 * h),3);
            if((x >= baseptk - h) && (x <= baseptk))
                return constant * (Math.Pow(h, 3) + (3*Math.Pow(h, 2)*(x - (baseptk - h))) + (3*h*Math.Pow(x - (baseptk - h), 2)) - (3*Math.Pow(x - (baseptk - h), 3)));
            if ((x <= baseptk + h) && (x >= baseptk))
                return constant * (Math.Pow(h, 3) + 3*Math.Pow(h, 2)*((baseptk + h) - x) + 3*h*Math.Pow((baseptk + h) - x, 2) - 3*Math.Pow((baseptk + h) - x, 3));
            if ((x <= baseptk + 2 * h) && (x >= baseptk + h))
                return constant * Math.Pow((baseptk + 2 * h) - x,3);
            else
                return 0;  
        }
    }


}
