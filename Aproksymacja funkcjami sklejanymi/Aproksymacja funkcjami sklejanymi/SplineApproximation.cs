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
        public double StartOfInterval { get; set; }
        public double EndOfInterval { get; set; }
        public int NumberOfIntervals { get; set; }

        public double H { get; }
        public double[] X { get; set; }
        public double[] Y { get; set; }


        public SplineApproximation()
        {

        }

        public SplineApproximation(double a, double b, int n, double[] x, double[] y)
        {
            StartOfInterval = a;
            EndOfInterval = b;
            NumberOfIntervals = n;
            H = (EndOfInterval - StartOfInterval) / NumberOfIntervals;
            X = x;
            Y = y;
        }

        private double[,] GetFactors(double[] array, double[] BasePtks, double h)
        {
            double sum = 0; 
            double[,] factors = new double[NumberOfIntervals + 3, NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
                sum = 0;
                for (int j = 0; j < NumberOfIntervals + 3; j++)
                {
                    sum = 0;
                    foreach (var ptk in X)
                        sum += BaseSpline(BasePtks[i], h, ptk) * BaseSpline(BasePtks[j], h, ptk);
                    factors[i, j] = sum;
                }
            }
            return factors;
        }

        private double[] GetConstantTerms(double[] BasePtks, double h)
        {
            double sum = 0;
            double[] constantTerms = new double[NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
                sum = 0;
                for (int j = 0; j < X.Length; j++)
                    sum += Y[j] * BaseSpline(BasePtks[i], h, X[j]);
                constantTerms[i] = sum;
            }
            return constantTerms;
        }

        public double[] Approximate()
        {
            double h = (EndOfInterval - StartOfInterval) / NumberOfIntervals;
            double[] BasePtks = new double[NumberOfIntervals + 3];
            BasePtks[0] = StartOfInterval - h; BasePtks[NumberOfIntervals + 2] = EndOfInterval + h;
            for (int i = 1; i < NumberOfIntervals + 2; i++)
                BasePtks[i] = BasePtks[i - 1] + h;

            double[,] factors = GetFactors(X, BasePtks, h);
            double[] constantTerms = GetConstantTerms(BasePtks,h);

            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            return A.SolveIterative(B, new MlkBiCgStab()).ToArray(); 
        }

        private static double BaseSpline(double baseptk, double h, double x)
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
        private double FunctionValue(double[] result, double x)
        {
            double value = 0;
            for (int i = 0; i < result.Length; i++)
                value += result[i] * BaseSpline(StartOfInterval + (H * (-1 + i)), H, x);
            return value;
        }
        public double[,] TestSplineSinus()
        {
            Random r = new Random();
            double range = Math.PI / 2;
            double rDouble = r.NextDouble() * range; //for doubles
            double[] x = new double[100];
            for(int i = 0; i < 100; i++)
            {
                x[i] = r.NextDouble() * range;
            }
            double[] y = new double[100];
            for (int i = 0; i < 100; i++)
            {
                y[i] = Math.Sin(x[i]);
            }
            SplineApproximation sp = new SplineApproximation(0, Math.PI / 2, 6, x, y);
            double h = (sp.EndOfInterval - sp.StartOfInterval) / sp.NumberOfIntervals;
            double[] result = sp.Approximate();
            for (int i = 0; i < result.Length; i++)
                Console.WriteLine("Współczynnik: "+ i +" " +result[i] + " Ilosc wspolczynników: " + result.Length);
            double[,] sinsprox = new double[100, 2];
            for(int i = 0; i<100; i++)
            {
                Console.WriteLine("Punkt:" + x[i]);
                //sinsprox[i, 0] = result[0] * BaseSpline(sp.StartOfInterval - h, h, x[i]) + result[1] * BaseSpline(sp.StartOfInterval, h, x[i]) + result[2] * BaseSpline(sp.StartOfInterval + h, h, x[i]) + result[3] * BaseSpline(sp.StartOfInterval + 2 * h, h, x[i]) + result[4] * BaseSpline(sp.StartOfInterval + 3 * h, h, x[i]) + result[5] * BaseSpline(sp.StartOfInterval + 4 * h, h, x[i]) + result[6] * BaseSpline(sp.StartOfInterval + 5 * h, h, x[i]) + result[7] * BaseSpline(sp.StartOfInterval + 6 * h, h, x[i]);
                sinsprox[i, 0] = sp.FunctionValue(result,x[i]);
                sinsprox[i, 1] = Math.Sin(x[i]);
            }
            return sinsprox;
        }


    }


}
