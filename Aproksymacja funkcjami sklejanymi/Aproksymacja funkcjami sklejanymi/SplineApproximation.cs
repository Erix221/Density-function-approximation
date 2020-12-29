using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using System.Diagnostics;

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

        private double[,] GetFactors1(double[] array, double[] BasePtks, double h)
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
                        sum += BaseSpline1(BasePtks[i], h, ptk) * BaseSpline1(BasePtks[j], h, ptk);
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

        private double[] GetConstantTerms1(double[] BasePtks, double h)
        {
            double sum = 0;
            double[] constantTerms = new double[NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
                sum = 0;
                for (int j = 0; j < X.Length; j++)
                    sum += Y[j] * BaseSpline1(BasePtks[i], h, X[j]);
                constantTerms[i] = sum;
            }
            return constantTerms;
        }

        public double[] Approximate1()
        {
            double h = (EndOfInterval - StartOfInterval) / NumberOfIntervals;
            double[] BasePtks = new double[NumberOfIntervals + 3];
            BasePtks[0] = StartOfInterval - h; BasePtks[NumberOfIntervals + 2] = EndOfInterval + h;
            for (int i = 1; i < NumberOfIntervals + 2; i++)
                BasePtks[i] = BasePtks[i - 1] + h;

            double[,] factors = GetFactors1(X, BasePtks, h);
            double[] constantTerms = GetConstantTerms1(BasePtks,h);

            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            return A.SolveIterative(B, new MlkBiCgStab()).ToArray(); 
        }

        public double[] Approximate()
        {
            double h = (EndOfInterval - StartOfInterval) / NumberOfIntervals;
            double[] BasePtks = new double[NumberOfIntervals + 3];
            BasePtks[0] = StartOfInterval - h; BasePtks[NumberOfIntervals + 2] = EndOfInterval + h;
            for (int i = 1; i < NumberOfIntervals + 2; i++)
                BasePtks[i] = BasePtks[i - 1] + h;

            double[,] factors = GetFactors(X, BasePtks, h);
            double[] constantTerms = GetConstantTerms(BasePtks, h);

            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            return A.SolveIterative(B, new MlkBiCgStab()).ToArray();
        }

        private static double BaseSpline1(double baseptk, double h, double x)
        {
            double h2 = h * h;
            double h3 = h2 * h;
            double constant = (1 / h3);
            double interval1 = x - (baseptk - 2*h);
            double interval2 = x - (baseptk - h);
            double interval3 = (baseptk + h) - x;
            double interval4 = (baseptk + 2*h) - x;
            if ((x >= baseptk - 2 * h) && (x <= baseptk - h))
                return constant * (interval1* interval1* interval1);
            if ((x >= baseptk - h) && (x <= baseptk))
                return constant * (h3 + 3*((h2 * (interval2)) + (h * (interval2* interval2)) - (interval2* interval2* interval2)));
            if ((x <= baseptk + h) && (x >= baseptk))
                return constant * (h3 + 3 * (h2 * (interval3) + h *(interval3* interval3) - (interval3* interval3* interval3)));
            if ((x <= baseptk + 2 * h) && (x >= baseptk + h))
                return constant * (interval4* interval4* interval4);
            else
                return 0;

        }

        private static double BaseSpline(double baseptk, double h, double x)
        {
            double constant = (1 / Math.Pow(h, 3));
            if ((x >= baseptk - 2 * h) && (x <= baseptk - h))
                return constant * Math.Pow(x - (baseptk - 2 * h), 3);
            if ((x >= baseptk - h) && (x <= baseptk))
                return constant * (Math.Pow(h, 3) + (3 * Math.Pow(h, 2) * (x - (baseptk - h))) + (3 * h * Math.Pow(x - (baseptk - h), 2)) - (3 * Math.Pow(x - (baseptk - h), 3)));
            if ((x <= baseptk + h) && (x >= baseptk))
                return constant * (Math.Pow(h, 3) + 3 * Math.Pow(h, 2) * ((baseptk + h) - x) + 3 * h * Math.Pow((baseptk + h) - x, 2) - 3 * Math.Pow((baseptk + h) - x, 3));
            if ((x <= baseptk + 2 * h) && (x >= baseptk + h))
                return constant * Math.Pow((baseptk + 2 * h) - x, 3);
            else
                return 0;
        }


        public double FunctionValue(double[] result, double x)
        {
            double value = 0;
            for (int i = 0; i < result.Length; i++)
                value += result[i] * BaseSpline(StartOfInterval + (H * (-1 + i)), H, x);
            return value;
        }

        public double FunctionValue1(double[] result, double x)
        {
            double value = 0;
            for (int i = 0; i < result.Length; i++)
                value += result[i] * BaseSpline1(StartOfInterval + (H * (-1 + i)), H, x);
            return value;
        }

        public double[] GetPointsforTest() 
        {
            Random r = new Random();
            double range = 2*Math.PI;
            double rDouble = r.NextDouble() * range; //for doubles
            double[] x = new double[1000000];
            for (int i = 0; i < 1000000; i++)
            {
                x[i] =(r.NextDouble() * range);
            }
            return x;
        }

    }


}
