using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using System;

namespace Spline_Approximation_dll
{
    public class SplineApproximation
    {
        public double StartOfInterval { get; set; }
        public double EndOfInterval { get; set; }
        public int NumberOfIntervals { get; set; }
        public double[] CoefficientForApproximatedFunction { get; set; }
        public double H { get; set; }
        public double[] X { get; set; }
        public double[] Y { get; set; }
        public double[] BasePtks { get; set; }



        public SplineApproximation()
        {

        }

        public double[] GetBasePtks()
        {
            return BasePtks;
        }

        public double GetH()
        {
            return H;
        }

        public int GetNumberOfIntervals()
        {
            return NumberOfIntervals;
        }

        public double GetStartOfInterval()
        {
            return StartOfInterval;
        }

        public double GetEndOfInterval()
        {
            return EndOfInterval;
        }

        public SplineApproximation(double a, double b, int n, double[] x, double[] y)
        {
            StartOfInterval = a;
            EndOfInterval = b;
            NumberOfIntervals = n;
            X = x;
            Y = y;
        }


        private double[,] GetFactors(double[] BasePtks, double h)
        {
            double sum;
            double[,] factors = new double[NumberOfIntervals + 3, NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
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
            decimal decimalH = (decimal)((EndOfInterval - StartOfInterval) / NumberOfIntervals);
            double h = Convert.ToDouble(decimalH);
            H = h;
            decimal[] BasePtksdecimal = new decimal[NumberOfIntervals + 3];
            double[] BasePtk = new double[NumberOfIntervals + 3];
            BasePtksdecimal[0] = (decimal)StartOfInterval - decimalH; BasePtksdecimal[NumberOfIntervals + 2] = (decimal)EndOfInterval + decimalH;
            for (int i = 1; i < NumberOfIntervals + 2; i++)
                BasePtksdecimal[i] = BasePtksdecimal[i - 1] + decimalH;
            BasePtk[0] = Convert.ToDouble(BasePtksdecimal[0]); BasePtk[NumberOfIntervals + 2] = Convert.ToDouble(BasePtksdecimal[1]);
            for(int i = 1; i < NumberOfIntervals + 2; i++)
                BasePtk[i] = Convert.ToDouble(BasePtksdecimal[i]);
            BasePtks = BasePtk;


            double[,] factors = GetFactors(BasePtks, h);
            double[] constantTerms = GetConstantTerms(BasePtks, h);

            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            CoefficientForApproximatedFunction = A.SolveIterative(B, new MlkBiCgStab()).ToArray();
            return CoefficientForApproximatedFunction;
        }


        private static double BaseSpline(double baseptk, double h, double x)
        {
            double h2 = h * h;
            double h3 = h2 * h;
            double constant = (1 / h3);
            double interval1 = x - (baseptk - 2 * h);
            double interval2 = x - (baseptk - h);
            double interval3 = (baseptk + h) - x;
            double interval4 = (baseptk + 2 * h) - x;
            if ((x >= baseptk - 2 * h) && (x <= baseptk - h))
                return constant * (interval1 * interval1 * interval1);
            if ((x >= baseptk - h) && (x <= baseptk))
                return constant * (h3 + 3 * ((h2 * interval2) + (h * (interval2 * interval2)) - (interval2 * interval2 * interval2)));
            if ((x <= baseptk + h) && (x >= baseptk))
                return constant * (h3 + 3 * (h2 * (interval3) + h * (interval3 * interval3) - (interval3 * interval3 * interval3)));
            if ((x <= baseptk + 2 * h) && (x >= baseptk + h))
                return constant * (interval4 * interval4 * interval4);
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

        public double[] BaseSplineCoefficients(double start, double end, double baseptk)
        {
            double h2 = H * H;
            double h3 = h2 * H;
            double constant = (1 / h3);
            decimal decimalstart = (decimal)(start);
            decimal decimalend = (decimal)(end);
            decimal decimalbaseptk = (decimal)(baseptk);
            decimal baseptk2m = (decimal)(baseptk - H - H);
            decimal baseptk1m = (decimal)(baseptk - H);
            decimal baseptk1p = (decimal)(baseptk + H );
            decimal baseptk2p = (decimal)(baseptk + H + H);
            if ((decimalstart >= baseptk2m) && (decimalend <= baseptk1m)) 
            {
                double[] coefficients = { constant, constant * (-3) * (baseptk - 2 * H), constant * 3 * (baseptk - 2 * H) * (baseptk - 2 * H),
                                          constant * (-1) * (baseptk - 2 * H) * (baseptk - 2 * H) * (baseptk - 2 * H) };
                return coefficients;
            }
            if ((decimalstart >= baseptk1m) && (decimalend <= decimalbaseptk))
            {
                double[] coefficients = { constant * (-3), constant * (3 * H + 9 * (baseptk - H)), 
                                          constant * (3 * (h2) - 6 * (baseptk - H) * H - 9 * ((baseptk - H) * (baseptk - H))),
                                          constant * (3 * (baseptk - H) * (baseptk - H) * (baseptk - H) + 3 * (baseptk - H) 
                                          * (baseptk - H) * H - 3 * (baseptk - H) * h2 + h3) };
                return coefficients;
            }
            if ((decimalend <= baseptk1p) && (decimalstart >= decimalbaseptk)) 
            {
                double[] coefficients = { constant * 3, constant * (3 * H - 9 * (baseptk + H)), 
                                          constant * ((-3) * (h2) - 6 * (baseptk + H) * H + 9 * ((baseptk + H) * (baseptk + H))), 
                                          constant * ((-3) * (baseptk + H) * (baseptk + H) * (baseptk + H) + 3 
                                          * (baseptk + H) * (baseptk + H) * H + 3 * (baseptk + H) * h2 + h3) };
                return coefficients;
            }
            if ((decimalend <= baseptk2p) && (decimalstart >= baseptk1p)) 
            {
                double[] coefficients = { -constant, constant * 3 * (baseptk + 2 * H), 
                                          -constant * 3 * (baseptk + 2 * H) * (baseptk + 2 * H),
                                           constant * (baseptk + 2 * H) * (baseptk + 2 * H) * (baseptk + 2 * H) };
                return coefficients;
            }
            double[] coe = { 0, 0, 0, 0 };
            return coe;
        }


        public double[] CoefficientsforEquation(double[] result, double start, double end)
        {
            double[] coefficients = new double[4];
            double[] coe;
            for (int i = 0; i < result.Length; i++)
            {
                coe = BaseSplineCoefficients(start, end, StartOfInterval + (H * (-1 + i)));
                coefficients[0] += (result[i] * coe[0]);
                coefficients[1] += (result[i] * coe[1]);
                coefficients[2] += (result[i] * coe[2]);
                coefficients[3] += (result[i] * coe[3]);
            }
            return coefficients;
        }

    }


}
