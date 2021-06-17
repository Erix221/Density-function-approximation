using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using System;
using System.Collections.Generic;

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
        public double[] ApproximationKnots { get; set; }
        public ThirdDegreeSpline[] BaseSplines { get; set; }



        public SplineApproximation()
        {

        }

        public double[] GetBasePtks()
        {
            return ApproximationKnots;
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
            H = Math.Round((b - a) / n, 15);
            X = x;
            Y = y;
        }


        private double[,] GetFactors()
        {
            double sum;
            double[,] factors = new double[NumberOfIntervals + 3, NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
                for (int j = 0; j < NumberOfIntervals + 3; j++)
                {
                    sum = 0;
                    foreach (var ptk in X)
                        sum += BaseSplines[i].Value(ptk) * BaseSplines[j].Value(ptk);
                    factors[i, j] = sum;
                }
            }
            return factors;
        }


        private double[] GetConstantTerms()
        {
            double sum = 0;
            double[] constantTerms = new double[NumberOfIntervals + 3];
            for (int i = 0; i < NumberOfIntervals + 3; i++)
            {
                sum = 0;
                for (int j = 0; j < X.Length; j++)
                    sum += Y[j] * BaseSplines[i].Value(X[j]);
                constantTerms[i] = sum;
            }
            return constantTerms;
        }

        private double[] CreateKnots()
        {
            double[] knots = new double[NumberOfIntervals + 3];
            knots[0] = Math.Round(StartOfInterval - H,15); 
            knots[NumberOfIntervals + 2] = Math.Round(EndOfInterval + H,15);
            for (int i = 0; i < NumberOfIntervals + 2; i++)
            {
                knots[i + 1] = Math.Round(StartOfInterval + H * i, 15);
            }  
            return knots;
        }

       private ThirdDegreeSpline[] CreateBaseSplines()
        {
            List<ThirdDegreeSpline> splines = new List<ThirdDegreeSpline>();
            List<CubicFunction> cubics = new List<CubicFunction>();
            double[] coe;
            foreach(var knot in ApproximationKnots)
            {
                coe = BaseSplineCoefficientsReal(knot, 1);
                cubics.Add( new CubicFunction(coe[0], coe[1], coe[2], coe[3]));
                coe = BaseSplineCoefficientsReal(knot, 2);
                cubics.Add(new CubicFunction(coe[0], coe[1], coe[2], coe[3]));
                coe = BaseSplineCoefficientsReal(knot, 3);
                cubics.Add(new CubicFunction(coe[0], coe[1], coe[2], coe[3]));
                coe = BaseSplineCoefficientsReal(knot, 4);
                cubics.Add(new CubicFunction(coe[0], coe[1], coe[2], coe[3]));
                splines.Add(new ThirdDegreeSpline(Math.Round(knot - 2 * H,15), Math.Round(knot + 2 * H, 15), 4, cubics.ToArray()));
                cubics.Clear();
            }
            return splines.ToArray();    
        }

        private ThirdDegreeSpline CreateApproximatedSpline(double[] coe)
        {
            List<CubicFunction> cubics = new List<CubicFunction>();
            for(int i = 0; i<ApproximationKnots.Length - 3; i++)
            {
                cubics.Add(BaseSplines[i].CubicFunctions[CheckRangeOfBaseSpline(i,i)]*coe[i] 
                    + BaseSplines[i + 1].CubicFunctions[CheckRangeOfBaseSpline(i + 1,i)] * coe[i + 1]
                    + BaseSplines[i + 2].CubicFunctions[CheckRangeOfBaseSpline(i + 2,i)] * coe[i + 2]
                    + BaseSplines[i + 3].CubicFunctions[CheckRangeOfBaseSpline(i + 3,i)] * coe[i + 3]);
            }
            return new ThirdDegreeSpline(StartOfInterval, EndOfInterval, NumberOfIntervals, cubics.ToArray());
        }

        private int CheckRangeOfBaseSpline(int spline, int i)
        {
            return BaseSplines[spline].CheckRange((ApproximationKnots[i + 2] + ApproximationKnots[i + 1]) / 2);
        }

       public ThirdDegreeSpline Approximate()
        {
            ApproximationKnots = CreateKnots();
            BaseSplines = CreateBaseSplines();

            double[,] factors = GetFactors();
            double[] constantTerms = GetConstantTerms();

            var A = Matrix<double>.Build.DenseOfArray(factors);
            var B = Vector<double>.Build.Dense(constantTerms);
            CoefficientForApproximatedFunction = A.SolveIterative(B, new MlkBiCgStab()).ToArray();

            return CreateApproximatedSpline(CoefficientForApproximatedFunction);
        }


        public double[] BaseSplineCoefficientsReal(double knot, int interval)
        {
            double h2 = H * H;
            double h3 = h2 * H;
            double c = (1 / h3);
            double[] coefficients;
            switch (interval)
            {
                case 1:
                    double knot2m = Math.Round(knot - 2 * H, 15);
                    coefficients = new[]{ c, c * (-3) * knot2m,
                                            c * 3 * knot2m * knot2m,
                                           -c * knot2m * knot2m * knot2m };
                    return coefficients;
                case 2:
                    double knot1m = Math.Round(knot - H, 15);
                    coefficients = new[]{ c * (-3), c * 3 * (H + 3 * knot1m),
                                          c * 3 * (h2 - knot1m *(2 * H - 3 * knot1m)),
                                          c * (3 * (knot1m * (knot1m * knot1m + knot1m * H - h2)) + h3) };
                    return coefficients;
                case 3:
                    double knot1p = Math.Round(knot + H, 15);
                    coefficients = new[]{ c * 3, c * 3 * (H - 3 * knot1p),
                                          c * 3 * (-h2 - knot1p * (2 * H + 3 * knot1p)),
                                          c * (3 * (knot1p*( -knot1p * knot1p  + knot1p * H +  h2)) + h3) };
                    return coefficients;
                case 4:
                    double knot2p = Math.Round(knot + 2*H, 15);
                    coefficients = new[]{ -c, c * 3 * knot2p,
                                          -c * 3 * knot2p * knot2p,
                                           c * knot2p * knot2p * knot2p };
                    return coefficients;
            }
            coefficients = new[]{ 0.0, 0.0, 0.0, 0.0 };
            return coefficients;
        }

    }


}
