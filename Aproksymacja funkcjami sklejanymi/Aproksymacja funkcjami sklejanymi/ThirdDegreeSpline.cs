using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Spline_Approximation_dll
{
    public class ThirdDegreeSpline
    {
        public double StartOfInterval { get; set; }
        public double EndOfInterval { get; set; }
        public int NumberOfIntervals { get; set; }
        public double H { get; set; }
        public double[] Knots { get; set; }
        public CubicFunction[] CubicFunctions { get; set; }
        public CubicFunction ZeroFunction { get; set; }

        public ThirdDegreeSpline(double start, double end, int n, CubicFunction[] functions)
        {
            StartOfInterval = start;
            EndOfInterval = end;
            NumberOfIntervals = n;
            H = Math.Round((EndOfInterval - StartOfInterval) / NumberOfIntervals,15);
            CubicFunctions = functions;
            ZeroFunction = new CubicFunction(0, 0, 0, 0);
            Knots = CalculateKnots();
        }


        public double[] CalculateRanges(ThirdDegreeSpline spline)
        {
            double[] FirstBasePtks = Knots;
            double[] SecondBasePtks = spline.Knots;
            List<double> result = FirstBasePtks.Concat(SecondBasePtks).ToList();
            result.Sort();

            double FirstPoint = FirstBasePtks[0] >= SecondBasePtks[0] ? FirstBasePtks[0] : SecondBasePtks[0];
            double LastPoint = FirstBasePtks[^1] <= SecondBasePtks[^1] ? FirstBasePtks[^1] : SecondBasePtks[^1];
            int startIndex = result.BinarySearch(FirstPoint);
            int endIndex = result.BinarySearch(LastPoint);
            List<double> final = result.GetRange(startIndex, endIndex - startIndex);
            return final.Distinct().ToArray();
        }

        public double[] CalculateIntersectionPoints(ThirdDegreeSpline spline)
        {
            if (CheckforOverlap(spline))
            {
                double[] Ranges = CalculateRanges(spline);
                CubicFunction[] functions = new CubicFunction[Ranges.Length - 1];
                for (int i = 0; i < Ranges.Length - 1; i++)
                {
                    functions[i] = GetFunction((Ranges[i + 1] + Ranges[i])/2) 
                                   - spline.GetFunction((Ranges[i + 1] + Ranges[i]) / 2);
                }                    
                    
                ThirdDegreeSpline newSpline = new ThirdDegreeSpline(Ranges[0], Ranges[^1],
                                                                    Ranges.Length, functions);
                return newSpline.CalculateIntersectionPointsOX();
            }
            double[] roots = new double[0];
            return roots;
        }

        public double[] CalculateIntersectionPointsOX()
        {
            List<double> finalroots = new List<double>();
            foreach (var func in CubicFunctions)
            {
                double[] roots = func.CalculateRealRoots();
                for (int i = 0; i < roots.Length; i++)
                {
                    if((roots[i] >= StartOfInterval) && (roots[i] <= EndOfInterval))
                    finalroots.Add(roots[i]);
                }
                    
            }                  
                return finalroots.ToArray();
        }

        public bool CheckforOverlap(ThirdDegreeSpline spline)
        {
            return StartOfInterval < spline.EndOfInterval && spline.StartOfInterval < EndOfInterval;
        }

        private double[] CalculateKnots()
        {
            double[] knots = new double[NumberOfIntervals + 1];
            for (int i = 0; i < NumberOfIntervals + 1; i++)
                knots[i] = Math.Round(StartOfInterval + H * (i), 15);
            return knots;
        }

        public double Value(double x)
        {
            int Check = CheckRange(x);
            if ((Check >= 0) && (Check < CubicFunctions.Length))
                return CubicFunctions[Check].Value(x);
            else
                return 0;
        }

        public double FirstDerivative(double x)
        {
            int Check = CheckRange(x);
            if ((Check >= 0) && (Check < CubicFunctions.Length))
                return CubicFunctions[Check].FirstDerivative(x);
            else
                return 0;
        }

        public int CheckRange(double x)
        {
            double result = ((x - StartOfInterval) / H);
            if (result < 0)
                return -1;
            return (int)result;
        }

        public CubicFunction GetFunction(double x)
        {
            int Check = CheckRange(x);
            if ((Check >= 0) && (Check < CubicFunctions.Length))
                return CubicFunctions[Check];
            else
                return ZeroFunction;

        }
    }
}
