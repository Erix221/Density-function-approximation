using System;
using System.Collections.Generic;
using System.Linq;

namespace Spline_Approximation_dll
{
    public class IntersectionPoints
    {
        public SplineApproximation FirstFunction { get; set; }
        public SplineApproximation SecondFunction { get; set; }
        public double[] Ranges { get; set; }

        public IntersectionPoints()
        {

        }

        public IntersectionPoints(SplineApproximation First, SplineApproximation Second)
        {
            FirstFunction = First;
            SecondFunction = Second;
        }

        public bool CheckforOverlap()
        {
            bool overlap = FirstFunction.GetStartOfInterval() < SecondFunction.GetEndOfInterval() && SecondFunction.GetStartOfInterval() < FirstFunction.GetEndOfInterval();
            return overlap;
        }

        public double[] CalculateRanges()
        {
            double[] FirstBasePtks = FirstFunction.GetBasePtks();
            double[] SecondBasePtks = SecondFunction.GetBasePtks();
            double[] result = FirstBasePtks.Concat(SecondBasePtks).ToArray();
            Array.Sort(result);

            double FirstPoint = FirstBasePtks[1] <= SecondBasePtks[1] ? FirstBasePtks[1] : SecondBasePtks[1];
            double LastPoint = FirstBasePtks[FirstBasePtks.Length - 2] >= SecondBasePtks[SecondBasePtks.Length - 2] 
                               ?  FirstBasePtks[FirstBasePtks.Length - 2] : SecondBasePtks[SecondBasePtks.Length - 2];
            int Startindex = Array.BinarySearch(result, FirstPoint);
            int Endindex = Array.BinarySearch(result, LastPoint);
            double[] Endresult = new double[Endindex - Startindex + 1];
            Array.Copy(result, Startindex, Endresult, 0, Endindex - Startindex + 1);
            double[] final = Endresult.Distinct().ToArray();

            return final;
        }

        public double[] CalculateIntersectionPoints(double[] ApproximatedCoefficientsFirst, double[] ApproximatedCoefficientsSecond)
        {
            if (CheckforOverlap())
            {
                Ranges = CalculateRanges();
                List<double> points = new List<double>();
                for (int i = 0; i < Ranges.Length - 1; i++)
                {
                    double[] coe1 = FirstFunction.CoefficientsforEquation(ApproximatedCoefficientsFirst, Ranges[i], Ranges[i + 1]);
                    double[] coe2 = SecondFunction.CoefficientsforEquation(ApproximatedCoefficientsSecond, Ranges[i], Ranges[i + 1]);
                    double[] coefinal = new double[4];
                    for (int j = 0; j < 4; j++)
                        coefinal[j] = coe1[j] - coe2[j];
                    Cubic_equation cb = new Cubic_equation(coefinal);
                    double[] roots = cb.CalculateRoots();
                    for (int k = 0; k < roots.Length; k++)
                        if (roots[k] >= Ranges[i] && roots[k] <= Ranges[i + 1])
                            points.Add(roots[k]);
                }
                return points.ToArray();
            }
            double[] roots1 = new double[0];
            roots1[0] = double.NaN;
            return roots1;
        }

        public double[] GetRanges()
        {
            return Ranges;
        }

    }

}
