using System;
using System.Collections.Generic;

namespace Spline_Approximation_dll
{
    public class DensityFunction
    {
        public double StartOfInterval { get; set; }
        public double EndOfInterval { get; set; }
        public double Neighbourhood { get; set; }
        public double StepSize { get; set; }
        public double[] X { get; set; }

        public DensityFunction()
        {

        }

        public DensityFunction(double a, double b, double prze, double ptk, double[] x)
        {
            StartOfInterval = a;
            EndOfInterval = b;
            Neighbourhood = prze;
            StepSize = ptk;
            X = x;
        }

        private void Getlowest(double current, double[] points, ref int lIndex)
        {
            for (int i = lIndex; i < points.Length; i++)
            {
                if (points[i] >= current - Neighbourhood && points[i] <= current + Neighbourhood)
                {
                    lIndex = i;
                    break;
                }
                else
                    lIndex = 0;
            }
        }

        private void Gethighest(double current, double[] points, ref int hIndex)
        {
            for (int i = hIndex; i < points.Length - 1; i++)
            {
                if (points[i] <= current + Neighbourhood && points[i] >= current - Neighbourhood)
                {
                    if (!(points[i + 1] <= current + Neighbourhood && points[i] >= current - Neighbourhood))
                    {
                        hIndex = i;
                        break;
                    }
                    else
                        hIndex = i + 1;
                }
            }
        }

  
        public double[][] CalculateDensityFunctionPoints()
        {
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            double[] points = X;
            double current = StartOfInterval, count = 0;
            int lIndex = 0, hIndex = 0;
            double[][] result = new double[2][];
            Array.Sort(points);
            while (current <= EndOfInterval)
            {
                if (current - Neighbourhood <= points[points.Length - 1] && points[0] <= current + Neighbourhood) 
                {
                    Getlowest(current, points, ref lIndex);
                    Gethighest(current, points, ref hIndex);
                    if (lIndex == hIndex)
                    {
                        if (!(points[lIndex] <= current + Neighbourhood && points[lIndex] >= current - Neighbourhood))
                            count = 0;
                        else
                            count = 1;
                    }
                    else
                        count = hIndex - lIndex + 1;
                    x.Add(current);
                    y.Add(count);
                }
                else
                {
                    x.Add(current);
                    y.Add(0);
                }
                current = Math.Round(current + StepSize,15);

            }
            result[0] = x.ToArray();
            result[1] = y.ToArray();
            for (int i = 0; i < result[1].Length; i++)
                result[1][i] /= X.Length; 

            return result;
        }
    }
}
