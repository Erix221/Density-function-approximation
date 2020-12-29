using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace Aproksymacja_funkcjami_sklejanymi
{
    public class DensitiyFunction
    {
        public double StartOfInterval { get; set; }
        public double EndOfInterval { get; set; }
        public double Neighbourhood { get; set; }
        public double StepSize { get; set; }
        public double[] X { get; set; }

        public DensitiyFunction()
        {

        }

        public DensitiyFunction(double a, double b, double prze, double ptk, double[] x)
        {
            StartOfInterval = a;
            EndOfInterval = b;
            Neighbourhood = prze;
            StepSize = ptk;
            X = x;
        }

  public double[][] CalculateBinary()
        {
            double current = StartOfInterval;
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            double[] points = X;
            double[][] result = new double[2][];
            Array.Sort(points);
            double count = GFG.CountInRange(points, points.Length, StartOfInterval- Neighbourhood, StartOfInterval+ Neighbourhood);
            double newcount;
            x.Add(StartOfInterval);
            y.Add(count);
            current += StepSize;
            int counter = 0;
            while (current <= EndOfInterval)
            {
                newcount = GFG.CountInRange(points, points.Length, current - Neighbourhood, current + Neighbourhood);
                if (newcount != count || counter == 3)
                {
                    counter = 0;
                    x.Add(current);
                    y.Add(newcount);
                    count = newcount;
                }
                else
                    counter += 1;
                current += StepSize;

            }
            result[0] = x.ToArray();
            result[1] = y.ToArray();
            return result;
        }

        public double[][] CalculateSmart()
        {
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            double[] points = X;
            double current = StartOfInterval, count;
            int lIndex = 0, hIndex = 0;
            double[][] result = new double[2][];
            Array.Sort(points);
            while(current <= EndOfInterval)
            {
                for (int i = lIndex; i < points.Length; i++)
                {
                    if (points[i] >= current - Neighbourhood)
                    {
                        lIndex = i;
                        break;
                    }             
                }
                for (int i = hIndex; i < points.Length; i++)
                {
                    if (points[i] <= current - Neighbourhood)
                    {
                        hIndex = i;
                        break;
                    }
                }
                count = hIndex - lIndex + 1;
                x.Add(current);
                y.Add(count);
                current += StepSize;

            }    
            result[0] = x.ToArray();
            result[1] = y.ToArray();
            return result;
        }
    }
}
