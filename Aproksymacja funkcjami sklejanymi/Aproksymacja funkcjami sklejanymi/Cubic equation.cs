using System;

namespace Spline_Approximation_dll
{
    public class Cubic_equation
    {
        public double[] Factors { get; set; }


        public Cubic_equation(double[] coe)
        {
            Factors = new double[4];
            Factors[0] = coe[0];
            Factors[1] = coe[1];
            Factors[2] = coe[2];
            Factors[3] = coe[3];
        }

        public double[] CalculateRoots()
        {
            double f, g, h;
            double[] roots;
            f = (Factors[2] / Factors[0]) - ((Factors[1] * Factors[1]) / (3 * (Factors[0] * Factors[0])));
            g = ((2 * Factors[1] * Factors[1] * Factors[1]) / (27 * Factors[0] * Factors[0] * Factors[0])) 
                - ((Factors[1] * Factors[2]) / (3 * (Factors[0] * Factors[0]))) + (Factors[3] / Factors[0]);
            h = ((g * g) / 4) + ((f * f * f) / 27);
            if (h > 0)
            {
                double x1 = -(g / 2) + Math.Sqrt(h);
                double x2 = -(g / 2) - Math.Sqrt(h);
                double x3 = Factors[1] / (3 * Factors[0]);
                double root1 = Math.Cbrt(x1) + Math.Cbrt(x2) - x3;
                roots = new double[1];
                roots[0] = root1;
                return roots;

            }
            else if (f == 0 && g == 0)
            {
                double root1 = Math.Cbrt(Factors[3] / Factors[0]);
                roots = new double[3];
                for (int i = 0; i < 3; i++)
                    roots[i] = root1;
                return roots;
            }
            else
            {
                double i = Math.Sqrt(((g * g) / 4) - h);
                double j = Math.Cbrt(i);
                double k = Math.Acos(-g / (2 * i));
                double m = Math.Cos(k / 3);
                double n = Math.Sqrt(3) * Math.Sin(k / 3);
                double p = (-Factors[1]) / (3 * Factors[0]);
                double root1 = (2 * j * m) + p;
                double root2 = (-j * (m + n)) + p;
                double root3 = (-j * (m - n)) + p;
                roots = new double[3];
                roots[0] = root1;
                roots[1] = root2;
                roots[2] = root3;
                return roots;
            }
        }
    }
}
