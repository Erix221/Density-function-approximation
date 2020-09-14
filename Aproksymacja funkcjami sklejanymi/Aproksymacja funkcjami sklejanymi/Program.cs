using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Numerics;


namespace Aproksymacja_funkcjami_sklejanymi
{
    class Program
    {
        static void Main(string[] args)
        {
            SplineApproximation sp = new SplineApproximation();
            double h = 1;
            double baseptk = 0;
            double x = 1.1;
            int n = 6;
            Dictionary<double, double> dic = new Dictionary<double, double>();
            dic.Add(0, Math.Sin(0));
            dic.Add(0.523, Math.Sin(0.523));
            dic.Add(1.57, Math.Sin(1.57));
            dic.Add(0.785, Math.Sin(0.785));
            dic.Add(1.046, Math.Sin(1.046));
            dic.Add(1.3089, Math.Sin(1.3089));
            dic.Add(0.2617, Math.Sin(0.2617));
            dic.Add(0.139, Math.Sin(0.139));
            //Console.WriteLine(sp.BaseSpline(baseptk,h,x));
            double[] xi = sp.Approximate(0,2*Math.PI,n,dic);
            //   for (int i = 0; i < n+3; i++)
            //  {
            //    for (int j = 0; j < n+3; j++)
            //     {
            //         Console.Write(string.Format("{0} ", xi[i, j]));
            //     }
            //     Console.Write(Environment.NewLine + Environment.NewLine);
            // }
            //  Console.ReadLine();
            Console.WriteLine("Współczynniki");
            for (int i = 0; i < xi.Length; i++)
                Console.WriteLine(xi[i]);



            double[] sinx = new double[6];
            double wartosc1 = Math.Sin(Math.PI / 2);
            double wartosc2 = xi[0]*sp.BaseSpline(0 - h, h, Math.PI / 2) + xi[1] * sp.BaseSpline(0, h, Math.PI / 2) + xi[2] * sp.BaseSpline(0 + h, h, Math.PI / 2) + xi[3] * sp.BaseSpline(0 + 2*h, h, Math.PI / 2) + xi[4] * sp.BaseSpline(0 + 3*h, h, Math.PI / 2);
            double wartosc3 = Math.Sin(0);
            double wartosc4 = xi[0] * sp.BaseSpline(0 - h, h, 0) + xi[1] * sp.BaseSpline(0, h, 0) + xi[2] * sp.BaseSpline(0 + h, h, 0) + xi[3] * sp.BaseSpline(0 + 2 * h, h, 0) + xi[4] * sp.BaseSpline(0 + 3 * h, h, 0);
            Console.WriteLine("Porównanie wartości");
            Console.WriteLine(wartosc1);
            Console.WriteLine(wartosc2);
            Console.WriteLine(wartosc3);
            Console.WriteLine(wartosc4);





        }
    }
}
