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
           // SplineApproximation sp = new SplineApproximation();
           // double[] points = sp.GetPointsforTest();
           // double[,] Testresult = sp.TestSplineSinus(points);
          //  double[,] Testresult1 = sp.TestSplineSinus1(points);

          //  Console.WriteLine("Porównanie wartości");
          //  for (int i = 0; i < 100; i++)
         //   {
           //     Console.WriteLine("Wartość aproksymowana(stara metoda): " + Testresult[i, 0] + "| Wartość dokładna: " + Testresult[i, 1]);
         //       Console.WriteLine("Wartość aproksymowana(nowa metoda): " + Testresult1[i, 0] + "| Wartość dokładna: " + Testresult1[i, 1]);
          //      Console.WriteLine("Różnica(stara): " + (Testresult[i, 1] - Testresult[i, 0]));
          //      Console.WriteLine("Różnica(nowa): " + (Testresult1[i, 1] - Testresult1[i, 0]));
          //      Console.WriteLine("Różnica między metodami(nowa - stara): " + (Testresult1[i, 0] - Testresult[i, 0]));
           // }






            /*
            DensitiyFunction df = new DensitiyFunction();
          double[][] result = df.Test();
            SplineApproximation sp = new SplineApproximation(4,6,4,result[0],result[1]);
            double[] wsp = sp.Approximate1();
            Console.WriteLine("Aproksymacja");
            double ptk = 4;
            for(int i=0; i<150; i++)
            {
                Console.WriteLine("Wartosc funkcji w " + ptk + " to: " + sp.FunctionValue(wsp,ptk)/ result[1].Length);
                ptk += 0.1;
            }
            for (int i = 0; i < result[1].Length; i++)
            {
                Console.WriteLine("Prawdziwa wartosc "  + result[1][i]/ result[1].Length);
            }
            // funkcjonalność programu: punkty*/

        }
    }
}
