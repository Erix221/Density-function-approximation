using System;
using System.Collections.Generic;
using System.Text;

namespace Aproksymacja_funkcjami_sklejanymi
{
    class GFG
    {

        static int LowerIndex(double[] arr, int n,
                              double x)
        {

            int l = 0, h = n - 1;
            while (l <= h)
            {
                int mid = (l + h) / 2;
                if (arr[mid] >= x)
                    h = mid - 1;
                else
                    l = mid + 1;
            }
            return l;
        }

        static int UpperIndex(double[] arr, int n,
                              double y)
        {
            int l = 0, h = n - 1;
            while (l <= h)
            {
                int mid = (l + h) / 2;
                if (arr[mid] <= y)
                    l = mid + 1;
                else
                    h = mid - 1;
            }
            return h;
        }

        public static int CountInRange(double[] arr, int n,
                                double x, double y)
        {
            int count = UpperIndex(arr, n, y) - LowerIndex(arr, n, x) + 1;
            return count;
        }

    }
}
