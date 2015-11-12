using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Structures
{
    public class GeneralHelper
    {
        public const double EPSILON = 2.2204460492503131E-16;
        /// <summary>
        /// Vergleicht zwei fließkomma Zahlen bis auf Maschinengenaugikeit.
        /// </summary>
        /// <param name="x">Erste zu vergleichende Zahl</param>
        /// <param name="y">Zweite zu vergleichende Zahl</param>
        public static bool isXAlmostEqualToY(double x, double y)
        {
            
            bool almostEqual = false;
            if (x == 0 || y == 0)
            {
                if (Math.Abs(x - y) <= 2 * EPSILON)
                    almostEqual = true;
            }
            else
            {
                if (Math.Abs(x - y) <= EPSILON * Math.Abs(x) && Math.Abs(x - y) <= EPSILON * Math.Abs(y))
                    almostEqual = true;
            }
            return almostEqual;
        }


        public static void WriteOutputText(String path, String text)
        {
            if (File.Exists(path))
                File.Delete(path);

            using (TextWriter writer = File.CreateText(path))
            {
                writer.Write(text);
            }
        }
    }
}
