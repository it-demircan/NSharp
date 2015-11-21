using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskManagement.SecondProject;
using TaskManagement.Kunoth;

namespace TaskManagement
{
    class Program
    {
        static void Main(string[] args)
        {
            //TaskTwo tt = new TaskTwo();
            //tt.evaluate();
            ThetaSolverProject tsp = new ThetaSolverProject();
            tsp.evaluate();
            Console.ReadKey();
        }
    }
}
