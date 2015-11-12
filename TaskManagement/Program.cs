using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskManagement.FirstProjekt;

namespace TaskManagement
{
    class Program
    {
        static void Main(string[] args)
        {
            TaskTwo two = new TaskTwo();
            two.evaluate();

            TaskThree three = new TaskThree();
            three.evaluate();

            TaskFour tt = new TaskFour();
            tt.evaluate();

            Console.ReadKey();
        }
    }
}
