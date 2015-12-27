using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskManagement.SecondProject;

namespace TaskManagement
{
    class Program
    {
        static void Main(string[] args)
        {
            TaskThree three = new TaskThree();
            three.evaluate();
            Console.ReadKey();
        }
    }
}
