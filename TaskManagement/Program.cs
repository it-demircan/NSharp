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
            try
            {
                FourthProject.TaskTwo.Test();
                Console.Write("Calculated");
            }
            catch (Exception err)
            {
                Console.Write(err.StackTrace);
            }
            
            Console.ReadKey();
        }
    }
}
