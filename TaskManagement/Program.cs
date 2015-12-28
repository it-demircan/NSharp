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
                ThirdProject.TaskOne one = new ThirdProject.TaskOne();
                one.TestSystemDG();
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
