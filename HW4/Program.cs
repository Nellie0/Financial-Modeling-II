using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

namespace MyApp
{
    class Program
    {
        static void Main(string[] args)
        {
            // User will enter a correlation value so that we can generate pairs of joint normally distributed random values
            Console.WriteLine("Please enter a value for the correlation (between -1 and 1)");
            double rho = Convert.ToDouble(Console.ReadLine());

            while ((rho < -1) || (rho > 1)) // This checks to see if the correlation is in the appropriate interval
            {
                Console.WriteLine("Your correlation value is not between -1 and 1. Please try again");
                rho = Convert.ToDouble(Console.ReadLine());
            }
            
            // Call each function and print each of their respective values
            // Each function (except sum_twelve) returns a tuple. We can call each individual item from the tuples to look cleaner in the console.
            double s1 = sum_twelve();
            double s2 = sum_twelve();
            Console.WriteLine($"Sum Twelve values: {s1}, {s2}");
            var b = box_muller();
            Console.WriteLine($"Box-Muller values: {b.Item1}, {b.Item2}");
            var p = polar_rejection();
            Console.WriteLine($"Polar Rejection values: {p.Item1}, {p.Item2}");
            var j = joint(rho);
            Console.WriteLine($"Correlated normal values: {j.Item1}, {j.Item2}");
        }
        
        // Sum Twelve method
        public static double sum_twelve()
        {
            Random rnd = new Random(); // Random class
            List<double> rand = new List<double>(); // Define a list of doubles
            for (int i=1; i<=12; i++)
            {
                rand.Add(rnd.NextDouble()); // Add 12 random numbers to our list
            }
            double sum = rand.Sum(); // Sum the numbers in the list
            sum = sum - 6; // Subtract the sum by six

            return sum;
        }

        // Box-Muller method
        public static (double, double) box_muller()
        {
            Random rnd = new Random(); // Random class
            // Define two random uniform variables. Use the Box-Muller equation to find z1 and z2 which are normally distributed
            double x1 = rnd.NextDouble();
            double x2 = rnd.NextDouble();
            double z1 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Cos(2 * Math.PI * x2);
            double z2 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Sin(2 * Math.PI * x2);

            return (z1,z2);
        }

        // Polar rejection method
        public static (double, double) polar_rejection()
        {
            Random rnd = new Random(); // Random class
            double x1, x2, w; // Initialize variables
            do // Loop creates x1 and x2 uniform variables with meana 0, then finds w. This repeats until w<=1
            {
                // NextDouble() is in the interval [0,1). We want to double this initial interval to [0,2), then subtract by 1 so that it's an interval of [-1,1). On a uniform distribution, this gives a mean of 0 which is what we want (everything will be centered on 0). Then we can actually calculate the Gaussian numbers.
                // This was ultimately inferred from calculating 1000 values, putting them into an Excel sheet, reading it into the test.ipynb file, and looking at histograms and a shapiro p-test.
                x1 = 2 * rnd.NextDouble() - 1;
                x2 = 2 * rnd.NextDouble() - 1;
                w = x1*x1 + x2*x2;
            } while(w > 1);

            // Find c, then multiply x1 and x2 by c. Return z1 and z2
            double c = Math.Sqrt((-2 * Math.Log(w)) / w);
            double z1 = c * x1;
            double z2 = c * x2;

            return (z1, z2);
            
        }
        public static (double, double) joint(double rho)
        {
            // I'll use normal variables from the Box-Muller method we created
            var b = box_muller();
            double x1 = b.Item1; // Calls the first variable from Box-Muller
            double x2 = b.Item2; // Calls the second variable
            double z1 = x1; //Defining z1 and z2, which will be correlated normal variables
            double z2 = rho * x1 + Math.Sqrt(1 - Math.Pow(rho, 2)) * x2;

            return (z1, z2);

        }
    }
}
