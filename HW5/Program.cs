using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

// Calculate the price of European calls and puts using Monte Carlo
// Calculate the Greeks associated with the calls and puts
// Calculate the standard error of the simulation
// Permit the user a choice in number of steps and number of simulations

namespace MyApp
{
    // Class of Normalized randoms that is inheriting from System's Random class
    public class NormalRandom : Random
    {
        // RetainedRandoms will be the array of normalized randoms
        public double[] RetainedRandoms { get; set; }
        
        // Box-Muller method to return a normalized random number
        public double NextGaussDouble()
        {
            double x1 = NextDouble();
            double x2 = NextDouble();
            double z1 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Cos(2 * Math.PI * x2);
            double z2 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Sin(2 * Math.PI * x2);
            return z1;
        }
        // Generates a certain amount of normal random numbers
        // The amount will be the number of simulations
        public double[] GenerateRetainedGaussDoubles(int quantity)
        {
            RetainedRandoms = new double[quantity];

            for(int  i = 0; i < quantity; i++)
            {
                RetainedRandoms[i] = NextGaussDouble();
            }
            return RetainedRandoms;
        }
    }
    // This class designates each of our variables that will later be printed in the future
    public class EvaluationResult
    {
        // Option price
        public double Price { get; set; }
        // Greeks of the option
        public double Delta { get; set; }
        public double Gamma { get; set; }
        public double Vega { get; set; }
        public double Theta { get; set; }
        public double Rho { get; set; }
        // Standard error of the monte carlo simulator
        public double StandardError { get; set; }

        // Overrides ToString to print the above variables
        public override string ToString()
        {
            return "Price:" + Price + "\nDelta:" + Delta + "\nGamma:" + Gamma + "\nVega:" + Vega + "\nTheta:" + Theta + "\nRho:" + Rho + "\nStandard Error:" + StandardError;
        }
    }
    // Parameters specifically for the simulation
    public class SimulationParameters
    {
        // This is our array of randoms
        public double[] Randoms { get; set; }
        // The number of time steps (prompted by user)
        public int Steps { get; set; }
        // The number of monte carlo simulations (prompted by user)
        public int Simulations { get; set; }
    }
    // Parameters of an option
    public abstract class Option
    {
        // The underlying is a parameter of the option, of class underlying
        public Underlying Underlying { get; set; }
        // The expiration date is a parameter of the option, of class datetime
        public DateTime ExpirationDate { get; set; }
        // Function to get the price (to be overridden)
        public abstract double GetPrice();
        // Function to get delta (to be overridden)
        public abstract double GetDelta();
        // Function to get gamma (to be overridden)
        public abstract double GetGamma();
        // Function to get vega (to be overridden)
        public abstract double GetVega();
        // Function to get theta (to be overridden)
        public abstract double GetTheta();
        // Function to get rho (to be overridden)
        public abstract double GetRho();
        // Function to get the standard error (to be overridden)
        public abstract double GetStandardError();
    }
    // Class for specifically European options (which this project considers)
    public class European : Option
    {
        // Strike is of European options
        public double Strike { get; set; }
        // A bool for whether the European option is a call or not
        public bool IsCall { get; set; }
        // Instance of EvaluationResult. Has variables Price, Delta, Gamma, and so on.
        EvaluationResult e = new EvaluationResult();
        // Each of these functions fills out the variables that we labeled in EvaluationResult.
        // We are overriding the functions from class Option
        public override double GetPrice()
        {
            return e.Price;
        }
        public override double GetDelta()
        {
            return e.Delta;
        }
        public override double GetGamma()
        {
            return e.Gamma;
        }
        public override double GetVega()
        {
            return e.Vega;
        }
        public override double GetTheta()
        {
            return e.Theta;
        }
        public override double GetRho()
        {
            return e.Rho;
        }
        public override double GetStandardError()
        {
            return e.StandardError;
        }
    }
    // Volatility changes over time, hence the name VolatilitySurface
    public class VolatilitySurface
    {
        // Returning volatility = 0.5 for simplicity
        public double Volatility { get { return 0.5;} }
    }
    // Properties of the Underlying
    public class Underlying
    {
        // We only consider last price for S
        public double LastPrice { get; set; }
        // Volatility Surface is a property of the underlying
        public VolatilitySurface Surface { get; set; }
    }
    // There is usually a yield curve, which is often used for risk-free rate
    public class YieldCurve
    {
        // At a certain time tenor, there is a rate.
        // 0.5 is returned for simplicity
        public double GetRateForTenor(double TenorInYears)
        {
            return 0.05;
        }
    }
    // This class is the actual simulator
    public static class Simulator
    {
        // This function is of class EvaluationResult because we're going to fill all of its variables using this simulator
        // Variables are the option, simulation parameters, and the yield curve
        public static EvaluationResult Evaluate(Option o, SimulationParameters p, YieldCurve c)
        {
            // black scholes inputs
            // S for underlying last price
            var S = o.Underlying.LastPrice;
            // K for European option strike price
            var K = ((European)o).Strike;
            // T for expiration date in terms of days divided by the number of days in a year
            var T = (o.ExpirationDate - DateTime.Today).Days / 252d;
            // r for risk-free rate (from yield curve)
            var r = c.GetRateForTenor(T);
            // vol for volatility from our volatility surface
            var vol = o.Underlying.Surface.Volatility;

            // The Monte Carlo simulator
            // Variables are S, vol, T, r, and StandardError
            // S, vol, T, and r will be changed later when determining greeks
            // StandardError is referenced
            double MonteCarlo(double S, double vol, double T, double r, out double StandardError)
            {
                // Initialize variables
                double Payoff; // Payoff function (call of put dependent)
                double[] PayoffArray = new double[p.Simulations]; // Array of payoffs of size simulations
                double AveragePayoff = 0.0; // AveragePayoff will be the option price
                double FuturePrice = S; // Future price is an estimated time price one step in the future
                double dt = T / p.Steps; // time over the number of time steps
                double SquaredDifference = 0; // squared difference used for standard error later

                // Grabbing random normal numbers
                NormalRandom rnd = new NormalRandom();
                double[] Randoms = rnd.GenerateRetainedGaussDoubles(p.Simulations);

                // Initiailizing payoff function
                double PayoffFunction(double FuturePrice)
                {
                    if (((European)o).IsCall == true)
                    {
                    // Call Payoff
                    Payoff = Math.Max(0, FuturePrice- K);
                    }
                    else
                    {
                    // Put Payoff
                    Payoff = Math.Max(0, K - FuturePrice);
                    }
                    return Payoff;
                }

                // Loop over our number of simulations
                for (int i = 0; i < p.Simulations; i++)
                {
                    // Loop over our number of time steps
                    for (int j = 0; j < p.Steps; j++)
                    {
                        // Estimating future price using geometric brownian motion
                        FuturePrice = S * Math.Exp((r - 0.5 * vol * vol) * dt + vol * Math.Sqrt(dt) * Randoms[i]);
                    }
                    // Payoff amount is determined by our estimated future price
                    Payoff = PayoffFunction(FuturePrice);
                    // Multiply this payoff by e^(-rT) and put it into the array
                    PayoffArray[i] = Payoff * Math.Exp(-r * T);
                }
                
                // Sum of squared differences to get standard error
                foreach (double payoff in PayoffArray)
                {   
                    // Sum up our found values
                    AveragePayoff += payoff;
                }
                // Divide by number of simulations
                AveragePayoff /= p.Simulations;

                // Now calculate the standard error using sum of squared errors
                foreach (double payoff in PayoffArray)
                {
                    double deviation = payoff - AveragePayoff;
                    SquaredDifference += Math.Pow(deviation, 2);
                }
                double StandardDev = Math.Sqrt(SquaredDifference / (p.Simulations - -1));
                // Computing standard error using standard deviation
                StandardError = StandardDev /  Math.Sqrt(p.Simulations);

                return AveragePayoff; // This returns our monte carlo option price
            }            
            // Greeks and errors
            // Delta is the change in stock price
            double CalculateDelta(double epsilon = 0.01)
            {
                // This uses the default underlying price
                double PV = MonteCarlo(S, vol, T, r, out _);
                // This uses a small chance in the underlying price
                double PV_up = MonteCarlo(S * (1 + epsilon), vol, T, r, out _);
                return (PV_up - PV) / (2 * S * epsilon); // Using formulas provided in class
            }
            // Gamma is the change in delta
            double CalculateGamma(double epsilon = 0.01)
            {
                // Can also probably write this in terms of delta values
                // Normal value
                double PV = MonteCarlo(S, vol, T, r, out _);
                // Small movement up with the stock price
                double PV_up = MonteCarlo(S * (1 + epsilon), vol, T, r, out _);
                // Small movement down
                double PV_down = MonteCarlo(S * (1 - epsilon), vol, T, r, out _);
                return (PV_up - 2 * PV + PV_down) / Math.Pow(S * (1 + epsilon), 2);
            }
            // Vega is the change in volatility
            double CalculateVega(double epsilon = 0.01)
            {
                // Small movement up in volatility
                double PV_up = MonteCarlo(S, vol + epsilon, T, r, out _);
                // Small movement down in volatility
                double PV_down = MonteCarlo(S, vol - epsilon, T, r, out _);
                return (PV_up - PV_down) / (2 * (vol + epsilon));
            }
            // Theta is the change in price over time
            double CalculateTheta(double epsilon = 0.01)
            {
                // Normal value
                double PV = MonteCarlo(S, vol, T, r, out _);
                // Small positive change in the time
                double PV_up = MonteCarlo(S, vol, T + epsilon, r, out _);
                return (PV_up - PV) / (T + epsilon);
            }
            // Rho is the change in rate over time
            double CalculateRho(double epsilon = 0.01)
            {
                // Small positive change in the risk-free rate
                double PV_up = MonteCarlo(S, vol, T, r + epsilon, out _);
                // Small negative change in the risk-free rate
                double PV_down = MonteCarlo(S, vol, T, r - epsilon, out _);
                return (PV_up - PV_down) / (2 * (r + epsilon));
            }
            
            
            // returning results
            // EvaluationResult contains all of our desired values that we'll be printing
            // So the rest of this is filling out those values from EvaluationResult
            EvaluationResult result = new EvaluationResult();
            double StandardError;
            result.Price = MonteCarlo(S, vol, T, r, out StandardError);
            result.Delta = CalculateDelta();
            result.Gamma = CalculateGamma();
            result.Vega = CalculateVega();
            result.Theta = CalculateTheta();
            result.Rho = CalculateRho();
            result.StandardError = StandardError;

            return result; // returns created instance of EvaluationResult
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            // setup your vol
            VolatilitySurface v = new VolatilitySurface();

            // setup your underlying
            Underlying u = new Underlying();
            u.Surface = v;
            u.LastPrice = 50; // This will be S

            // setup your derivative
            European e = new European();
            e.Underlying = u; // from instance of our underlying class
            e.Strike = 50;
            e.ExpirationDate = new DateTime(2024, 09, 10);
            e.IsCall = true; // We're letting this be a call option

            // define simulation parameters
            SimulationParameters p = new SimulationParameters();
            NormalRandom rnd = new NormalRandom();
            // Ask for number of steps
            Console.WriteLine("How many steps do you want?");
            var steps = Convert.ToInt32(Console.ReadLine());
            // Ask for number of simulations
            Console.WriteLine("How many simulations do you want?");
            var simulations = Convert.ToInt32(Console.ReadLine());
            p.Randoms = rnd.GenerateRetainedGaussDoubles(simulations);
            p.Steps = steps;
            p.Simulations = simulations;

            // define your yield curve
            YieldCurve yc = new YieldCurve();

            // price your derivative (one option price)
            var result = Simulator.Evaluate(e, p, yc);

            // A list of option prices
            /*
            List<Option> portfolio = new List<Option>();
            foreach(var o in portfolio)
            {
                var result = Simulator.Evaluate(o, p, yc);
            }
            */

            // Prints out to the user the option price, greeks, and standard error
            Console.WriteLine(result.ToString());
        }
    }
}