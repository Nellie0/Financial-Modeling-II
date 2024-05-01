using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

// Implement antithetic variance reduction and Van der Corput sequences
// Antithetic should exhibit the following features:
    // Allow the user to apply antithetic, or not
    // Apply antithetic when calculating option price
    // Apply antithetic when calculation option Greeks
    // Produce a correct standard error under both antithetic and non-antithetic circumstances
// The Van der Corput functionality will be special. 
    // It need only exist in one dimension (one step), and it may use a fixed number of simulations
    // The user must be able to enter the base values for the simulation. 
    // Option price and Greeks - but not standard error - needs to be returned to the user.

namespace MyApp
{
    // Class of Normalized randoms that is inheriting from System's Random class
    public class NormalRandom : Random
    {
        // RetainedRandoms will be the array of normalized randoms
        public double[] RetainedRandoms { get; set; }
        
        // Box-Muller method to return a normalized random number
        // This takes some doubles x1 and x2 so that we can either choose to use
        // NextDouble() or VanDerCorput doubles
        public double NextGaussDouble(double x1, double x2)
        {
            double z1 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Cos(2 * Math.PI * x2);
            double z2 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Sin(2 * Math.PI * x2);
            return z1;
        }
        // Generates a certain amount of normal random numbers
        // The amount will be the number of simulations
        public double[] GenerateRetainedGaussDoubles(int quantity)
        {
            RetainedRandoms = new double[quantity];
            double x1 = NextDouble();
            double x2 = NextDouble();

            for(int  i = 0; i < quantity; i++)
            {
                RetainedRandoms[i] = NextGaussDouble(x1, x2);
            }
            return RetainedRandoms;
        }
    }
    public class VanDerCorput
    {
        // Base to be used in Van der Corput sequence
        public int BaseValue { get; set; }
        // Starting index
        public int Index { get; set; }
        // Checks if the user wants to use Van der Corput or not
        public bool IsVanDerCorput { get; set; }
        // The array of quasi randoms
        public double[] QuasiRandoms { get; set; }
        NormalRandom rnd = new NormalRandom();

        // Generates the Van der Corput numbers
        public double GenerateVanDerCorput(int BaseValue, int Index)
            {
                // Initialize variable
                double result = 0.0;
                double f = 1.0 / BaseValue;

                while (Index > 0)
                {
                    result += (Index % 2) * f;
                    Index /= BaseValue;
                    f /= BaseValue;
                }

                return result;
            }
        // Puts the Van der Corput numbers into Box Muller
        public double[] VanDerCorputSequence(int BaseValue, int Index, int quantity)
        {
            double[] QuasiRandoms = new double[quantity];
            for (int i = 0; i < quantity; i++)
            {
                QuasiRandoms[i] = rnd.NextGaussDouble(GenerateVanDerCorput(BaseValue, i + 1), GenerateVanDerCorput(BaseValue + 1, i + 1));
            }
            return QuasiRandoms;
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
            return "Price: " + Price + "\nDelta: " + Delta + "\nGamma: " + Gamma + "\nVega: " + Vega + "\nTheta: " + Theta + "\nRho: " + Rho + "\nStandard Error: " + StandardError;
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
    public class AntitheticReduction : SimulationParameters
    {
        // The user will decide if the simulation applies antithetic or not
        public bool IsAntithetic { get; set; }
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
        // 0.05 is returned for simplicity
        public double GetRateForTenor(double TenorInYears)
        {
            return 0.05;
        }
    }
    // This class is the actual simulator
    public static class Simulator
    {
        // This function is of class EvaluationResult because we're going to fill all of its variables using this simulator
        // Variables are the option, simulation parameters, yield curve, and VanDerCorput
        public static EvaluationResult Evaluate(Option o, SimulationParameters p, YieldCurve c, VanDerCorput van)
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
            // Grabbing random normal numbers
            NormalRandom rnd = new NormalRandom();

            // The Monte Carlo simulator
            // Variables are S, vol, T, r, and StandardError
            // S, vol, T, and r will be changed later when determining greeks
            // StandardError is referenced
            double MonteCarlo(double S, double vol, double T, double r, double[] Randoms, out double StandardError)
            {
                // Initialize variables
                double Payoff; // Payoff function (call of put dependent)
                double[] PayoffArray = new double[p.Simulations]; // Array of payoffs of size simulations
                double AveragePayoff = 0.0; // AveragePayoff will be the option price
                double FuturePrice = S; // Future price is an estimated time price one step in the future
                double dt = T / p.Steps; // time over the number of time steps
                double SquaredDifference = 0; // squared difference used for standard error later

                // Initiailizing payoff function
                double PayoffFunction(double FuturePrice)
                {
                    if (((European)o).IsCall == true)
                    {
                    // Call Payoff
                    Payoff = Math.Max(0, FuturePrice- K);
                    }
                    else if (((European)o).IsCall == false)
                    {
                    // Put Payoff
                    Payoff = Math.Max(0, K - FuturePrice);
                    }
                    else
                    {
                        throw new Exception("Something went wrong with checking if the option is a call or put.");
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
                
                // if the user wants to use antithetic
                if (((AntitheticReduction)p).IsAntithetic == true)
                {
                    double SumPayoffs = 0.0;
                    double SquaredPayoffs = 0.0;

                    for (int i = 0; i < p.Simulations; i+=2)
                    {
                        // antithetic uses two payoffs
                        double Payoff1 = PayoffArray[i];
                        double Payoff2 = PayoffArray[i+1];

                        double AntitheticEst = 0.5 * (Payoff1 + Payoff2);
                        SumPayoffs += AntitheticEst;
                        SquaredPayoffs += AntitheticEst * AntitheticEst;
                    }

                    double OptionPrice = Math.Exp(-r * T) * (SumPayoffs / (p.Simulations / 2)); 

                    // N is the number of pairs
                    int N = p.Simulations / 2;
                    // Calculate the standard error
                    AveragePayoff = SumPayoffs / (2 * N);
                    // definition of variance
                    double Variance = SquaredPayoffs / N - (AveragePayoff * AveragePayoff);
                    // taking the square root of the variance to produce the standard error
                    StandardError = Math.Sqrt(Variance) / Math.Sqrt(N);

                    return OptionPrice;

                }
                // if the user does not want to use antithetic
                else if (((AntitheticReduction)p).IsAntithetic == false)
                {
                    // This returns our monte carlo option price (non-antithetic)

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

                    return AveragePayoff;
                }
                // Lets us know something went wrong
                else 
                {
                    StandardError = 0.0;
                    return 0.0;
                }
            }
            // Greeks and errors
            // These functions will automatically use antithetic or non-antithetic
            // depending on what the user chose
            // since they call the above monte carlo function

            // Delta is the change in stock price
            double CalculateDelta(double[] Randoms, double epsilon = 0.01)
            {
                // This uses the default underlying price
                double PV = MonteCarlo(S, vol, T, r, Randoms, out _);
                // This uses a small chance in the underlying price
                double PV_up = MonteCarlo(S * (1 + epsilon), vol, T, r, Randoms, out _);
                return (PV_up - PV) / (2 * S * epsilon); // Using formulas provided in class
            }
            // Gamma is the change in delta
            double CalculateGamma(double[] Randoms, double epsilon = 0.01)
            {
                // Can also probably write this in terms of delta values
                // Normal value
                double PV = MonteCarlo(S, vol, T, r, Randoms, out _);
                // Small movement up with the stock price
                double PV_up = MonteCarlo(S * (1 + epsilon), vol, T, r, Randoms, out _);
                // Small movement down
                double PV_down = MonteCarlo(S * (1 - epsilon), vol, T, r, Randoms, out _);
                return (PV_up - 2 * PV + PV_down) / Math.Pow(S * (1 + epsilon), 2);
            }
            // Vega is the change in volatility
            double CalculateVega(double[] Randoms, double epsilon = 0.01)
            {
                // Small movement up in volatility
                double PV_up = MonteCarlo(S, vol + epsilon, T, r, Randoms, out _);
                // Small movement down in volatility
                double PV_down = MonteCarlo(S, vol - epsilon, T, r, Randoms, out _);
                return (PV_up - PV_down) / (2 * (vol + epsilon));
            }
            // Theta is the change in price over time
            double CalculateTheta(double[] Randoms, double epsilon = 0.01)
            {
                // Normal value
                double PV = MonteCarlo(S, vol, T, r, Randoms, out _);
                // Small positive change in the time
                double PV_up = MonteCarlo(S, vol, T + epsilon, r, Randoms, out _);
                return (PV_up - PV) / (T + epsilon);
            }
            // Rho is the change in rate over time
            double CalculateRho(double[] Randoms, double epsilon = 0.01)
            {
                // Small positive change in the risk-free rate
                double PV_up = MonteCarlo(S, vol, T, r + epsilon, Randoms, out _);
                // Small negative change in the risk-free rate
                double PV_down = MonteCarlo(S, vol, T, r - epsilon, Randoms, out _);
                return (PV_up - PV_down) / (2 * (r + epsilon));
            }
            
            // returning results
            // EvaluationResult contains all of our desired values that we'll be printing
            // So the rest of this is filling out those values from EvaluationResult
            
            // if the user chose to use Van der Corput quasi randoms
            if (van.IsVanDerCorput == true)
            {
                EvaluationResult result2 = new EvaluationResult();
                double[] QuasiRandoms = van.VanDerCorputSequence(van.BaseValue, van.Index, p.Simulations);
                double StandardError;
                result2.Price = MonteCarlo(S, vol, T, r, QuasiRandoms, out StandardError);
                result2.Delta = CalculateDelta(QuasiRandoms);
                result2.Gamma = CalculateGamma(QuasiRandoms);
                result2.Vega = CalculateVega(QuasiRandoms);
                result2.Theta = CalculateTheta(QuasiRandoms);
                result2.Rho = CalculateRho(QuasiRandoms);
                result2.StandardError = 0.0; // Van der Corput standard error is 0 for this project
                return result2;
            }
            // if the user chose to not use Van der Corput
            else if (van.IsVanDerCorput == false)
            {
                EvaluationResult result1 = new EvaluationResult();
                double[] Randoms = rnd.GenerateRetainedGaussDoubles(p.Simulations);
                double StandardError;
                result1.Price = MonteCarlo(S, vol, T, r, Randoms, out StandardError);
                result1.Delta = CalculateDelta(Randoms);
                result1.Gamma = CalculateGamma(Randoms);
                result1.Vega = CalculateVega(Randoms);
                result1.Theta = CalculateTheta(Randoms);
                result1.Rho = CalculateRho(Randoms);
                result1.StandardError = StandardError;
                return result1;
            }
            // If something went wrong, this will throw an error and end the simulation
            else
            {
                throw new Exception("Oops! Something went wrong.");
            }
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
            AntitheticReduction a = new AntitheticReduction();
            NormalRandom rnd = new NormalRandom();
            // Ask for number of steps
            Console.WriteLine("How many steps do you want for the Monte Carlo simulator?");
            var steps = Convert.ToInt32(Console.ReadLine());
            // Ask for number of simulations
            // Bug the user if they input an odd number of simulations
            Console.WriteLine("How many simulations do you want for the Monte Carlo simulator? Please enter an even number.");
            var simulations = Convert.ToInt32(Console.ReadLine());
            // We only want an even number of simulations since antithetic requires pairs of numbers
            if (simulations % 2 != 0)
            {
                throw new ArgumentOutOfRangeException(nameof(simulations), "Something went wrong. Please enter an even number.");
            }
            // Ask the user if they want to use antithetic variance reduction
            Console.WriteLine("Do you want to use antithetic variance reduction?");
            string anti = Console.ReadLine();
            // Setting IsAntithetic to true or false
            if ((anti.ToLower() == "yes") || (anti.ToLower() == "y"))
            {
                a.IsAntithetic = true;
            }
            else if ((anti.ToLower() == "no") || (anti.ToLower() == "n"))
            {
                a.IsAntithetic = false;
            }
            else
            {
                throw new ArgumentOutOfRangeException(nameof(anti), "Answer yes or no for antithetic variance reduction.");
            }
            a.Randoms = rnd.GenerateRetainedGaussDoubles(simulations);
            a.Steps = steps;
            a.Simulations = simulations;

            // define your yield curve
            YieldCurve yc = new YieldCurve();

            // define parameters for Van der Corput Sequence
            VanDerCorput van = new VanDerCorput();
            // Asking the user if they want to use Van der Corput for quasi randoms
            Console.WriteLine("Do you want to use a Van der Corput sequence?");
            var vanbool = Console.ReadLine();
            // setting IsVanDerCorput true or false
            if ((vanbool.ToLower() == "yes") || (vanbool.ToLower() == "y"))
            {
                van.IsVanDerCorput = true;
                // We only need to ask for base and index values if IsVanDerCorput is true
                Console.WriteLine("Please input the base value for the Van der Corput sequence.");
                van.BaseValue = Convert.ToInt32(Console.ReadLine());
                // The base value should be a prime number
                if (van.BaseValue == 1)
                {
                    throw new Exception("The base value must be a prime number.");
                }
                for (int i = 2; i<= Math.Ceiling(Math.Sqrt(van.BaseValue)); i++)
                {
                    if (van.BaseValue == 2)
                    {
                        break;
                    }
                    if (van.BaseValue % i == 0)
                    {
                        throw new Exception("The base value must be a prime number.");
                    }
                }
                // Now ask for index value
                Console.WriteLine("Please input the index value for the Van der Corput sequence.");
                van.Index = Convert.ToInt32(Console.ReadLine());
            }
            else if ((vanbool.ToLower() == "no") || (vanbool.ToLower() == "n"))
            {
                van.IsVanDerCorput = false;
            }
            // throws an error and ends the simulation if something went wrong
            else
            {
                throw new ArgumentOutOfRangeException(nameof(vanbool), "Answer yes or no for Van der Corput sequence.");
            }

            // price your derivative (one option price)
            var result = Simulator.Evaluate(e, a, yc, van);

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
