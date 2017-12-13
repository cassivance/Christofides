using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;
using Priority_Queue;
using System.Linq;

namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        //creates distance matrix based on the cities provided
        public double[][] makeDistanceMatrix()
        {
            double[][] distanceMatrix = new double[Cities.Length][];
            for (int i = 0;  i < Cities.Length; i++)
            {
                distanceMatrix[i] = new double[Cities.Length];
                for (int j = 0; j < Cities.Length; j++)
                {
                    double cost = Cities[i].costToGetTo(Cities[j]);
                    if (cost != 0)
                    {
                        distanceMatrix[i][j] = cost;
                    }
                    else
                    {
                        distanceMatrix[i][j] = double.PositiveInfinity;
                    }
                }
            }
            return distanceMatrix;
        }

        //copy distance matrix to avoid passing by reference
        public double[][] copyDistanceMatrix(double[][] currentDistanceMatrix)
        {
            double[][] distanceMatrix = new double[Cities.Length][];
            for (int i = 0; i < Cities.Length; i++)
            {
                distanceMatrix[i] = new double[Cities.Length];
                for (int j = 0; j < Cities.Length; j++)
                {
                    distanceMatrix[i][j] = currentDistanceMatrix[i][j];
                }
            }
            return distanceMatrix;
        }

        //Updates distance matrix when going from city A to city B
        public double[][] goToNode(double[][] distanceMatrix, int cityA, int cityB)
        {
            for (int i = 0; i < Cities.Length; i++)
            {
                distanceMatrix[cityA][i] = double.PositiveInfinity;
                distanceMatrix[i][cityB] = double.PositiveInfinity;
            }
            distanceMatrix[cityB][cityA] = double.PositiveInfinity;
            return distanceMatrix;
        }

        //calculates priority based on depth and lower bound
        public float getPriority(PathNode currentNode)
        {
            return (float)(currentNode.getLowerBound() / (currentNode.getDepth() + 1));
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];
            Stopwatch timer = new Stopwatch();
            List<City> pathSoFar = new List<City>();
            SimplePriorityQueue<PathNode> priorityQueue = new SimplePriorityQueue<PathNode>();
            int numStates = 0;
            int depth = 1, cityIndex = 0;
            int totalStates = 0;
            int BSSFUpdates = 0;
            int prunedStates = 0; 

            timer.Start();
            //set initial BSSF by finding random path
            defaultSolveProblem();

            //Create first state
            double[][] distanceMatrix = makeDistanceMatrix();
            pathSoFar.Add(Cities[0]);
            double[][] startDistanceMatrix = copyDistanceMatrix(distanceMatrix);
            PathNode parentNode = new PathNode(pathSoFar, startDistanceMatrix, Cities.Length, depth, cityIndex);
            numStates++;
            totalStates++;
            parentNode.reducedCostMatrix();
            priorityQueue.Enqueue(parentNode, getPriority(parentNode));

            //expand nodes
            while (priorityQueue.Count != 0)
                //&& timer.ElapsedMilliseconds < (TIME_LIMIT * 1000))
            {
                PathNode currentNode = priorityQueue.Dequeue();
                if(currentNode.getDepth() != Cities.Length && currentNode.getLowerBound() < costOfBssf()) {
                    //explore each path in the current node
                    depth = currentNode.getDepth() + 1;
                    for (int i = 1; i < Cities.Length; i++)
                    {
                        distanceMatrix = currentNode.getDistancesMatrix();
                        double costOfMoving = distanceMatrix[currentNode.getCityPosition()][i];
                        if (costOfMoving < double.PositiveInfinity)
                        {
                            //make new PathNode
                            List<City> newPathSoFar = new List<City>();
                            for(int j = 0; j < currentNode.getPartialPath().Count; j++)
                            {
                                newPathSoFar.Add(currentNode.getPartialPath()[j]);
                            }
                            newPathSoFar.Add(Cities[i]);
                            double lowerBound = currentNode.getLowerBound() + costOfMoving;
                            double[][] newDistanceMatrix = copyDistanceMatrix(distanceMatrix);
                            newDistanceMatrix = goToNode(newDistanceMatrix, currentNode.getCityPosition(), i);
                            PathNode childNode = new PathNode(newPathSoFar, newDistanceMatrix, lowerBound, Cities.Length, depth, i);
                            totalStates++;
                            childNode.reducedCostMatrix();
                            lowerBound = childNode.getLowerBound();

                            //If this path could potentially be better than best path, add it to the queue
                            if (lowerBound < costOfBssf())
                            {
                                priorityQueue.Enqueue(childNode, getPriority(currentNode));
                                if (priorityQueue.Count > numStates)
                                {
                                    numStates = priorityQueue.Count;
                                }
                            }
                            else
                            {
                                prunedStates++;
                            }
                        }
                    }
                }
                else
                {
                    if(currentNode.getLowerBound() < costOfBssf() && currentNode.getPartialPath().Count == Cities.Length)
                    {
                        Route.Clear();
                        List<City> currentPath = currentNode.getPartialPath();
                        for (int i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                        {
                            Route.Add(currentPath[i]);
                        }
                        bssf = new TSPSolution(Route);
                        BSSFUpdates++;
                    }
                    else
                    {
                        prunedStates++;
                    }
                }
            }
            
            timer.Stop();

            Console.WriteLine("BSSF UPDATES: " + BSSFUpdates);
            Console.WriteLine("STATES AT A GIVEN TIME: " + numStates);
            Console.WriteLine("TOTAL STATES: " + totalStates);
            Console.WriteLine("PRUNED STATES: " + prunedStates);
            results[COST] = costOfBssf().ToString();    // load results into array here
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = numStates.ToString();

            return results;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            Stopwatch timer = new Stopwatch();
            timer.Start();


            ArrayList bestRoute = new ArrayList(); // the best route we have found so far
            ArrayList route = new ArrayList(); // represents to temp route
            List<int> visitedCities = new List<int>();
            double cost = 0;
            double bestCost = Double.PositiveInfinity;
            int bestCostStartIndex = -1;
            int solutionCount = 0;


            // Try every starting point for the route
            for (int beginIndex = 0; beginIndex < Cities.Length; beginIndex++) {

                int sourceIndex = beginIndex; // This will be updated every time a new city is added to 'route' 
                bool noGreedyPathExists = false;
                route.Clear();
                visitedCities.Clear();
                cost = 0;

                while (visitedCities.Count != Cities.Length && !noGreedyPathExists)
                {
                    // look for the nearest city to the sourceIndex that hasn't been visited and isn't infinity
                    City cityA = Cities[sourceIndex];
                    double shortestSubRoute = Double.PositiveInfinity; // represents the shortest path from source City to any neihbor
                    int shortestDestIndex = -1;

                    for (int destIndex = 0; destIndex < Cities.Length; destIndex++)
                    {
                        // Ignore the city if it is itself or already visited
                        if (destIndex == sourceIndex || visitedCities.Contains(destIndex))
                        {
                            continue;
                        }

                        City cityB = Cities[destIndex];
                        double path = cityA.costToGetTo(cityB);
                        if (path < shortestSubRoute)
                        {
                            shortestSubRoute = path;
                            shortestDestIndex = destIndex;
                        }
                    }

                    if (!Double.IsInfinity(shortestSubRoute))
                    {
                        route.Add(Cities[shortestDestIndex]);
                        visitedCities.Add(shortestDestIndex);
                        cost += shortestSubRoute;
                        sourceIndex = shortestDestIndex;
                    } else 
                    {
                        noGreedyPathExists = true;
                    }
                }

                // Check if we can complete the route by going back to the source node
                if (visitedCities.Count == Cities.Length)
                {
                    City startCity = Cities[beginIndex];
                    City lastCity = Cities[sourceIndex];
                    double lastSubRouteCost = lastCity.costToGetTo(startCity);

                    if (!Double.IsPositiveInfinity(lastSubRouteCost))
                    {
                        cost += lastSubRouteCost;
                        solutionCount++;

                        if (cost < bestCost)
                        {
                            bestCost = cost;
                            bestCostStartIndex = beginIndex;
                            bestRoute = route;
                        }
                    }
                }
            }

            timer.Stop();

            string resultsCost = "";

            if (bestCostStartIndex != -1) {
                bssf = new TSPSolution(bestRoute);
                resultsCost = costOfBssf().ToString();
            } else {
                resultsCost = "No Solution";
            }

            results[COST] = resultsCost;
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = solutionCount.ToString();

            return results;
        }

        public void nodesMatched(int node1, int node2, double[][] oddVerticiesMatrix)
        {
            for (int i = 0; i < Cities.Count(); i++)
            {
                oddVerticiesMatrix[node1][i] = double.PositiveInfinity;
                oddVerticiesMatrix[node2][i] = double.PositiveInfinity;
                oddVerticiesMatrix[i][node1] = double.PositiveInfinity;
                oddVerticiesMatrix[i][node2] = double.PositiveInfinity;
            }
        }

        //Returns a pretty good perfect match set from the odd edge vertices
        public Dictionary<int, ArrayList> getPerfectMatches(double[][] oddVerticiesMatrix)
        {
            Dictionary<int, ArrayList> matches = new Dictionary<int, ArrayList>();
            int i;
            bool added = true;
            while (added == true)
            {
                bool changed = true;
                while (changed == true)
                {
                    changed = false;
                    for (i = 0; i < Cities.Count(); i++)
                    {
                        int numTraversed = oddVerticiesMatrix[i].Where(x => x == double.PositiveInfinity).Count();
                        if (numTraversed == (oddVerticiesMatrix.Length - 1))
                        {
                            ArrayList destinationNodes = new ArrayList();
                            int destinationIndex = Array.IndexOf(oddVerticiesMatrix[i], oddVerticiesMatrix[i].Min());
                            destinationNodes.Add(destinationIndex);
                            matches.Add(i, destinationNodes);

                            ArrayList reversePath = new ArrayList();
                            reversePath.Add(i);
                            matches.Add(destinationIndex, reversePath);

                            nodesMatched(i, destinationIndex, oddVerticiesMatrix);
                            changed = true;
                        }
                    }
                }

                added = false;
                i = 0;
                while (added == false && i < Cities.Count())
                {
                    double min = oddVerticiesMatrix[i].Min();
                    if (min < double.PositiveInfinity)
                    {
                        ArrayList destinationNodes = new ArrayList();
                        int destinationIndex = Array.IndexOf(oddVerticiesMatrix[i], oddVerticiesMatrix[i].Min());
                        destinationNodes.Add(destinationIndex);
                        matches.Add(i, destinationNodes);

                        ArrayList reversePath = new ArrayList();
                        reversePath.Add(i);
                        matches.Add(destinationIndex, reversePath);
                        
                        nodesMatched(i, destinationIndex, oddVerticiesMatrix);
                        added = true;
                    }
                    i++;
                }
            }
            return matches;
        }

        public Dictionary<int, ArrayList> combineMSTandOdds(double[][] mstEdgeMatrix, Dictionary<int, ArrayList> matches)
        {
            for (int i = 0; i < Cities.Count(); i++)
            {
                for (int j = 0; j < Cities.Count(); j++)
                {
                    if (mstEdgeMatrix[i][j] != double.PositiveInfinity)
                    {
                        if (matches.ContainsKey(i))
                        {
                            if (!matches[i].Contains((object)j))
                            {
                                matches[i].Add(j);
                            }
                        }
                        else
                        {
                            ArrayList destinationNodes = new ArrayList();
                            destinationNodes.Add(j);
                            matches.Add(i, destinationNodes);
                        }
                    }
                }
            }
            return matches;
        }



        public ArrayList findPath(ArrayList curElement, Dictionary<int, ArrayList> matches) {
            ArrayList finalPath = new ArrayList();
            while (curElement.Count > 0) {
                finalPath.Add(curElement[0]);
                object temp = curElement[0];
                curElement.RemoveAt(0);
                curElement = matches[(int)temp];
            }
            return finalPath;
        }

        public ArrayList turnHamiltonian(ArrayList finalPath)
        {
            ArrayList hamiltonian = new ArrayList();
            for (int i = 0; i < finalPath.Count; i++)
            {
                if (!hamiltonian.Contains(finalPath[i]))
                {
                    hamiltonian.Add(finalPath[i]);
                }
            }
            return hamiltonian;
        }

        public ArrayList turnHamiltonianReverse(ArrayList finalPath)
        {
            ArrayList hamiltonian = new ArrayList();
            for (int i = finalPath.Count - 1; i > 0; i--)
            {
                if (!hamiltonian.Contains(finalPath[i]))
                {
                    hamiltonian.Add(finalPath[i]);
                }
            }
            return hamiltonian;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];
            Stopwatch timer = new Stopwatch();

            timer.Start();

            MSTWizard mstWizard = new MSTWizard(Cities);
            mstWizard.CreateMST();


            /*
             * This is a matrix of values and Infinities. mstEdgeMatrix[1][4] would be the weight at
             * the 1st row and the 4th column, which represents the weight for an edge from city 1 to 4.
             * If the value is not infinity, there is an edge there. If the value is Double.PositiveInfinity, there is no edge there. 
             */
            double[][] mstMatrix = mstWizard.GetMSTMatrix();  

            /*
             * Similar to the mstEdgeMatrix in its format. If there is Double.PositiveInfinity in a single matrix,
             * then there is no edge there. If the whole row is Double.PositiveInfinity, then that vertex (row) is not of odd degree.
             * Only vertecies with an odd number of connected edges should have values in their respective row. 
             */
            double[][] oddVerticiesMatrix = mstWizard.GetOddVerticiesMatrix();

            Dictionary<int, ArrayList> matches = getPerfectMatches(oddVerticiesMatrix);

            matches = combineMSTandOdds(mstMatrix, matches);



            ArrayList curElement = matches[0];
            ArrayList finalPath = findPath(curElement, matches);
            bool finished = false;
            while (finished == false) {
                ArrayList possibleStarts = new ArrayList();
                for (int i = 0; i < finalPath.Count; i++)
                {
                    if (matches[(int)finalPath[i]].Count != 0)
                    {
                        possibleStarts.Add((int)finalPath[i]);
                    }
                }

                if (possibleStarts.Count > 0)
                {
                    ArrayList tempPath = new ArrayList();
                    for (int i = 0; i < finalPath.Count; i++)
                    {
                        int startIndex = finalPath.IndexOf(possibleStarts[0]) + 1;
                        tempPath.Add(finalPath[(i + startIndex) % finalPath.Count]);
                    }
                    finalPath = tempPath;

                    ArrayList addToFinal = findPath(matches[(int)possibleStarts[0]], matches);
                    for (int i = 0; i < addToFinal.Count; i++)
                    {
                        finalPath.Add(addToFinal[i]);
                    }
                }

                else
                {
                    finished = true;
                }
              }



            ArrayList tempList = new ArrayList();
            tempList = turnHamiltonian(finalPath);

            ArrayList cityList = new ArrayList();
            for (int i = 0; i < tempList.Count; i++)
            {
                cityList.Add(Cities[(int)tempList[i]]);
            }
            bssf = new TSPSolution(cityList);
            double firstBssfCost = costOfBssf();

            ArrayList secondTempList = new ArrayList();

            secondTempList = turnHamiltonianReverse(finalPath);

            cityList.Clear();
            for (int i = 0; i < tempList.Count; i++)
            {
                cityList.Add(Cities[(int)tempList[i]]);
            }
            bssf = new TSPSolution(cityList);
            double secondBssfCost = costOfBssf();

            if (firstBssfCost < secondBssfCost)
            {
                for (int i = 0; i < tempList.Count; i++)
                {
                    cityList.Add(Cities[(int)tempList[i]]);
                }
                bssf = new TSPSolution(cityList);
            }


            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = "null";

            return results;
        }
        #endregion
    }

}