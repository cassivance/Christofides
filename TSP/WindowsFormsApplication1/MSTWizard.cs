using System;
using Priority_Queue;
using System.Collections.Generic;
using System.Collections;

namespace TSP
{
    struct Edge {
        public int Start, End;

        public Edge(int s, int e){
            Start = s;
            End = e;
        }
    }

    class MSTWizard
    {
        const double NOT_VISITED = 1;
        const double VISITED = Double.PositiveInfinity;
        private City[] _cities;
        private double[][] _edgeWeightsMatrix;
        private double[][] _mstMatrix;
        private double[][] _oddVerticiesMatrix;
        private double[][] _visitedEdges;



        public MSTWizard(City[] cities)
        {
            _cities = cities;
            _edgeWeightsMatrix = new double[cities.Length][];
            for(int i = 0; i < cities.Length; i++)
            {
                _edgeWeightsMatrix[i] = new double[cities.Length];
            }

            _oddVerticiesMatrix = new double[cities.Length][];
            for(int i = 0; i < cities.Length; i++)
            {
                _oddVerticiesMatrix[i] = new double[cities.Length];

            }

            _visitedEdges = new double[cities.Length][];
            for (int i = 0; i < cities.Length; i++)
            {
                _visitedEdges[i] = new double[cities.Length];
            }

            _mstMatrix = new double[cities.Length][];
            for (int i = 0; i < cities.Length; i++)
            {
                _mstMatrix[i] = new double[cities.Length];
            }


        }

        public void CreateMST()
        {
            SimplePriorityQueue<Edge> priorityQueue = new SimplePriorityQueue<Edge>();
            InitializeEdgeCostMatrixFromCities();
            InitializeVisitedEdgesMatrix();


            InitializeMSTMatrix();  
            List<int> visitedVertices = new List<int>();

            // put every edge from the first row into the priority queue as a starting point
            int firstRow = 0;
            for (int col = 0; col < _cities.Length; col++)
            {
                double edgeWeight = _edgeWeightsMatrix[firstRow][col];
                if (!Double.IsPositiveInfinity(edgeWeight)) // nothing has been visited yet so we don't have to check against that
                {
                    Edge edge = new Edge(firstRow, col);
                    priorityQueue.Enqueue(edge, (float) edgeWeight);
                }
            }
            visitedVertices.Add(firstRow);

            while(priorityQueue.Count > 0 && visitedVertices.Count < _cities.Length)
            {
                Edge bestEdge = priorityQueue.Dequeue();

                bool isEdgeVisited = _visitedEdges[bestEdge.Start][bestEdge.End].Equals(VISITED);
                bool areVerticiesVisited = visitedVertices.Contains(bestEdge.Start) && visitedVertices.Contains(bestEdge.End);
                if (!isEdgeVisited && !areVerticiesVisited)
                {
                    // mark the edge as visited
                    _visitedEdges[bestEdge.Start][bestEdge.End] = VISITED;
                    _visitedEdges[bestEdge.End][bestEdge.Start] = VISITED;

                    // mark the new vertex as visited
                    visitedVertices.Add(bestEdge.End);

                    // add the edge to the mstMatrix
                    double forwardEdgeWeight = _edgeWeightsMatrix[bestEdge.Start][bestEdge.End];
                    double backwardEdgeWeight = _edgeWeightsMatrix[bestEdge.End][bestEdge.Start];
                    _mstMatrix[bestEdge.Start][bestEdge.End] = forwardEdgeWeight;
                    _mstMatrix[bestEdge.End][bestEdge.Start] = backwardEdgeWeight;

                    // add the new possible edges to the priority queue
                    int newEdgeStart = bestEdge.End;
                    for (int col = 0; col < _cities.Length; col++)
                    {
                        double newEdgeWeight = _edgeWeightsMatrix[newEdgeStart][col];
                        bool isNewEdgeVisited = _visitedEdges[newEdgeStart][col].Equals(VISITED);
                        if (!Double.IsPositiveInfinity(newEdgeWeight) && !isNewEdgeVisited)
                        {
                            Edge newEdge = new Edge(newEdgeStart, col);
                            priorityQueue.Enqueue(newEdge, (float) newEdgeWeight);
                        }
                    }
                }
            }
        }

        /**
         * initializes _edgeWeightsMatrix with the edge weights from city to city
         * Double.PositiveInfinity represents no edge...
         */
        private void InitializeEdgeCostMatrixFromCities(){
            for (int row = 0; row < _cities.Length; row++) 
            {
                for (int col = 0; col < _cities.Length; col++)
                {
                    if (row == col) {
                        _edgeWeightsMatrix[row][col] = Double.PositiveInfinity;
                        continue;
                    }

                    City cityA = _cities[row];
                    City cityB = _cities[col];
                    double cost = cityA.costToGetTo(cityB);
                    _edgeWeightsMatrix[row][col] = cost;
                }
            }


            // TEST MATRIX
            //_edgeWeightsMatrix[0] = new double[] { Double.PositiveInfinity, 5, 10, Double.PositiveInfinity };
            //_edgeWeightsMatrix[1] = new double[] { 5, Double.PositiveInfinity, 4, 11 };
            //_edgeWeightsMatrix[2] = new double[] { 10, 4, Double.PositiveInfinity, 5 };
            //_edgeWeightsMatrix[3] = new double[] { Double.PositiveInfinity, 11, 5, Double.PositiveInfinity };
        }

        /**
         * WARNING: THIS MUST BE CALLED AFTER InitializeEdgeCostMatrixFromCities HAS BEEN CALLED
         * 
         * Initialize the _visitedEdges matrix. 
         * If the edge exits, mark it with a 1. If it doesn't exist mark it with an infinity.
         */
        private void InitializeVisitedEdgesMatrix(){
            for (int row = 0; row < _cities.Length; row++)
            {
                for (int col = 0; col < _cities.Length; col++)
                {
                    if (!Double.IsPositiveInfinity(_edgeWeightsMatrix[row][col])){
                        _visitedEdges[row][col] = NOT_VISITED;
                    } else
                    {
                        _visitedEdges[row][col] = VISITED;
                    }
                }
            }


        }

        private void InitializeMSTMatrix()
        {

            for (int row = 0; row < _cities.Length; row++)
            {
                for (int col = 0; col < _cities.Length; col++)
                {
                    _mstMatrix[row][col] = Double.PositiveInfinity;
                }
            }
        }

        public double[][] GetMSTMatrix(){
            return _mstMatrix;
        }

        /*
         * Returns a matrix with the mstRows with odd degree infinitied out
         */

        public double[][] GetOddVerticiesMatrix()
        {

            double[][] oddVerticiesMatrix = new double[_cities.Length][];
            double[] infinityRow = new double[_cities.Length];

            for (int i = 0; i < _cities.Length; i++)
            {
                infinityRow[i] = Double.PositiveInfinity;
            }

            for (int row = 0; row < _cities.Length; row++)
            {
                int edgeCount = 0;
                for (int col = 0; col < _cities.Length; col++)
                {
                    double mstEdgeWeight = _mstMatrix[row][col];
                    if (!Double.IsPositiveInfinity(mstEdgeWeight))
                    {
                        edgeCount++;
                    }
                }

                if (edgeCount % 2 == 0)
                {
                    oddVerticiesMatrix[row] = infinityRow;
                }
                else
                {
                    oddVerticiesMatrix[row] = _edgeWeightsMatrix[row];
                }
            }

            return oddVerticiesMatrix;


        }


    }
}
