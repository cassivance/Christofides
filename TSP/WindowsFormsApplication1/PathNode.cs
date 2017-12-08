using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    class PathNode
    {
        public PathNode(List<City> partialPath, double[][] distanceMatrix, int numCities, int depth, int cityPosition)
        {
            _PartialPath = partialPath;
            _DistanceMatrix = distanceMatrix;
            _LowerBound = 0.0;
            _NumCities = numCities;
            _Depth = depth;
            _CityPosition = cityPosition;
        }
        public PathNode(List<City> partialPath, double[][] distanceMatrix, double lowerBound, int numCities, int depth, int cityPosition)
        {
            _PartialPath = partialPath;
            _DistanceMatrix = distanceMatrix;
            _LowerBound = lowerBound;
            _NumCities = numCities;
            _Depth = depth;
            _CityPosition = cityPosition;
        }
        private List<City> _PartialPath;
        private double[][] _DistanceMatrix;
        private double _LowerBound;
        private int _NumCities;
        private int _Depth;
        private int _CityPosition;

        public int getDepth()
        {
            return _Depth;
        }

        public int getCityPosition()
        {
            return _CityPosition;
        }
        public List<City> getPartialPath()
        {
            return _PartialPath;
        }

        public void setPartialPath(List<City> value)
        {
           _PartialPath = value;
        }

        public double[][] getDistancesMatrix()
        {
            return _DistanceMatrix;
        }
        
        public void setDistanceMatrix(double[][] value)
        { 
            _DistanceMatrix = value;
        }

        public double getLowerBound()
        {
            return _LowerBound;
        }

        public void setLowerBound(double value)
        { 
            _LowerBound = value;
        }

        public void reducedCostMatrix()
        {
            //Reduce rows
            for (int i = 0; i < _NumCities; i++)
            {
                double minValue = _DistanceMatrix[i].Where(x => x >= 0).Min();
                if (minValue != 0 && minValue != double.PositiveInfinity)
                { 
                    _LowerBound = _LowerBound + minValue;
                    for (int j = 0; j < _NumCities; j++)
                    {
                        _DistanceMatrix[i][j] = _DistanceMatrix[i][j] - minValue;
                    }
                }
            }
            //Reduce columns
            for (int j = 0; j < _NumCities; j++)
            {
                double minColumn = double.PositiveInfinity;
                for (int i = 0; i < _NumCities; i++)
                {
                    if(_DistanceMatrix[i][j] < minColumn)
                    {
                        minColumn = _DistanceMatrix[i][j];
                    }
                }
                if (minColumn != 0 && minColumn != double.PositiveInfinity)
                {
                    _LowerBound = _LowerBound + minColumn;
                    for (int i = 0; i < _NumCities; i++)
                    {
                        _DistanceMatrix[i][j] = _DistanceMatrix[i][j] - minColumn;
                    }
                }
            }
        }
    }
}
