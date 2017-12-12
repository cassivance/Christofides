using System;
namespace TSP
{
    class MSTWizard
    {
        private City[] _cities;
        private double[][] _edgesMatrix;
        private double[][] _oddVerteciesMatrix;

        public MSTWizard(City[] cities)
        {
            _cities = cities;
            _edgesMatrix = new double[cities.Length][];
            for(int i = 0; i < cities.Length; i++)
            {
                _edgesMatrix[i] = new double[cities.Length];
            }

            _oddVerteciesMatrix = new double[cities.Length][];
            for(int i = 0; i < cities.Length; i++)
            {
                _oddVerteciesMatrix[i] = new double[cities.Length];
            }

        }

        public void CreateMST()
        {

        }

        public double[][] GetEdgesMatrix(){
            return _edgesMatrix;
        }

        public double[][] GetOddVerticiesMatrix(){
            return _oddVerteciesMatrix;
        }


    }
}
