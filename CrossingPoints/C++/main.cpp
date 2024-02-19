#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <fstream>
#include "Iir.h"

using namespace std;


class CrossingPointsAnalyser {
private:

    Eigen::ArrayXd separations;
    Eigen::ArrayXd crossingPoints;
    Eigen::ArrayXd Y;
    Eigen::ArrayXd X;
    double const samplingFrequency, reference;

    static vector<Eigen::ArrayXd> readDataFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("could not open file: " + filename);
        }

        vector<vector<double>> columns;
        string line;

        if (getline(file, line)) {
            istringstream iss(line);
            string token;
            while (getline(iss, token, ' ')) {
                columns.emplace_back();
            }
        } else {
            throw runtime_error("file is empty: " + filename);
        }

        file.clear();
        file.seekg(0);

        while (getline(file, line)) {
            istringstream iss(line);
            string token;
            size_t columnIndex = 0;
            while (getline(iss, token, ' ')) {
                if (columnIndex < columns.size()) {
                    columns[columnIndex].push_back(stod(token));
                    columnIndex++;
                }
            }
        }

        vector<Eigen::ArrayXd> arrays;
        for (const auto& col : columns) {
            Eigen::ArrayXd arr = Eigen::Map<const Eigen::ArrayXd>(col.data(), col.size());
            arrays.push_back(arr);
        }

        return arrays;
    }

    [[nodiscard]] vector<Eigen::ArrayXd> applyFilter(Eigen::ArrayXd x, Eigen::ArrayXd y) const{
        const int order = 2;
        const double percent = 0.05;

        Iir::Butterworth::HighPass<order> filter;
        filter.setup(samplingFrequency, 1);



        Eigen::ArrayXd filtered_y(y.size());

        for (Eigen::Index i = 0; i < y.size(); ++i){
            filtered_y(i) = filter.filter(y(i));
        }

        int lowerBound = (int)  filtered_y.size()*percent;
        int upperBound = (int) filtered_y.size()*(1-percent);




        filtered_y = filtered_y.segment(lowerBound, (upperBound -lowerBound)+1);


        x = x.segment(lowerBound, (upperBound -lowerBound)+1) - x.minCoeff();

        vector<Eigen::ArrayXd> result = {x, filtered_y};

        return result;

    }

    Eigen::ArrayXd calculateCrossingPointsSeparation(Eigen::ArrayXd x, Eigen::ArrayXd y){
        double m, c, crossingX;
        // Removed the local vector<double> crossingPointsContainer;
        crossingPoints = Eigen::ArrayXd::Zero(x.size()); // Adjust the size as necessary.
        int numCrossings = 0; // Keep track of the number of crossings.

        for (Eigen::Index i = 0; i < y.size() -1 ; i++){
            if ((y[i] <= 0 && y[i+1] >=0) || (y[i] >= 0 && y[i+1] <=0)){
                m = (y[i+1] - y[i])/(x[i+1] - x[i]);
                c = y[i] - m*x[i];
                if (m != 0){
                    crossingX = -c/m;
                    crossingPoints[numCrossings++] = crossingX; // Assign to the class member.
                }
            }
        }
        crossingPoints.conservativeResize(numCrossings); // Resize to the actual number of crossings.

        separations = (crossingPoints.tail(numCrossings - 1)
                       - crossingPoints.head(numCrossings - 1)).abs();
        return separations;
    }


    int writeXY(const string& filename) {
        ofstream outFile(filename);
        if (!outFile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return -1;
        }
        
        for (Eigen::Index i = 0; i < X.size(); ++i) {
            outFile << X[i] << ", " << Y[i] << endl;
        }
        
        outFile.close();

        cout << "Arrays have been written to " << filename << endl;

        return 0;
    }

    [[nodiscard]] static int writeOther(const string& filename, Eigen::ArrayXd data) {
        ofstream outFile(filename);
        if (!outFile.is_open()) {
            cerr << "Failed to open the file for writing." << endl;
            return -1;
        }
        
        for (Eigen::Index i = 0; i < data.size(); ++i) {
            outFile << data[i] << endl;
        }
        
        outFile.close();

        cout << "Arrays have been written to " << filename << endl;

        return 0;
    }

    static double calculateStdDev(const Eigen::ArrayXd& data) {
        double mean = data.mean(); // Compute the mean
        Eigen::ArrayXd diff = data - mean; // Compute the difference from the mean
        double sqSum = diff.square().sum(); // Sum of squared differences
        double variance = sqSum / (data.size()-1); // Variance
        return sqrt(variance); // Standard deviation is the square root of variance
    }


public:

    CrossingPointsAnalyser(const string& datafile, double samplingFrequency, double reference)
    : samplingFrequency(samplingFrequency), reference(reference) {
        vector<Eigen::ArrayXd> data = readDataFile(datafile);
        X = data[5];
        if (data[0].isZero(0)) {
            Y = data[1];
        } else {
            Y = data[0];
        }

    }

    ~CrossingPointsAnalyser()= default;

    int run(){

        vector<Eigen::ArrayXd> XY = applyFilter(X, Y);


        X = XY[0];
        Y = XY[1];

        separations = calculateCrossingPointsSeparation(X, Y);

        writeXY("filtered_xy.txt");
        writeOther("separations.txt", separations);
        writeOther("crossing_points.txt", crossingPoints);

        cout << "mean difference between crossing points is " << fixed << setprecision(0)
        << separations.mean() << "\n" <<"the standard deviation is " << calculateStdDev(separations)
         <<endl;

        cout << "distance moved per microstep is " << std::scientific << std::setprecision(5)
        << reference/(2*separations.mean()) << endl;

        return 0;
    }





};


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    string dataFileName(argv[1]);

    CrossingPointsAnalyser analyser = CrossingPointsAnalyser(dataFileName, 50, 532e-9);
    analyser.run();

    return 0;
}