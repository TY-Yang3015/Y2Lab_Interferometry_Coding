#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <fstream>

using namespace std;


class InterferogramSimulator {
private:

    double N_SIGMA;
    int N_MODE;

    Eigen::ArrayXd x, y;



    [[nodiscard]] Eigen::ArrayXd gaussianSpectrum(const double& intensity) const{
        Eigen::ArrayXd spectrum = Eigen::ArrayXd::LinSpaced(N_MODE, -N_SIGMA, N_SIGMA);
        spectrum = (-spectrum.square()/2.).exp();
        spectrum /= spectrum.sum();
        spectrum *= intensity;

        return spectrum;
    }


public:

    InterferogramSimulator(double N_SIGMA, double N_MODE, Eigen::ArrayXd x, Eigen::ArrayXd y)
        : N_SIGMA(N_SIGMA), N_MODE(N_MODE), x(x), y(y){
    }

    ~InterferogramSimulator()= default;

    static std::pair<Eigen::ArrayXd, Eigen::ArrayXd> auto_mesh(const double& sampling_frequency=50, const double& motor_speed=30000,
                                                               const double& start_position=-1e7, const double& conversion_factor=1e-11) {
        double sampling_distance = motor_speed / sampling_frequency;
        int n_sample = static_cast<int>(2. * -start_position / sampling_distance);

        Eigen::ArrayXd x = Eigen::ArrayXd::LinSpaced(n_sample, start_position, -start_position) * 2. * conversion_factor;
        Eigen::ArrayXd y = Eigen::ArrayXd::Zero(n_sample);

        return make_pair(x, y);
    }





    Eigen::ArrayXd addGaussian(const float& wavelength, const float& sigma
                                , const float& intensity) {
        Eigen::ArrayXd wavelengths = Eigen::ArrayXd::LinSpaced(N_MODE, -N_SIGMA*sigma, N_SIGMA*sigma) + wavelength;
        Eigen::ArrayXd spectrum = gaussianSpectrum(intensity);

        Eigen::ArrayXd spatial = Eigen::ArrayXd::Zero(x.size());
        for (Eigen::Index i=0; i < wavelengths.size(); i++) {
            spatial += spectrum[i] * ((2. * M_PI * x / wavelengths[i]).sin()) + spectrum[i];
        }
        y += spatial;

        return spatial;
    }






    Eigen::ArrayXd addSquare(const double& start, const double& width, const double& intensity) {
        Eigen::ArrayXd wavelengths = Eigen::ArrayXd::LinSpaced(N_MODE, start, start+width);
        double spectrum = intensity/N_MODE;

        Eigen::ArrayXd spatial = Eigen::ArrayXd::Zero(x.size());
        for (Eigen::Index i = 0; i < wavelengths.size(); i++){
            spatial += spectrum * ((2. * M_PI * x / wavelengths[i]).sin()) + spectrum;
        }

        y += spatial;

        return spatial;
    }

    int writeTxt(const string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return -1;
        }

        // Assuming x and y are of the same size
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            outFile << x(i) << ", " << y(i) << std::endl;
        }

        // Close the file
        outFile.close();

        std::cout << "Arrays have been written to " << filename << std::endl;

        return 0;
    }



};

int main() {
    auto [x, y] = InterferogramSimulator::auto_mesh();
    InterferogramSimulator simulator(5, 50, x, y);

    //simulator.addGaussian(550e-9f, 0.1e-9f, 1.f);
    //simulator.addGaussian(560e-9f, 0.1e-9f, 1.f);

    simulator.addSquare(550e-9f,10e-9f, 1.f);
    simulator.addSquare(560e-9f, 10e-9f, 1.f);

    simulator.writeTxt("./cpp_simulation_result.txt");

    return 0;
}


