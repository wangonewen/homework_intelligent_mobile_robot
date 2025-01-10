#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <Eigen/Dense>


class TrajectoryGenerator {
public:
    TrajectoryGenerator() = default;
    ~TrajectoryGenerator() = default;

    Eigen::MatrixXd PolyQPGeneration(
        const int order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time);
    
    int Factorial(int x);
    Eigen::MatrixXd getM(int n_seg, int n_order, const Eigen::VectorXd& ts);
    Eigen::MatrixXd getCt(int n_seg, int n_order);
    Eigen::MatrixXd getQ(int n_seg, int n_order, const Eigen::VectorXd& ts);

private:
    double _qp_cost;
    Eigen::MatrixXd _Q;
    Eigen::VectorXd _Px, _Py, _Pz;
};

#endif // TRAJECTORY_H