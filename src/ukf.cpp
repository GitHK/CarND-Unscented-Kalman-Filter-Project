#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define NIS_RADAR_FILE "../NIS/NIS_radar.txt"
#define NIS_LIDAR_FILE "../NIS/NIS_lidar.txt"
#define KF_2PI 2.f * M_PI

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // state vector dimension
    x_length_ = 5;

    // initial state vector
    x_ = VectorXd(x_length_);

    // initial covariance matrix
    P_ = MatrixXd(x_length_, x_length_);
    P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1.5;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    /**
    TODO:
    Complete the initialization. See ukf.h for other member properties.
    Hint: one or more values initialized above might be wildly off...
    */

    lambda_ = 3 - x_length_;

    // augmented state dimension
    n_aug_length_ = x_length_ + 2;

    // sigma points dimension
    n_sig_length_ = 2 * n_aug_length_ + 1;

    // initialize weights.
    weights_ = VectorXd(n_sig_length_);
    weights_.fill(0.5 / (n_aug_length_ + lambda_));
    weights_(0) = lambda_ / (lambda_ + n_aug_length_);

    // initialize measurement noise covariance matrix
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    R_lidar_ = MatrixXd(2, 2);
    R_lidar_ << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;

    // open files for writing NIS values
    NIS_radar_file.open(NIS_RADAR_FILE);
    NIS_lidar_file.open(NIS_LIDAR_FILE);

    if (!NIS_radar_file.is_open()) {
        cout << "Error opening " << NIS_RADAR_FILE << endl;
        exit(1);
    }

    if (!NIS_lidar_file.is_open()) {
        cout << "Error opening " << NIS_LIDAR_FILE << endl;
        exit(1);
    }

}

UKF::~UKF() {
    // close files when memory gets deallocated
    NIS_radar_file.close();
    NIS_lidar_file.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
    /**
    TODO:
    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if (!is_initialized_) {
        if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
            cout << "EKF: first measurement RADAR" << endl;

            double rho = measurement_package.raw_measurements_[0];
            double phi = measurement_package.raw_measurements_[1];
            double rho_dot = measurement_package.raw_measurements_[2];
            double vx = rho_dot * cos(phi);
            double vy = rho_dot * sin(phi);

            double x = rho * cos(phi);
            double y = rho * sin(phi);
            double v = sqrt(vx * vx + vy * vy);

            // initialize with RADAR values
            x_ << x, y, v, 0.f, 0.f;
        } else {
            cout << "EKF: first measurement LIDAR" << endl;

            double x = measurement_package.raw_measurements_[0];
            double y = measurement_package.raw_measurements_[1];

            // initialize with LASER values
            x_ << x, y, 0.f, 0.f, 0.f;
        }

        // Saving first timestamp in seconds
        previous_timestamp_ = measurement_package.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }


    double delta_t = (measurement_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_package.timestamp_;

    // Prediction
    Prediction(delta_t);

    // Update
    if (measurement_package.sensor_type_ == MeasurementPackage::LASER)
        UpdateLidar(measurement_package);
    else if (measurement_package.sensor_type_ == MeasurementPackage::RADAR)
        UpdateRadar(measurement_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO:
    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    /*
     * Generate sigma points
     * */

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_length_);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_length_, n_aug_length_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(x_length_, x_length_) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    // Creating sigma points.
    //create sigma point matrix
    int x_aug_length = x_aug.size();
    MatrixXd Xsig_aug = MatrixXd(x_aug_length, n_sig_length_);

    MatrixXd P_squared = P_aug.llt().matrixL();

    double sqrt_lamnda_plus_x_aug_length = sqrt(lambda_ + x_aug_length);

    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < x_aug_length; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt_lamnda_plus_x_aug_length * P_squared.col(i);
        Xsig_aug.col(i + 1 + x_aug_length) = x_aug - sqrt_lamnda_plus_x_aug_length * P_squared.col(i);
    }


    /*
     * Predict Sigma Points
     * */

    Xsig_pred_ = MatrixXd(x_length_, n_sig_length_);
    //predict sigma points
    for (int i = 0; i < n_sig_length_; i++) {
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.0001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;


        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }


    /*
     * Predict Mean and Covariance
     * */

    x_ = Xsig_pred_ * weights_;

    P_.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++) {
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        AngleNormalization(x_diff, 3);
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_package) {
    /**
    TODO:
    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the lidar NIS.
    */

    /*
     * Predict
     * */

    int n_z = 2;
    MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_length_);

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++)
        z_pred = z_pred + weights_(i) * Zsig.col(i);


    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    S = S + R_lidar_;

    /*
     * Update
     * */

    VectorXd z = measurement_package.raw_measurements_;

    MatrixXd Tc = MatrixXd(x_length_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    //NIS Lidar Update
    double NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
    NIS_lidar_file << NIS_lidar_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_package) {
    /**
    TODO:
    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the radar NIS.
    */

    /*
     * Predict
     */

    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, n_sig_length_);
    //transform sigma points into measurement space
    for (int i = 0; i < n_sig_length_; i++) {
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        double v_x = cos(yaw) * v;
        double v_y = sin(yaw) * v;

        double r = sqrt(p_x * p_x + p_y * p_y);
        double phi = atan2(p_y, p_x);
        double r_dot = (p_x * v_x + p_y * v_y) / sqrt(p_x * p_x + p_y * p_y);

        Zsig(0, i) = r;
        Zsig(1, i) = phi;
        Zsig(2, i) = r_dot;
    }

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++)
        z_pred = z_pred + weights_(i) * Zsig.col(i);

    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        AngleNormalization(z_diff, 1);
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    S = S + R_radar_;

    /*
     * Update
     */

    VectorXd z = measurement_package.raw_measurements_;
    MatrixXd Tc = MatrixXd(x_length_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < n_sig_length_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        AngleNormalization(z_diff, 1);
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        AngleNormalization(x_diff, 3);
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;
    AngleNormalization(z_diff, 1);


    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    //NIS Update
    double NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
    NIS_radar_file << NIS_radar_ << endl;
}


void UKF::AngleNormalization(VectorXd y, int i) {
    while (y(i) > M_PI) y(i) -= KF_2PI;
    while (y(i) < -M_PI) y(i) += KF_2PI;
}
