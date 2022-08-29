#include "ukf.h"
#include <iostream>

#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3.0 - n_x_;

  // create vector for weights
  weights_radar_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  weights_radar_.fill(0.5 / (n_aug_ + lambda_));
  weights_radar_(0) = lambda_ / (lambda_ + n_aug_);
  std::cout << "Weights: " << weights_radar_ << std::endl;

  x_.fill(0.0);
  P_ = MatrixXd::Identity(n_x_, n_x_);

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
    {
      x_.head(2) << meas_package.raw_measurements_;
    }
    else
    {
      x_[0] = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      x_[1] = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }
  else
  {
    double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
    {
      UpdateLidar(meas_package);
    }
    else
    {
      UpdateRadar(meas_package);
    }
  }
}

void UKF::Prediction(double delta_t)
{
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

  /* Generate Sigma Points */
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // augmented state covariance
  P_aug.fill(0.0);
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_aug_ + 1); // sigma point matrix
  Xsig.fill(0.0);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); // augmented mean vector
  Xsig_aug.fill(0.0);
  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  GenerateSigmaPoints(P_, Xsig);
  AugmentedSigmaPoints(P_, P_aug, Xsig_aug);
  SigmaPointPrediction(Xsig_aug, Xsig_pred_, delta_t);
  PredictMeanAndCovariance(Xsig_pred_, x_, P_);
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function!
   * Use lidar data to update the belief about the object's position.
   * Modify the state vector, x_, and covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred = VectorXd(n_z_lidar_);
  z_pred.fill(0.0);
  MatrixXd Zsig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1); // matrix for sigma points in measurement space
  Zsig.fill(0.0);
  MatrixXd S = MatrixXd(n_z_lidar_, n_z_lidar_);
  S.fill(0.0);

  PredictLidarMeasurement(Xsig_pred_, Zsig, z_pred, S);
  UpdateStateByLidar(Xsig_pred_, z, z_pred, Zsig, S, x_, P_);
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * Use radar data to update the belief about the object's position.
   * Modify the state vector, x_, and covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred = VectorXd(n_z_radar_);
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1); // matrix for sigma points in measurement space
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  PredictRadarMeasurement(Xsig_pred_, Zsig, z_pred, S);
  UpdateState(Xsig_pred_, z, z_pred, Zsig, S, x_, P_);
}

void UKF::GenerateSigmaPoints(MatrixXd &P_in, MatrixXd &Xsig_out)
{

  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  Xsig.fill(0.0);

  // calculate square root of P
  MatrixXd A = P_in.llt().matrixL();

  std::cout<< "GenerateSigmaPoints" << std::endl;
  std::cout << "P :" << P_in << std::endl;
  std::cout << "A :" << A << std::endl;

  // calculate sigma points ...
  Xsig.col(0) = x_;
  // set sigma points as columns of matrix Xsig
  double coeff = sqrt(lambda_ + n_x_);
  for (int i = 0; i < 5; ++i)
  {
    Xsig.col(i + 1) = x_ + coeff * A.col(i);
    Xsig.col(i + 1 + n_x_) = x_ - coeff * A.col(i);
  }
  std::cout << "Xsig : \n"
            << Xsig << std::endl;
  Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd &P_in, MatrixXd &P_aug, MatrixXd &Xsig_aug_out)
{

  // create sigma point matrix
  // create augmented mean vector
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  std::cout<< "AugmentedSigmaPoints" << std::endl;
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.segment<5>(0) = x_;
  std::cout << "x_aug: " << x_aug << std::endl;

  // create augmented covariance matrix
  P_aug.block<5, 5>(0, 0) = P_in;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  std::cout << "Xsig_aug : \n"
            << Xsig_aug << std::endl;
  std::cout << "P_aug : \n"
            << P_aug << std::endl;
  Xsig_aug_out = Xsig_aug;
  P_aug = P_aug;
}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, MatrixXd &Xsig_out, double delta_t)
{

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0);
  double dt = delta_t;

  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_acc = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yaw_acc = Xsig_aug(6, i);

    // avoid division by zero
    if (yaw_acc < std::numeric_limits<double>::epsilon())
    {
      Xsig_pred(0, i) = p_x + v * cos(yaw) * dt;
      Xsig_pred(1, i) = p_y + v * sin(yaw) * dt;
    }
    else
    {
      Xsig_pred(0, i) = p_x + v / yaw_acc * (sin(yaw + yaw_acc * dt) - sin(yaw));
      Xsig_pred(1, i) = p_y + v / yaw_acc * (-cos(yaw + yaw_acc * dt) + cos(yaw));
    }
    Xsig_pred(0, i) += 0.5 * pow(dt, 2) * cos(yaw) * nu_a;
    Xsig_pred(1, i) += 0.5 * pow(dt, 2) * sin(yaw) * nu_a;
    Xsig_pred(2, i) = v + dt * nu_a;
    Xsig_pred(3, i) = yaw + yaw_acc * dt + 0.5 * pow(dt, 2) * nu_yaw_acc;
    Xsig_pred(4, i) = yaw_acc + dt * nu_yaw_acc;
    // write predicted sigma points into right column
  }

  std::cout<< "SigmaPointPrediction" << std::endl;
  std::cout << "Xsig_out: \n"
            << Xsig_pred << std::endl;
  // write result
  Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred_in, VectorXd &x_out, MatrixXd &P_out)
{

  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  x = Xsig_pred_in * weights_radar_;

  // predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    MatrixXd A = Xsig_pred_in.col(i) - x_;
    A = A * A.transpose();
    P += weights_radar_(i) * A;
  }
  std::cout<< "PredictMeanAndCovariance" << std::endl;
  std::cout << "x_pred: \n"
            << x << std::endl;
  std::cout << "P_pred: \n"
            << P << std::endl;
  // write result
  x_out = x;
  P_out = P;
}

void UKF::PredictRadarMeasurement(MatrixXd &Xsig_pred_in,
                                  MatrixXd &Zsig_in,
                                  VectorXd &z_out,
                                  MatrixXd &S_out)
{

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // transform sigma points into measurement space
    double p_x = Xsig_pred_in(0, i);
    double p_y = Xsig_pred_in(1, i);
    double v = Xsig_pred_in(2, i);
    double yaw = Xsig_pred_in(3, i);
    double yaw_acc = Xsig_pred_in(4, i);

    double rho = sqrt(pow(p_x, 2) + pow(p_y, 2));
    double phi = atan(p_y / p_x);
    double rho_acc = (p_x * cos(yaw) * v + p_y * sin(yaw) * v) / rho;

    // calculate mean predicted measurement
    Zsig_in(0, i) = rho;
    Zsig_in(1, i) = phi;
    Zsig_in(2, i) = rho_acc;
    z_pred(0) += weights_radar_(i) * Zsig_in(0, i);
    z_pred(1) += weights_radar_(i) * Zsig_in(1, i);
    z_pred(2) += weights_radar_(i) * Zsig_in(2, i);
  }

  MatrixXd R(3, 3);
  R(0, 0) = pow(std_radr_, 2);
  R(1, 1) = pow(std_radphi_, 2);
  R(2, 2) = pow(std_radrd_, 2);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // calculate innovation covariance matrix S
    MatrixXd A = (Zsig_in.col(i) - z_pred);
    A = A * A.transpose();
    S += weights_radar_(i) * A;
  }
  S += R;

  std::cout << "PredictRadarMeasurement" << std::endl;
  std::cout << "Z_pred : \n" << z_pred << std::endl;
  std::cout << "S : \n" << S << std::endl;
  // write result
  z_out = z_pred;
  S_out = S;
}

void UKF::PredictLidarMeasurement(MatrixXd &Xsig_pred_in,
                                  MatrixXd &Zsig_in,
                                  VectorXd &z_out,
                                  MatrixXd &S_out)
{
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_lidar_);
  z_pred.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // transform sigma points into measurement space
    Zsig_in(0, i) = Xsig_pred_in(0, i);
    Zsig_in(1, i) = Xsig_pred_in(1, i);
  }

  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_radar_(i) * Zsig_in.col(i);
  }

  MatrixXd R(n_z_lidar_, n_z_lidar_);
  R(0, 0) = pow(std_laspx_, 2);
  R(1, 1) = pow(std_laspy_, 2);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_lidar_, n_z_lidar_);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // calculate innovation covariance matrix S
    MatrixXd A = (Zsig_in.col(i) - z_pred);
    A = A * A.transpose();
    S += weights_radar_(i) * A;
  }
  S += R;

  std::cout << "PredictLidarMeasurement" << std::endl;
  std::cout << "Xsig_pred: \n" << Xsig_pred_ << std::endl;
  std::cout << "Zsig : \n" << Zsig_in << std::endl;
  std::cout << "Z_pred : \n" << z_pred << std::endl;
  std::cout << "S : \n" << S << std::endl;
  // write result
  z_out = z_pred;
  S_out = S;
}

void UKF::UpdateState(MatrixXd &Xsig_pred_in,
                      VectorXd &z,
                      VectorXd &z_pred_in,
                      MatrixXd &Zsig_in,
                      MatrixXd &S_in,
                      VectorXd &x_out,
                      MatrixXd &P_out)
{

  // create example vector for predicted state mean
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  // create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  Tc.fill(0.0);
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
  S.fill(0.0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    Tc += weights_radar_(i) * ((Xsig_pred_in.col(i) - x_) * (Zsig_in.col(i) - z_pred_in).transpose());
    S = S_in + weights_radar_(i) * (Zsig_in.col(i) - z_pred_in) * (Zsig_in.col(i) - z_pred_in).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S_in.inverse();

  // update state mean and covariance matrix
  x = x_ + K * (z - z_pred_in);
  P = P_ - K * S_in * K.transpose();

  std::cout << "UpdateState" << std::endl;
  std::cout << "X : \n" << x << std::endl;
  std::cout << "P : \n" << P << std::endl;
  // write result
  x_out = x;
  P_out = P;
}

void UKF::UpdateStateByLidar(MatrixXd &Xsig_pred_in,
                             VectorXd &z,
                             VectorXd &z_pred_in,
                             MatrixXd &Zsig_in,
                             MatrixXd &S_in,
                             VectorXd &x_out,
                             MatrixXd &P_out)
{

  // create example vector for predicted state mean
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  // create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);
  Tc.fill(0.0);
  MatrixXd S = MatrixXd(n_z_lidar_, n_z_lidar_);
  S.fill(0.0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    Tc += weights_radar_(i) * ((Xsig_pred_in.col(i) - x_) * (Zsig_in.col(i) - z_pred_in).transpose());
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S_in.inverse();

  // update state mean and covariance matrix
  x = x_ + K * (z - z_pred_in);
  P = P_ - K * S_in * K.transpose();

  std::cout << "UpdateStateByLidar" << std::endl;
  std::cout << "X : \n" << x << std::endl;
  std::cout << "P : \n" << P << std::endl;
  // write result
  x_out = x;
  P_out = P;
}
