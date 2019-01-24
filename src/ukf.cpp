#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#define PI 3.14159265

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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // P_ should be initialized as a symetric 
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  

  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(5, 2*n_aug_+1);

  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to CTRV, lazy initialization. we put psi, psi_dot and v to 0
      Because we have no idea how the car is oriented. QUESTION : Can we do better ? Assuming phi = psi ??
      */
      double p_x = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      double p_y = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_ << p_x, p_y, 0, 0, 0;
    } 
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0,0;
    }

    if (fabs(x_(0)) < 0.0001 and fabs(x_(1)) < 0.0001){
      x_(0) = 0.0001;
      x_(1) = 0.0001;
    }

    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /****************************************************************************
   *  Prediction                                                              *
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);

  /****************************************************************************
   *  Update                                                                  *
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }




}

MatrixXd UKF::GenerateSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug << x_,
           0,
           0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(7, 15);

  //create augmented covariance matrix
  MatrixXd Fill_down = Eigen::MatrixXd::Zero(2, P_.cols());
  MatrixXd Fill_right = Eigen::MatrixXd::Zero(P_.rows(), 2);
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_, 0,
       0 , std_yawdd_*std_yawdd_;
  P_aug << P_, Fill_right,
          Fill_down, Q;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  MatrixXd R = sqrt(lambda_ +n_aug_)*A;
  MatrixXd Xsig_1 =  R.colwise() + x_aug;
  MatrixXd Xsig_2 =  - R;
  Xsig_2 = Xsig_2.colwise() + x_aug;
  
  Xsig_aug << x_aug, Xsig_1, Xsig_2;

  return Xsig_aug;

}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /****************************************************************************
   *  Sigma points generation                                                 *
   ****************************************************************************/
  std::cout << "Prediction!!!" << std::endl;
  MatrixXd Xsig_aug = MatrixXd(7, 15);
  // MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug = GenerateSigmaPoints();


  /****************************************************************************
   *  Sigma points prediction                                                 *
   ****************************************************************************/

  //predict sigma points
  for (int i = 0; i< 15; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    // while (yaw_p> PI) {
    //   yaw_p-=2.*PI;
    //   std::cout << yaw_p << std::endl;
    // }
    // while (yaw_p<-PI) yaw_p+=2.*PI;

    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }


   /****************************************************************************
   *  Predict meand and covariance matrix from the predcition ofsigma points   *
   ****************************************************************************/
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_+ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> PI) x_diff(3)-=2.*PI;
    while (x_diff(3)<-PI) x_diff(3)+=2.*PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  MatrixXd Zsig = Xsig_pred_.block(0,0,2,2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S and Cross relation matrix
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);

  MatrixXd Tc = MatrixXd(5, 2);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    S = S + weights_(i) * z_diff * z_diff.transpose();
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(2,2);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_meas = VectorXd(2);
  z_meas << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z_meas - z_pred;

  //NIS computation
  double NIS = z_diff.transpose()*S.inverse()*z_diff;
  ofstream myfile;
  myfile.open ("NIS.txt", ios_base::app);
  myfile << NIS << endl;
  myfile.close();

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot

  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S and Cross relation matrix
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);

  MatrixXd Tc = MatrixXd(5, 3);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> PI) z_diff(1)-=2.*PI;
    while (z_diff(1)<-PI) z_diff(1)+=2.*PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> PI) x_diff(3)-=2.*PI;
    while (x_diff(3)<-PI) x_diff(3)+=2.*PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(3,3);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_meas = VectorXd(3);
  z_meas << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  VectorXd z_diff = z_meas - z_pred;

  //angle normalization
  while (z_diff(1)> PI) z_diff(1)-=2.*PI;
  while (z_diff(1)<-PI) z_diff(1)+=2.*PI;

  double NIS = z_diff.transpose()*S.inverse()*z_diff;
  ofstream myfile;
  myfile.open ("NIS.txt", ios_base::app);
  myfile << NIS << endl;
  myfile.close();

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}
