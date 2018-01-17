#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;

  n_aug_ = 7;

  n_radar_ = 3;

  n_laser_ = 2;

  lambda_ = 3 - n_aug_;

  //x_aug_ = VectorXd( n_aug_ );

  //deltax_ = VectorXd( n_aug_ );

  Xsig_aug_ = MatrixXd( n_aug_, 2*n_aug_+1 );

  Xsig_pred_ = MatrixXd( n_x_, 2*n_aug_+1 );

  deltax_ = VectorXd( n_aug_ );
  x_diff = VectorXd( n_aug_ );
  //P_aug_ = MatrixXd( n_aug_, n_aug_ );

  //L_ = MatrixXd( n_aug_, n_aug_ );

  /*
  weights_ = VectorXd( 2*n_aug_+1 );
  weights_(0) = lambda_/( lambda_ + n_aug_ );
  for( int i=1; i<2*n_aug_+1; i++ )
     weights_(i) = 0.5/( n_aug_ + lambda_ );
  */
  // Variables used for radar update
  z_pred_radar_ = VectorXd( n_radar_ );
  deltaz_radar_ = VectorXd( n_radar_ );
  Zsig_radar_ = MatrixXd( n_radar_, 2*n_aug_+1 );
  // Radar measurement noise covariance matrix is constant/persistent
  R_radar_ = MatrixXd( n_radar_, n_radar_ );
  R_radar_.fill(0.);
  R_radar_(0,0) = std_radr_*std_radr_;
  R_radar_(1,1) = std_radphi_*std_radphi_;
  R_radar_(2,2) = std_radrd_*std_radrd_;
  S_radar_ = MatrixXd( n_radar_, n_radar_ );
  Tc_radar_ = MatrixXd( n_x_, n_radar_ );
  K_radar_ = MatrixXd( n_x_, n_radar_ );

  // Variables used for laser update
  z_pred_laser_ = VectorXd( n_laser_ );
  deltaz_laser_ = VectorXd( n_laser_ );
  Zsig_laser_ = MatrixXd( n_laser_, 2*n_aug_+1 );
  // Laser measurement noise covariance matrix is constant/persistent
  R_laser_ = MatrixXd( n_laser_, n_laser_ );
  R_laser_.fill(0.);
  R_laser_(0,0) = std_laspx_*std_laspx_;
  R_laser_(1,1) = std_laspy_*std_laspy_;
  S_laser_ = MatrixXd( n_laser_, n_laser_ );
  Tc_laser_ = MatrixXd( n_x_, n_laser_ );
  K_laser_ = MatrixXd( n_x_, n_laser_ );

  //set vector for weights
  weights = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  is_initialized_ = false;
  previous_timestamp_ = 0;
}

UKF::~UKF() {}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = n_x_;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P_
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x_ - sqrt(lambda_+n_x) * A.col(i);
  }

  //print result
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  //write result
  *Xsig_out = Xsig;

}

/*
 * input :
 *      x_ , state vector
 *      P_ , state covariance matrix
 * output :
 *      Xsig_aug_ , augumented sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  int used_test_data = false;
  if (used_test_data == true){
        assignedTestValues(1);
  }
  else{
        //debugInputValues(1);
  }
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_+n_aug) * L.col(i);
  }

  //print result
  std::cout << std::fixed;
  std::cout << std::endl<<"1. Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;

}

/*
 * input :
 *      Xsig_aug_ , augumented sigma points
 * output :
 *      Xsig_pred_ , predicted sigma points matrix
 */
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out,double delta_t) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

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
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //print result
  std::cout << std::endl<< "2. Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}

/*
 * input :
 *      Xsig_pred_ , predicted sigma points matrix
 * output :
 *      x_pred_ , predicted state mean
 *      P_pred_ , predicted state covariance matrix
 */
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

  //predicted state mean
  x.fill(0.0);
  VectorXd Xsig_pred_col = VectorXd( n_x_ );
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
    Xsig_pred_col = weights(i) * Xsig_pred_.col(i);
    x = x+ Xsig_pred_col;
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  std::cout << std::endl<< "3.1 Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << std::endl<< "3.2 Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

/*
 * RADAR
 *
 */

/*
 * input :
 *      Xsig_pred_ , predicted sigma points matrix
 * output :
 *      z_pred_ , mean predicted measurement
 *      S_ , predicted measurement covariance matrix
 *      Zsig_radar_ , measurement sigma point matrix
 */
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = n_radar_;

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_radar_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_radar_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred_radar = VectorXd(n_z);
  z_pred_radar.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred_radar = z_pred_radar + weights(i) * Zsig_radar_.col(i); //??? PROBLEM
  }

  //innovation covariance matrix S
  MatrixXd S_radar = MatrixXd(n_z,n_z);
  S_radar.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd deltaz_radar_ = Zsig_radar_.col(i) - z_pred_radar_;

    //angle normalization
    while (deltaz_radar_(1)> M_PI) deltaz_radar_(1)-=2.*M_PI;
    while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1)+=2.*M_PI;

    S_radar = S_radar + weights(i) * deltaz_radar_ * deltaz_radar_.transpose();
  }

  S_radar = S_radar + R_radar_;

  //print result
  std::cout << std::endl<< "4.1 z_pred_radar: " << std::endl << z_pred_radar << std::endl;
  std::cout << std::endl<< "4.2 S_radar: " << std::endl << S_radar << std::endl;

  //write result
  *z_out = z_pred_radar_;
  *S_out = S_radar;
}

void UKF::debugInputValues(int testcase)
{
  std::cout << std::endl<<"into debugInputValues" << std::endl;
  if (testcase == 1)
        std::cout <<std::endl << "x_: " << x_ << std::endl << "P_: " << P_<< std::endl <<std::endl;
  else if (testcase == 5)
        std::cout <<std::endl<< "Xsig_pred_: " << Xsig_pred_<<std::endl << "x_: " << x_ << std::endl << "P_: " << P_<< std::endl << "Zsig_radar_: " << Zsig_radar_ << std::endl<< "z_pred_radar_: "<< z_pred_radar_ << std::endl << "S_radar_" << S_radar_ <<std::endl;
  std::cout << std::endl<<"exit debugInputValues" << std::endl;
}

void UKF::assignedTestValues(int testcase)
{

  if (testcase == 1){
        //set example state
        x_ <<   5.7441,
                 1.3800,
                2.2049,
                0.5015,
                0.3528;

        //create example covariance matrix
        P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
                -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
                0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
                -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
                -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
  }
  else if (testcase == 5){
	  Xsig_pred_ <<
		 5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
		   1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
		  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
		 0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
		  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

	  x_ <<
	     5.93637,
	     1.49035,
	     2.20528,
	    0.536853,
	    0.353577;

	  P_ <<
	  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
	  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
	  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
	 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
	 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

	  Zsig_radar_ <<
	      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
	     0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
	      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

	  z_pred_radar_ <<
	      6.12155,
	     0.245993,
	      2.10313;

	  S_radar_ <<
	      0.0946171, -0.000139448,   0.00407016,
	   -0.000139448,  0.000617548, -0.000770652,
	     0.00407016, -0.000770652,    0.0180917;
	  }
}

/*
 * input :
 *      z , incoming radar measurement
 *      z_pred_ , mean predicted measurement
 *      Zsig_ , measurement sigma point matrix
 *      Xsig_pred_ , predicted sigma points matrix
 *      x_pred_ , predicted state mean
 * output :
 *      x_ , state vector
 *      P_ , state covariance matrix
 */
void UKF::UpdateRadarState(VectorXd* x_out, MatrixXd* P_out,MeasurementPackage meas_package) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = n_radar_;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  int used_test_data = false;
  if (used_test_data == true){
        assignedTestValues(5);

        z <<
        5.9214,
        0.2187,
        2.0062;
  }
  else{
        //debugInputValues(5);
  }
  std::cout << "UpdateRadarState raw measurement: " << std::endl << z << std::endl;

  //calculate cross correlation matrix
  Tc_radar_.fill(0.0);
  std::cout << std::endl<<"UpdateRadarState 0" << std::endl;
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    deltaz_radar_ = Zsig_radar_.col(i) - z_pred_radar_;
    //angle normalization
    std::cout << std::endl<<"UpdateRadarState 0.1" << std::endl;
    while (deltaz_radar_(1)> M_PI) deltaz_radar_(1)-=2.*M_PI;
    while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1)+=2.*M_PI;

    std::cout << std::endl<<"UpdateRadarState 0.2" << std::endl;
    // state difference
    x_diff = Xsig_pred_.col(i) - x_;
    std::cout << std::endl<<"UpdateRadarState 0.3" << std::endl << "x_diff : " << x_diff << std::endl;
    //angle normalization
    uint count = 0;
    uint max_count = 1000000; //PROBLEM
    while ((x_diff(3)> M_PI)) {x_diff(1)-=2.*M_PI; count++; if(count < max_count){count=0;break;}}
    while ((x_diff(3)<-M_PI)) {x_diff(1)+=2.*M_PI; count++; if(count < max_count){count=0;break;}}

    std::cout << std::endl<<"UpdateRadarState 0.4" << std::endl;
    Tc_radar_ = Tc_radar_ + weights(i) * x_diff * deltaz_radar_.transpose();
    std::cout << std::endl<<"UpdateRadarState 0.5" << std::endl;
  }

  std::cout << std::endl<<"UpdateRadarState 1" << std::endl;
  //std::cout << "Tc_radar_ : " << std::endl << Tc_radar_ << std::endl<<"S_radar_ : "<<std::endl<<S_radar_ << std::endl<<"S_radar_.inverse : "<<std::endl<<S_radar_.inverse() << std::endl;
  //Kalman gain K;
  MatrixXd K_radar_ = Tc_radar_ * S_radar_.inverse();

  //residual
  VectorXd deltaz_radar_ = z - z_pred_radar_;

  //std::cout << "z: " << std::endl << z <<std::endl<< "z_pred_radar_:" << z_pred_radar_<< std::endl << " deltaz_radar_ :" << deltaz_radar_ <<std::endl;
  //angle normalization
  while (deltaz_radar_(1)> M_PI) deltaz_radar_(1)-=2.*M_PI;
  while (deltaz_radar_(1)<-M_PI) deltaz_radar_(1)+=2.*M_PI;

  //update state mean and covariance matrix
  VectorXd x = VectorXd(n_x);
  MatrixXd P = MatrixXd(n_x,n_x);
  x = x_ + K_radar_ * deltaz_radar_;
  P = P_ - K_radar_*S_radar_*K_radar_.transpose();

  //std::cout << "x_: " << std::endl << x_ << std::endl << "K_radar_ :" << K_radar_ << std::endl << "deltaz_radar_ : " <<deltaz_radar_<<std::endl;
  //print result
  std::cout << std::endl<< "5.1 Updated state x: " << std::endl << x << std::endl;
  std::cout << std::endl<< "5.2 Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;

}

/*
 * LIDAR
 *
 */

/*
 * input :
 *      Xsig_pred_ , predicted sigma points matrix
 * output :
 *      z_pred_ , mean predicted measurement
 *      S_ , predicted measurement covariance matrix
 *      Zsig_ , measurement sigma point matrix
 */
void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = n_laser_;

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    Zsig_laser_(0,i) = Xsig_pred_(0,i);
    Zsig_laser_(1,i) = Xsig_pred_(1,i);
  }

  //mean predicted measurement
  z_pred_laser_.fill(0.);
  for (int i=0; i < 2*n_aug+1; i++) {
      //z_pred = z_pred + weights(i) * Zsig_.col(i);
      z_pred_laser_ = z_pred_laser_ + weights(i)*Zsig_laser_.col(i);
  }

  //innovation covariance matrix S
  S_laser_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    deltaz_laser_ = Zsig_laser_.col(i) - z_pred_laser_;
    S_laser_ = S_laser_ + weights(i)*deltaz_laser_*deltaz_laser_.transpose();
  }

  //add measurement noise covariance matrix
  S_laser_ = S_laser_ + R_laser_;
  //print result
  std::cout << "4.1 z_pred_laser: " << std::endl << z_pred_laser_ << std::endl;
  std::cout << "4.2 S_laser_: " << std::endl << S_laser_ << std::endl;

  //write result
  *z_out = z_pred_laser_;
  *S_out = S_laser_;
}

/*
 * input :
 *      z , incoming radar measurement
 *      z_pred_ , mean predicted measurement
 *      Zsig_ , measurement sigma point matrix
 *      Xsig_pred_ , predicted sigma points matrix
 *      x_pred_ , predicted state mean
 * output :
 *      x_ , state vector
 *      P_ , state covariance matrix
 */
void UKF::UpdateLidarState(VectorXd* x_out, MatrixXd* P_out,MeasurementPackage meas_package) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = n_laser_;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;
  std::cout << "raw measurement: " << std::endl << z << std::endl;
  //float rho = measurement_pack.raw_measurements_[0];
  //float phi = measurement_pack.raw_measurements_[1];
  ////x << rho*cos(phi), rho*sin(phi), 0.f, 0.f;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc_laser_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    deltax_ = Xsig_pred_.col(i) - x_;
    deltaz_laser_ = Zsig_laser_.col(i) - z_pred_laser_;
    Tc_laser_ = Tc_laser_ + weights(i)*deltax_*deltaz_laser_.transpose();
  }

  std::cout << "Tc_laser_ : " << std::endl << Tc_laser_ << "S_laser_ : "<<std::endl<<S_laser_ << std::endl;
  //Kalman gain K;
  MatrixXd K_laser_ = Tc_laser_ * S_laser_.inverse();

  //residual
  deltaz_laser_ = z - z_pred_laser_;

  //update state mean and covariance matrix
  VectorXd x = VectorXd(n_x);
  MatrixXd P = MatrixXd(n_x,n_x);
  x = x_ + K_laser_ * deltaz_laser_;
  P = P_ - K_laser_*S_laser_*K_laser_.transpose();

  std::cout << "x_: " << std::endl << x_ << std::endl << "K_laser_ :" << K_laser_ << std::endl << "deltaz_laser_ : " << deltaz_laser_<<std::endl;
  //print result
  std::cout << "5.1 Updated state x: " << std::endl << x << std::endl;
  std::cout << "5.2 Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this
fy the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  std::cout << "UKF::Prediction " << std::endl;

  //MatrixXd Xsig = MatrixXd(5, 11);
  //GenerateSigmaPoints(&Xsig);

  AugmentedSigmaPoints(&Xsig_aug_);
  SigmaPointPrediction(&Xsig_pred_,delta_t);
  PredictMeanAndCovariance(&x_pred_, &P_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout << "UKF UpdateLidar " << endl;
  PredictLidarMeasurement(&z_pred_laser_, &S_laser_);
  UpdateLidarState(&x_, &P_,meas_package);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  cout << "UKF UpdateRadar " << endl;
  PredictRadarMeasurement(&z_pred_radar_, &S_radar_);
  UpdateRadarState(&x_, &P_,meas_package);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  cout << endl <<"UKF ProcessMeasurement  sensor type: " << measurement_pack.sensor_type_<<endl<<endl;
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
   TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF init: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout<<"into radar measurement pack " <<std::endl;
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      x_ << rho*cos(phi), rho*sin(phi), 0.f, 0.f, 0.f;
      cout <<" Init RADAR measurement : "<< measurement_pack.raw_measurements_[0] << " "<<measurement_pack.raw_measurements_[1]<<endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
      cout <<" Init LIDAR measurement : "<< measurement_pack.raw_measurements_[0] << " "<<measurement_pack.raw_measurements_[1]<<endl;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    cout << "before init state: " << endl;
    int n_x = n_x_;
    //set augmented dimension
    int n_aug = n_aug_;
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    P_.fill(0.);
    P_(0,0) = 1.;
    P_(1,1) = 1.;
    P_(2,2) = 1.;
    P_(3,3) = 1.;
    P_(4,4) = 1.;
    //Xsig_aug_ = MatrixXd(7, 15);

    //Xsig_pred_ = MatrixXd(15, 5);

    x_pred_ = VectorXd(5);
    P_pred_ = MatrixXd(5, 5);

    Zsig_ = MatrixXd(n_z, 2 * n_aug + 1);
    z_pred_ = VectorXd(3);
    S_ = MatrixXd(3, 3);

    is_initialized_ = true;
    cout << "after init: " << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   cout << "UKF: start prediction" << endl;
	//compute the time elapsed between the current and previous measurements
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
   previous_timestamp_ = measurement_pack.timestamp_;
   cout << "UKF: prediction " << dt<< endl;
   if (dt > 0.0001 )
   {
   	Prediction(dt);
   }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //UpdateRadar( measurement_pack.raw_measurements_ );
    UpdateRadar( measurement_pack);
  } else {
    // Laser updates
    UpdateLidar(measurement_pack);
  }
}
