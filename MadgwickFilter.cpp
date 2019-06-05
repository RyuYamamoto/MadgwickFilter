/*
 * Software License Agreement (BSD License)
 *
 * Copyright (c) 2019, Ryu Yamamoto.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "MadgwickFilter.h"

MadgwickFilter::MadgwickFilter(double gyro_error)
{
	beta = std::sqrt(3/4) * (pi * (gyro_error/180.0));

	euler = Eigen::Vector3d::Zero();
	quaternion = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0);
}

MadgwickFilter::~MadgwickFilter()
{

}

void MadgwickFilter::reset()
{
	euler = Eigen::Vector3d::Zero();
	quaternion = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0);
}

void MadgwickFilter::estimate_pose(Eigen::Vector3d acc, Eigen::Vector3d gyro)
{
	try
	{
		std::chrono::system_clock::time_point time_stamp = std::chrono::system_clock::now(); // get time stamp

		Eigen::Vector4d d_quaternion;
		d_quaternion[0] = 0.5*(-quaternion[1]*gyro.x()-quaternion[2]*gyro.y()-quaternion[3]*gyro.z());
		d_quaternion[1] = 0.5*( quaternion[0]*gyro.x()+quaternion[2]*gyro.z()+quaternion[3]*gyro.y());
		d_quaternion[2] = 0.5*( quaternion[0]*gyro.y()-quaternion[1]*gyro.z()+quaternion[3]*gyro.x());
		d_quaternion[3] = 0.5*( quaternion[0]*gyro.z()+quaternion[1]*gyro.y()-quaternion[2]*gyro.x());

		acc = acc.normalized();

		Eigen::MatrixXd J; J.resize(3,4);
		J << -2*quaternion[2], 2*quaternion[3], -2*quaternion[0], 2*quaternion[1],
				  2*quaternion[1], 2*quaternion[0],  2*quaternion[3], 2*quaternion[2],
					0, -4*quaternion[1], -4*quaternion[2], 0;

		Eigen::Vector3d f;
		f.x() = 2*(quaternion[1]*quaternion[3] - quaternion[0]*quaternion[2]) - acc.x();
		f.y() = 2*(quaternion[0]*quaternion[1] - quaternion[2]*quaternion[3]) - acc.y();
		f.z() = 2*(0.5 - std::pow(quaternion[1], 2) - std::pow(quaternion[2], 2)) - acc.z();

		Eigen::VectorXd f_nable; f_nable.resize(4,1);
		f_nable = J.transpose() * f;
		f_nable = f_nable.normalized();

		double delta = std::chrono::duration_cast<std::chrono::milliseconds>(time_stamp-old_time_stamp).count() / 1000.0;
		quaternion += (d_quaternion-(beta*f_nable)) * delta;
		quaternion = quaternion.normalized();

		euler.x() = std::atan2(2*(quaternion[0]*quaternion[1] + quaternion[2]*quaternion[3]), -2*(std::pow(quaternion[1], 2)) + std::pow(quaternion[2], 2) + 1);
		euler.y() = -std::asin(-2*(quaternion[1]*quaternion[3] + quaternion[0]*quaternion[2]));
		euler.z() = std::atan2(2*(quaternion[0]*quaternion[3] + quaternion[1]*quaternion[2]), -2*(std::pow(quaternion[2], 2)) + std::pow(quaternion[0], 2) + 1);

		old_time_stamp = time_stamp;
	}
	catch(std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}

