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

#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/math/constants/constants.hpp>

const double pi		= boost::math::constants::pi<double>();
const double eps	= 1e-8;

class MadgwickFilter
{
	public:
		MadgwickFilter(double gyro_error);
		~MadgwickFilter();
		void estimate_pose(Eigen::Vector3d acc, Eigen::Vector3d gyro);
		void reset();
		Eigen::Vector3d get_euler(){return euler;}
		Eigen::VectorXd get_quaternion(){return quaternion;};
	private:
		double beta;
		std::chrono::system_clock::time_point old_time_stamp;
		Eigen::Vector3d euler;
		Eigen::Vector4d quaternion;
};
