/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	num_particles = 201;

	 std_x = std[0];
	 std_y = std[1];
	 std_theta = std[2];


// distributioon for x,y and theta

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);


//Creating a number of particle

	for (int i = 0; i < num_particles; ++i) {
            Particle p;
            p.x = dist_x(gen);
            p.y = dist_y(gen);
            p.theta = dist_theta(gen);
            p.weight = 1.0;
            p.id = i;
            weights.push_back(p.weight);
            particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
    // generate random Gaussian noise
    normal_distribution<double> N_x(0, std_pos[0]);
    normal_distribution<double> N_y(0, std_pos[1]);
    normal_distribution<double> N_theta(0, std_pos[2]);
    for(auto& p: particles){
        if( fabs(yaw_rate) < 0.0001){
            //if velocity is too close to zero
            p.x += velocity * delta_t * cos(p.theta);
            p.y += velocity * delta_t * sin(p.theta);
        } else{
            p.x += velocity / yaw_rate * ( sin( p.theta + yaw_rate*delta_t ) - sin(p.theta) );
            p.y += velocity / yaw_rate * ( cos( p.theta ) - cos( p.theta + yaw_rate*delta_t ) );
            p.theta += yaw_rate * delta_t;
        }
    // predicted particles with added sensor noise
        p.x += N_x(gen);
        p.y += N_y(gen);
        p.theta += N_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	double closest_distance, distance;
	for(int i = 0; i < observations.size(); i++) {
        //for every observation in observations
		LandmarkObs observation = observations[i];
		for(int j = 0; j < predicted.size(); j++){
		    //for every prediction in predictions
			LandmarkObs prediction = predicted[j];
            distance = dist(prediction.x,prediction.y,observation.x,observation.y);
			if(j==0 || distance < closest_distance){
			    //for initial or if the distance is lesser
				closest_distance = distance;
				observations[i].id = j;
    // map observations with the predicted landmark
			}
		}
	//cout<<observations[i].x << ',' << observations[i].y << ',' << predicted[observations[i].id].x << ',' << predicted[observations[i].id].y << endl;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	default_random_engine gen;
	for (int i = 0; i< num_particles;i++){
            double pred_x = particles[i].x;
            double pred_y = particles[i].y;
            double pred_theta = particles[i].theta;
            vector<LandmarkObs> transformed_obs;
            for (int j=0;j<observations.size();j++){
                double x = observations[j].x*cos(pred_theta) - observations[j].y *sin(pred_theta) + pred_x;
                double y = observations[j].y*cos(pred_theta) + observations[j].x *sin(pred_theta) + pred_y;
                //Transform the observation into the global coordinate.
                transformed_obs.push_back(LandmarkObs{ observations[j].id, x, y });
            }
            vector<LandmarkObs> predictions;
            for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
                float lm_x = map_landmarks.landmark_list[j].x_f;
                float lm_y = map_landmarks.landmark_list[j].y_f;
                int lm_id = map_landmarks.landmark_list[j].id_i;
                double distance2 = dist(lm_x,lm_y,pred_x,pred_y);
                if (distance2 <= sensor_range){
                        //Predictioned landmarks which are close to the car
                    predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
                }
            }
            //Mapping the data and associated value
            dataAssociation(predictions,transformed_obs);

            particles[i].weight = 1.0;
            for(auto& obs_m: transformed_obs){
                auto landmark = predictions[obs_m.id];
                double x_term = pow(obs_m.x - landmark.x, 2) / (2 * pow(std_landmark[0], 2));
                double y_term = pow(obs_m.y - landmark.y, 2) / (2 * pow(std_landmark[1], 2));
                double w = exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
                particles[i].weight *=  w;
                //update weights
    }
    weights[i] = particles[i].weight;
}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
	vector<Particle> resampled;
	discrete_distribution<int> dist_w(weights.begin(), weights.end());
	for (auto &particle : particles){
		resampled.push_back(particles[dist_w(gen)]);
	}
	particles = resampled;


}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

