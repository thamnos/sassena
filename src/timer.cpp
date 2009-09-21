/*
 *  timer.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "timer.hpp"

#include <deque>

#include "log.hpp"
#include "parameters.hpp"

using namespace std;

double Timer::t_diff(timeval start, timeval end) {
	int seconds = (end.tv_sec-start.tv_sec);
	int microseconds = (end.tv_usec-start.tv_usec);
	seconds = (microseconds < 0) ? seconds-1 : seconds;
	microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
	return seconds + (microseconds/1000000.0);
}

double Timer::t_diff(Timer_timeval start, timeval end) {
	int seconds = (end.tv_sec-start.tv_sec);
	int microseconds = (end.tv_usec-start.tv_usec);
	seconds = (microseconds < 0) ? seconds-1 : seconds;
	microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
	return seconds + (microseconds/1000000.0);
}

void Timer::start(std::string tk)  { 
	states[tk]=true;

	timeval t;
	gettimeofday(&t, 0);
	starttimes[tk] = t;
	
	if (Params::Inst()->debug.timer) Info::Inst()->write(string("Starting timer for <")+tk+string(">"));
}

void Timer::stop(std::string tk)   { 
	if (states.find(tk)==states.end()) return; 
	if (states[tk]) { 
		states[tk]=false;
		timeval t;
		gettimeofday(&t, 0);
//		times[tk](t_diff(starttimes[tk],t));
		times[tk].push_back(t_diff(starttimes[tk],t));
	}
	else {
		cerr << "ERROR>> " << " Timer not started, but stopped: " << tk << endl;
		throw;
	}
	if (Params::Inst()->debug.timer) Info::Inst()->write(string("Stopping timer for <")+tk+string(">"));
}

void Timer::clear() {
	times.clear();
	starttimes.clear();
	states.clear();
}

Timer& Timer::operator+=(Timer& other) {
	map<string,vector<double> >::iterator i = other.times.begin();
	for(; i != other.times.end(); i++)
	{
		for(vector<double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
		{
			this->times[i->first].push_back(*j);			
		}
	}
	return *this;
}

Timer Timer::operator+( Timer& other) {
	Timer result;
	result = *this;
	result += other;
	return result;
}

vector<string> Timer::keys() {
	vector<string> ret;
//	for (map<string,times_type>::iterator ti=times.begin();ti!=times.end();ti++) {
	for (map<string,vector<double> >::iterator ti=times.begin();ti!=times.end();ti++) {
		ret.push_back(ti->first);
	}
	return ret;
}

double Timer::mean(std::string tk) {
	double sum=0.0;
	vector<double>& ts = times[tk];
	for(size_t i = 0; i < ts.size(); ++i)
	{
		sum+=ts[i];
	}
	return sum/ts.size();
//	return boost::accumulators::mean(times[tk]);
}

double Timer::variance(std::string tk) {
	double sum=0;
	double sum2=0;
	vector<double>& ts = times[tk];
	for(size_t i = 0; i < ts.size(); ++i)
	{
		sum+=ts[i];
		sum2+=ts[i]*ts[i];
	}
	return ((sum2/ts.size()) - (sum/ts.size())*(sum/ts.size()) );
//	return boost::accumulators::variance(times[tk]);
}

double Timer::sum(std::string tk) {
	double sum=0;
	vector<double>& ts = times[tk];
	for(size_t i = 0; i < ts.size(); ++i)
	{
		sum+=ts[i];
	}
	return sum;
//	return boost::accumulators::sum(times[tk]);
}
double Timer::min(std::string tk) {
	vector<double>& ts = times[tk];
	if (ts.size()<1) return 0;

	double min=ts[0];
	
	for(size_t i = 1; i < ts.size(); ++i)
	{
		if (ts[i]<min) min = ts[i];
	}	
	return min;
//	return boost::accumulators::min(times[tk]);
}
double Timer::max(std::string tk) {
	vector<double>& ts = times[tk];
	if (ts.size()<1) return 0;

	double max=ts[0];
	
	for(size_t i = 1; i < ts.size(); ++i)
	{
		if (ts[i]>max) max = ts[i];
	}	
	return max;
//	return boost::accumulators::max(times[tk]);
}

double Timer::count(std::string tk) {
	return times[tk].size();
//	return boost::accumulators::count(times[tk]);
}

// end of file