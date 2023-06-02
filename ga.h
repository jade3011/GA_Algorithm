#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <functional>
#include <cmath>
#include <fstream>
#include<string>


class GA
{
public:
	int pop_size_;
	int generation_;
	double pc_;
	double pm_;

	double x1_min_;
	double x1_max_;
	double x2_min_;
	double x2_max_;

	std::vector<std::vector<int>> population_;
	std::vector<double> best_fitness_;

	GA();
	void setParameter();
	void setPopulation();
	int calcDecimal(const std::vector<int> chromosome, int start_n, int finsh_n);

	void calcFitness(std::vector<double> & fitness);
	double calcTotalFitness(const std::vector<double>& fitness);
	void calcProbability(const std::vector<double>& fitness, std::vector<double> & probability_fitness);
	void calcCumulative(const std::vector<double>& probability_fitness, std::vector<double>cumulative_probability);
	void selectChromosome(const std::vector<double>& rand2select, const std::vector<double> cumulative_probability);
	void setReproduction();


	void paringChromosome(const std::vector<double> &rand2select, std::vector<std::pair<int, int>> &pair);
	void cutReplace(const std::vector<int> & rand2pair, const std::vector<std::pair<int, int>>&pair);
	void setCrossover();

	void setMuation();

	void print();

	void process();

	template<typename T>
	T randGenerate(T start_range, T finsh_range);
};