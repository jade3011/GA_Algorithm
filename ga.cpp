#include "ga.h"

GA::GA()
{
	setParameter();
	setPopulation();
	x1_min_ = -3.0;
	x1_max_ = 12.1;
	x2_min_ = 4.1;
	x2_max_ = 5.8;
}

void GA::setParameter()
{
	std::cout << "input pop_size" << std::endl;
	std::cin >> pop_size_;

	std::cout << "probability crossover" << std::endl;
	std::cin >> pc_;

	std::cout << "probability mutation" << std::endl;
	std::cin >> pm_;

	std::cout << "generation" << std::endl;
	std::cin >> generation_;
}

void GA::setPopulation()
{
	std::mt19937 engine((unsigned int)time(NULL));
	std::uniform_int_distribution<int> distribution(0, 1);
	auto generator = bind(distribution, engine);

	for (int j = 0; j < pop_size_; j++)
	{
		std::vector<int> chromosome;
		for (int i = 0; i < 33; i++)
		{
			chromosome.push_back(generator());
		}
		population_.push_back(chromosome);

	}

}

int GA::calcDecimal(const std::vector<int> chromosome, int start_n, int finsh_n)
{
	int decimal_value = 0;
	int degree = 0;
	for (int i = finsh_n-1; i >= start_n; i--)
	{
		decimal_value = chromosome.at(i)*pow(2, degree) + decimal_value;
		degree++;
	}
	return decimal_value;
}

void GA::calcFitness(std::vector<double> & fitness)
{
	double x1;
	double x2;
	double fitness_value;

	for (int i = 0; i < population_.size(); i++)
	{
		x1 = x1_min_ + calcDecimal(population_.at(i),0,18) * ((x1_max_ - x1_min_) / (pow(2, 18) - 1));
		x2 = x2_min_ + calcDecimal(population_.at(i), 18, 33) * ((x2_max_ - x2_min_) / (pow(2, 15) - 1));
		
		fitness_value = 21.5 + x1 * sin(4 * M_PI * x1) + x2 * sin(20 * M_PI * x2);
		fitness.push_back(fitness_value);
	}
}

double GA::calcTotalFitness(const std::vector<double>& fitness)
{
	double total_fitness = 0.0;
	for (int i = 0; i < fitness.size(); i++)
	{
		total_fitness = fitness.at(i) + total_fitness;
	}
	return total_fitness;
}

void GA::calcProbability(const std::vector<double> & fitness, std::vector<double>& probability_fitness)
{
	double total_fitness;
	total_fitness = calcTotalFitness(fitness);
	for (int i = 0; i < fitness.size(); i++)
	{
		double probability;
		probability = fitness.at(i) / total_fitness;
		probability_fitness.push_back(probability);
	}
}

void GA::calcCumulative(const std::vector<double>& probability_fitness, std::vector<double>cumulative_probability)
{
	double cumulative = 0.0;
	for (int i = 0; i < probability_fitness.size(); i++)
	{
		cumulative = probability_fitness.at(i) + cumulative;
		cumulative_probability.push_back(cumulative);
	}
}

template<typename T>
T GA::randGenerate(T start_range, T finsh_range)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	if constexpr (std::is_same_v<T, int>)
	{
		std::uniform_int_distribution<int>dist(start_range, finsh_range);
		int random_num = dist(gen);
		return random_num;
	}
	else
	{
		std::uniform_real_distribution<T> dist(start_range, finsh_range);
		T random_num = dist(gen);
		return random_num;
	}


}

void GA::selectChromosome(const std::vector<double>& rand2select,const std::vector<double> cumulative_probability)
{
	std::vector<std::vector<int>> temp_population;
	for (int i = 0; i < population_.size(); i++)
	{
		temp_population.push_back(population_.at(i));
	}

	for (int i = 0; i < rand2select.size(); i++)
	{
		for (int j = 0; j < cumulative_probability.size(); j++)
		{
			if (rand2select.at(i) < cumulative_probability.at(j))
			{
				population_.at(i) = temp_population.at(j);
				break;
			}
		}
		
	}
}

void GA::setReproduction()
{
	std::vector<double> fitness;
	std::vector<double> probability_fitness;
	std::vector<double> cumulative_probability;
	std::vector<double> rand2select;

	calcFitness(fitness);
	calcProbability(fitness, probability_fitness);
	calcCumulative(probability_fitness,cumulative_probability);
	for (int i = 0; i < pop_size_; i++)
	{
		rand2select.push_back(randGenerate(0.0, 1.0));
	}
	selectChromosome(rand2select, cumulative_probability);
}

void GA::paringChromosome(const std::vector<double>& rand2select, std::vector<std::pair<int, int>>& pair)
{
	std::vector<int> select_pair;// 확률이 맞아서 선택된 갯수가 홀수인지 짝수인지 판단해서 홀수일경우 맨앞에것 하나 더 넣기위해 사용
	
	for (int i = 0; i < rand2select.size(); i++)
	{
		if (rand2select.at(i) < pc_)
		{
			select_pair.push_back(i);
		}
	}
	if (select_pair.size() % 2 == 0) //짝수일 때
	{
		if (select_pair.size() == 0)
		{ }
		else
		{
			for (int i = 0; i < select_pair.size() - 1; i = i + 2)
			{
				pair.push_back(std::make_pair(select_pair.at(i), select_pair.at(i + 1)));
			}
		}
	}
	else //홀수일 때
	{
		if (select_pair.size() == 1)
		{ }
		else
		{
			for (int i = 1; i < select_pair.size() - 1; i = i + 2)
			{
				pair.push_back(std::make_pair(select_pair.at(i-1), select_pair.at(i)));
			}
			pair.push_back(std::make_pair(select_pair.at(select_pair.size() - 1), select_pair.at(0)));
		}
	}
}

void GA::cutReplace(const std::vector<int>& rand2pair, const std::vector<std::pair<int,int>>&pair)
{
	int value1, value2;
	std::vector<std::vector<int>> temp_population;
	for (int i = 0; i < population_.size(); i++)
	{
		temp_population.push_back(population_.at(i));
	}

	for (int i = 0; i < rand2pair.size(); i++)
	{
		value1 = pair.at(i).first;
		value2 = pair.at(i).second;

		for (int j = rand2pair.at(i); j < 33; j++)
		{
			population_.at(value1).at(j) = population_.at(value2).at(j);
			population_.at(value2).at(j) = temp_population.at(value1).at(i);
		}

	}
}

void GA::setCrossover()
{
	std::vector<double> rand2select_crossover;
	std::vector<std::pair<int, int>>pair;
	std::vector<int>rand2pair;
	for (int i = 0; i < pop_size_; i++)
	{
		rand2select_crossover.push_back(randGenerate(0.0, 1.0));
	}
	paringChromosome(rand2select_crossover, pair);
	for (int i = 0; i < pair.size(); i++)
	{
		rand2pair.push_back(randGenerate(1, 32));
	}
	cutReplace(rand2pair,pair);
}

void GA::setMuation()
{
	std::vector<double> rand2mutation;
	std::vector<int> select_mutation;
	for (int i = 0; i < pop_size_*33; i++)
	{
		rand2mutation.push_back(randGenerate(0.0, 1.0));
	}

	for (int i = 0; i < rand2mutation.size(); i++)
	{
		if (rand2mutation.at(i) < pm_)
		{
			select_mutation.push_back(i);
		}
	}

	for (int i = 0; i < select_mutation.size(); i++)
	{
		int population_num;
		int chromosome_num;
		population_num = select_mutation.at(i) / 33;
		chromosome_num = select_mutation.at(i) % 33;
		if (population_.at(population_num).at(chromosome_num) == 0)
		{
			population_.at(population_num).at(chromosome_num) = 1;
		}
		else
		{
			population_.at(population_num).at(chromosome_num) = 0;
		}
	}
	
}

void GA::print()
{
	std::vector<double> fitness;
	double total_fitness;
	double max_fitness = DBL_MIN;
	double mean_fitness;
	calcFitness(fitness);
	total_fitness = calcTotalFitness(fitness);

	for (int i = 0; i < fitness.size(); i++)
	{
		if (max_fitness < fitness.at(i))
		{
			max_fitness = fitness.at(i);
		}
	}
	mean_fitness = total_fitness / pop_size_;
	std::cout << "max fitness" << max_fitness << std::endl << "mean fitness" << mean_fitness << std::endl;
	best_fitness_.push_back(max_fitness);
	std::fstream fp;
    fp.open("fitness.csv",std::ios::app);

   if (fp.is_open())
   {
       fp << std::to_string(max_fitness) << "," <<std::to_string(mean_fitness) << "\n"; 
	   fp.close();
   }
}

void GA::process()
{
	double best = DBL_MIN;
	for (int i = 0; i < generation_; i++)
	{
		setReproduction();
		setCrossover();
		setMuation();
		print();
	}
	for (int i = 0; i < best_fitness_.size(); i++)
	{
		if (best < best_fitness_.at(i))
		{
			best = best_fitness_.at(i);
		}
	}
	std::cout << "best" <<best<< std::endl;
}