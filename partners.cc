/*

  For each agent on each iteration there are four events which should be
  executed in this order:

  breakupPartnership:

     if num_partners > 0:
         if uniform rng < (relationship_stickiness_attribute / num_partners)
             break up last formed partnership

  formNewPartnership:

     if uniform rng  < partner_forming_attribute /  (num_partners + 1)
        find a new matching partner and stick at back of partners queue

     MORE COMPLEX VERSION

     if num_partners == 0:
         uniform rng < partner_forming_attribute
     else:
         uniform rng < concurrency_attribute

  haseSex:
     if num_partners > 0 and uniform rng < sexual_drive_attribute:
         partner = min(geometric_distribution(preference_fifs_attribute),
                                              num_partners)
	 determine hiv transmission risk

  So we have the following attributes per agent:

  - relationship_stickiness_attribute: higher implies stays in relationships
                                       longer.
    Initialize to 1 - beta distribution(ALPHA_STICKINESS, BETA_STICKINESS)

  - partner_forming_attribute: higher value implies more likely to form
                               new relationship if no partners

    Initialize to beta distribution(ALPHA_PARTNER_FORM, BETA_PARTNER_FORM)

  - concurrency_attribute: higher value implies more likely to form new
                           partners if in partnership

    Initialize to beta distribution(ALPHA_CONCURRENCY, BETA_CONCURRENCY)

  - sexual_drive_attribute: higher value implies more likely to have sex

    Initialize to beta distribution(ALPHA_DRIVE, BETA_DRIVE)

  - preference_fifs_attribute: preference to have sex with least recently
                               formed partner (fifs = first in, first sex)

    Initialize to beta distribution(ALPHA_FIFS, BETA_FIFS)

    The above parameters can also be differentiated by sex. E.g. ALPHA_MALE_FIFS
    and ALPHA_FEMALE_FIFS.



*/

#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "stats.hh"

thread_local std::mt19937 rng;

const double YEAR_IN_DAYS = 365.25;
const double YEAR = 1.0;
const double MONTH = 1.0 / 12.0;
const double WEEK = 1.0 / 52.0;
const double DAY = 1.0 / YEAR_IN_DAYS;
const double HOUR = DAY / 24.0;

class Agent;

typedef std::unordered_map<const char *, double> ParameterMap;
typedef std::vector<Agent *> AgentVector;

enum Sex {
  MALE = 0,
  FEMALE = 1
};

struct Agent {
  unsigned id;
  Sex sex;
  double age;
  /* 0=HIV-
     1=HIV+ primary infection
     2=HIV+ CDC stage 1
     ...
     5=HIV+ CDC stage 4
   */
  unsigned hiv;
  bool alive;
  AgentVector partners;

  /* Attributes */
  double relationship_stickiness_attribute;
  double partner_forming_attribute;
  double concurrency_attribute;
  double sexual_drive_attribute;
  double preference_fifs_attribute;
  double force_infection_attribute;

  void init(unsigned i, const ParameterMap& parameters)
  {
    id = i;
    sex = std::bernoulli_distribution(0.5)(rng) == 0 ? MALE : FEMALE;
    age = std::uniform_real_distribution<double>(15.0, 20.0)(rng);
    hiv = std::min(std::geometric_distribution<int>(0.9)(rng), 5);
    alive = true;
    relationship_stickiness_attribute =
      sim::beta_distribution<>(2.0, parameters.at("MEAN_PARTNERSHIP_TIME") /
			       parameters.at("TIME_STEP") * 2.0) (rng);
    partner_forming_attribute =
      sim::beta_distribution<>(2.0, parameters.at("MEAN_TIME_UNTIL_PARTNER") /
			       parameters.at("TIME_STEP") * 2.0) (rng);
    concurrency_attribute =
      sim::beta_distribution<>(2.0, parameters.at("MEAN_TIME_CONCURRENT") /
			       parameters.at("TIME_STEP") * 2.0) (rng);
    sexual_drive_attribute =
      sim::beta_distribution<>(2.0, parameters.at("MEAN_TIME_SEX") /
			       parameters.at("TIME_STEP") * 2.0) (rng);
    preference_fifs_attribute =
      sim::beta_distribution<>(2.0, 2.0 / parameters.at("PREFERENCE_FIFS")
			       - 2.0)(rng);
    force_infection_attribute = sex == MALE ?
      sim::beta_distribution<>(2.0, 2.0 / parameters.at("MEAN_RISK_HET_MALE_SEX")
			       - 2.0)(rng) :
      sim::beta_distribution<>(2.0,
			       2.0 / parameters.at("MEAN_RISK_HET_FEMALE_SEX")
			       - 2.0) (rng);
  }

  // EVENTS

  void simple_infection_event(const double prevalence_males,
			      const double prevalence_females)
  {
    if (hiv == 0) {
      double prevalence = sex == MALE ? prevalence_females : prevalence_males;
      double risk_infection = force_infection_attribute *
	partner_forming_attribute * prevalence;
      if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < risk_infection)
	hiv = 1;
    }
  }

  void stage_advance_event(const double prob_leave_acute_infection)
  {
    if (hiv == 1 && std::uniform_real_distribution<double>(0.0, 1.0)(rng) <
	prob_leave_acute_infection)
      ++hiv;
  }

  // Every agent has to age on each iteration of the simulation
  void age_event(const double time_elapsed)
  {
    age += time_elapsed;
  }


};


void
initialize_agents(AgentVector& agents, const ParameterMap parameters)
{
  unsigned i = 0;
  for (auto &agent : agents) {
    agent = new Agent();
    agent->init(i++, parameters);
  }
}

void destroy_agents(AgentVector& agents)
{
  for (auto &agent: agents)
    delete agent;
}


struct Prevalence {
  unsigned males_alive = 0;
  unsigned females_alive = 0;
  unsigned males_infected = 0;
  unsigned females_infected = 0;
  double male_prevalence, female_prevalence;
};

static Prevalence calc_prevalence(const AgentVector& agents)
{
  Prevalence p;
  for (auto& agent: agents) {
    if (agent->alive) {
      if (agent->sex == MALE) {
	++p.males_alive;
	if (agent->hiv > 0)
	  ++p.males_infected;
      } else {
	++p.females_alive;
	if (agent->hiv > 0)
	  ++p.females_infected;
      }
    }
  }
  p.male_prevalence = (double) p.males_infected / p.males_alive;
  p.female_prevalence = (double) p.females_infected / p.females_alive;
  return p;
}

// On each step of the iteration we write out CSV data
void report(double date,  const AgentVector& agents)
{
  // date, num agents, num alive, num infected, num alive infected
  Prevalence p = calc_prevalence(agents);

  unsigned hiv[6] = {0,0,0,0,0,0};
  for (auto & agent: agents)
    ++hiv[agent->hiv];

  std::cout << date << ", "
	    << agents.size() << ", "
	    << p.males_alive + p.females_alive << ", "
	    << p.males_infected + p.females_infected << ", "
	    << (double) (p.males_infected + p.females_infected) /
    (p.males_alive + p.females_alive) << ", "
	    << p.males_alive << ", "
	    << p.males_infected << ", "
	    << p.male_prevalence << ", "
	    << p.females_alive << ", "
	    << p.females_infected << ", "
	    << p.female_prevalence << ", "
	    << hiv[0] << ", " << hiv[1] << ", " << hiv[2] << ", "
	    << hiv[3] << ", " << hiv[4] << ", " << hiv[5]
	    << std::endl;
}

void summary(const unsigned sim_no, const char* description,
	     const AgentVector& agents, ParameterMap &outputs)
{
  unsigned males = 0, females, hiv_males = 0, hiv_females = 0;
  std::vector<unsigned> hiv(6);
  double avg_age = 0.0;
  double youngest = agents[0]->age;
  double oldest = agents[0]->age;
  for (auto& agent: agents) {
    ++hiv[agent->hiv];
    if (agent->sex == MALE) {
      ++males;
      if (agent->hiv > 0) ++hiv_males;
    } else {
      if (agent->hiv > 0) ++hiv_females;
    }
    avg_age += agent->age;
    if (agent->age > oldest) oldest = agent->age;
    if (agent->age < youngest) youngest = agent->age;
  }
  std::ostringstream prefix_stream;
  prefix_stream << "summary," << sim_no << "," << description << ",";
  std::string prefix = prefix_stream.str();
  std::cout << prefix
	    << "males: " << males << std::endl;
  females = agents.size() - males;
  std::cout << prefix
	    << "females," << agents.size() - males << std::endl;
  std::cout << prefix
	    << "youngest," << youngest << std::endl;
  std::cout << prefix
	    << "oldest," << oldest << std::endl;
  std::cout << prefix
	    << "Average age," << avg_age / agents.size() << std::endl;
  for (size_t i = 0; i < 6; ++i)
    std::cout << prefix << "HIV " << i << " " << hiv[i] << std::endl;
  std::cout << prefix
	    << "Male prevalence: " << (double) hiv_males / males << std::endl;
  std::cout << prefix
	    << "Female prevalence: " << (double) hiv_females / females
	    << std::endl;
  // Incidence
  if(outputs.find("HIV_MALES") != outputs.end()) {
    unsigned diff_hiv_males = hiv_males - outputs.at("HIV_MALES");
    unsigned diff_hiv_females = hiv_females - outputs.at("HIV_FEMALES");
    std::cout << prefix << "Male incidence: " << (double) diff_hiv_males / males
	      << std::endl;
    std::cout << prefix
	      << "Female incidence: " << (double) diff_hiv_females / females
	      << std::endl;
    std::cout << prefix
	      << "Incidence: " << (double) (diff_hiv_males + diff_hiv_females) /
      agents.size() << std::endl;
  }
  outputs["HIV_MALES"] = hiv_males;
  outputs["HIV_FEMALES"] = hiv_females;
}

void simulate(AgentVector& agents,
	      const ParameterMap& parameters)
{
  double num_years = parameters.at("NUM_YEARS");
  double time_step = parameters.at("TIME_STEP");
  double start_date = parameters.at("START_DATE");
  double prob_leave_acute_infection = parameters.at("LEAVE_ACUTE_INFECTION");
  unsigned num_iterations = num_years / time_step;

  for (unsigned i = 0; i < num_iterations; ++i) {
    shuffle(agents.begin(), agents.end(), rng);

    Prevalence p = calc_prevalence(agents);

    for (auto & agent: agents) {
      if (agent->alive) {
	agent->simple_infection_event(p.male_prevalence, p.female_prevalence);
	agent->stage_advance_event(prob_leave_acute_infection);
	agent->age_event(time_step);
      }
    }
    report(start_date + time_step * i, agents);
  }
}


int main(int argc, char *argv[])
{
  // Set our parameters
  ParameterMap parameters, outputs;

  parameters["NUM_YEARS"] = 2.0;
  parameters["TIME_STEP"] = DAY;
  parameters["START_DATE"] = 2015.0;

  /* Parameters to estimate */
  parameters["MEAN_TIME_UNTIL_PARTNER"] = YEAR / 4.0;
  parameters["MEAN_PARTNERSHIP_TIME"] = YEAR / 4.0;
  parameters["MEAN_TIME_CONCURRENT"] = YEAR;
  parameters["MEAN_TIME_SEX"] = DAY;
  parameters["PREFERENCE_FIFS"] = 0.5;
  parameters["MEAN_RISK_HET_MALE_SEX"] = 0.01;
  parameters["MEAN_RISK_HET_FEMALE_SEX"] = 0.02;
  parameters["LEAVE_ACUTE_INFECTION"] = 0.0238095238;

  // Seed our Mersenne Twister to some arbitrarily chosen number
  rng.seed(23);
  AgentVector agents(10000);
  initialize_agents(agents, parameters);
  summary(0, "begin", agents, outputs);
  std::cout << "year, agents, alive, infected, prevalence, males_alive, "
    "males_infected, male_prevalence, females_alive, females_infected, "
    "female_prevalence, hiv_neg, hiv_p, cdc1, cdc2, cdc3, cdc4"
	    << std::endl;
  report(parameters["START_DATE"], agents);
  simulate(agents, parameters);
  summary(0, "end", agents, outputs);
  destroy_agents(agents);
}
