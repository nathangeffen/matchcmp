# An algorithm for modelling sex among youths

## Notes by Nathan Geffen

### 18 June 2016

We are interested in modelling sexual networks among youths for the purpose of
understanding the transmission of sexually transmitted infections, especially
HIV.

Most models have very simple assumptions about the formation of sexual
relationships. Consider the South African HIV epidemic. The most widely cited
paper, Granich et al. 2009, barely considers heterogeneity, i.e. the differences
in behaviour across individuals.

More sophisticated models, such as those of the Actuarial Society of South
Africa, typically divide adults into four categories of varying sexual
behaviour: NOT (people having no sexual contact or in long-term monogamous
relationships who are not at risk of HIV), RSK (people in stable relationships
but at some risk either because they or their partners have more than one sexual
relationship), STD (people with high levels of sexually transmitted infections)
and PRO (described by Doyle as "people categorised by sexual mobility or
promiscuity, e.g. prostitutes and their frequent clients, etc.").

This was for deterministic models. For agent-based models we can introduce much
more sophisticated heterogeneity. For example, instead of placing agents in four
compartments, they can be modelled on a real-valued continuum that determines
their risk.

We want to capture these characteristics of sexual behaviour: youths form
partnerships, and have sex with their partners. They break up with their
partners and then find new partner with whom they have sex. Some youths will
form concurrent partnerships, i.e. have form a new partnership while already in
one.

While the vast majority of youths are heterosexual, many will experiment
with homosexuality, some will be bisexual and others exclusively
homosexual. But as a first version of this algorithm, only heterosexuality is
considered. However, we do not envisage that it will be hard to extend it to
include sexuality modelled on a continuum from exclusive heterosexuality,
through to exact bisexuality, through to exclusive homosexuality.

There is additional complexity our algorithm does not model, such as a youth
having casual sex with multiple people at the same time.

A typical agent based model is a simulation divided into multiple time
steps. Here the time step is one day, and it is executed 365 times, to model a
year's worth of events.

The algorithm has multiple parameters. We have no idea what
their values are, but given known outputs, such as HIV prevalence, and HIV
incidence over a one-year period, we can estimate or fit the model parameters.

## Description of Algorithm

For each agent on each iteration (or time step) of the simulation there are
three events which should be executed in this order: (1) breaking an existing
partnership, (2) forming a new partnership, (3) having sex with a partner.

- breakupPartnership

    Each agent has an attribute determining how likely they are to stay in a
    relationship. We call this relationship stickiness. Each agent has a
    relationship stickiness attribute. To determine if an agent breaks a
    relationship we generate a uniform random number. If it is less than the
    relationship stickiness attribute divided by the number of partners, we
    remove the partner with the lowest relationship stickiness value.

        if num_partners > 0:
            if uniform rng < (relationship_stickiness_attribute / num_partners)
                remove partner with lowest relationship_stickiness

- formNewPartnership

    Each agent has a propensity to form new partnerships. In the simple version
    of this algorithm if we generate a uniform random number less than the
    partner's propensity to form new partnerships value (divided by the number
    of partners plus one so as to make new partnerships increasingly unlikely
    as the number of existing partnerships for this agent increases) we find a
    new partner for this agent and add to the list of partners.

    In the more complex version. there's a separate attribute determining how
    likely agents already in relationships are to form a concurrent relationship.

        if uniform rng  < partner_forming_attribute /  (num_partners + 1)
            find a new matching partner and stick at back of partners queue

     MORE COMPLEX VERSION

        if num_partners == 0:
            if uniform rng < partner_forming_attribute
			find and add new partner to agent's partner list
        else:
            if uniform rng < concurrency_attribute / num_partners
			find and append new partner to agent's partners list

- haseSex

	Each agent has a sexual drive attribute. If the agent is in a partnership
    and a uniform random number indicates is less than the sexual drive
    attribute, the agent has sex with one of his or her partners randomly
    selected according to a geometric distribution (where the parameter
    indicates the preference for the earliest formed partnerships).

        if num_partners > 0 and uniform rng < sexual_drive_attribute:
            partner = min(geometric_distribution(preference_fifs_attribute),
                                                 num_partners)
	        determine hiv transmission risk


So we have the following attributes per agent that have to be estimated:

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

The above parameters should also be differentiated by sex.
E.g. ALPHA_MALE_FIFS and ALPHA_FEMALE_FIFS.

## Discussion

The algorithm described here is a sophisticated modelling of sexual
heterogeneity. But does such a complex algorithm necessarily capture reality any
better than simpler algorithms? I don't know the answer to this question. But
perhaps what we can answer is whether with different partner matching algorithms
we get vastly different estimates for the model parameters. If yes, then it does
suggest a limitation of this method, based as it is on almost complete ignorance
of the parameters. If no, then perhaps we are a small step closer to
understanding the sexual dynamics driving the HIV epidemic.
