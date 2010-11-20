#ifndef EARLEY_STATE
#define EARLEY_STATE

#include "GrammarRule.h"
#include "BookKeeping.h"

template <typename T>
class EarleyState : public Stringable
{
    public:
        EarleyState();
        EarleyState(const GrammarRule<T> &rule, unsigned int start_position, unsigned int dot_position);
        EarleyState(const GrammarRule<T> &rule, unsigned int start_position, unsigned int dot_position, double forward_prob, double inner_prob);

        GrammarRule<T>& rule();
        unsigned int& start();
        unsigned int& dot();
        double& alpha();
        double& gamma();

        const GrammarRule<T>& rule() const;
        const unsigned int& start() const;
        const unsigned int& dot() const;
        const double& alpha() const;
        const double& gamma() const;

        bool operator==(const EarleyState &other) const;
        bool operator<(const EarleyState &other) const;

        const T& lastSymbol() const;
        const T& nextSymbol() const;
        bool isComplete() const;
        bool isFinished() const;

        std::vector<Progenitor>& progenitors();
        const std::vector<Progenitor>& progenitors() const;
        void addProgenitor(const Progenitor &progenitor);
        const Progenitor& viterbiProgenitor() const;
        const double& viterbiProbability() const;
        bool hasProgenitors() const;

        virtual std::string toString() const;

    private:
        GrammarRule<T> the_rule;
        unsigned int the_start, the_dot;
        double the_alpha, the_gamma;

        std::vector<Progenitor> the_progenitors;
        unsigned int the_viterbi_index;
};



// implementation

template <typename T>
EarleyState<T>::EarleyState()
: the_rule(), the_start(0), the_dot(0), the_alpha(0.0), the_gamma(0.0)
{}

template <typename T>
EarleyState<T>::EarleyState(const GrammarRule<T> &rule, unsigned int start_position, unsigned int dot_position)
: the_rule(rule), the_start(start_position), the_dot(dot_position), the_alpha(0.0), the_gamma(0.0)
{}

template <typename T>
EarleyState<T>::EarleyState(const GrammarRule<T> &rule, unsigned int start_position, unsigned int dot_position, double forward_prob, double inner_prob)
: the_rule(rule), the_start(start_position), the_dot(dot_position), the_alpha(forward_prob), the_gamma(inner_prob)
{}

template <typename T>
GrammarRule<T>& EarleyState<T>::rule()
{
    return the_rule;
}

template <typename T>
unsigned int& EarleyState<T>::start()
{
    return the_start;
}

template <typename T>
unsigned int& EarleyState<T>::dot()
{
    return the_dot;
}

template <typename T>
double& EarleyState<T>::alpha()
{
    return the_alpha;
}

template <typename T>
double& EarleyState<T>::gamma()
{
    return the_gamma;
}

template <typename T>
const GrammarRule<T>& EarleyState<T>::rule() const
{
    return the_rule;
}

template <typename T>
const unsigned int& EarleyState<T>::start() const
{
    return the_start;
}

template <typename T>
const unsigned int& EarleyState<T>::dot() const
{
    return the_dot;
}

template <typename T>
const double& EarleyState<T>::alpha() const
{
    return the_alpha;
}

template <typename T>
const double& EarleyState<T>::gamma() const
{
    return the_gamma;
}

template <typename T>
bool EarleyState<T>::operator==(const EarleyState<T> &other) const
{
    return (the_start == other.the_start) && (the_dot == other.the_dot) && (the_rule == other.the_rule);
}

template <typename T>
bool EarleyState<T>::operator<(const EarleyState &other) const
{
    if(the_start > other.the_start)
        return true;

    if(the_start == other.the_start)
    {
        if(the_dot > other.the_dot)
            return true;

        if(the_dot == other.the_dot)
            return the_rule.RHS() > other.the_rule.RHS();
    }

    return false;
}

template <typename T>
const T& EarleyState<T>::lastSymbol() const
{
	return the_rule.RHS()[the_dot-1];
}

template <typename T>
const T& EarleyState<T>::nextSymbol() const
{
	return the_rule.RHS()[the_dot];
}

template <typename T>
bool EarleyState<T>::isComplete() const
{
    return the_dot == the_rule.RHS().size();
}

template <typename T>
bool EarleyState<T>::isFinished() const
{
    return (the_start == 0) && isComplete();
}

template <typename T>
std::vector<Progenitor>& EarleyState<T>::progenitors()
{
    return the_progenitors;
}

template <typename T>
const std::vector<Progenitor>& EarleyState<T>::progenitors() const
{
    return the_progenitors;
}

template <typename T>
void EarleyState<T>::addProgenitor(const Progenitor &progenitor)
{
    the_progenitors.push_back(progenitor);

    if(the_progenitors.size() == 1)
        the_viterbi_index = 0;
    else
        if(progenitor.weight() > the_progenitors[the_viterbi_index].weight())
            the_viterbi_index = the_progenitors.size()-1;
}

template <typename T>
const Progenitor& EarleyState<T>::viterbiProgenitor() const
{
    return the_progenitors[the_viterbi_index];
}

template <typename T>
const double& EarleyState<T>::viterbiProbability() const
{
    if(the_viterbi_index < the_progenitors.size())
        return the_progenitors[the_viterbi_index].weight();
    else
        return the_gamma;   //viterbi is the same as gamma if it has never been completed, i.e. addCompletion was called
}

template <typename T>
bool EarleyState<T>::hasProgenitors() const
{
    return the_viterbi_index < the_progenitors.size();
}

template <typename T>
std::string EarleyState<T>::toString() const
{
    std::ostringstream sout;

    sout << the_start << " " << the_rule.LHS() << " -->";
    for(unsigned int i = 0; i < the_dot; i++)
        sout << " " << the_rule.RHS()[i];
    sout << ".";
    for(unsigned int i = the_dot; i < the_rule.RHS().size(); i++)
        sout << the_rule.RHS()[i] << " ";
    sout << "   [" << the_alpha << " " << the_gamma << "]   \t";
    sout << "(" << viterbiProbability() << ")";

    return sout.str();
}

#endif
