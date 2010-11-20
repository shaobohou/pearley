#ifndef GRAMMAR_RULE
#define GRAMMAR_RULE

#include <vector>
#include <string>
#include <sstream>

#include "Stringable.h"

template <typename T>
class GrammarRule : public Stringable
{
    public:
        GrammarRule()
        : the_weight(0.0)
        {}

        GrammarRule(const T &LHS, const std::vector<T> &RHS, double weight)
        :the_LHS(LHS), the_RHS(RHS), the_weight(weight)
        {}

        T& LHS()
        {
            return the_LHS;
        }

        const T& LHS() const
        {
            return the_LHS;
        }

        std::vector<T>& RHS()
        {
            return the_RHS;
        }

        const std::vector<T>& RHS() const
        {
            return the_RHS;
        }

        double& weight()
        {
            return the_weight;
        }

        const double& weight() const
        {
            return the_weight;
        }

        bool operator==(const GrammarRule<T> &other) const
        {
            return (the_LHS == other.the_LHS) && (the_RHS == other.the_RHS);
        }

        virtual std::string toString() const
        {
            std::ostringstream sout;

            sout << the_LHS << " -->";
            for(unsigned int i = 0; i < the_RHS.size(); i++)
                sout << " " << the_RHS[i];
            sout << "   [" << the_weight << "]";

            return sout.str();
        }

    private:
        T the_LHS;
        std::vector<T> the_RHS;
        double the_weight;
};



template <typename T>
class RulesList
{
    public:
        RulesList(const T &nonterminal)
        {
            the_nonterminal = nonterminal;
        }

        T& nonterminal()
        {
            return the_nonterminal;
        }

        std::vector<unsigned int>& rules()
        {
            return the_rules;
        }
        
        std::vector<unsigned int>& nonterminal_corners()
        {
            return the_nonterminal_corners;
        }
        
        std::vector<unsigned int>& terminal_corners()
        {
            return the_terminal_corners;
        }

        const T& nonterminal() const
        {
            return the_nonterminal;
        }

        const std::vector<unsigned int>& rules() const
        {
            return the_rules;
        }
        
        const std::vector<unsigned int>& nonterminal_corners() const
        {
            return the_nonterminal_corners;
        }
        
        const std::vector<unsigned int>& terminal_corners() const
        {
            return the_terminal_corners;
        }

    private:
        T the_nonterminal;
        std::vector<unsigned int> the_rules;
        std::vector<unsigned int> the_nonterminal_corners;
        std::vector<unsigned int> the_terminal_corners;
};

#endif
