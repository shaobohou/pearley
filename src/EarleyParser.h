#ifndef EARLEY_PARSER
#define EARLEY_PARSER

#include "tnt/inv.h"
#include "GrammarRule.h"
#include "EarleyState.h"
#include "ParseTree.h"

#include <algorithm>
#include <map>
#include <set>

using TNT::eye;


// Earley Parser template definition
template <typename T>
class EarleyParser : public Stringable
{
    public:
        typedef std::vector<EarleyState<T> > EarleySet;

        EarleyParser();
        EarleyParser(const std::vector<GrammarRule<T> > &the_rules, const T &start_symbol, const T &wild_card);

        std::vector<std::vector<unsigned int> > parse(const std::vector<T> &sentence, std::vector<EarleySet> &earley_chart) const;

        std::vector<unsigned int> parse_init(std::vector<EarleySet> &earley_chart) const;
        std::vector<unsigned int> parse_word(std::vector<EarleySet> &earley_chart, const T &in_word, bool do_prediction) const;

        ParseTree<T> getViterbiParse(const EarleyState<T> &end_state, const std::vector<EarleySet> &earley_chart) const;
        std::vector<std::pair<T, double> > getNextWordTransitions(const EarleySet &in_states) const;

        void printEarleySet(const EarleySet &early_set) const;
        virtual std::string toString() const;

    private:
        // given private variables
        std::vector<GrammarRule<T> > rules;
        const T the_start_symbol, the_wild_card, the_dummy_symbol;

        // relationships
        Array2D<double> left_nonterminal_corner;
        Array2D<double> transitive_left_nonterminal_corner;
        Array2D<double> left_terminal_corner;
        Array2D<double> transitive_left_terminal_corner;

        // constructed private variables
        std::vector<T> terminals;
		std::map<T, unsigned int> terminal2index;
        std::vector<RulesList<T> > rules_lists;
		std::map<T, unsigned int> nonterminal2rules_list;

        void checkInputs() const;
        void normaliseWeights();

        // in_parse_column refers to the column that in_states originated from in the parse_table
        EarleySet predict(const EarleySet &in_states, unsigned int in_parse_column) const;
        EarleySet scan(const EarleySet &in_states, unsigned int in_parse_column, const std::set<T> &constraint_terminals) const;
        EarleySet complete(const EarleySet &in_states, const std::vector<EarleySet> &in_chart) const;

        static void accumulateCompletionProbability(EarleySet &completed_states, ProgenitorsLists &progenitors_lists, unsigned int index, const std::vector<EarleySet> &in_chart);
};



// implementations

using std::vector;
using std::string;
using std::pair;
using std::map;
using std::set;

// public template functions

template <typename T>
EarleyParser<T>::EarleyParser()
{}

template <typename T>
EarleyParser<T>::EarleyParser(const vector<GrammarRule<T> > &the_rules, const T &start_symbol, const T &wild_card)
:rules(the_rules), the_start_symbol(start_symbol), the_wild_card(wild_card), the_dummy_symbol(" ")
{
    if(rules.size() == 0)
    {
        std::cout << "There should be at least one rule." << std::endl;
        exit(1);
    }

    // added rules into mapping from LHS to RHS
    for(unsigned int i = 0; i < rules.size(); i++)
    {
        const T &LHS = rules[i].LHS();
        if(nonterminal2rules_list.find(LHS) == nonterminal2rules_list.end())
        {
            nonterminal2rules_list[LHS] = rules_lists.size();
            rules_lists.push_back(RulesList<T>(LHS));
        }

        unsigned int nonterminal_index = nonterminal2rules_list[LHS];
        rules_lists[nonterminal_index].rules().push_back(i);
    }

    // add a symbol to the terminals iff it is not a nonterminal (not in the map lhs2rules) and is not already in the terminals
    for(unsigned int i = 0; i < rules.size(); i++)
        for(unsigned int j = 0; j < rules[i].RHS().size(); j++)
        {
            const T &test_symbol = rules[i].RHS()[j];
            if(nonterminal2rules_list.find(test_symbol) == nonterminal2rules_list.end())
                if(find(terminals.begin(), terminals.end(), test_symbol) == terminals.end())
                {
                    terminals.push_back(test_symbol);
                    terminal2index[test_symbol] = terminals.size()-1;
                }
        }

    // check inputs and normalise weights
    checkInputs();
    normaliseWeights();

    // compute left nonterminal corner relationship
    left_nonterminal_corner = Array2D<double>(rules_lists.size(), rules_lists.size(), 0.0);
    for(unsigned int i = 0; i < rules_lists.size(); i++)
    {
        const T &LHS = rules_lists[i].nonterminal();
        for(unsigned int j = 0; j < rules_lists.size(); j++)
        {
            const vector<unsigned int> &rules_index = rules_lists[j].rules();
            for(unsigned int k = 0; k < rules_index.size(); k++)
                if(LHS == rules[rules_index[k]].RHS()[0])
                    left_nonterminal_corner(j, i) += rules[rules_index[k]].weight();
        }
    }

    // compute transitive left nonterminal corner relationship
//    transitive_left_nonterminal_corner = diag(Array1D<double>(left_nonterminal_corner.dim1(), 1.0)) - left_nonterminal_corner;
    transitive_left_nonterminal_corner = eye<double>(left_nonterminal_corner.dim1()) - left_nonterminal_corner;
    transitive_left_nonterminal_corner = inv(transitive_left_nonterminal_corner);
    for(int i = 0; i < transitive_left_nonterminal_corner.dim1(); i++)
        for(int j = 0; j < transitive_left_nonterminal_corner.dim2(); j++)
        {
            if(transitive_left_nonterminal_corner(i, j) > 0.0)
                rules_lists[i].nonterminal_corners().push_back(j);

            if(transitive_left_nonterminal_corner(i, j) <= 0.0)
                transitive_left_nonterminal_corner(i, j) = 0.0;
        }

    // compute left terminal corner relationship
    left_terminal_corner = Array2D<double>(rules_lists.size(), terminals.size(), 0.0);
    for(unsigned int i = 0; i < rules_lists.size(); i++)
        for(unsigned int j = 0; j < terminals.size(); j++)
        {
            const vector<unsigned int> &rules_index = rules_lists[i].rules();
            for(unsigned int k = 0; k < rules_index.size(); k++)
                if(rules[rules_index[k]].RHS()[0] == terminals[j])
                    left_terminal_corner(i, j) += rules[rules_index[k]].weight();
        }

    // compute transitive left terminal corner relationship
    transitive_left_terminal_corner = transitive_left_nonterminal_corner * left_terminal_corner;
    for(unsigned int i = 0; i < rules_lists.size(); i++)
        for(unsigned int j = 0; j < terminals.size(); j++)
            if(transitive_left_terminal_corner(i, j) > 0.0)
                rules_lists[i].terminal_corners().push_back(j);

    // added rules for wild card state
    nonterminal2rules_list[the_wild_card] = rules_lists.size();
    rules_lists.push_back(RulesList<T>(the_wild_card));
    for(int i = 0; i < transitive_left_nonterminal_corner.dim2(); i++)
        rules_lists.back().nonterminal_corners().push_back(i);
    for(unsigned int i = 0; i < terminals.size(); i++)
        rules_lists.back().terminal_corners().push_back(i);
}

template <typename T>
vector<vector<unsigned int> > EarleyParser<T>::parse(const vector<T> &sentence, vector<EarleySet> &earley_chart) const
{
    for(unsigned int i = 0; i < sentence.size(); i++)
        std::cout << sentence[i] << " ";
    std::cout << std::endl << std::endl;

    vector<vector<unsigned int> > states_counts;
    earley_chart.clear();

    // initialise parse
    states_counts.push_back(parse_init(earley_chart));
    if(earley_chart.back().size() == 0)
    {
        std::cout << "No new states predicted." << std::endl;
        return states_counts;
    }

    // parse each word
    for(unsigned int i = 0; i < sentence.size(); i++)
        states_counts.push_back(parse_word(earley_chart, sentence[i], (i < (sentence.size()-1))));

    return states_counts;
}

template <typename T>
vector<unsigned int> EarleyParser<T>::parse_init(vector<EarleySet> &earley_chart) const
{
    vector<unsigned int> states_count(3, static_cast<unsigned int>(0));

    // initialise the first state
    EarleyState<T> start_state(GrammarRule<T>(the_dummy_symbol, std::vector<T>(1, the_start_symbol), 1.0), 0, 0, 1.0, 1.0);
    earley_chart = vector<EarleySet>(1, EarleySet(1, start_state));
    states_count[1] = earley_chart.back().size();
    std::cout << "-----------------------------------------------------" << std::endl << std::endl;
    std::cout << earley_chart[0][0] << std::endl << std::endl;

    // perform initial prediction
    std::cout << "predicted" << std::endl;
    EarleySet new_set = predict(earley_chart[0], 0);
    earley_chart[0].insert(earley_chart[0].end(), new_set.begin(), new_set.end());
    states_count[2] = earley_chart.back().size();
    std::cout << std::endl << std::endl << std::endl;

    return states_count;
}

template <typename T>
vector<unsigned int> EarleyParser<T>::parse_word(vector<EarleySet> &earley_chart, const T &in_word, bool do_prediction) const
{
    vector<unsigned int> states_count(3, static_cast<unsigned int>(0));

    std::cout << in_word << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl << std::endl;

    std::cout << "scanned" << std::endl;
    set<T> constraint_terminals;
    constraint_terminals.insert(in_word);
    // constraint_terminals.insert(terminals.begin(), terminals.end());
    EarleySet scanned_states = scan(earley_chart.back(), earley_chart.size()-1, constraint_terminals);
    earley_chart.push_back(scanned_states);
    states_count[0] = earley_chart.back().size();
    std::cout << std::endl;
    if(scanned_states.size() == 0)
    {
        std::cout << "No new states scanned." << std::endl;
        return states_count;
    }

    std::cout << "completed" << std::endl;
    EarleySet completed_states = complete(earley_chart.back(), earley_chart);
    earley_chart.back().insert(earley_chart.back().end(), completed_states.begin(), completed_states.end());
    states_count[1] = earley_chart.back().size();
    std::cout << std::endl;

    if(!do_prediction)
    {
        states_count[2] = earley_chart.back().size();
        return states_count;
    }

    std::cout << "predicted" << std::endl;
    EarleySet predicted_states = predict(earley_chart.back(), earley_chart.size()-1);
    earley_chart.back().insert(earley_chart.back().end(), predicted_states.begin(), predicted_states.end());
    states_count[2] = earley_chart.back().size();
    std::cout << std::endl << std::endl << std::endl;

    return states_count;
}

template <typename T>
ParseTree<T> EarleyParser<T>::getViterbiParse(const EarleyState<T> &end_state, const std::vector<EarleySet> &earley_chart) const
{
    if(end_state.dot() != 0)
    {
        assert(end_state.hasProgenitors());

        const Progenitor &progenitor = end_state.viterbiProgenitor();
        const ParseEntry &incompleted_index = progenitor.incompleted_index();
        const EarleyState<T> &incompleted_state = earley_chart[incompleted_index.col()][incompleted_index.row()];

        // recursively call parse tree building on the first predecessor state
        ParseTree<T> left_tree = getViterbiParse(incompleted_state, earley_chart);

        const T &end_symbol = end_state.lastSymbol();
        typename map<T, unsigned int>::const_iterator nonterminal2rules_list_it = nonterminal2rules_list.find(end_symbol);
        if(nonterminal2rules_list_it == nonterminal2rules_list.end())
            left_tree.attach(1, ParseTree<T>(std::vector<T>(1, end_symbol)));
        else
        {
            const ParseEntry &complete_index = progenitor.complete_index();
            const EarleyState<T> &complete_state = earley_chart[complete_index.col()][complete_index.row()];

            // recursively call parse tree building on the second predecessor state
            ParseTree<T> right_tree = getViterbiParse(complete_state, earley_chart);
            left_tree.attach(1, right_tree);
        }

        return left_tree;
    }
    else
        return ParseTree<T>(std::vector<T>(1, end_state.rule().LHS()));
}

template <typename T>
std::vector<std::pair<T, double> > EarleyParser<T>::getNextWordTransitions(const EarleySet &in_states) const
{
    double prob_sum = 0.0;
    map<T, double> symbol2prob;
    for(unsigned int i = 0; i < in_states.size(); i++)
    {
        if(in_states[i].isComplete())
            continue;

        const T &symbol_after_dot = in_states[i].nextSymbol();
        typename map<T, unsigned int>::const_iterator nonterminal2rules_list_it = nonterminal2rules_list.find(symbol_after_dot);
        if(nonterminal2rules_list_it != nonterminal2rules_list.end())
            continue;

        prob_sum += in_states[i].alpha();
        typename map<T, double>::iterator symbol2prob_it = symbol2prob.find(symbol_after_dot);
        if(symbol2prob_it == symbol2prob.end())
            symbol2prob[symbol_after_dot] = in_states[i].alpha();
        else
            symbol2prob_it->second += in_states[i].alpha();
    }

    for(typename map<T, double>::iterator it = symbol2prob.begin(); it != symbol2prob.end(); it++)
        it->second /= prob_sum;

    return std::vector<std::pair<T, double> >(symbol2prob.begin(), symbol2prob.end());
}

template <typename T>
void EarleyParser<T>::printEarleySet(const EarleySet &early_set) const
{
    for(unsigned int i = 0; i < early_set.size(); i++)
        std::cout << early_set[i] << std::endl;
}

template <typename T>
std::string EarleyParser<T>::toString() const
{
    std::ostringstream sout;

    sout << "Terminals" << std::endl;
    for(unsigned int i = 0; i < terminals.size(); i++)
        sout << i << " = " << terminals[i] << std::endl;
    sout << std::endl;

    sout << "Grammar Rules:" << std::endl;
    for(unsigned int i = 0; i < rules.size(); i++)
        sout << i << " = " << rules[i] << std::endl;

    return sout.str();
}



// private template functions

template <typename T>
void EarleyParser<T>::checkInputs() const
{
    if(nonterminal2rules_list.find(the_wild_card) != nonterminal2rules_list.end())
    {
        std::cout << "The wild card should not be one of the nonterminals i.e. appear on the RHS of any rules." << std::endl;
        exit(1);
    }

    if(find(terminals.begin(), terminals.end(), the_wild_card) != terminals.end())
    {
        std::cout << "The wild card should not be one of the terminals i.e. appear on the LHS of any rules." << std::endl;
        exit(1);
    }

    if(nonterminal2rules_list.find(the_dummy_symbol) != nonterminal2rules_list.end())
    {
        std::cout << "The dummy symbol should not be one of the nonterminals i.e. appear on the RHS of any rules." << std::endl;
        exit(1);
    }

    if(find(terminals.begin(), terminals.end(), the_dummy_symbol) != terminals.end())
    {
        std::cout << "The dummy symbol should not be one of the terminals i.e. appear on the LHS of any rules." << std::endl;
        exit(1);
    }
}

template <typename T>
void EarleyParser<T>::normaliseWeights()
{
    for(unsigned int i = 0; i < rules_lists.size(); i++)
    {
        const vector<unsigned int> &rules_index = rules_lists[i].rules();

        double total_weight = 0.0;
        for(unsigned int j = 0; j < rules_index.size(); j++)
            total_weight += rules[rules_index[j]].weight();

        for(unsigned int j = 0; j < rules_index.size(); j++)
            rules[rules_index[j]].weight() = rules[rules_index[j]].weight()/total_weight;
    }
}

template <typename T>
typename EarleyParser<T>::EarleySet EarleyParser<T>::predict(const EarleySet &in_states, unsigned int in_parse_column) const
{
    //assert(in_states.size() > 0);

    // accumulate alphas of inputs states with the same nonterminal Z to the right of the dot
    map<unsigned int, double> Z2alphas;
    for(unsigned int i = 0; i < in_states.size(); i++)
    {
        // if the current_state is already complete, proceed to the next state
        const EarleyState<T> &current_state = in_states[i];
        if(current_state.isComplete())
            continue;

        // if the next symbol after the dot is a terminal, proceed to the next state
        const T& symbol_after_dot = current_state.nextSymbol();
        typename map<T, unsigned int>::const_iterator nonterminal2rules_list_it = nonterminal2rules_list.find(symbol_after_dot);
        if(nonterminal2rules_list_it == nonterminal2rules_list.end())
            continue;

        // accumulate alpha for state with the same unexpanded nonterminal to the right of the dot (Z)
        unsigned int list_index = nonterminal2rules_list_it->second;
        map<unsigned int, double>::iterator Z2alphas_it = Z2alphas.find(list_index);
        if(Z2alphas_it == Z2alphas.end())
            Z2alphas[list_index]  = in_states[i].alpha();
        else
            Z2alphas_it->second  += in_states[i].alpha();
    }

    // accumulate alphas of the nonterminal Y to the right of the dot through transitive left corner relationship, Z ==> Y
    map<unsigned int, double> Y2alphas;
    for(map<unsigned int, double>::const_iterator it = Z2alphas.begin(); it != Z2alphas.end(); it++)
    {
        const unsigned int Z_index = it->first;
        const double Z_alpha = it->second;

        // accumulate alpha for each possible production rule
        const std::vector<unsigned int> &left_corners = rules_lists[Z_index].nonterminal_corners();
        for(unsigned int j = 0; j < left_corners.size(); j++)
        {
            unsigned int left_corner_index = left_corners[j];
            map<unsigned int, double>::iterator Y2alphas_it = Y2alphas.find(left_corner_index);
            if(Y2alphas_it == Y2alphas.end())
                Y2alphas[left_corner_index]  = Z_alpha*transitive_left_nonterminal_corner(Z_index, left_corner_index);
            else
                Y2alphas[left_corner_index] += Z_alpha*transitive_left_nonterminal_corner(Z_index, left_corner_index);
        }
    }

    // add each state Y --> v, update the alphas and gammas
    EarleySet predicted_states;
    for(map<unsigned int, double>::const_iterator it = Y2alphas.begin(); it != Y2alphas.end(); it++)
    {
        const unsigned int Y_index = it->first;
        const double Y_alpha = it->second;
        const vector<unsigned int> rules_index = rules_lists[Y_index].rules();
        for(unsigned int j = 0; j < rules_index.size(); j++)
        {
            unsigned int temp_rule_index = rules_index[j];
            predicted_states.push_back(EarleyState<T>(rules[temp_rule_index], in_parse_column, 0));
            predicted_states.back().alpha() = Y_alpha*rules[temp_rule_index].weight();
            predicted_states.back().gamma() = rules[temp_rule_index].weight();
        }
    }

    // print the new states
    printEarleySet(predicted_states);

    return predicted_states;
}

template <typename T>
typename EarleyParser<T>::EarleySet EarleyParser<T>::scan(const EarleySet &in_states, unsigned int in_parse_column, const set<T> &constraint_terminals) const
{
    //assert(in_states.size() > 0);

    EarleySet scanned_states;
    for(unsigned int i = 0; i < in_states.size(); i++)
    {
        // if current_state is already complete, proceed to the next state
        const EarleyState<T> &current_state = in_states[i];
        if(current_state.isComplete())
            continue;

        // get the next symbol after the dot and test if it is a terminal by looking for it in the RHS of the rules
        const T& symbol_after_dot = current_state.nextSymbol();
        typename map<T, unsigned int>::const_iterator nonterminal2rules_list_it = nonterminal2rules_list.find(symbol_after_dot);

        // if symbol is a nonterminal then proceed to the next state
        if(nonterminal2rules_list_it != nonterminal2rules_list.end())
            continue;

        // if the terminal symbol is not one of the constraint terminal then proceed to the next state
        if(constraint_terminals.find(symbol_after_dot) == constraint_terminals.end())
            continue;

        // add scanned state and move the dot one position to the right
        scanned_states.push_back(EarleyState<T>(current_state.rule(), current_state.start(), current_state.dot()+1, current_state.alpha(), current_state.gamma()));
        scanned_states.back().addProgenitor(Progenitor(ParseEntry(in_parse_column, i), scanned_states.back().gamma()));    // assume some knowledge about the structure of the parse table
    }

    // print the new states
    printEarleySet(scanned_states);

    return scanned_states;
}

template <typename T>
typename EarleyParser<T>::EarleySet EarleyParser<T>::complete(const EarleySet &in_states, const vector<EarleySet> &in_chart) const
{
    //assert(in_states.size() > 0);

    // initialise the completed states with the input states, erase them later
    EarleySet completed_states = in_states;
    unsigned int number_of_in_states = in_states.size();

    // for each completed state contains list of indices of states that it depends on
    ProgenitorsLists progenitors_lists(number_of_in_states, vector<Progenitor>());
    // maps from a completed state to its index in the return list
    map<EarleyState<T>, unsigned int> completed2index;
    for(unsigned int i = 0; i < completed_states.size(); i++)
    {
        // if the completed_state is complete, then use it to get more completed states
        if(!completed_states[i].isComplete())
            continue;

        // get the first unused complete state in the list and use it to complete more states
        const T complete_nonterminal = completed_states[i].rule().LHS();
        const unsigned int start_position = completed_states[i].start();

        for(unsigned int j = 0; j < in_chart[start_position].size(); j++)
        {
            // if current_state is already complete, then proceed to the next state
            const EarleyState<T> &current_in_state = in_chart[start_position][j];
            if(current_in_state.isComplete())
                continue;

            // if the nonterminal symbol to the right of the dot does not match the complete LHS, then proceed to the next state
            const T &symbol_after_dot = current_in_state.nextSymbol();
            if(symbol_after_dot != complete_nonterminal)
                continue;

            // maps from index into the all input states to index into the completed states
            EarleyState<T> new_completed_state(current_in_state.rule(), current_in_state.start(), current_in_state.dot()+1);
            typename map<EarleyState<T>, unsigned int>::iterator completed2index_it = completed2index.find(new_completed_state);
            // if the input state has not already been completed, do so
            if(completed2index_it == completed2index.end())
            {
                // add a new completed state
                completed_states.push_back(new_completed_state);
                // add the mapping from all input states and initialise the dependent states list
                completed2index[new_completed_state] = completed_states.size()-1;
                progenitors_lists.push_back(vector<Progenitor>());
            }

            // update dependency list, the new completed state depends on the current complete state used
            unsigned int completed_out_index = completed2index[new_completed_state];
            progenitors_lists[completed_out_index].push_back(Progenitor(ParseEntry(start_position, j), ParseEntry(in_chart.size()-1, i)));
        }
    }

    // call recursive function to compute the probabilities
    for(unsigned int i = number_of_in_states; i < completed_states.size(); i++)
        accumulateCompletionProbability(completed_states, progenitors_lists, i, in_chart);

    // sort the states and correct the row number of the complete state used, only necessary if want the print result to look nice
    if(true)
    {
        // sort and record mapping
        sort(completed_states.begin()+number_of_in_states, completed_states.end());
        std::vector<unsigned int> unsorted2sorted(completed_states.size(), 0);
        for(unsigned int i = 0; i < number_of_in_states; i++)
            unsorted2sorted[i] = i;
        for(unsigned int i = number_of_in_states; i < completed_states.size(); i++)
            unsorted2sorted[completed2index[completed_states[i]]] = i;

        // correct the row index
        for(unsigned int i = number_of_in_states; i < completed_states.size(); i++)
            for(unsigned int j = 0; j < completed_states[i].progenitors().size(); j++)
            {
                unsigned int old_index = completed_states[i].progenitors()[j].complete_index().row();
                completed_states[i].progenitors()[j].complete_index().row() = unsorted2sorted[old_index];
            }
    }

    // erase complete states from the scanning step
    completed_states.erase(completed_states.begin(), completed_states.begin()+number_of_in_states);

    // print the new states
    printEarleySet(completed_states);

    return completed_states;
}

template <typename T>
void EarleyParser<T>::accumulateCompletionProbability(EarleySet &completed_states, ProgenitorsLists &progenitors_lists, unsigned int index, const vector<EarleySet> &earley_chart)
{
    if(progenitors_lists[index].size() > 0)
    {
        completed_states[index].alpha() = 0.0;
        completed_states[index].gamma() = 0.0;

        // recursively call the function on any state it depends on and accumulate the probabilities
        for(unsigned int i = 0; i < progenitors_lists[index].size(); i++)
        {
            const Progenitor &progenitor = progenitors_lists[index][i];
            const ParseEntry &incompleted_index = progenitors_lists[index][i].incompleted_index();
            const ParseEntry &complete_index = progenitors_lists[index][i].complete_index();
            const EarleyState<T> &incompleted_state = earley_chart[incompleted_index.col()][incompleted_index.row()];
            const EarleyState<T> &complete_state = completed_states[complete_index.row()];

            // recursively update dependent completed states
            accumulateCompletionProbability(completed_states, progenitors_lists, complete_index.row(), earley_chart);
            completed_states[index].alpha() += incompleted_state.alpha()*complete_state.gamma();
            completed_states[index].gamma() += incompleted_state.gamma()*complete_state.gamma();
            const double new_progenitor_weight = incompleted_state.viterbiProbability()*complete_state.viterbiProbability();
            completed_states[index].addProgenitor(Progenitor(progenitor.incompleted_index(), progenitor.complete_index(), new_progenitor_weight));
        }

        progenitors_lists[index].clear();
    }
}

#endif
