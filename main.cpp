#include <iostream>
#include <fstream>

#include "EarleyParser.h"

using std::cout;
using std::endl;


std::vector<GrammarRule<std::string> > readRules(std::istream &in)
{
    std::vector<GrammarRule<std::string> > rules;
    unsigned int line_count = 0;
    while(!in.eof())
    {
        line_count++;

        std::string line;
        getline(in, line);
        std::stringstream ss(line);

        if(line.size() == 0)
            continue;

        double weight;
        ss >> weight;

        std::string LHS;
        ss >> LHS;

        std::string temp;
        ss >> temp;

        if(temp != "-->")
        {
            std::cout << "line " << line_count << ": [" << line << "] is not a rule." << std::endl;
            exit(1);
        }

        std::vector<std::string> RHS;
        while(!ss.eof())
        {
            std::string symbol;
            ss >> symbol;
            RHS.push_back(symbol);
        }

        rules.push_back(GrammarRule<std::string>(LHS, RHS, weight));
    }

    std::cout << rules.size() << " rules loaded." << std::endl;

    return rules;
}

void readGrammarFromFile(const std::string &filename, std::vector<GrammarRule<std::string> > &rules, std::string &start_symbol, std::string &wild_card)
{
    std::ifstream in(filename.c_str(), std::ios::in);
    if (!in.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    if(!in.eof())
    {
        std::string line;
        getline(in, line);
        std::stringstream ss(line);

        ss >> start_symbol;
        ss >> wild_card;
    }
    else
    {
        std::cout << "Unable to load start symbol and/or dummy symbol from first line in " << filename << std::endl;
        exit(1);
    }

    rules = readRules(in);
}

std::vector<std::vector<std::string> > readSentencesFromFile(const std::string &filename)
{
    std::ifstream in(filename.c_str(), std::ios::in);
    if (!in.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    std::vector<std::vector<std::string> > sentences;
    unsigned int line_count = 0;
    while(!in.eof())
    {
        line_count++;

        std::string line;
        getline(in, line);
        std::stringstream ss(line);

        if(line.size() == 0)
            continue;

        std::vector<std::string> sentence;
        while(!ss.eof())
        {
            std::string symbol;
            ss >> symbol;
            sentence.push_back(symbol);
        }

        sentences.push_back(sentence);
    }

    std::cout << sentences.size() << " sentence loaded." << std::endl;

    return sentences;
}

int main(int argc, char *argv[])
{
    if(argc < 3)
    {
        cout << "Usage:" << endl;
        cout << "EarleyParser <rules_filename> <sentences_filename>" << endl;
        exit(1);
    }

    std::string grammar_filename(argv[1]);
    std::string start_symbol, wild_card;
    std::vector<GrammarRule<std::string> > rules;
    readGrammarFromFile(grammar_filename, rules, start_symbol, wild_card);
    std::cout << std::endl;
    EarleyParser<std::string> parser(rules, start_symbol, wild_card);
    std::cout << parser << std::endl;

    std::cout << std::endl;
    std::string sentences_filename(argv[2]);
    std::vector<std::vector<std::string> > sentences = readSentencesFromFile(sentences_filename);
    std::cout << std::endl << std::endl;

    std::vector<EarleyParser<std::string>::EarleySet> earley_chart;
    std::vector<std::vector<unsigned int> > states_counts = parser.parse(sentences[0], earley_chart);
    std::cout << std::endl << "Finished parsing." << std::endl << std::endl;

    ParseTree<std::string> viterbi_parse_tree;
    for(unsigned int i = 0; i < earley_chart.back().size(); i++)
        if(earley_chart.back()[i].isFinished() && (earley_chart.back()[i].rule().LHS() == start_symbol))
        {
            std::cout << "Computing parse tree from " << i << "th state in the final set:   " << earley_chart.back()[i] << std::endl;
            viterbi_parse_tree = parser.getViterbiParse(earley_chart.back()[i], earley_chart);
        }
    viterbi_parse_tree.print(0, 0);
    
    for(unsigned int i = 0; i < states_counts.size(); i++)
    {
        for(unsigned int j = 0; j < states_counts[i].size(); j++)
            std::cout << states_counts[i][j] << " ";
        std::cout << std::endl;
    }
    
    for(unsigned int i = 0; i < earley_chart.size(); i++)
    {
        std::vector<std::pair<std::string, double> > transition_probs = parser.getNextWordTransitions(earley_chart[i]);
        for(unsigned int j = 0; j < transition_probs.size(); j++)
            std::cout << transition_probs[j].first << " --> " << transition_probs[j].second << "    ";
        std::cout << std::endl;
    }
}
