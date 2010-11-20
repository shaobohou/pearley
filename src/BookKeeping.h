#ifndef BOOKKEEPING
#define BOOKKEEPING

class ParseEntry
{
    public:
        ParseEntry()
        : the_col(static_cast<unsigned int>(-1)), the_row(static_cast<unsigned int>(-1))
        {}

        ParseEntry(unsigned int col, unsigned int row)
        : the_col(col), the_row(row)
        {}

        unsigned int& row()
        {
            return the_row;
        }

        const unsigned int& col() const
        {
            return the_col;
        }

        const unsigned int& row() const
        {
            return the_row;
        }

    private:
        unsigned int the_col, the_row;
};

class Progenitor
{
    public:
        Progenitor(const ParseEntry &incompleted_index, double weight)
        : the_incompleted_index(incompleted_index), the_weight(weight)
        {}

        Progenitor(const ParseEntry &incompleted_index, const ParseEntry &complete_index)
        : the_incompleted_index(incompleted_index), the_complete_index(complete_index), the_weight(0.0)
        {}

        Progenitor(const ParseEntry &incompleted_index, const ParseEntry &complete_index, double weight)
        : the_incompleted_index(incompleted_index), the_complete_index(complete_index), the_weight(weight)
        {}

        ParseEntry& complete_index()
        {
            return the_complete_index;
        }

        const ParseEntry& incompleted_index() const
        {
            return the_incompleted_index;
        }

        const ParseEntry& complete_index() const
        {
            return the_complete_index;
        }

        const double& weight() const
        {
            return the_weight;
        }

    private:
        ParseEntry the_incompleted_index, the_complete_index;
        double the_weight;
};

typedef std::vector<std::vector<Progenitor> > ProgenitorsLists;

#endif
