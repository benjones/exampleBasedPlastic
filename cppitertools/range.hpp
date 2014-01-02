#ifndef __RANGE__H__
#define __RANGE__H__

// range() for range-based loops with start, stop, and step.
//
// Acceptable forms are:
// for (auto i : range(stop)) { ... } // start = 0, step = 1
// for (auto i : range(start, stop)) { ... } // step = 1
// for (auto i : range(start, stop, step)) { ... }
//
// The start may be greater than the stop if the range is negative
// The range will effectively be empty if:
//   1) step is positive and start > stop
//   2) step is negative and start < stop
//
// If a step of 0 is provided, a RangeException will be thrown

#include <exception>

namespace iter {

    // Thrown when step 0 occurs
    class RangeException : public std::exception {
        virtual const char *what() const throw() { 
            return "Range() step must be non-zero";
        }
    };

    //Forward declarations of Enumerable and enumerate
    template <typename T>
    class Range;

    template <typename T>
    Range<T> range(T);
    template <typename T>
    Range<T> range(T, T);
    template <typename T>
    Range<T> range(T, T, T);

    template <typename T>
    class Range {
        friend Range range<T>(T);
        friend Range range<T>(T, T);
        friend Range range<T>(T, T, T);
        private:
            const T start; 
            const T stop;
            const T step;

            Range(T stop) :
                start(0),
                stop(stop),
                step(1)
            { }

            Range(T start, T stop, T step=1) :
                start(start),
                stop(stop),
                step(step)
            { }

        public:
            Range() = delete;
            Range(const Range &) = default;
            class Iterator {
                private:
                    T value;
                    const T step;
                public:
                    Iterator(T val, T step) :
                        value(val),
                        step(step)
                    { }

                    T operator*() const {
                        return this->value;
                    }

                    Iterator & operator++() {
                        this->value += this->step;
                        return *this;
                    }

                    // This operator would more accurately read as "in bounds"
                    // or "incomplete" because exact comparison with the end
                    // isn't good enough for the purposes of this Iterator.
                    // There are two odd cases that need to be handled
                    //
                    // 1) The Range is infinite, such as
                    // Range (-1, 0, -1) which would go forever down toward
                    // infinitely (theoretically).  If this occurs, the Range
                    // will instead effectively be empty
                    //
                    // 2) (stop - start) % step != 0.  For
                    // example Range(1, 10, 2).  The iterator will never be
                    // exactly equal to the stop value.
                    bool operator!=(const Range::Iterator & other) const {
                        return !(this->step > 0 && this->value >= other.value) 
                            && !(this->step < 0 && this->value <= other.value);
                    }
            };

            Iterator begin() const {
                return Iterator(start, step);
            }

            Iterator end() const { 
                return Iterator(stop, step);
            }
    };


    template <typename T>
    Range<T> range(T stop) {
        return Range<T>(stop);
    }

    template <typename T>
    Range<T> range(T start, T stop) {
        return Range<T>(start, stop);
    }

    template <typename T>
    Range<T> range(T start, T stop, T step) {
        if (step == 0) {
            throw RangeException();
        }
        return Range<T>(start, stop, step);
    }
}

#endif //ifndef __RANGE__H__
