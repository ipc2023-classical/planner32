#ifndef NUMERIC_DOMINANCE_INT_EPSILON_H
#define NUMERIC_DOMINANCE_INT_EPSILON_H


#include <cassert>

#include <iostream>

class IntEpsilon {
    int value; 
    int epsilon; 

public: 
    IntEpsilon() : value (0), epsilon(0) {}

    IntEpsilon(int value_) :  value(value_), epsilon (0) {   
    }

    IntEpsilon(int value_, int epsilon_) : value(value_), epsilon(epsilon_) {} 

    IntEpsilon (const IntEpsilon &) = default;
    IntEpsilon & operator= ( const IntEpsilon & ) = default;	
    IntEpsilon (IntEpsilon &&) = default;
    IntEpsilon & operator= ( IntEpsilon && ) = default;	
    ~IntEpsilon () = default;

    bool operator< (const IntEpsilon & other) const ;
    bool operator<=(const IntEpsilon & other) const ;
    bool operator==(const IntEpsilon & other) const ;
    bool operator!=(const IntEpsilon & other) const ;
    bool operator> (const IntEpsilon & other) const ;
    bool operator>=(const IntEpsilon & other) const ;

    IntEpsilon & operator+=(const IntEpsilon & other) ;
    IntEpsilon & operator+=(const int & other) ;
    IntEpsilon & operator-=(const int & other) ;
    IntEpsilon & operator-=(const IntEpsilon & other) ;
    IntEpsilon operator -() const;	

    template <typename T> 
    const IntEpsilon operator-(const T &other) const {
	return IntEpsilon(*this) -= other;
    }

    template <typename T> 
    const IntEpsilon operator+(const T &other) const {
	return IntEpsilon(*this) += other;
    }

    int get_value() const {
	return value;
    }

    int get_epsilon() const {
	assert(epsilon >= -1 && epsilon <= 1);
	return epsilon;
    } 

    friend std::ostream &operator<<(std::ostream &os, const IntEpsilon & other);
};


template <typename T> inline void set_epsilon(T & t) {
    t = 0;
}

template <> inline void set_epsilon(IntEpsilon & t) {
    t = IntEpsilon(0, 1);
}


#endif
