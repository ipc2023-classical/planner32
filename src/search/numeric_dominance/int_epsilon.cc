#include "int_epsilon.h"

#include <cassert>
using namespace std;


bool IntEpsilon::operator< (const IntEpsilon & other) const {
    return value < other.value || 
		   (value == other.value &&
		    epsilon < other.epsilon);
}

bool IntEpsilon::operator<=(const IntEpsilon & other) const {
    return value < other.value || 
		   (value == other.value &&
		    epsilon <= other.epsilon);
}


bool IntEpsilon::operator> (const IntEpsilon & other) const {
    return value > other.value || 
		   (value == other.value &&
		    epsilon > other.epsilon);
}

bool IntEpsilon::operator>=(const IntEpsilon & other) const {
    return value > other.value || 
		   (value == other.value &&
		    epsilon >= other.epsilon);
}

bool IntEpsilon::operator==(const IntEpsilon & other) const {
    return value == other.value &&
	epsilon == other.epsilon;
}

bool IntEpsilon::operator!=(const IntEpsilon & other) const {
    return value != other.value || epsilon != other.epsilon;
}

IntEpsilon & IntEpsilon::operator+=(const IntEpsilon & other) {
    assert (epsilon*other.epsilon >= 0); //Check that both have the same sign    
    value += other.value;
    if (epsilon == 0) { //Epsilon cannot be cancelled out
	epsilon += other.epsilon;
    }
    assert(epsilon >= -1 && epsilon <= 1);
    return *this;
}

IntEpsilon & IntEpsilon::operator-=(const IntEpsilon & other) {
    assert (epsilon*other.epsilon >= 0); //Check that both have the same sign    
    value -= other.value;
    if (epsilon == 0) { //Epsilon cannot be cancelled out
	epsilon -= other.epsilon;
    }
    assert(epsilon >= -1 && epsilon <= 1);
    return *this;
}


IntEpsilon & IntEpsilon::operator+=(const int & other_value) {
    value += other_value;
    return *this;
}

IntEpsilon & IntEpsilon::operator-=(const int & other_value) {
    value -= other_value;
    return *this;
}
  
IntEpsilon IntEpsilon::operator -() const {
    return IntEpsilon(-value, -epsilon);
}

ostream &operator<<(ostream &os, const IntEpsilon & other) {
    os << other.value; 
    if(other.epsilon > 0) {
	assert(other.epsilon == 1);
	os << "(e)";	
    }else if(other.epsilon < 0) {
	assert(other.epsilon == -1);
	os << "(-e)";	
    }
    return os; 
}
