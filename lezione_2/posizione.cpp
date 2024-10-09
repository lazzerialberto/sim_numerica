#include "posizione.h"

#include <iostream>
#include <cmath>

using namespace std;

posizione::posizione(){
	m_x=0;
	m_y=0;
	m_z=0;
}

posizione::posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
}


double posizione::GetX() const {
	return m_x;
}

double posizione::GetY() const{
	return m_y;
}

double posizione::GetZ() const{
	return m_z;
}

double posizione::GetR() const {
	return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}

double posizione::GetPhi() const {
	return atan2(m_y,m_x);
}

double posizione::GetTheta() const {
	return acos(m_z/GetR());
}

double posizione::GetRho() const {
	return sqrt(m_x*m_x+m_y*m_y);
}

double posizione::Distanza(const posizione& a) const {
	return sqrt(pow(GetX()-a.GetX(),2) + pow(GetY()-a.GetY(),2) + pow(GetZ()-a.GetZ(),2));
}
