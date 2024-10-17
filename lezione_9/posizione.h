#ifndef __posizione_h__
#define __posizione_h__

class posizione{

	public:

        posizione();
        posizione(double x, double y, double z);

        //methods

        void AddPos(double,double,double);
        void SetPos(posizione &);
        void SetPos(double,double,double);
        void SetPosPolar(double,double,double);
        void SetX(double);
        void SetY(double);
        void SetZ(double);

        double GetX() const;
        double GetY() const;
        double GetZ() const;
        double GetR() const;
        double GetPhi() const;
        double GetTheta() const;
        double GetRho() const;
        double Distanza(const posizione&) const;

	private:

	    double m_x, m_y, m_z;


};



#endif