// GPU N-body utilities, 2010/7/3
// Keigo Nitadori, keigo@riken.jp

struct DblFloat{
	float hi, lo;

	__device__ __host__
	DblFloat(){}

	__host__
	DblFloat(const double x){
		hi = float(x);
		lo = float(x - double(hi));
	}

	__device__
	DblFloat(const float _hi, const float _lo) : hi(_hi), lo(_lo) {}

	__host__
	operator double() const{
		return double(hi) + double(lo);
	}

	__device__
	operator float() const{
		return hi + lo;
	}

	__device__
	DblFloat operator - (const DblFloat &rhs) const {
		return DblFloat(hi-rhs.hi, lo-rhs.lo);
	}

	__device__
	void operator += (const float f){
		const float tmp = hi + f;
		lo += (f - (tmp - hi)); // tmp == hi + f + eps, 
		                        // eps = (tmp -hi) - f
		                        // lo -= eps;
		hi = tmp;
	}

	__device__
	void operator += (const DblFloat &d){
		operator += (d.hi);
		operator += (d.lo);
	}

	__device__
	void regularize(){
		const float tmp = hi + lo;
		lo =  lo - (tmp - hi); // tmp = hi + lo + eps, 
		                       // eps = (tmp - hi) - lo, 
							   // lo = -eps.
		hi = tmp;
	}
};

template <typename real>
struct Gvec3{
	real x, y, z;

	__device__ __host__
	Gvec3(){}

	__device__ __host__
	Gvec3(const real &r) : x(r), y(r), z(r) {}

	__device__ __host__
	Gvec3(const real &_x, const real &_y, const real &_z) : 
		x(_x), y(_y), z(_z) {}

	template <typename Real>
	__host__
	Gvec3(const Real &_x, const Real &_y, const Real &_z) : 
		x(real(_x)), y(real(_y)), z(real(_z)) {}

	template <typename Real>
	__host__
	Gvec3(const Real *p) : 
		x(real(p[0])), y(real(p[1])), z(real(p[2])) {}

	template <typename Real>
	__host__
	void write(Real *p) const
	{
		p[0] = Real(x);
		p[1] = Real(y);
		p[2] = Real(z);
	}

	template <typename Real> 
	__host__ __device__
	operator Gvec3<Real> () const {
		return Gvec3<Real> (Real(x), Real(y), Real(z));
	}

	__device__
	Gvec3 operator + (const Gvec3 &rhs) const {
		return Gvec3(x+rhs.x, y+rhs.y, z+rhs.z);
	}

	__device__
	Gvec3 operator - (const Gvec3 &rhs) const {
		return Gvec3(x-rhs.x, y-rhs.y, z-rhs.z);
	}

	__device__
	real operator * (const Gvec3 &rhs) const {
		return (x*rhs.x + y*rhs.y + z*rhs.z);
	}

	__device__
	Gvec3 operator * (const real &s) const {
		return Gvec3(s*x, s*y, s*z);
	}

	__device__
	friend Gvec3 operator * (const real &s, const Gvec3 &vec){
		// return s * vec;
		return Gvec3(s*vec.x, s*vec.y, s*vec.z);
	}

	template <typename Real> 
	__device__
	void operator += (const Gvec3<Real> &rhs){
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
	}

#if 0
	template <typename Real> 
	__device__
	static Gvec3<Real> diff(const Gvec3 &lhs, const Gvec3 &rhs) 
	{
		return Gvec3<Real> (lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
	}
#endif

};

#if 0
__device__ inline
Gvec3<float> operator - (
		const Gvec3<DblFloat> &lhs, 
		const Gvec3<DblFloat> &rhs){
	return Gvec3<float> (lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}	
#endif
