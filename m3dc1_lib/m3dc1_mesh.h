#ifndef M3DC1_MESH_H
#define M3DC1_MESH_H

#include <math.h>
#include <string>

#define TOL 1e-4

//Foward class declarations
class m3dc1_stellarator_mesh;
class m3dc1_3d_mesh;

class m3dc1_mesh {
 private:
  int memory_depth;

  int hits, misses;
  
  int last_elm; 
  int** next_elm;

  void clear_memory();

 protected:
  int* nneighbors;
  int** neighbor;

  virtual int max_neighbors() const {return 3;}

  //  static const double tol = 1e-1;

  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    double t = (a[i] + b[i] + c[i])*TOL;
    if(eta + t < 0.) return false;
    if(eta - t > c[i]) return false;
    double x = 1. - eta/c[i];
    if(xi + t < -b[i]*x) return false;
    if(xi - t > a[i]*x) return false;
    return true;
  }
  virtual void global_to_local(const int i, 
			       const double X,const double Phi,const double Z,
			       double* xi, double* zi, double* eta) const
  {
    *xi  =  (X - x[i])*co[i] + (Z - z[i])*sn[i] - b[i];
    *eta = -(X - x[i])*sn[i] + (Z - z[i])*co[i];
  }
  int shared_nodes(const int i, const int j);
  bool elements_are_neighbors(const int i, const int j);

  int in_element_threadsafe(double X, double Phi, double Z, 
			    double* xi=0, double* zi=0, double* eta=0, 
			    int guess=-1);
  int in_element_memory(double X, double Phi, double Z, 
			double* xi=0, double* zi=0, double* eta=0, 
			int guess=-1);

 public:
  int nelms;
  double* a;
  double* b;
  double* c;
  double* co;
  double* sn;
  double* x;
  double* z;
  int* bound;
  int* region;
  double period;
  bool toroidal;
  int nplanes;

 public:
  m3dc1_mesh(int n);
  m3dc1_mesh(const m3dc1_mesh&);
  m3dc1_mesh(const m3dc1_3d_mesh&);
  m3dc1_mesh(const m3dc1_stellarator_mesh&);
  virtual ~m3dc1_mesh(); 
  m3dc1_mesh& operator=(const m3dc1_mesh&);

  bool set_memory_depth(int d);
  inline int get_memory_depth() const {
    return memory_depth;
  }
  inline int get_hits() const {
    return hits;
  }
  inline int get_misses() const {
    return misses;
  }
  inline int get_last_elm() const {
    return last_elm;
  }
  inline int** get_next_elm() const {
    return next_elm;
  }
  inline int get_next_elm(const int i, const int j) const {
    return next_elm[i][j];
  }

  virtual void find_neighbors();

  bool is_in_element(const int i, 
		     const double X, const double Phi, const double Z,
		     double* xi_out=0, double* zi_out=0, double* eta_out=0)
    const
  {
    double xi, zi, eta;
    global_to_local(i, X, Phi, Z, &xi, &zi, &eta);
    if(is_in_element_local(i, xi, zi, eta)) {
      if(xi_out) *xi_out = xi;
      if(zi_out) *zi_out = zi;
      if(eta_out) *eta_out = eta;
      return true;
    }
    return false;
  }
  virtual int in_element(double X, double Phi, double Z, 
			 double* xi=0, double* zi=0, double* eta=0, 
			 int guess=-1);

  virtual void extent(double *X0, double* X1,
		      double *Phi0, double* Phi1,
		      double *Z0, double* Z1) const;

};

class m3dc1_3d_mesh : public m3dc1_mesh {
 protected:
  virtual bool is_in_element_local(const int i, const double xi, 
				   const double zi, const double eta) const
  {
    if(!m3dc1_mesh::is_in_element_local(i, xi, zi, eta)) return false;
    if(zi + TOL*d[i] < 0.) return false;
    if(zi - TOL*d[i] > d[i]) return false;
    return true;
  }
  virtual void global_to_local(const int i, 
			       const double X,const double Phi,const double Z,
			       double* xi, double* zi, double* eta) const
  {
    m3dc1_mesh::global_to_local(i, X, Phi, Z, xi, zi, eta);
    *zi = Phi - phi[i];
  }

  int shared_nodes(const int i, const int j);
  bool elements_are_neighbors(const int i, const int j);
  virtual int max_neighbors() const {return 5;}

 public:
  double *phi;
  double *d;

  virtual void find_neighbors();

  virtual int in_element(double X, double Phi, double Z, 
			 double* xi=0, double* zi=0, double* eta=0, 
			 int guess=-1)
  {
    while(Phi < 0      ) Phi += period;
    while(Phi >= period) Phi -= period;

    return m3dc1_mesh::in_element(X, Phi, Z, xi, zi, eta, guess);
  }


  m3dc1_3d_mesh(int n);
  m3dc1_3d_mesh(m3dc1_mesh*);
  m3dc1_3d_mesh(m3dc1_3d_mesh*);
  m3dc1_3d_mesh(m3dc1_stellarator_mesh*);
  m3dc1_3d_mesh(const m3dc1_mesh&);
  m3dc1_3d_mesh(const m3dc1_3d_mesh&);
  m3dc1_3d_mesh(const m3dc1_stellarator_mesh&);
  virtual ~m3dc1_3d_mesh();
};

#endif
