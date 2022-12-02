#include "m3dc1_stellarator.h"

/*
extern "C" {
    void LogicalFromCylindrical(double*, double*, std::string);
    void CylindricalFromLogical(double*, double*, std::string);
    void VmecFromCylindrical(double *, double *, std::string);
    void inverseTransform(double*, double*, std::string, char*);
}
*/

m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(const int n,
                        m3dc1_stellarator_mesh::eq_type_enum eq_type,
                        std::string eq_file) :
                        m3dc1_3d_mesh(n), eq_type(eq_type), eq_file(eq_file)
{
  rst = new double[n];
  zst = new double[n];
}

m3dc1_mesh::m3dc1_mesh(m3dc1_stellarator_mesh* m)
{
    m3dc1_mesh((m3dc1_mesh*)m);
}

m3dc1_mesh::m3dc1_mesh(const m3dc1_stellarator_mesh& m)
{
    nelms = m.nelms;

    allocate_memory(nelms);
    set_memory_depth(m.get_memory_depth());

    std::copy(m.a, m.a + nelms, a);
    std::copy(m.b, m.b + nelms, b);
    std::copy(m.c, m.c + nelms, c);
    std::copy(m.co, m.co + nelms, co);
    std::copy(m.sn, m.sn + nelms, sn);
    std::copy(m.x, m.x + nelms, x);
    std::copy(m.z, m.z + nelms, z);
    std::copy(m.bound, m.bound + nelms, bound);
    std::copy(m.region, m.region + nelms, region);
    period = m.period;
    toroidal = m.toroidal;
    nplanes = m.nplanes;
    last_elm = m.get_last_elm();
    hits = m.get_hits();
    misses = m.get_misses();
    std::copy(m.nneighbors, m.nneighbors + nelms, nneighbors);
    for(int i=0; i<nelms; i++)
    {
        neighbor[i] = new int[max_neighbors()];
        std::copy(m.neighbor[i], m.neighbor[i] + max_neighbors(), neighbor[i]);
    }
    for(int d=0; d<max_neighbors(); d++)
        for(int i=0; i<nelms; i++)
            next_elm[d][i] = m.get_next_elm(d, i);
    
}

m3dc1_3d_mesh::m3dc1_3d_mesh(const m3dc1_stellarator_mesh& m) : m3dc1_mesh(m)
{
    d = new double[nelms];
    phi = new double[nelms];

    std::copy(m.d, m.d + nelms, d);
    std::copy(m.phi, m.phi + nelms, phi);
}

m3dc1_3d_mesh::m3dc1_3d_mesh(m3dc1_stellarator_mesh* m) : m3dc1_mesh(m)
{
    d = m->d;
    phi = m->phi;
}


m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(const m3dc1_mesh& m,
                        m3dc1_stellarator_mesh::eq_type_enum eq_type,
                        std::string eq_file) : 
                        m3dc1_3d_mesh(m), eq_type(eq_type), eq_file(eq_file)
{
    rst = new double[nelms];
    zst = new double[nelms];
}


m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(m3dc1_mesh* m,
                        m3dc1_stellarator_mesh::eq_type_enum eq_type,
                        std::string eq_file) :
                        m3dc1_3d_mesh(m), eq_type(eq_type), eq_file(eq_file)
{
    rst = new double[nelms];
    zst = new double[nelms];
}

m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(const m3dc1_3d_mesh& m,
                        m3dc1_stellarator_mesh::eq_type_enum eq_type,
                        std::string eq_file) : 
                        m3dc1_3d_mesh(m), eq_type(eq_type), eq_file(eq_file)
{
    rst = new double[nelms];
    zst = new double[nelms];
}

m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(m3dc1_3d_mesh* m,
                        m3dc1_stellarator_mesh::eq_type_enum eq_type,
                        std::string eq_file) : 
                        m3dc1_3d_mesh(m), eq_type(eq_type), eq_file(eq_file)
{
    rst = new double[nelms];
    zst = new double[nelms];
}

m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(m3dc1_stellarator_mesh* m) : m3dc1_3d_mesh(m)
{
    rst = m->rst;
    zst = m->zst;
    eq_type = m->eq_type;
    eq_file = m->eq_file;
}

m3dc1_stellarator_mesh::m3dc1_stellarator_mesh(const m3dc1_stellarator_mesh& m) : m3dc1_3d_mesh(m)
{
    rst = new double[nelms];
    zst = new double[nelms];

    std::copy(m.rst, m.rst + nelms, rst);
    std::copy(m.zst, m.zst + nelms, zst);
    eq_type = m.eq_type;
    eq_file = m.eq_file;
}

m3dc1_stellarator_mesh::~m3dc1_stellarator_mesh()
{
    delete[] rst;
    delete[] zst;
}

void m3dc1_stellarator_mesh::global_to_logical(const double Rst, 
                             const double Phi, const double Zst,
                             double* xl, double* zl) const
{
    double logical_coords[3] = {0.0};
    double cyl_coords [3] = {Rst, Phi, Zst};
    //LogicalFromCylindrical(logical_coords, cyl_coords, eq_file);
    *xl = logical_coords[0];
    *zl = logical_coords[2];
}

void m3dc1_stellarator_mesh::logical_to_global(const double xl,
                             const double phi, const double zl,
                             double* Rst, double* Zst) const
{
    double s, theta;
    double cyl_coords [3] = {0.0};

    theta = std::atan2(zl, xl);
    s = std::pow(xl, 2.0) + std::pow(zl, 2.0);
    double logical_coords [3] = {s, theta, phi};
    //CylindricalFromLogical(cyl_coords, logical_coords, eq_file);
    *Rst = cyl_coords[0];
    *Zst = cyl_coords[2];
}

void m3dc1_stellarator_mesh::logical_to_local(int i, const double xl,
                             const double phi, const double zl,
                             double* xi, double* zi, double* eta) const
{
    m3dc1_3d_mesh::global_to_local(i, xl, phi, zl, xi, zi, eta);
}

void m3dc1_stellarator_mesh::global_to_local(const double Rst,
                             const double Phi, const double Zst,
                             double* xi, double* zi, double* eta)
{
    double xl, zl;
    m3dc1_stellarator_mesh::global_to_logical(Rst, Phi, Zst, &xl, &zl);
    int i = m3dc1_3d_mesh::in_element(xl, Phi, zl);
    m3dc1_3d_mesh::global_to_local(i, xl, Phi, zl, xi, zi, eta);
}

m3dc1_3d_field::m3dc1_3d_field(m3dc1_stellarator_mesh* m)
{ 
    mesh = (m3dc1_3d_mesh*)m;
    data = new double[mesh->nelms*nbasis*tbasis];
}

m3dc1_3d_field::m3dc1_3d_field(const m3dc1_stellarator_mesh& m)
{
    mesh = new m3dc1_3d_mesh(m);
    data = new double[mesh->nelms*nbasis*tbasis];
}

m3dc1_stellarator_field::m3dc1_stellarator_field(const m3dc1_mesh& m,
                                                 m3dc1_stellarator_mesh::eq_type_enum eq_type,
                                                 std::string eq_file) :
                                                 m3dc1_3d_field()
{
    stell_mesh = new m3dc1_stellarator_mesh(m, eq_type, eq_file);
    data = new double[stell_mesh->nelms*nbasis*tbasis];
}


m3dc1_stellarator_field::m3dc1_stellarator_field(m3dc1_stellarator_mesh* m) : m3dc1_3d_field()
{
    stell_mesh = m;
    data = new double[stell_mesh->nelms*nbasis*tbasis];
}

m3dc1_stellarator_field::~m3dc1_stellarator_field()
{
    delete stell_mesh;
    delete[] data;
}

/*
bool m3dc1_stellarator_field::eval(const double rst, const double phi, const double zst,
			                  const m3dc1_get_op op, double* val, int* element)
{

  for(int q=0; q<tbasis; q++) {
    for(int p=0; p<nbasis; p++) {      
      if(getval) {
	temp = data[j]*pow(xi,mi[p])*pow(eta,ni[p]);
	v[OP_1] = temp;
      }
      
      if((op & GET_DVAL) == GET_DVAL) {
	temp = data[j]*mi[p]*pow(xi,mi[p]-1)*pow(eta,ni[p]);
	v[OP_DR] = mesh->co[e]*temp;
	v[OP_DZ] = mesh->sn[e]*temp;
	
	temp = data[j]*ni[p]*pow(xi,mi[p])*pow(eta,ni[p]-1);
	v[OP_DR] -= mesh->sn[e]*temp;
	v[OP_DZ] += mesh->co[e]*temp;	
      }

      if((op & GET_DDVAL) == GET_DDVAL) {
	temp = data[j]*(mi[p]-1)*mi[p]*pow(xi,mi[p]-2)*pow(eta,ni[p]);
	v[OP_DRR] =  co2*temp;
	v[OP_DRZ] = cosn*temp;
	v[OP_DZZ] =  sn2*temp;
	
	temp = data[j]*mi[p]*ni[p]*pow(xi,mi[p]-1)*pow(eta,ni[p]-1);
	v[OP_DRR] -=     2.*cosn*temp;
	v[OP_DRZ] += (co2 - sn2)*temp;
	v[OP_DZZ] +=     2.*cosn*temp;
	
	temp = data[j]*(ni[p]-1)*ni[p]*pow(xi,mi[p])*pow(eta,ni[p]-2);
	v[OP_DRR] +=  sn2*temp;
	v[OP_DRZ] -= cosn*temp;
	v[OP_DZZ] +=  co2*temp;
      }

      temp = pow(zi, li[q]);
      if((op & GET_VAL) == GET_VAL) {
	val[OP_1  ] += v[OP_1  ]*temp;
      }
      if((op & GET_DVAL) == GET_DVAL) {
	val[OP_DR ] += v[OP_DR ]*temp;
	val[OP_DZ ] += v[OP_DZ ]*temp;
      }
      if((op & GET_DDVAL) == GET_DDVAL) {
	val[OP_DRR] += v[OP_DRR]*temp;
	val[OP_DRZ] += v[OP_DRZ]*temp;
	val[OP_DZZ] += v[OP_DZZ]*temp;
      }

      if((op & GET_PVAL) == GET_PVAL) {
	temp = pow(zi, li[q]-1)*li[q];
	val[OP_DP  ] += v[OP_1  ]*temp;
	if((op & GET_DVAL) == GET_DVAL) {
	  val[OP_DRP ] += v[OP_DR ]*temp;
	  val[OP_DZP ] += v[OP_DZ ]*temp;
	}
	if((op & GET_DDVAL) == GET_DDVAL) {
	  val[OP_DRRP] += v[OP_DRR]*temp;
	  val[OP_DRZP] += v[OP_DRZ]*temp;
	  val[OP_DZZP] += v[OP_DZZ]*temp;
	}
      }

      if((op & GET_PPVAL) == GET_PPVAL) {
	temp = pow(zi, li[q]-2)*li[q]*(li[q]-1);
	val[OP_DPP  ] += v[OP_1  ]*temp;
	if((op & GET_DVAL) == GET_DVAL) {
	  val[OP_DRPP ] += v[OP_DR ]*temp;
	  val[OP_DZPP ] += v[OP_DZ ]*temp;
	}
	if((op & GET_DDVAL) == GET_DDVAL) {
	  val[OP_DRRPP] += v[OP_DRR]*temp;
	  val[OP_DRZPP] += v[OP_DRZ]*temp;
	  val[OP_DZZPP] += v[OP_DZZ]*temp;
	}
      }

      j++;
    }
  }
  return true;
}

void m3dc1_stellarator_field::load_field(char* field_id)
{
  double local_data[OP_NUM];
  for(int i=0; i<OP_NUM; i++) val[i] = 0.;

  double s, theta, dsdx, dsdy, dtdx, dtdy;
  double dRds, dRdt, dRdz, dZds, dZdt, dZdz;
  double internal_coords[3], cyl_coords[3], logical_coords[3];
  char no_op[] = ":none";
  char deriv_op[] = ":ds";
  // Evaluate chain rules for equilibrium to logical to global values
  int j = e*nbasis*tbasis;
  for(int i=0; i<mesh->nelms; i++) {
    logical_coords[0] = stell_mesh->x[i];
    logical_coords[1] = stell_mesh->phi[i];
    logical_coords[2] = stell_mesh->z[i];
    cyl_coords[0] = stell_mesh->rst[i];
    cyl_coords[1]
    VmecFromCylindrical(internal_coords, cyl_coords, stell_mesh->eq_file);

    val[OP_1] = logical_val[OP_1];
    deriv_op[2] = 's';
    inverseTransform(&dRds, internal_coords, stell_mesh->eq_file, deriv_op);
  }
}
*/