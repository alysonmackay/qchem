#include "orb_transform.h"

using namespace std;

int spatial(int k) {return k / 2; }
bool same_spin(int k, int l) {return k%2 == l%2; } // spin parity

vector<double> ao2mo(const vector<double>& TEI_AO, const Matrix& coeff, int nao) {
	// transform TEI_AO to MO basis 
	int M = nao*(nao+1)/2;
	Matrix TMP(M,M); // buffer for half-transformed integrals
	Matrix dense(nao,nao);
	int size = TEI_AO.size();
	vector<double> TEI_MO(size);
	
	for(int p=0, pq=0; p < nao; p++) {
		for(int q=0; q <= p; q++, pq++) {
			// unpack rs into full dense matrix 
			for(int r=0; r < nao; r++) 
				for(int s=0; s < nao; s++)
					dense(r,s) = TEI_AO[compound(pq, compound(r,s))];
			Matrix X = coeff.transpose() * dense * coeff; // r->k and s->l
			
			// repack: only lower tri of X is unique
			for(int k=0, kl=0; k < nao; k++)
				for(int l=0; l <= k; l++, kl++)
					TMP(pq,kl) = X(k,l);
		}
	}

	for(int k=0, kl=0; k < nao; k++) {
		for(int l=0; l <= k; l++, kl++) {
			// unpack pq into full dense matrix 
			for(int p=0; p < nao; p++) 
				for(int q=0; q < nao; q++) 
					dense(p,q) = TMP(compound(p,q),kl);
			Matrix X = coeff.transpose() * dense * coeff; // p->i and q->j
			
			for(int i=0, ij=0; i < nao && ij <= kl; i++) // repack
				for(int j=0; j <= i && ij <= kl; j++, ij++) {
					int ijkl = compound(ij,kl);
					TEI_MO[ijkl] = X(i,j);
				}
		}
	}
	
	return TEI_MO;
}

// brute force N^8 version of algorithm; keep for testing
vector<double> ao2mo_brute(const vector<double>& TEI_AO, const Matrix& coeff, int nao) {
	int size = TEI_AO.size();
	vector<double> TEI_MO(size);
	int i, j, k, l, ijkl; // MO indices 
	int p, q, r, s, pq, rs, pqrs; // AO indices
	for(i=0, ijkl=0; i < nao; i++) {
		for(j=0; j <= i; j++) {
			for(k=0; k <= i; k++) {
				for(l=0; l <= (i==k ? j : k); l++, ijkl++) { 

					for(p=0; p < nao; p++) {
						for(q=0; q < nao; q++) { 
							pq = compound(p,q);
							for(r=0; r < nao; r++) {
								for(s=0; s < nao; s++) {
									rs = compound(r,s); 
									pqrs = compound(pq, rs); 
									
									TEI_MO[ijkl] += coeff(p,i) * coeff(q,j) * coeff(r,k) * coeff(s,l) * TEI_AO[pqrs];
								}
							}
						}
					}
				}
			}
		}
	}
	return TEI_MO;
}

Tensor4 mo2so(const vector<double>& TEI_MO, int nao) {
	int nso = 2 * nao; // spin orbitals 
	Tensor4 TEI_SO(nso, Vec3(nso, Vec2(nso, Vec1(nso, 0.0))));

	for(int p=0; p < nso; p++)
		for(int q=0; q < nso; q++)
			for(int r=0; r < nso; r++) 
				for(int s=0; s < nso; s++) {
					int pr = compound(spatial(p),spatial(r));
					int qs = compound(spatial(q),spatial(s));
					int pqrs = compound(pr,qs);
					double direct = TEI_MO[pqrs] * same_spin(p,r) * same_spin(q,s); // Dirac <pq|rs>
					int ps = compound(spatial(p),spatial(s));
					int qr = compound(spatial(q),spatial(r));
					int psqr = compound(ps,qr);
					double exchange = TEI_MO[psqr] * same_spin(p,s) * same_spin(q,r); // <pq|sr>
					TEI_SO[p][q][r][s] = direct - exchange; // antisymmetrized <pq||rs> 
				}
	return TEI_SO;
}
