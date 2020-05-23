#include "volume_helper.hpp"


Polytope* Polytope_new_box(int n, FT r) {
   Polytope* p = Polytope_new(n, 2*n);

   for(int i=0; i<n; i++) {// for each dim
      Polytope_set_b(p, i,   r);
      Polytope_set_b(p, i+n, r);
      for(int x=0; x<n; x++) {
         Polytope_set_a(p, i,   x, (x==i)?1:0);
         Polytope_set_a(p, i+n, x, (x==i)?-1:0);
      }
   }

   return p;
}

PolytopeT* PolytopeT_new_box(int n, FT r) {
   PolytopeT* p = PolytopeT_new(n, 2*n);

   for(int i=0; i<n; i++) {// for each dim
      PolytopeT_set_b(p, i,   r);
      PolytopeT_set_b(p, i+n, r);
      for(int x=0; x<n; x++) {
         PolytopeT_set_a(p, i,   x, (x==i)?1:0);
         PolytopeT_set_a(p, i+n, x, (x==i)?-1:0);
      }
   }
   PolytopeT_fix_inv(p);
   return p;
}

void make_random_poly(const std::vector<double> &ell, int m, Polytope **ret){

  // create n random hyperplanes around ellipsoid ell
  int n = ell.size();
  *ret = Polytope_new(n, m);
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0,1);

  
  for (int i = 0; i < m; i++){

    // sample uar on unit dim-sphere
    std::vector<double> pnt(n);
    double rad = 0;

    for (int j = 0; j < n; j++) {
      double number = distribution(generator);
      rad += number * number;
      pnt[j] = number;
    }
    rad = sqrt(rad);

    
    // normalize point (on surface) and map point to ellipsoid
    // (note it's no longer uniform at that point but maybe good enough, we could make it uniform by discarding each point with a probability proportional to the distortion factor of the matrix inducing the ellipsoid at that point)
    for (int j = 0; j < n; j++) {
        pnt[j] *= ell[j]/rad;
    }

    // add the tangent plane of ellispoid at pnt to poly
    // note that gradient of the ellipsoid at x is [2*ell_i x_i]_{i = 1 to dim}
    double sum = 0;
    for (int j = 0; j < n; j++){
        Polytope_set_a(*ret, i, j, 2*pnt[j] / (ell[j] * ell[j]));
        sum += Polytope_get_a(*ret, i, j) * pnt[j];
    }
    Polytope_set_b(*ret, i, sum);    
  }

}



void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


int read_polyvest_p(string filename, Polytope **P){

    ifstream file;
    file.open(filename);


    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    int n, m;
    file >> m >> n;

    *P = Polytope_new(n, m);

    FT num;
    for (int i = 0; i < m; i++){
        file >> num;
        Polytope_set_b(*P, i, num);
        for (int j = 0; j < n; j++){
            file >> num;
            Polytope_set_a(*P, i, j, num);
        }
    }

    return 0;
}


// copied from vinci
FT sread_rational_value (char *s);
FT sread_rational_value (char *s) {
    
   char *numerator_s, *denominator_s = NULL, *position, token;
   int sign = 1, i;
   FT numerator, denominator;

   /* determine the sign of the number */
   numerator_s = s;
   if (s [0] == '-')
   {  sign = -1;
      numerator_s++;
   }
   else if (s [0] == '+')
      numerator_s++;

   /* look for a sign '/' and in this case split the number in numerator and denominator */
   position = strchr (numerator_s, '/');
   if (position != NULL)
   {  *position = '\0'; /* terminates the numerator */
      denominator_s = position + 1;
   };

   /* determine the floating point values of numerator and denominator */
   numerator = 0;
   for (i = 0; i < strlen (numerator_s); i++)
   {  token = numerator_s [i];
      if (strchr ("0123456789", token)) /* token is a cypher */
         numerator = 10 * numerator + (int) token - 48;
   }

   if (position != NULL)
   {  denominator = 0;
      for (i = 0; i < strlen (denominator_s); i++)
      {  token = denominator_s [i];
         if (strchr ("0123456789", token)) /* token is a cypher */
            denominator = 10 * denominator + (int) token - 48;
      }
   }
   else denominator = 1;

   return sign * numerator / denominator;
}


FT read_vinci_nr(string in, string type){

    istringstream instr(in);
    if (!type.compare("integer")){
        FT d;
        instr >> d;
        return d;
    }
    else if (!type.compare("real")){
        FT f;
        instr >> f;
        return f;
    }
    else if (!type.compare("rational")){
        
        std::regex rgx("(.*)/(.*)");
        std::smatch matches;

        if(std::regex_search(in, matches, rgx)) {

            FT num, den;
            istringstream ns(matches[1]);
            istringstream ds(matches[2]);
            ns >> num;
            ds >> den;
            //cout << "print " << num/den << " for " << in << "\n";
            return num/den;
        }
        else {
            FT d;
            istringstream s(in);
            s >> d;
            return d;
        }
    }

    assert(0 && "not a valid number");
    return 0.0;

    

}


FT read_vinci(string filename, Polytope **P, FT *vol){
    ifstream file;
    file.open(filename);

    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    std::string line, type, in;
    std::string b = "begin";
    bool found = false;
    while (std::getline(file, line)){
        if (!b.compare(line)) {
            found = true;
            break;
        }
    }

    if (!found){
        printf("no begin found in vinci file\n");
        return 1;
    }


    int m, n;
    file >> m >> n >> type;
    n--;

    *P = Polytope_new(n,m);
                      
    for (int i = 0; i < m; i++){
        file >> in;
        Polytope_set_b(*P, i, read_vinci_nr(in, type));
        for (int j = 0; j < n; j++){
            file >> in;
            Polytope_set_a(*P, i, j, read_vinci_nr(in, type));            
        }
    }

    // extract solved value
    while (std::getline(file, line)){
        std::regex rgx("Volume:(.*)");
        std::smatch matches;

        if(std::regex_search(line, matches, rgx)) {

            istringstream os(matches[0].str());
            string s;
            os >> s >> *vol;
            return 0;
        }
    }

    *vol = 0;
    return 0;
    
}

int polytope_count_pairs(Polytope* p) {
    const int n= p->n;
    const int m= p->m;

    int count = 0;
    int nz = 0;
    for(int i=0;i<n;i++) {
       double last = 0;
       for(int j=0;j<m;j++) {
          double aij = Polytope_get_a(p,j,i);
	  if(aij!=0) {
             nz++;
	     if(last!=0 && aij*last > 0) {count++;}
	  }
	  if(last!=0) {last =0;} else {last = aij;}
       }
    }
    std::cout << "NZ: " << nz << ", pairs: "<< count << "\n";
    return count;
}

struct OElement {
   int u,v;
   int overlap;
};

Polytope* optimize_polytope(Polytope *p) {
   polytope_count_pairs(p);
   Polytope* q = (Polytope*)Polytope_T.clone(p);
   
   // generator element for each pair of constraints:
   const int n= p->n;
   const int m= p->m;
   std::vector<OElement> oel;
   for(int u=0;u<m;u++) {
      for(int v=0;v<u;v++) {
         int overlap = 0;
	 for(int i=0;i<n;i++) {
	    FT ui = Polytope_get_a(p,u,i);
	    FT vi = Polytope_get_a(p,v,i);
	    if(ui*vi>0) {overlap++;}
	 }
	 if(overlap > 0) {oel.push_back({u,v,overlap});}
      }
   }
   std::sort(oel.begin(),oel.end(), [](const OElement &a,const OElement &b){return a.overlap > b.overlap;});
   //std::cout << "oel " << oel.size() <<  "\n";
   
   std::vector<int> neighbor_1(m,-1);
   std::vector<int> neighbor_2(m,-1);
   evp::Union_find uf(m);
   for(const OElement &el : oel) {
      // check if even possible:
      if(neighbor_2[el.u]!=-1 || neighbor_2[el.v]!=-1) {continue;}
      // check for cycles:
      if(uf.find(el.u) == uf.find(el.v)) {continue;}
      uf.Union(el.u,el.v);

      // allocate u:
      if(neighbor_1[el.u]==-1) {
         neighbor_1[el.u] = el.v;
      } else {
         neighbor_2[el.u] = el.v;
      }
      // allocate v:
      if(neighbor_1[el.v]==-1) {
         neighbor_1[el.v] = el.u;
      } else {
         neighbor_2[el.v] = el.u;
      }
   }
   
   ///  // print check:
   ///  for(int i=0;i<m;i++) {
   ///     std::cout << "n: " << i << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
   ///  }

   // go traverse:
   std::vector<bool> marked(m,false);
   std::vector<int> permutation;
   for(int i=0;i<m;i++) {
      //std::cout << "n: " << i << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
      if(!marked[i] && neighbor_2[i]==-1) {
         // start here!
         int current = i;
	 int last = -1;
	 while(true) {
            permutation.push_back(current);
	    assert(!marked[current]);
	    marked[current] = true;
	    int next = (neighbor_2[current]==last)?neighbor_1[current]:neighbor_2[current];
	    //std::cout << "push " << current << ", next: "<< next<< "\n";
	    if(next==-1) {break;}
	    last = current;
	    current = next;
	 }
      }
   }

   ///  // print check:
   ///  for(int i=0;i<m;i++) {
   ///     std::cout << "fin: " << i << " " << marked[i] << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
   ///  }
   //std::cout << "perm " << permutation.size() << " " << m << "\n";
   assert(permutation.size()==m);
   
   // copy rows according to permutation:
   for(int i=0;i<n;i++) {
      double last = 0;
      for(int j=0;j<m;j++) { 
	 double aij = Polytope_get_a(p,permutation[j],i);
         Polytope_set_a(q,j,i,aij);
      }
   }

   polytope_count_pairs(q);
   return q;
}




