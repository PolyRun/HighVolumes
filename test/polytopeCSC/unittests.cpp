#include "../../src/volume/volume_helper.hpp"





int main(){

    { 
        Solved_Body *p = solved_body_generator()->get("cube_r1.0_10", false);
        p->polytopeCSC();

        Solved_Body *p_clone = p->clone();
        
        p->print();
        p_clone->print();

        assert(p != p_clone);

        int n = p->n;
        FT *x = (FT *) malloc(n * sizeof(FT));
        
        {
            for (int i = 0; i < n; i++) { x[i] = 0; }
            for (int b = 0; b < p->bcount; b++){
                assert(p->type[b]->inside(p->body[b], x));
            }

            FT *dir = (FT *) calloc(n, sizeof(FT));
            dir[0] = 1;
            FT t0, t1;
            p->type[0]->intersect(p->body[0], x, dir, &t0, &t1);
            std::cout << "t0 = " << t0 << " t1 = " << t1 << std::endl;

            free(dir);
        }
        {
            for (int i = 0; i < n; i++) { x[i] = 1; }
            for (int b = 0; b < p->bcount; b++){
                assert(p->type[b]->inside(p->body[b], x));
            }
        }
        {
            for (int i = 0; i < n; i++) { x[i] = 2; }
            for (int b = 0; b < p->bcount; b++){
                assert(!p->type[b]->inside(p->body[b], x));
            }
        }


        
        free(x);
    }

    std::cout << "tests passed!\n";
    
    return 0;
    
    
}
