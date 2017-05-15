/NO_K/ || /no_k/ { ni = nk; nk  = 1; }
($1 == "NODESI") { $3 = ni ; 
		 $6 = nj ;
		 $9 = nk ;
               }
{print $0;}
