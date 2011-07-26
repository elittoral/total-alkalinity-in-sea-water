function [ K1, K2, KW ] = ConstsNaCl( CNaCl, T )

TK = 273.15 + T;

     if (abs(CNaCl - 0.7) > 0.05)
         WRITE (6,'(" ConstsNaCl: C(NaCl) .ne. 0.7")')
     end
     if (abs(T - 25) > 0.1)
         WRITE (6,'(" ConstsNaCl: T .ne. 25 ")')
     end

     if ((abs(CNaCl - 0.7) > 0.05) || (abs(T - 25) > 0.1)) STOP


     K1 = exp(-13.82);
     K2 = exp(-21.97) ;
     KW = exp(-31.71);


end

