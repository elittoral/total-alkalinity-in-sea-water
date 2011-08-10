function [ output_args ] = DensNaCl( CNaCl, T )

  mNaCl = CNaCl/(1 - 0.058443*CNaCl) ;


      DH2O = 999.842594 + 6.793952e-2 * T - 9.095290e-3 * T^2 ...
     +     1.001685e-4  * T^3 - 1.120083e-6 * T^4 - 6.536332e-9 * T^5;

      DNaCl = DH2O + mNaCl*(46.5655 - 0.2341*t + 3.4128e-3 * T^2 ...
     +                     -  2.7030e-5 * T^3  + 1.4037e-7 * T^4)   ...
     + mNaCl^1.5  *  (-1.8527 + 5.3956e-2  * T - 6.2635e-4   * T^2) ...
     +  mNaCl^2  *  (-1.6368 - 9.5653e-4  * T + 5.2829e-5 *   T^2) ...
     +          0.2274  * mNaCl^2.5;

      output_args = 1e-3 *   DNaCl;

end

