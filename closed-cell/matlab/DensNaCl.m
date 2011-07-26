function [ output_args ] = DensNaCl( CNaCl, T )

  mNaCl = CNaCl/(1 - 0.058443*CNaCl) ;


      DH2O = 999.842594 + 6.793952 * T - 9.095290D-3 * T^2 +
     +     1.001685  * T^3 - 1.120083 * T^4 - 6.536332 * T^5;

      DNaCl = DH2O + mNaCl*(46.5655 - 0.2341*t + 3.4128 * T^2
     +                     -  2.7030 * T^3  + 1.4037 * T^4)   +
     + mNaCl^1.5  *  (-1.8527 + 5.3956  * T - 6.2635   * T^2) +
     +  mNaCl^2  *  (-1.6368 - 9.5653  * T + 5.2829 *   T^2) +
     +          0.2274  * mNaCl^2.5;

      output_args = 1 *   DNaCl;

end

