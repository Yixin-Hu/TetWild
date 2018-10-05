// Input: 
// Each V is a vector of size count.
// count is the number of tetrahedra to process
// Output:
// E is a vector of size count, it must be preallocated

export void energy_ispc(
  uniform double V1_x[], 
  uniform double V1_y[], 
  uniform double V1_z[], 
  uniform double V2_x[], 
  uniform double V2_y[], 
  uniform double V2_z[], 
  uniform double V3_x[], 
  uniform double V3_y[], 
  uniform double V3_z[], 
  uniform double V4_x[], 
  uniform double V4_y[], 
  uniform double V4_z[], 
  uniform double E[], 
  uniform int count) 
  {
    foreach (index = 0 ... count) 
    {
      // Load the appropriate tetrahedron vertices and compute the energy
      double helper_0[12];
      helper_0[0] = V1_x[index];
      helper_0[1] = V1_y[index];
      helper_0[2] = V1_z[index];
      helper_0[3] = V2_x[index];
      helper_0[4] = V2_y[index];
      helper_0[5] = V2_z[index];
      helper_0[6] = V3_x[index];
      helper_0[7] = V3_y[index];
      helper_0[8] = V3_z[index];
      helper_0[9] = V4_x[index];
      helper_0[10] = V4_y[index];
      helper_0[11] = V4_z[index];
      double helper_1 = helper_0[2];
      double helper_2 = helper_0[11];
      double helper_3 = helper_0[0];
      double helper_4 = helper_0[3];
      double helper_5 = helper_0[9];
      double helper_6 = 0.577350269189626*helper_3 - 1.15470053837925*helper_4 + 0.577350269189626*helper_5;
      double helper_7 = helper_0[1];
      double helper_8 = helper_0[4];
      double helper_9 = helper_0[7];
      double helper_10 = helper_0[10];
      double helper_11 = 0.408248290463863*helper_10 + 0.408248290463863*helper_7 + 0.408248290463863*helper_8 - 1.22474487139159*helper_9;
      double helper_12 = 0.577350269189626*helper_10 + 0.577350269189626*helper_7 - 1.15470053837925*helper_8;
      double helper_13 = helper_0[6];
      double helper_14 = -1.22474487139159*helper_13 + 0.408248290463863*helper_3 + 0.408248290463863*helper_4 + 0.408248290463863*helper_5;
      double helper_15 = helper_0[5];
      double helper_16 = helper_0[8];
      double helper_17 = 0.408248290463863*helper_1 + 0.408248290463863*helper_15 - 1.22474487139159*helper_16 + 0.408248290463863*helper_2;
      double helper_18 = 0.577350269189626*helper_1 - 1.15470053837925*helper_15 + 0.577350269189626*helper_2;
      double helper_19 = 0.5*helper_13 + 0.5*helper_4;
      double helper_20 = 0.5*helper_8 + 0.5*helper_9;
      double helper_21 = 0.5*helper_15 + 0.5*helper_16;
      double result_0 = -(helper_1*(-1.5*helper_1 + 0.5*helper_2 + helper_21) + helper_10*(-1.5*helper_10 + helper_20 + 0.5*helper_7) + helper_13*(-1.5*helper_13 + 0.5*helper_3 + 0.5*helper_4 + 0.5*helper_5) + helper_15*(0.5*helper_1 - 1.5*helper_15 + 0.5*helper_16 + 0.5*helper_2) + helper_16*(0.5*helper_1 + 0.5*helper_15 - 1.5*helper_16 + 0.5*helper_2) + helper_2*(0.5*helper_1 - 1.5*helper_2 + helper_21) + helper_3*(helper_19 - 1.5*helper_3 + 0.5*helper_5) + helper_4*(0.5*helper_13 + 0.5*helper_3 - 1.5*helper_4 + 0.5*helper_5) + helper_5*(helper_19 + 0.5*helper_3 - 1.5*helper_5) + helper_7*(0.5*helper_10 + helper_20 - 1.5*helper_7) + helper_8*(0.5*helper_10 + 0.5*helper_7 - 1.5*helper_8 + 0.5*helper_9) + helper_9*(0.5*helper_10 + 0.5*helper_7 + 0.5*helper_8 - 1.5*helper_9))*pow(pow((helper_1 - helper_2)*(helper_11*helper_6 - helper_12*helper_14) - (-helper_10 + helper_7)*(-helper_14*helper_18 + helper_17*helper_6) + (helper_3 - helper_5)*(-helper_11*helper_18 + helper_12*helper_17), (double)2), (double)-0.333333333333333);

      // Write the result back
      E[index] = result_0;
    }
}
