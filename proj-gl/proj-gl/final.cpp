void rasterize_ellipse(int a, int b) {

  i_start, j_start;

  for(int quadrant = 1; quadrant < 4; ++quadrant) {
    if(quadrant == 1) {
      i_start = a;
      j_start = b;
    } else if (quadrant == 2) {
      i_start = -a;
      j_start = b;
    } else if (quadrant == 3) {
      i_start = -a;
      j_start = -b;
    } else {
      i_start = a;
      j_start = -b;
    }

    //LOOP A
    if(quadrant == 1 ||  quadrant == 4) {
      for(int i = 0; j = j_start; (b*b * i_ab) <= (a*a * j_ab); ++i) {
        i_ab = std::abs(i);
        j_ab = std::abs(j);
        draw(i, j);
        if(!(update(i+1, j))) {
          if(quadrant == 1) {
            j--;
          } else {
            j++;
          }
        }
      }

    } else { //quadrant 2 & quadrant 3
      for(int i = 0; j = j_start; (b*b * i_ab) <= (a*a * j_ab); ++j) {
        i_ab = std::abs(i);
        j_ab = std::abs(j);
        draw(i, j);
        if(!(update(i - 1, j))) {
          if(quadrant == 2) {
            j--;
          } else {
            j++;
          }
        }
      }
    }

    //LOOP B
    if(quadrant == 1 || quadrant == 2) {
      for(int i = i_start; j = 0; (b*b * i_ab) >= (a*a * j_ab); ++j) {
        i_ab = std::abs(i);
        j_ab = std::abs(j);
        draw(i, j);
        if(!(update(i, j+1))) {
          if(quadrant == 1) {
            i--;
          } else {
            i++;
          }
        }
      }
    } else { //quadrant 3 and quadrant 4
      for(int i = i_start; j = 0; (b*b * i_ab) >= (a*a * j_ab); --j) {
        i_ab = std::abs(i);
        j_ab = std::abs(j);
        draw(i, j);
        if(!(update(i, j - 1))) {
          if(quadrant == 4) {
            i--
          } else {
            i++;
          }
        }
      }
    }
  }
}
