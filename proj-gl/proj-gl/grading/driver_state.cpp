#include "driver_state.h"
#include <cstring>
#include <algorithm>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;

    unsigned int p = width * height;
    state.image_color = new pixel[p];
    state.image_depth = new float[p];
    for(unsigned int i = 0; i < p; ++i) {
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = 2;
    }

    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
      data_vertex vert;
      int f_p_v = state.floats_per_vertex; //amount of data per vertex
      data_geometry triangle[state.num_vertices];
      switch(type) {
        case render_type::triangle: {

          //the array of vertex_data can be divided into num_vertices sub-arrays
          //iterate to the 1st element of each sub-array by adding multiples of floats_per_vertex
          //put 1st element of each sub-array into a data_vertex variable and a data_geometry variable
          for(int i = 0; i < state.num_vertices; ++i) { //iterate through each vertex
            vert.data = state.vertex_data + (i * f_p_v);
            triangle[i].data = vert.data;
            state.vertex_shader(vert, triangle[i], state.uniform_data);
          }

          //every 3 vertices is 1 triangle
          //call clip_triangle on each set of 3 triangles
          for(int k = 0; k < state.num_vertices; k += 3) {
            clip_triangle(state, triangle[k], triangle[k + 1], triangle[k + 2], 0);
          }

          break;
        } case render_type::indexed: {
          int tri = state.num_triangles * 3; //indices per triangle
          int index; //index of index_data array
          data_geometry triangle[tri];

          //the index_data array refers to a vertex
          //iterate through index_data to find its corresponding vertex and vertex_data entry
          //put each element found in vertex_data into a data_vertex variable and a data_geometry variable
          for(int i = 0; i < tri; ++i) {
            index = *(state.index_data + i);
            vert.data = state.vertex_data + (index * f_p_v);
            triangle[i].data = vert.data;
            state.vertex_shader(vert, triangle[i], state.uniform_data);
          }

          for(int k = 0; k < tri; k += 3) {
            clip_triangle(state, triangle[k], triangle[k + 1], triangle[k + 2], 0);
          }

          break;
        } case render_type::fan: {
          for(int i = 0; i < state.num_vertices; ++i) { //iterate through each vertex
            vert.data = state.vertex_data + (i * f_p_v);
            triangle[i].data = vert.data;
            state.vertex_shader(vert, triangle[i], state.uniform_data);
          }

          // the 1st vertex is always triangle[0]
          //call clip_triangle on each set of 3 triangles
          for(int k = 1; k < state.num_vertices - 1; k++) {
            clip_triangle(state, triangle[0], triangle[k], triangle[k + 1], 0);
          }
          break;
        } case render_type::strip: {
          for(int i = 0; i < state.num_vertices; ++i) { //iterate through each vertex
            vert.data = state.vertex_data + (i * f_p_v);
            triangle[i].data = vert.data;
            state.vertex_shader(vert, triangle[i], state.uniform_data);
          }

          for(int k = 0; k < state.num_vertices - 2; k++) {
            clip_triangle(state, triangle[k], triangle[k + 1], triangle[k + 2], 0);
          }
          break;
        }
        default:
          break;
        }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    //clipping is done by cutting the triangle along each edge of the image-space
    //since image-space is a cube, there are 6 edges(cases)
    //the image-space ranges from -w to w, check if vertices A, B, and C are within the image-space
    bool pos, a = true, b = true, c = true;
    int index = 0;
    switch(face) {
      case 0: {//x = w
        pos = true;
        index = 0;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } case 1: {//x = -w
        pos = false;
        index = 0;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } case 2: {//y = w
        pos = true;
        index = 1;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } case 3: {//y = -w
        pos = false;
        index = 1;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } case 4: {//z = w
        pos = true;
        index = 2;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } case 5: {//z = -w
        pos = false;
        index = 2;
        validity_check(v0, v1, v2, pos, a, b, c, index);
        break;
      } default: {
        break;
      }
    }

    //there are 8 possible cases when checking if the vertices are within range
    //if out-of-range, create new vertices
    data_geometry p, q, r;
    if(a && !b && !c) {
      //A is valid
      new_vertex(state, v0, v1, p, pos, index);
      new_vertex(state, v0, v2, q, pos, index);
      r.gl_Position = v0.gl_Position;
      r.data = v0.data;
      clip_triangle(state, r, p, q, face+1); //triangle APQ
    } else if (!a && b && !c) {
      //B is valid
      new_vertex(state, v1, v0, p, pos, index);
      new_vertex(state, v1, v2, q, pos, index);
      r.gl_Position = v1.gl_Position;
      r.data = v1.data;
      clip_triangle(state, r, q, p, face+1); //triangle BQP
    } else if (!a && !b && c) {
      //C is valid
      new_vertex(state, v2, v0, p, pos, index);
      new_vertex(state, v2, v1, q, pos, index);
      r.gl_Position = v2.gl_Position;
      r.data = v2.data;
      clip_triangle(state, r, p, q, face+1); //triangle CPQ
    } else if (a && b && !c) {
      //A and B are valid
      new_vertex(state, v0, v2, p, pos, index);
      r.gl_Position = v0.gl_Position;
      r.data = v0.data;
      q.gl_Position = v1.gl_Position;
      q.data = v1.data;
      clip_triangle(state, r, q, p, face+1); //triangle ABP

      new_vertex(state, v1, v2, q, pos, index);
      r.gl_Position = v1.gl_Position;
      r.data = v1.data;
      clip_triangle(state, r, q, p, face+1);  //triangle BQP
    } else if (a && !b && c) {
      //A and C are valid
      new_vertex(state, v2, v1, p, pos, index);
      r.gl_Position = v0.gl_Position;
      r.data = v0.data;
      q.gl_Position = v2.gl_Position;
      q.data = v2.data;
      clip_triangle(state, q, r, p, face+1); //triangle CAP

      new_vertex(state, v0, v1, q, pos, index);
      r.gl_Position = v0.gl_Position;
      r.data = v0.data;
      clip_triangle(state, r, q, p, face+1);  //triangle AQP
    } else if (!a && b && c) {
      //B and C are valid
      new_vertex(state, v1, v0, p, pos, index);
      r.gl_Position = v1.gl_Position;
      r.data = v1.data;
      q.gl_Position = v2.gl_Position;
      q.data = v2.data;
      clip_triangle(state, r, q, p, face+1); //triangle BCP

      new_vertex(state, v2, v0, q, pos, index);
      r.gl_Position = v2.gl_Position;
      r.data = v2.data;
      clip_triangle(state, r, q, p, face+1);  //triangle CQP
    } else if (a && b && c) {
      //triangle is fully inside
      clip_triangle(state, v0, v1, v2, face+1);
    } else {
      //triangle is fully outside
    }

}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    float ax, ay, bx, by, cx, cy, t_area, alpha, a_area, beta, b_area, gamma, g_area;
    size_t width = state.image_width;
    size_t height = state.image_height;

    //calculate the x and y coordinates of each vertex
    ax = (v0.gl_Position[0]/v0.gl_Position[3] + 1) * 0.5 * width;
    ay = (v0.gl_Position[1]/v0.gl_Position[3] + 1) * 0.5 * height;
    bx = (v1.gl_Position[0]/v1.gl_Position[3] + 1) * 0.5 * width;
    by = (v1.gl_Position[1]/v1.gl_Position[3] + 1) * 0.5 * height;
    cx = (v2.gl_Position[0]/v2.gl_Position[3] + 1) * 0.5 * width;
    cy = (v2.gl_Position[1]/v2.gl_Position[3] + 1) * 0.5 * height;

    // std::cout << "(" << ax << ", " << ay << ")" << std::endl;
    // std::cout << "(" << bx << ", " << by << ")" << std::endl;
    // std::cout << "(" << cx << ", " << cy << ")" << std::endl << std::endl;

    //minimize iterations by finding the upper and lower bounds of vertices
    size_t h_max = std::max({ay, by, cy});
    size_t h_min = std::min({ay, by, cy});
    size_t w_max = std::max({ax, bx, cx});
    size_t w_min = std::min({ax, bx, cx});

    t_area = 0.5 * ((bx*cy - cx*by) + (cx*ay - ax*cy) + (ax*by - bx*ay));

    //iterate through bounds to create triangle
    for(size_t i = w_min; i <= w_max; ++i) {
      for(size_t j = h_min; j <= h_max; ++j) {
        float px =  i;
        float py =  j;

        //calculating barycentric coordinates
        a_area = 0.5 * ((bx*cy - cx*by) + (cx*py - px*cy) + (px*by - bx*py));
        alpha = a_area / t_area;

        b_area = 0.5 * ((px*cy - cx*py) + (cx*ay - ax*cy) + (ax*py - px*ay));
        beta = b_area / t_area;

        g_area = 0.5 * ((bx*py - px*by) + (px*ay - ax*py) + (ax*by - bx*ay));
        gamma = g_area / t_area;

        //checking if current position is within triangle space
        if( (alpha >= 0)  && (beta >= 0) && (gamma >= 0) ) {

          //find depth of z-position for z-buffering
          float z_position = (alpha * v0.gl_Position[2] / v0.gl_Position[3]) + (beta * v1.gl_Position[2] / v1.gl_Position[3]) + (gamma * v2.gl_Position[2] / v2.gl_Position[3]);

          //if z-postion depth < current depth, change color/data and update image_depth
          if(z_position < state.image_depth[(j * width) + i]) {

            state.image_depth[(j * width) + i] = z_position;

            //create data_fragment variable for determining color at current position
            data_fragment frag;
            frag.data = new float[MAX_FLOATS_PER_VERTEX];
            data_output out;

            //iterate through interp_rules to determine how to interpolate vertex_data
            for(int k = 0; k < state.floats_per_vertex; ++k) {
              switch(state.interp_rules[k]) {
                case interp_type::flat: {
                  frag.data[k] = v0.data[k];
                  break;
                } case interp_type::noperspective: {
                  frag.data[k] = (alpha * v0.data[k]) + (beta * v1.data[k]) + (gamma * v2.data[k]);
                  break;
                } case interp_type::smooth: {
                    float divisor, r_alpha, r_beta, r_gamma, p_alpha, p_beta, p_gamma;

                    p_alpha = alpha / v0.gl_Position[3];
                    p_beta =  beta / v1.gl_Position[3];
                    p_gamma = gamma / v2.gl_Position[3];
                    divisor = p_alpha + p_beta + p_gamma;
                    r_alpha = p_alpha / divisor;
                    r_beta = p_beta / divisor;
                    r_gamma = p_gamma / divisor;

                    frag.data[k] = (r_alpha * v0.data[k]) + (r_beta * v1.data[k]) + (r_gamma * v2.data[k]);
                    break;
                  }
                default: {
                  break;
                }
              }
            }
            state.fragment_shader(frag, out, state.uniform_data);
            state.image_color[(j * width) + i] = make_pixel(out.output_color[0] * 255, out.output_color[1] * 255, out.output_color[2] * 255);
          }
        }
      }
    }
}

//helper function for clipping used to determine whether the vertex is outside the image-space
void validity_check(const data_geometry& v0, const data_geometry& v1,
  const data_geometry& v2, bool& pos, bool& a, bool& b, bool& c, int& index) {
  if(pos) { // x/y/z = w
    a = (v0.gl_Position[index] <= v0.gl_Position[3]) ? true : false;
    b = (v1.gl_Position[index] <= v1.gl_Position[3]) ? true : false;
    c = (v2.gl_Position[index] <= v2.gl_Position[3]) ? true : false;
  } else {// x/y/z = -w
    a = (v0.gl_Position[index] >= -v0.gl_Position[3]) ? true : false;
    b = (v1.gl_Position[index] >= -v1.gl_Position[3]) ? true : false;
    c = (v2.gl_Position[index] >= -v2.gl_Position[3]) ? true : false;
  }
}

//helper function for clipping used to create a new vertex within the image-space
void new_vertex(driver_state& state, const data_geometry& in,
  const data_geometry& out, data_geometry& point, bool& pos, int& coord) {
    float lambda, gamma;

    //find intersection point
    if(pos) {// x/y/z = w
      lambda = (out.gl_Position[3] - out.gl_Position[coord]) / (in.gl_Position[coord] - in.gl_Position[3] + out.gl_Position[3] - out.gl_Position[coord]);
    } else {// x/y/z = -w
      lambda = (-out.gl_Position[3] - out.gl_Position[coord]) / (in.gl_Position[coord] + in.gl_Position[3] - out.gl_Position[3] - out.gl_Position[coord]);
    }

    //create new vertex at intersection point
    point.gl_Position = (lambda * in.gl_Position) + ((1 - lambda) * out.gl_Position);
    point.data = new float[MAX_FLOATS_PER_VERTEX];

    //create data values for new vertex
    gamma = (lambda * in.gl_Position[3]) / (lambda * in.gl_Position[3] + (1 - lambda) * out.gl_Position[3]);
    for(int k = 0; k < state.floats_per_vertex; ++k) {
      switch(state.interp_rules[k]) {
        case interp_type::flat: {
          point.data[k] = in.data[k];
          break;
        } case interp_type::noperspective: {
          point.data[k] = gamma *  in.data[k] + (1.0 - gamma) * out.data[k];
          break;
        } case interp_type::smooth: {
          point.data[k] = lambda * in.data[k] + (1.0 - lambda) * out.data[k];
          break;
          }
        default: {
          break;
        }
      }
    }
  }
