#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <typeinfo>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom,
                       int _right, double graylevel, int width) {
    for (int w = -width / 2; w <= width / 2; w++) {
        int top = _top + w, left = _left + w, right = _right + w, bottom =
                _bottom + w;

        // if any of the coordinates are out-of-bounds, truncate them
        top = min(max(top, 0), input.rows() - 1);
        bottom = min(max(bottom, 0), input.rows() - 1);
        left = min(max(left, 0), input.cols() - 1);
        right = min(max(right, 0), input.cols() - 1);

        // draw top and bottom lines
        for (int j = left; j <= right; j++)
            input[top][j] = input[bottom][j] = graylevel;
        // draw left and right lines
        for (int i = top; i <= bottom; i++)
            input[i][left] = input[i][right] = graylevel;
    }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {
    NOTEHEAD = 0, QUARTERREST = 1, EIGHTHREST = 2
} Type;

class DetectedSymbol {
public:
    int row, col, width, height;
    Type type;
    char pitch;
    double confidence;
};

// Function that outputs the ascii detection output file
void write_detection_txt(const string &filename,
                         const vector<struct DetectedSymbol> &symbols) {
    ofstream ofs(filename.c_str());

    for (int i = 0; i < symbols.size(); i++) {
        const DetectedSymbol &s = symbols[i];
        ofs << s.row << " " << s.col << " " << s.width << " " << s.height
        << " ";
        if (s.type == NOTEHEAD)
            ofs << "filled_note " << s.pitch;
        else if (s.type == EIGHTHREST)
            ofs << "eighth_rest _";
        else
            ofs << "quarter_rest _";
        ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void write_detection_image(const string &filename,
                           const vector <DetectedSymbol> &symbols, const SDoublePlane &input) {
    SDoublePlane output_planes[3];
    for (int i = 0; i < 3; i++)
        output_planes[i] = input;

    for (int i = 0; i < symbols.size(); i++) {
        const DetectedSymbol &s = symbols[i];

        overlay_rectangle(output_planes[s.type], s.row, s.col,
                          s.row + s.height - 1, s.col + s.width - 1, 255, 2);
        overlay_rectangle(output_planes[(s.type + 1) % 3], s.row, s.col,
                          s.row + s.height - 1, s.col + s.width - 1, 0, 2);
        overlay_rectangle(output_planes[(s.type + 2) % 3], s.row, s.col,
                          s.row + s.height - 1, s.col + s.width - 1, 0, 2);

        if (s.type == NOTEHEAD) {
            char str[] = {s.pitch, 0};
            draw_text(output_planes[0], str, s.row, s.col + s.width + 1, 0, 2);
            draw_text(output_planes[1], str, s.row, s.col + s.width + 1, 0, 2);
            draw_text(output_planes[2], str, s.row, s.col + s.width + 1, 0, 2);
        }
    }

    SImageIO::write_png_file(filename.c_str(), output_planes[0],
                             output_planes[1], output_planes[2]);
}

int reflect(int pixel, int bound) { //need to consider filter of length more than 3
    //cout<<"reflect ("<<pixel<<","<<bound<<")";
    if (pixel < 0)
        pixel = -pixel - 1;
    else if (pixel > bound - 1)
        pixel = pixel - 1;
    //cout<<"return ("<<pixel<<","<<bound<<")\n\n";
    return pixel;
}


// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
//SDoublePlane convolve_separable(const SDoublePlane &input,
 //                               const SDoublePlane &row_filter, const SDoublePlane &col_filter) 
 
SDoublePlane convolve_separable(const SDoublePlane &input,
                                const SDoublePlane &row_filter) 
  {
    SDoublePlane output(input.rows(), input.cols());


    //splitting the matrix
    int rows = input.rows();
    int cols = input.cols();
    int rowHalf = rows / 2;
    double sum = 0;

    //first filter
    int **h1 = new int *[rowHalf];
    for (int i = 0; i < rowHalf; ++i)
        h1[i] = new int[cols];

    //second filter
    int **h2 = new int *[rows - rowHalf];
    for (int i = 0; i < rowHalf; ++i)
        h2[i] = new int[cols];

    for (int i = 0; i < input.rows(); i++) {
        for (int j = 0; j < input.cols(); j++) {
            sum = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
					if (i-m >= 0 && j-n >= 0 && i-m < input.rows() && j-n < input.cols()){
//						cout << "i:" << i << "  j:" << j << "  m:" << m << "  n:" << n << "  rows:" << input.rows() << "  cols:" << input.cols() << "\n";
//						continue;
						sum = sum + h1[m + 1][n + 1] * input[i - m][j - n];
					}						
                }
            }
            output[i][j] = sum;
        }

    }

    for (int i = 0; i < input.rows(); i++) {
        for (int j = 0; j < input.cols(); j++) {
            sum = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
					if (i-m >= 0 && j-n >= 0 && i-m < input.rows() && j-n < input.cols()){
//						cout << "i:" << i << "  j:" << j << "  m:" << m << "  n:" << n << "\n";
//						continue;
						sum = sum + h2[m + 1][+!n] * output[i - m][j - n];
					}
					
                }
            }
            output[i][j] = sum;
        }
    }

    return output;
}

// Convolve an image with a general convolution kernel

SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter) {
    SDoublePlane output(input.rows(), input.cols());

    // Convolution code here

    //filling pixels before boundaries


    double sum = 0.0;
    int r = input.rows();
    int c = input.cols();
    for (int i = 1; i < r - 1; i++) {
        for (int j = 1; j < c - 1; j++) {
            sum = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    sum = sum + filter[m + 1][n + 1] * input[i - m][j - n];
                }
            }
            output[i][j] = sum;
        }

    }

    //boundaries

    //top - input[0][j], bottom - input[input.rows()-1][j]
    for (int j = 0; j < c; j++) {
        sum = 0;
        for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
                //cout<<"Top: send ("<<0-m<<","<<r<<") and ("<<j-n<<","<<c<<")\n";
                sum = sum + filter[m + 1][n + 1] * input[reflect(0 - m, r)][reflect(j - n, c)];
            }
            output[0][j] = sum;
        }
        sum = 0;
        for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
                //cout<<"Bottom: send ("<<0-m<<","<<r<<") and ("<<j-n<<","<<c<<")\n";
                sum = sum + filter[m + 1][n + 1] * input[reflect(r - 1 - m, r)][reflect(j - n, c)];
            }
            output[r - 1][j] = sum;
        }
    }
    //left - input[i][0], right - input[i][input.cols()-1]
    for (int i = 0; i < r; i++) {
        sum = 0;
        for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
                sum = sum + filter[m + 1][n + 1] * input[reflect(i - m, r)][reflect(0 - n, c)];
            }
            output[i][0] = sum;
        }
        sum = 0;
        for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
                sum = sum + filter[m + 1][n + 1] * input[reflect(i - m, r)][reflect(c - 1 - n, c)];
            }
            output[i][c - 1] = sum;
        }
    }

    return output;
}

// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx) {
    SDoublePlane output(input.rows(), input.cols());

    // Implement a sobel gradient estimation filter with 1-d filters

    return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh = 0) {
    SDoublePlane output(input.rows(), input.cols());

    // Implement an edge detector of your choice, e.g.
    // use your sobel gradient operator to compute the gradient magnitude and threshold

    return output;
}

//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[]) {
    if (!(argc == 2)) {
        cerr << "usage: " << argv[0] << " input_image" << endl;
        return 1;
    }

    string input_filename(argv[1]);
    SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());

    // test step 2 by applying mean filters to the input image
    SDoublePlane mean_filter(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            mean_filter[i][j] = 1 / 9.0;
//    SDoublePlane output_image = convolve_general(input_image, mean_filter);
		 SDoublePlane output_image = convolve_separable(input_image, mean_filter);

    // randomly generate some detected symbols -- you'll want to replace this
    //  with your symbol detection code obviously!
    vector <DetectedSymbol> symbols;
    for (int i = 0; i < 10; i++) {
        DetectedSymbol s;
        s.row = rand() % input_image.rows();
        s.col = rand() % input_image.cols();
        s.width = 20;
        s.height = 20;
        s.type = (Type) (rand() % 3);
        s.confidence = rand();
        s.pitch = (rand() % 7) + 'A';
        symbols.push_back(s);
    }

    write_detection_txt("detected.txt", symbols);
    write_detection_image("detected.png", symbols, input_image);
    write_detection_image("detected2.png", symbols, output_image);
}
