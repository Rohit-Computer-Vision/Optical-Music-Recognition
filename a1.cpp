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
SDoublePlane convolve_separable(const SDoublePlane &input,
                               const SDoublePlane &row_filter, const SDoublePlane &col_filter)  
  {
    SDoublePlane output(input.rows(), input.cols());

	int k = row_filter.cols();
	int sr = k/2;
	int rows = row_filter.rows();
	double sum =0.0;
    for(int i =sr;i<input.rows() -sr;++i)
		for(int j = sr;j<input.cols() -sr; ++j){
			sum = 0.0;			
                for (int n = 0; n < row_filter.cols(); n++) {
					sum = sum +row_filter[rows-1][n] * input[i- (rows-1)][j-n+1];
				}

			output[i][j]=sum;
		}

	k = col_filter.rows();
	int ic = k/2;
	int cols = col_filter.cols();
	for(int i = ic;i<input.rows() - ic; ++i)
		for(int j = ic;j<input.cols() -ic;++j){
			sum=0.0;
			for (int m = 0; m < col_filter.rows(); m++) {

				 	sum = sum + col_filter[m][cols-1] * output[i-m+1][j - (cols-1)]; 

			 }
			 output[i][j] = sum;
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

// Detect symbols in the given input_image
vector <DetectedSymbol> detectSymbols(SDoublePlane input_image, SDoublePlane template_image ) {

	vector <DetectedSymbol> symbols;

	double sum = 0.0;
    int r = input_image.rows();
    int c = input_image.cols();

	/*
	cout << "Input Image: \n";
    for (int i = 0; i < r - 1; i++) {
        for (int j = 0; j < c - 1; j++) {
			cout << input_image[i][j] << "  " ;
             sum = 0;
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    sum = sum + filter[m + 1][n + 1] * input[i - m][j - n];
                }
            }
            output[i][j] = sum;
        }
		cout << "\n";
    }
	*/

    for (int i = 0; i < 1; i++) {
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


	return symbols;
}

// Print an image to a file
void printImg2File(string filename, SDoublePlane img){
	ofstream outFile;
	outFile.open(filename.c_str());

	int r = img.rows();
    int c = img.cols();
    for (int i = 0; i < r ; i++) {
        for (int j = 0; j < c ; j++) {
			outFile << img[i][j] << ",";
		}
		outFile << "\n";
	}			
	outFile.close();	
}

//Converts a grey scale image to binary image
SDoublePlane convert_binary(SDoublePlane &input){

	int rows = input.rows();
	int cols = input.cols();
	SDoublePlane output(rows,cols);

	double threshold = 200.0;
	for(int i= 0;i<rows;++i)
		for(int j=0;j<cols;++j){
			if(input[i][j] >= threshold)
				output[i][j] = 1;
			else
				output[i][j] =0;
		}
	return output;
}


SDoublePlane find_hamming_distance(SDoublePlane &img_input, SDoublePlane &img_template){
	
	int input_rows = img_input.rows();
	int input_cols = img_input.cols();
	int template_rows = img_template.rows();
	int template_cols = img_template.cols();
		
	double sum = 0.0;
	SDoublePlane output(input_rows, input_cols);
	for(int i =0;i<input_rows - template_rows ; ++i){
		for(int j =0;j<input_cols - template_cols; ++j){
			sum=0.0;
			for(int k =0;k<template_rows;++k){
				for(int l=0;l<template_cols;++l){
					sum = sum + (img_input[i+k][j+l] * img_template[k][l] + (1 - img_input[i+k][j+l])*(1 - img_template[k][l]));	
				}
			}
			output[i][j] = sum;
		}
	}
	printImg2File("input_img_file_Output.txt", img_input);
	printImg2File("img2fileOutput.txt", output);
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

	string TEMPLATE_NOTEHEAD = "template1.png";
	string TEMPLATE_QUARTERREST = "template2.png";
	string TEMPLATE_EIGHTHREST = "template3.png";

    SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());


    // test step 2 by applying mean filters to the input image
    SDoublePlane mean_filter(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            mean_filter[i][j] = 1 / 9.0;
  // SDoublePlane output_image = convolve_general(input_image, mean_filter);

	SDoublePlane row_filter(1,3);
	SDoublePlane col_filter(3,1);

	//1 row and three columns
	for (int i = 0; i < 1; i++)
        for (int j = 0; j < 3; j++)
            row_filter[i][j] = 1 / 3.0;

	for(int i =0; i<3;i++)
		for (int j=0;j<1;j++)
			col_filter[i][j]= 1/3.0;


		// Convolve General 2D Kernel	
	// SDoublePlane output_image = convolve_general(input_image, mean_filter); // Uncomment Later

	// Convolve Separable Kernel			 
	 SDoublePlane output_image = convolve_separable(input_image, row_filter, col_filter); // Uncomment Later

	// Read NOTEHEAD - Template 1
	SDoublePlane template_notehead = SImageIO::read_png_file(TEMPLATE_NOTEHEAD.c_str());
	SDoublePlane template_notehead_grey_scale = convert_binary(template_notehead);

	// Read QUARTERREST
	SDoublePlane template_quarterrest = SImageIO::read_png_file(TEMPLATE_QUARTERREST.c_str());
	SDoublePlane template_quarterrest_grey_scale = convert_binary(template_quarterrest);
	
 
 // Read EIGHTHREST
	SDoublePlane template_eighthrest = SImageIO::read_png_file(TEMPLATE_EIGHTHREST.c_str());
	SDoublePlane template_eighthrest_grey_scale = convert_binary(template_eighthrest);
	
//	vector <DetectedSymbol> symbols = detectSymbols(input_image, template_img_notehead);




	SDoublePlane convoluted_image = convolve_separable(input_image, row_filter, col_filter);
	SDoublePlane convoluted_template = convolve_separable(template_notehead, row_filter, col_filter);
	
	SDoublePlane binary_image = convert_binary(convoluted_image);
	SDoublePlane binary_template = convert_binary(convoluted_template);
	

	find_hamming_distance(binary_image, binary_template);


	/*
    // randomly generate some detected symbols -- you'll want to replace this
    //  with your symbol detection code obviously!
    vector <DetectedSymbol> symbols;
    for (int i = 0; i < 1; i++) { 
        DetectedSymbol s;
        s.row = rand() % input_image.rows();
        s.col = rand() % input_image.cols();
        s.width = 20;
        s.height = 20;
        s.type = (Type) (rand() % 3);
        s.confidence = rand();
        s.pitch = (rand() % 7) + 'A';
        symbols.push_back(s);
    } */

    // write_detection_txt("detected.txt", symbols);
    // write_detection_image("detected.png", symbols, input_image);
    // write_detection_image("detected2.png", symbols, output_image);

	//printImg2File("img2fileOutput.txt", template_img_notehead );
}
																																																																																																																																																																																																																																																																																																																											