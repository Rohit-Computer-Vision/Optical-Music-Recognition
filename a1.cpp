
#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <typeinfo>
#include <string>
#include <set>

using namespace std;

#define DEG2RAD(d) (d * M_PI / 180.0)

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

typedef enum {
	NOTEHEAD = 0, QUARTERREST = 1, EIGHTHREST = 2
} Type;

SDoublePlane make_mean_filter(int n) {
	SDoublePlane mean_filter(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			mean_filter[i][j] = 1.0 / (n*n);
	return mean_filter;
}

SDoublePlane make_row_filter(int n) {
	SDoublePlane row_filter(1, n);
	for (int j = 0; j < n; j++)
		row_filter[0][j] = 1.0 / n;
	return row_filter;
}

SDoublePlane make_col_filter(int n) {
	SDoublePlane col_filter(n, 1);
	for (int i = 0; i < n; i++)
		col_filter[i][0] = 1.0 / n;
	return col_filter;
}

class DetectedSymbol {
public:
	int row, col, width, height;
	Type type;
	char pitch;
	double confidence;
};

class LineLocation{
public:
		int row, col;
		char above_marker;
		char line_marker;
		char below_marker;
};

class HammingDistances{
public: 
	SDoublePlane hamming_matrix;
	double max_hamming_distance;
};


class Line {
public:
	int x1,y1,x2,y2; 
};

// Function that outputs the ascii detection output file
void write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols) {
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


void write_detection_image(const string &filename, const SDoublePlane &input) {
	SDoublePlane output_planes[3];
	for (int i = 0; i < 3; i++)
		output_planes[i] = input;
		
	SImageIO::write_png_file(filename.c_str(), output_planes[0],
			output_planes[1], output_planes[2]);
}

// Function that outputs a visualization of detected symbols
void write_detection_image(const string &filename,
		const vector<DetectedSymbol> &symbols, const SDoublePlane &input) {
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
			char str[] = { s.pitch, 0 };
			draw_text(output_planes[0], str, s.row, s.col + s.width + 1, 0, 2);
			draw_text(output_planes[1], str, s.row, s.col + s.width + 1, 0, 2);
			draw_text(output_planes[2], str, s.row, s.col + s.width + 1, 0, 2);
		}
	}

	SImageIO::write_png_file(filename.c_str(), output_planes[0],
			output_planes[1], output_planes[2]);
}

int reflect(int pixel, int bound) {
	if (pixel < 0)
		pixel = -pixel - 1;
	else if (pixel > bound - 1)
		pixel = pixel - 1;
	return pixel;
}

// Convolve an image with a separable convolution kernel
SDoublePlane convolve_separable(const SDoublePlane &input,
		const SDoublePlane &row_filter, const SDoublePlane &col_filter) {
	SDoublePlane output(input.rows(), input.cols());

	int k = row_filter.cols();
	int sr = k / 2;
	int rows = row_filter.rows();
	double sum = 0.0;
	for (int i = sr; i < input.rows() - sr; ++i)
		for (int j = sr; j < input.cols() - sr; ++j) {
			sum = 0.0;
			for (int n = 0; n < row_filter.cols(); n++) {
				sum = sum
						+ row_filter[rows - 1][n]
								* input[i - (rows - 1)][j - n + 1];
			}

			output[i][j] = sum;
		}

	k = col_filter.rows();
	int ic = k / 2;
	int cols = col_filter.cols();
	for (int i = ic; i < input.rows() - ic; ++i)
		for (int j = ic; j < input.cols() - ic; ++j) {
			sum = 0.0;
			for (int m = 0; m < col_filter.rows(); m++) {

				sum = sum
						+ col_filter[m][cols - 1]
								* output[i - m + 1][j - (cols - 1)];

			}
			output[i][j] = sum;
		}

	return output;

}

// Convolve an image with a general convolution kernel
SDoublePlane convolve_general(const SDoublePlane &input,
		const SDoublePlane &filter) {
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
				sum = sum
						+ filter[m + 1][n + 1]
								* input[reflect(0 - m, r)][reflect(j - n, c)];
			}
			output[0][j] = sum;
		}
		sum = 0;
		for (int m = -1; m < 2; m++) {
			for (int n = -1; n < 2; n++) {
				sum =
						sum
								+ filter[m + 1][n + 1]
										* input[reflect(r - 1 - m, r)][reflect(
												j - n, c)];
			}
			output[r - 1][j] = sum;
		}
	}
	
	for (int i = 0; i < r; i++) {
		sum = 0;
		for (int m = -1; m < 2; m++) {
			for (int n = -1; n < 2; n++) {
				sum = sum
						+ filter[m + 1][n + 1]
								* input[reflect(i - m, r)][reflect(0 - n, c)];
			}
			output[i][0] = sum;
		}
		sum = 0;
		for (int m = -1; m < 2; m++) {
			for (int n = -1; n < 2; n++) {
				sum =
						sum
								+ filter[m + 1][n + 1]
										* input[reflect(i - m, r)][reflect(
												c - 1 - n, c)];
			}
			output[i][c - 1] = sum;
		}
	}

	return output;
}

// Apply a sobel operator to an image
SDoublePlane sobel_filter(const SDoublePlane &input) {
	SDoublePlane output1(input.rows(), input.cols());
	SDoublePlane output2(input.rows(), input.cols());
	SDoublePlane output(input.rows(), input.cols());

	// Implement a sobel gradient estimation filter with 1-d filters
	SDoublePlane sobelHorFilter(3, 3);
	SDoublePlane sobelVerFilter(3, 3);
	
	//initialize
	sobelHorFilter[0][0] = -1;
	sobelHorFilter[0][1] = 0;
	sobelHorFilter[0][2] = 1;
	sobelHorFilter[1][0] = -2;
	sobelHorFilter[1][1] = 0;
	sobelHorFilter[1][2] = 2;
	sobelHorFilter[2][0] = -1;
	sobelHorFilter[2][1] = 0;
	sobelHorFilter[2][2] = 1;
	
	sobelVerFilter[0][0] = -1;
	sobelVerFilter[0][1] = -2;
	sobelVerFilter[0][2] = -1;
	sobelVerFilter[1][0] = 0;
	sobelVerFilter[1][1] = 0;
	sobelVerFilter[1][2] = 0;
	sobelVerFilter[2][0] = 1;
	sobelVerFilter[2][1] = 2;
	sobelVerFilter[2][2] = 1;
	
	output1 = convolve_general(input, sobelHorFilter);
	output2 = convolve_general(input, sobelVerFilter);

	for (int i = 0; i < input.rows(); i++) 
		for (int j = 0; j < input.cols(); j++){
			output[i][j] = sqrt(output1[i][j]*output1[i][j] + output2[i][j]*output2[i][j]);
		}
	
	return output;
}

SDoublePlane direct_sobel(const SDoublePlane &input) {
	SDoublePlane output1(input.rows(), input.cols());
	SDoublePlane output2(input.rows(), input.cols());
	SDoublePlane output(input.rows(), input.cols());

	// Implement a sobel gradient estimation filter with 1-d filters
	SDoublePlane sobelHorFilter(3, 3);
	SDoublePlane sobelVerFilter(3, 3);
	
	//initialize
	sobelHorFilter[0][0] = -1;
	sobelHorFilter[0][1] = 0;
	sobelHorFilter[0][2] = 1;
	sobelHorFilter[1][0] = -2;
	sobelHorFilter[1][1] = 0;
	sobelHorFilter[1][2] = 2;
	sobelHorFilter[2][0] = -1;
	sobelHorFilter[2][1] = 0;
	sobelHorFilter[2][2] = 1;
	
	sobelVerFilter[0][0] = -1;
	sobelVerFilter[0][1] = -2;
	sobelVerFilter[0][2] = -1;
	sobelVerFilter[1][0] = 0;
	sobelVerFilter[1][1] = 0;
	sobelVerFilter[1][2] = 0;
	sobelVerFilter[2][0] = 1;
	sobelVerFilter[2][1] = 2;
	sobelVerFilter[2][2] = 1;
	
	output1 = convolve_general(input, sobelHorFilter);
	output2 = convolve_general(input, sobelVerFilter);

	for (int i = 0; i < input.rows(); i++) 
		for (int j = 0; j < input.cols(); j++) {
			output[i][j] = sqrt(output1[i][j]*output1[i][j] + output2[i][j]*output2[i][j]);
			if(output[i][j] > 200)
				output[i][j] = 255;
			else
				output[i][j] = 0;
		}
				
	return output;
}

// Apply an edge detector to an image, returns the binary edge map
SDoublePlane find_edges(const SDoublePlane &input, double thresh = 0) {
	SDoublePlane output(input.rows(), input.cols());

	// Implement an edge detector of your choice, e.g.
	// use your sobel gradient operator to compute the gradient magnitude and threshold

	return output;
}

// Detect symbols in the given input_image

// Print an image to a file
void printImg2File(string filename, SDoublePlane img) {
	ofstream outFile;
	outFile.open(filename.c_str());

	for (int i = 0; i < img.rows(); i++) {
		for (int j = 0; j < img.cols(); j++) {
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

	double threshold = 100;
	for(int i= 0;i<rows;++i)
		for(int j=0;j<cols;++j){
			if(input[i][j] >= threshold)
				output[i][j] = 1;
			else
				output[i][j] = 0;
		}
	return output;
}

SDoublePlane convert_BW(SDoublePlane &input){

	int rows = input.rows();
	int cols = input.cols();
	SDoublePlane output(rows,cols);

	double threshold = 100.0;
	for(int i= 0;i<rows;++i)
		for(int j=0;j<cols;++j){
			if(input[i][j] >= threshold || input[i][j] == 1)
				output[i][j] = 255;
			else
				output[i][j] =0;
		}
	return output;
}

SDoublePlane convert_binary_to_BW(SDoublePlane &input){

	int rows = input.rows();
	int cols = input.cols();
	SDoublePlane output(rows,cols);

	for(int i= 0; i<rows; ++i)
		for(int j=0; j<cols; ++j){
			if(input[i][j] < 1)
				output[i][j] = 0;
			else
				output[i][j] = 255;
		}
	return output;
}

SDoublePlane display_binary(SDoublePlane &input){

	int rows = input.rows();
	int cols = input.cols();
	SDoublePlane output(rows,cols);

	for(int i= 0; i<rows; ++i)
		for(int j=0; j<cols; ++j){
			if(input[i][j] < 0.3)
				output[i][j] = 0;
			else
				output[i][j] = 255;
		}
	return output;
}

SDoublePlane convert_blur_to_binary(SDoublePlane &input){

	int rows = input.rows();
	int cols = input.cols();
	SDoublePlane output(rows,cols);

	for(int i= 0; i<rows; ++i)
		for(int j=0; j<cols; ++j){
			if(input[i][j] < 1)
				output[i][j] = 0;
			else
				output[i][j] = 1;
		}
	return output;
}

SDoublePlane central_edge(SDoublePlane &input){
	//input is a pure b&w image/template
	int rows = input.rows();
	int cols = input.cols();
	int count = 0, start = 0, mid = 0;
	SDoublePlane output(rows,cols);

	for(int i=0; i<rows; i++) {
		count = 0;
		for(int j=0; j<cols; j++){
			if(input[i][j] == 1) {
				count++;
				if(count == 1) {
					start = j;
					//printf("\nstart: %d", start);
				}
				//printf("\nif: (%d,%d,%d), ", i, j, count);
			}
			else if(input[i][j] == 0 or j == cols-1) {
				//printf("\nelse: (%d,%d,%d), ", i, j, count);
				if(count > 0) {
					mid = start + count/2;
					while(start < j) {
						if(start == mid) {
							output[i][start] = 255;
							//printf("\nOutput[%d][%d]=255", i, start);
						}
						start++;
					}
				}
				count = 0;
			}
		}
	}
	return output;
}

SDoublePlane alternate_edge(SDoublePlane &input){
	int rows = input.rows();
	int cols = input.cols();
	int count = 0, start = 0, mid = 0;
	SDoublePlane output(rows,cols);

	//making right and bottom boundary black
	for(int j=0; j<cols; j++)
		input[rows-1][j] = 0;
	for(int i=0; i<rows; i++)
		input[i][cols-1] = 0;
	
	for(int j=0; j<cols; j++){
		for(int i=0; i<rows; i++) {
			if(input[i][j] < input[i][j+1])
				input[i][j] = 0;
		}
	}
	return input;
}

int dmin(int a, int b, int c) {
	int min = 10000;
	if(a < min)	min = a;
	if(b < min)	min = b;
	if(c < min)	min = c;
	return min;
}

SDoublePlane calculate_D(SDoublePlane &binary_image_blur_sobel) {
	int i, j, input_rows = binary_image_blur_sobel.rows(), input_cols = binary_image_blur_sobel.cols();
	SDoublePlane D(input_rows, input_cols);

	//initialize edges as 0 and others as infinity(10000) in binary_image_blur_sobel
	for(i=0;i<input_rows;i++)
		for(j=0;j<input_cols;j++)
			if(binary_image_blur_sobel[i][j] == 1)
				D[i][j] = 0;
			else
				D[i][j] = 10000;
	
	//distance below and to the right of edges
	for(i=1;i<input_rows;i++)
		for(j=1;j<input_cols;j++)
			D[i][j] = dmin(D[i][j], D[i][j-1]+1, D[i-1][j]+1);
	
	//distance above and to the left of edges
	for(i=input_rows-2; i>=0; i--)
		for(j=input_cols-2; j>=0 ;j--)
			D[i][j] = dmin(D[i][j], D[i+1][j]+1, D[i][j+1]+1);
	
	//bottom row
	i = input_rows-1;	
	for(j=input_cols-2; j>=0 ;j--)
		D[i][j] = dmin(D[i][j], D[i-1][j]+1, D[i][j+1]+1);
	
	//rightmost column
	j = input_cols-1;
	for(i=input_rows-2; i>=0; i--)
		D[i][j] = dmin(D[i][j], D[i+1][j]+1, D[i][j-1]+1);
	
	return D;
}

SDoublePlane calculate_F(SDoublePlane D, SDoublePlane &binary_template, double threshold) {	
	double sum = 0.0;
	int m = binary_template.rows(), n = binary_template.cols();
	int i, j, k, l, input_rows = D.rows(), input_cols = D.cols();
	
	SDoublePlane output(input_rows, input_cols);
	for(i=0; i<input_rows - m; i++){
		for(j=0; j<input_cols - n; j++){
			sum=0.0;
			for(k =0; k<m; k++){
				for(l=0; l<n; l++){
					sum += binary_template[k][l] * D[i+k][j+l];
				}
			}
			output[i][j] = sum;
		}
	}
	//converting low scores to white n others black
	for(i=0; i<input_rows - m; i++)
		for(j=0; j<input_cols - n; j++)
			if(output[i][j] < threshold){
				output[i][j] = 255;
			}
			else
				output[i][j] = 0;
			
	return output;
}

HammingDistances find_hamming_distance(SDoublePlane &img_input, SDoublePlane &img_template){
	
	int input_rows = img_input.rows();
	int input_cols = img_input.cols();
	int template_rows = img_template.rows();
	int template_cols = img_template.cols();
	HammingDistances hm;
		
	double sum = 0.0;
	double max_hamming = 0.0;
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
			if(sum>max_hamming)
				max_hamming = sum;
		}
	}
	//printImg2File("input_img_file_Output.txt", img_input);
	//printImg2File("img2fileOutput.txt", output);
	hm.hamming_matrix = output;
	hm.max_hamming_distance= max_hamming;
	return hm;	
}

void set_line_tags(LineLocation &lineLoc, int count){
	
	switch(count){
		case 1:
		lineLoc.above_marker = 'G';
		lineLoc.line_marker = 'F';
		lineLoc.below_marker = 'E';
		break;
		
		case 2:
		lineLoc.line_marker = 'D';
		lineLoc.below_marker = 'C';
		break;
			
		case 3:
		lineLoc.line_marker = 'B';
		lineLoc.below_marker = 'A';
		break;
		
		case 4:
		lineLoc.line_marker = 'G';
		lineLoc.below_marker = 'F';
		break;

		case 5:
		lineLoc.line_marker = 'E';
		lineLoc.below_marker = 'D';
		break;

		case 6:
		lineLoc.above_marker = 'B';
		lineLoc.line_marker = 'A';
		lineLoc.below_marker = 'G';
		break;

		case 7:
		lineLoc.line_marker = 'F';
		lineLoc.below_marker = 'E';
		break;

		case 8:
		lineLoc.line_marker = 'D';
		lineLoc.below_marker = 'C';
		break;

		case 9:
		lineLoc.line_marker = 'B';
		lineLoc.below_marker = 'A';
		break;
		
		case 10:
		lineLoc.line_marker = 'G';
		lineLoc.below_marker = 'F';
		break;
		
	}
}

// Finding Line Location for Q4. As this was a open ended question, we found the number of lines using this approach.
vector<LineLocation> find_line_location(SDoublePlane &input){
	
	int rows = input.rows();
	int cols = input.cols();
	vector<LineLocation> allLinesLocVector;
	double sum;
	int lineCounter = 1;
	int j;
	
	for(int i=0;i<rows;i++)
	{	sum =0;
		for(j=cols - 45;j<cols -40;j++){
			sum += input[i][j]; 
		}
		if(sum/5 < 130){
			LineLocation lineLoc;
			lineLoc.row = i;
			lineLoc.col = j;
			set_line_tags(lineLoc, lineCounter);
			lineCounter++;
			allLinesLocVector.push_back(lineLoc);
		}
	
	}
	return allLinesLocVector;
}

void set_symbol_marker(DetectedSymbol &s, vector<LineLocation> allLinesLocVector){
	
	int symbol_start_row = s.row;
	int symbol_end_row = symbol_start_row + 7;
		
	// Finding avg distance between staff lines
	double avg_dist_staff_lines,tmpSum;
	
	for(int i=0; i<4; i++){
		tmpSum = allLinesLocVector[i+1].row - allLinesLocVector[i].row;
	}		
	avg_dist_staff_lines = tmpSum/4 + 6 ;
	
	// Marking the pitch
	for(int i=0; i<allLinesLocVector.size(); i++){
		
		// Boundary lines of both staves
		
		// Top line of High Pitch Staff
		if (i==0){
			if ( symbol_end_row < allLinesLocVector[i].row && symbol_start_row > allLinesLocVector[i].row - avg_dist_staff_lines){
				s.pitch = 'G';
				continue;
				}
		}
		
		// Bottom line of High Pitch Staff
		else if (i==4 && symbol_end_row < allLinesLocVector[5].row){
			if ( symbol_start_row >= allLinesLocVector[i].row && symbol_end_row <= allLinesLocVector[i].row + avg_dist_staff_lines  ){
				s.pitch = 'D';
				continue;							
			}
			
			else if ( symbol_start_row >= allLinesLocVector[i].row + avg_dist_staff_lines && symbol_end_row >= allLinesLocVector[i].row + s.height ){
				s.pitch = 'B';
				continue;
			}
			
			else if ( symbol_start_row >= allLinesLocVector[i].row && symbol_end_row >= allLinesLocVector[i].row + avg_dist_staff_lines ){
				s.pitch = 'C';
				continue;
			}

		}

		// Top line of Low Pitch Staff		
		else if (i==5 && symbol_start_row > allLinesLocVector[4].row + 2*s.height && symbol_end_row < allLinesLocVector[i].row){
			
			if ( symbol_end_row < allLinesLocVector[i].row && symbol_start_row <= allLinesLocVector[i].row - avg_dist_staff_lines - 5 ){
				s.pitch = 'C';
				continue;
			}
			
			else if ( symbol_end_row < allLinesLocVector[i].row && symbol_start_row > allLinesLocVector[i].row - avg_dist_staff_lines ){
				s.pitch = 'B';
				continue;
			}
						
		}
		
		// Bottom line of Low Pitch Staff
		else if (i==9 && symbol_start_row > allLinesLocVector[i].row){
			if ( symbol_start_row > allLinesLocVector[i].row && symbol_end_row <= allLinesLocVector[i].row + avg_dist_staff_lines  ){
				s.pitch = 'F';	
				continue;
			}			
			else if ( symbol_start_row > allLinesLocVector[i].row && symbol_end_row > allLinesLocVector[i].row + avg_dist_staff_lines ){
				s.pitch = 'E';
				continue;
			}				
		}
		
		else{
			
			if((allLinesLocVector[i].row < symbol_end_row) && (allLinesLocVector[i].row>symbol_start_row )){
				s.pitch = allLinesLocVector[i].line_marker;
				continue;
			}
		
		
			else if(symbol_start_row>=allLinesLocVector[i].row && symbol_end_row <= allLinesLocVector[i+1].row)
			{
				s.pitch = allLinesLocVector[i].below_marker;
				continue;
			}
		} 	
	}
}

// Find symbols in the given image for the given template
void find_symbols(HammingDistances hm, SDoublePlane input, SDoublePlane img_template, vector <DetectedSymbol> &symbols, Type symbol_type){
	
	int template_rows = img_template.rows();
	int template_cols = img_template.cols();
	
	double max_hamming_distance = hm.max_hamming_distance;
	SDoublePlane matrix = hm.hamming_matrix;
	
	vector<LineLocation> allLinesLocVector;
	double confidence_threshold;
	if (symbol_type == NOTEHEAD){
		confidence_threshold = 0.9;
		allLinesLocVector = find_line_location(input);
	}
	else
		confidence_threshold = 0.95;
		
	// Finding symbols
	for(int i =0;i<matrix.rows();i++){
		for(int j=0;j<matrix.cols();j++){
			double value = matrix [i][j];
			if( value >= confidence_threshold * max_hamming_distance) {
				DetectedSymbol s;
				s.row = i;
				s.col = j;
				s.width = template_cols;
				s.height = template_rows;
				s.type = symbol_type;
				s.confidence = value;
				
				s.pitch = ' ';
				
				if(symbol_type == NOTEHEAD)
					set_symbol_marker(s, allLinesLocVector);
				else
					s.pitch = ' ';
				
				symbols.push_back(s);
				// Marking the pixels of the template so that they are not detected again
				for (int x=i; x < i+s.height; x++){
					for (int y=j; y < j+s.width; y++){
						matrix[x][y] = -100;
					}	
				}
				
			}

		}
	}
}

// Hough Transform
// Inspired by Lecture from Prof. William Hoff, Colorado School of Mines, Engineering Division 
// https://www.youtube.com/watch?v=o-n7NoLArcs and http://goo.gl/TxbGWi
// 
SDoublePlane runHoughTransform(SDoublePlane &img){
	
	//printImg2File("sobelPNG.txt",img);

	int r = img.rows();
	int c = img.cols();
	double hough_height;

	// Initialize Accumulator matrix
	if (r>c)
		hough_height = r / sqrt(2);
	else
		hough_height = c / sqrt(2);
	
	int maxDist = round(hough_height * 2);
	int theta = 180;
	double rho,t;
	int iRho;
	
	SDoublePlane H(maxDist, theta);
	
	for (int i = 0; i < maxDist; i++)
	  for (int j = 0; j < theta; j++)
		H[i][j] = 0;
	
	
	double center_x = c/2;
	double center_y = r/2;
	
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (img[i][j] > 250){
				
				// Fill accumulator array
				for (int iTheta = 0 ; iTheta < theta ; iTheta++){
					
					// Getting angle in radians
					t = iTheta * M_PI / 180;
					
					// Calculate distance from origin for this angle				
					rho = ( j - center_x) * cos (t) + ( i - center_y) * sin(t);
					iRho = int(round(rho + hough_height));
					H[iRho][iTheta]++;
				}
			}
		}
	}
	return H;
}

// Find lines from Hough Accumulator matrix greater than threshold (percentage against max votes vaue)
vector<Line>  getLinesFromHoughSpace(SDoublePlane &accumulator, SDoublePlane &img, double threshold){
	
	int img_w = img.cols();
	int img_h = img.rows();
    int rows = accumulator.rows();
    int cols = accumulator.cols();
    int windowSize = 4;
    double maxValue = 0;
    vector<Line> houghLines;
	
	// Find max vote value
	for(int rho=0;rho<rows;rho++){
		for(int theta =0; theta<cols;theta++){
			if (accumulator[rho][theta] > maxValue)
				maxValue = accumulator[rho][theta] ;
		}
	}
    
	// Printing accumulator matrix
    //printImg2File("Matrix.txt", accumulator);

    for(int rho=0;rho<rows;rho++){
		for(int theta =0; theta<cols;theta++){
        
			if (accumulator[rho][theta] >= threshold * maxValue){
				
				int max = accumulator[rho][theta];
				for (int wx=-windowSize; wx<=windowSize; wx++){
					for (int wy=-windowSize; wy<=windowSize; wy++){
						if ((wx+rho >= 0 && wx+rho < rows) && (wy+theta>=0 && wy+theta<cols)){
							if ((int)accumulator[rho+wx][theta+wy] > max){
								max = accumulator[rho+wx][theta+wy];
								wx += 2;
								wy += 2;
							}
						}
					}
				}
				
				if (max > accumulator[rho][theta])
					continue;
				
				Line line;
				double degInRadian = DEG2RAD(theta);
				// Transform the selected lines back to Cartesian Space
				if (theta >= 45 && theta <= 135){
					// y = (rho - x cos(theta)) / sin(theta)
					line.x1 = 0;
					line.y1 = ((rho - rows/2) - ((line.x1 - img_w/2) * cos( degInRadian))) / sin(degInRadian) + img_h/2;					
					line.x2 = img_w - 0;
					line.y2 = ((rho - rows/2) - ((line.x2 - img_w/2) * cos(degInRadian))) / sin(degInRadian) + img_h/2;
			
				} else {
					
					// x = (rho - y sin(theta)) / cos(theta)
					line.y1 = 0;
					line.x1 = ((rho - rows/2) - ((line.y1 - img_h/2) * sin(degInRadian))) / cos(degInRadian) + img_w/2;
					line.y2 = img_h - 0;
					line.x2 = ((rho - rows/2) - ((line.y2 - img_h/2) * sin(degInRadian))) / cos(degInRadian) + img_w/2;
				}
				
				
				if (houghLines.empty()) houghLines.push_back(line);
				else{
					if (abs(line.y1 - (houghLines.back().y1)) > 5) houghLines.push_back(line);
				}
			}
        }
    }
    return houghLines;
}


double getAvgDistanceBetweenStaffLines(vector<Line> linesVector){
	
	int numberOfLines = linesVector.size();
	double allStaffHeightSum=0,avgDistanceBetweenStaffLines=0;
	int numberOfStaff = 0;
	if (numberOfLines % 5 == 0){
		for (int i=0; i<numberOfLines; i+=5){
			allStaffHeightSum += sqrt(pow(linesVector[i+4].x1-linesVector[i].x1,2)+pow(linesVector[i+4].y1-linesVector[i].y1,2));
			numberOfStaff++;
		}
		avgDistanceBetweenStaffLines = allStaffHeightSum/(double)numberOfStaff;
	}else{
		avgDistanceBetweenStaffLines = 0;
	}	
	
	if (avgDistanceBetweenStaffLines == 0)
		avgDistanceBetweenStaffLines = 50;
	
	// Avg distance between every staff line
	return avgDistanceBetweenStaffLines/4.0;
	
}

void write_staves_image(const string &filename, const SDoublePlane &img, vector<Line> linesVector){
	
	SDoublePlane output_planes[3];
	
	for (int i = 0; i < 3; i++)
		output_planes[i] = img;
					
	int r = img.rows();
	int c = img.cols();
	
	int x1, y1, x2, y2;
	for (int i = 0; i < linesVector.size(); i++) {
		x1 = linesVector[i].x1;
		y1 = linesVector[i].y1;
		x2 = linesVector[i].x2;
		y2 = linesVector[i].y2;
		//cout<<"(x1,y1): ("<<x1<<","<<y1<<"), (x2,y2):("<<x2<<","<<y2<<")\n"; 
		
		overlay_rectangle(output_planes[2], y1, x1, y2, x2, 255, 2);
		overlay_rectangle(output_planes[0], y1, x1, y2, x2, 0, 2);
		overlay_rectangle(output_planes[1], y1, x1, y2, x2, 0, 2);

	}
	
	SImageIO::write_png_file(filename.c_str(), output_planes[0],
			output_planes[1], output_planes[2]);
		
}


// Resize image
SDoublePlane resize_image(SDoublePlane &input, double newScaleRatio)
{
	if (newScaleRatio == 1) return input;
	
    int rows = input.rows();
    int cols = input.cols();
 
	int newWidth = newScaleRatio * cols;
	int newHeight = newScaleRatio * rows ;
	
    SDoublePlane output(newHeight, newWidth);
    
    set<int> mySetRows;
    
    while( mySetRows.size() < newHeight ){
        mySetRows.insert( rand()% rows +1 );
    }
    
    set<int> mySetCols;
    
    while( mySetCols.size() < newWidth ){
        mySetCols.insert( rand()% cols+1 );
    }
    
    int ii =0;
    for(int i = 1; i <rows;i++){
        if(mySetRows.find(i) == mySetRows.end()){
            continue;
        }
        int jj=0;
        for(int j = 1; j<cols;j++){
            if(mySetCols.find(j) == mySetCols.end())
            {
                continue;
            }
            
            
            if(i<rows-1 && j<cols-1){
                
                output[ii][jj] = (input[i][j] + input [i+1][j] + input[i][j+1] + input[i+1][j+1])/4;
            }
            jj++;
        }
        ii++;
        
    }
    
    //write_detection_image("Resized_pic.png", output);
    return output;
}

void detectSymbolsHammingDistance(const SDoublePlane& img,
		const SDoublePlane& col_filter, const SDoublePlane& row_filter,  
		const SDoublePlane& template_notehead,
		const SDoublePlane& template_quarterrest,
		const SDoublePlane& template_eighthrest, const string &filename){
		
	cout<<"Convoluting...\n";		
	SDoublePlane convoluted_image = convolve_separable(img, row_filter, col_filter);
	SDoublePlane convoluted_notehead_template = convolve_separable(template_notehead, row_filter, col_filter);
	SDoublePlane convoluted_quarterrest_template = convolve_separable(template_quarterrest, row_filter, col_filter);
	SDoublePlane convoluted_eighthrest_template = convolve_separable(template_eighthrest, row_filter, col_filter);
	
	
	// Convert to Binary
	cout<<"Converting to Binary...\n";		
	SDoublePlane binary_convoluted_image = convert_binary(convoluted_image);
	SDoublePlane binary_convoluted_notehead_template = convert_binary(convoluted_notehead_template);
	SDoublePlane binary_convoluted_quarterrest_template = convert_binary(convoluted_quarterrest_template);
	SDoublePlane binary_convoluted_eighthrest_template = convert_binary(convoluted_eighthrest_template);
	
	// Q4. Detecting symbols using Hamming Distances 
	cout<<"Applying HammingDistance approach...\n";
	HammingDistances hm_notehead = find_hamming_distance(binary_convoluted_image, binary_convoluted_notehead_template);
	HammingDistances hm_quarterrest = find_hamming_distance(binary_convoluted_image, binary_convoluted_quarterrest_template);
	HammingDistances hm_eighthrest = find_hamming_distance(binary_convoluted_image, binary_convoluted_eighthrest_template);
	
	cout<<"Detecting symbols...\n";
	vector <DetectedSymbol> symbols;
	find_symbols(hm_notehead, img, binary_convoluted_notehead_template, symbols, NOTEHEAD);
	find_symbols(hm_quarterrest, img, binary_convoluted_quarterrest_template, symbols, QUARTERREST);
	find_symbols(hm_eighthrest, img, binary_convoluted_eighthrest_template, symbols, EIGHTHREST);

	write_detection_image((filename+".png").c_str(), symbols, img);
	if(!filename.compare("detected7"))	write_detection_txt((filename+".txt").c_str(), symbols);
	cout<<"Output written: "<<filename<<".png\n";

}

void q5_score(SDoublePlane& input_image, SDoublePlane& mean_filter, SDoublePlane& template_notehead, SDoublePlane& template_quarterrest, SDoublePlane& template_eighthrest) {
	int rows = input_image.rows(), cols = input_image.cols();
	vector <DetectedSymbol> symbols;
	
	SDoublePlane sobelPNG = sobel_filter(input_image);
	SDoublePlane image_sobel = direct_sobel(input_image);
	write_detection_image("edges.png", symbols, convert_BW(image_sobel));
	
	SDoublePlane image_basic_sobel = sobel_filter(input_image);
	SDoublePlane image_basic_sobel_blur = sobel_filter(image_basic_sobel);
	SDoublePlane image_basic_sobel_blur_alternate = alternate_edge(image_basic_sobel_blur);	//direct sobel
	
	SDoublePlane image_blur = convolve_general(input_image, mean_filter);
	SDoublePlane image_blur_sobel = sobel_filter(image_blur);
	SDoublePlane binary_image_blur_sobel = convert_binary(image_blur_sobel);
	SDoublePlane image_edge_map = binary_image_blur_sobel;
	
	SDoublePlane template_notehead_blur = convolve_general(template_notehead, mean_filter);
	SDoublePlane template_notehead_blur_sobel = sobel_filter(template_notehead_blur);
	SDoublePlane binary_template_notehead_blur_sobel = convert_binary(template_notehead_blur_sobel);
	
	SDoublePlane template_quarterrest_blur = convolve_general(template_quarterrest, mean_filter);
	SDoublePlane template_quarterrest_blur_sobel = sobel_filter(template_quarterrest_blur);
	SDoublePlane binary_template_quarterrest_blur_sobel = convert_binary(template_quarterrest_blur_sobel);
	
	SDoublePlane template_eighthrest_blur = convolve_general(template_eighthrest, mean_filter);
	SDoublePlane template_eighthrest_blur_sobel = sobel_filter(template_eighthrest_blur);
	SDoublePlane binary_template_eighthrest_blur_sobel = convert_binary(template_eighthrest_blur_sobel);
	
	//Sobel on Image
	SDoublePlane binary_image = convert_binary(input_image);
	SDoublePlane binary_image_sobel = sobel_filter(binary_image);
	SDoublePlane binary_image_sobel_blur = convolve_general(binary_image_sobel, mean_filter);
	SDoublePlane binary_binary_image_sobel_blur = convert_blur_to_binary(binary_image_sobel_blur);
	SDoublePlane bw_binary_binary_image_sobel_blur = convert_binary_to_BW(binary_binary_image_sobel_blur);
	
	//Sobel on template_notehead
	SDoublePlane binary_template_notehead = convert_binary(template_notehead);
	SDoublePlane binary_template_notehead_sobel = sobel_filter(binary_template_notehead);
	SDoublePlane bw_binary_template_notehead_sobel = convert_binary_to_BW(binary_template_notehead_sobel);
	SDoublePlane binary_template_notehead_sobel_blur = convolve_general(binary_template_notehead_sobel, mean_filter);
	SDoublePlane binary_binary_template_notehead_sobel_blur = convert_blur_to_binary(binary_template_notehead_sobel_blur);
	
	//Sobel on template_quarterrest
	SDoublePlane binary_template_quarterrest = convert_binary(template_quarterrest);
	SDoublePlane binary_template_quarterrest_sobel = sobel_filter(binary_template_quarterrest);
	SDoublePlane bw_binary_template_quarterrest_sobel = convert_binary_to_BW(binary_template_quarterrest_sobel);
	SDoublePlane binary_template_quarterrest_sobel_blur = convolve_general(binary_template_quarterrest_sobel, mean_filter);
	SDoublePlane binary_binary_template_quarterrest_sobel_blur = convert_blur_to_binary(binary_template_quarterrest_sobel_blur);

	//Sobel on template_eighthrest
	SDoublePlane binary_template_eighthrest = convert_binary(template_eighthrest);
	SDoublePlane binary_template_eighthrest_sobel = sobel_filter(binary_template_eighthrest);
	SDoublePlane bw_binary_template_eighthrest_sobel = convert_binary_to_BW(binary_template_eighthrest_sobel);
	SDoublePlane binary_template_eighthrest_sobel_blur = convolve_general(binary_template_eighthrest_sobel, mean_filter);
	SDoublePlane binary_binary_template_eighthrest_sobel_blur = convert_blur_to_binary(binary_template_eighthrest_sobel_blur);
	
	//calculate D
	SDoublePlane D(rows, cols);
	D = calculate_D(binary_binary_image_sobel_blur);
	
	//calculate score function 'f'
	SDoublePlane F_template_notehead(rows, cols), F_template_quarterrest(rows, cols), F_template_eighthrest(rows, cols);
	F_template_notehead = calculate_F(D, binary_binary_template_notehead_sobel_blur, 3);
	F_template_quarterrest = calculate_F(D, binary_binary_template_quarterrest_sobel_blur, 5);
	F_template_eighthrest = calculate_F(D, binary_binary_template_eighthrest_sobel_blur, 2);
	
	//find symbols (using the class HammingDistances in order to match the fuctioning of find_symbols)
	HammingDistances hm_F_template_notehead, hm_F_template_quarterrest, hm_F_template_eighthrest;
	hm_F_template_notehead.hamming_matrix = F_template_notehead;
	hm_F_template_quarterrest.hamming_matrix = F_template_quarterrest;
	hm_F_template_eighthrest.hamming_matrix = F_template_eighthrest;
	hm_F_template_notehead.max_hamming_distance = hm_F_template_quarterrest.max_hamming_distance = hm_F_template_eighthrest.max_hamming_distance = 255;
	find_symbols(hm_F_template_notehead, input_image, binary_binary_template_notehead_sobel_blur, symbols, NOTEHEAD);
	find_symbols(hm_F_template_quarterrest, input_image, binary_binary_template_quarterrest_sobel_blur, symbols, QUARTERREST);
	find_symbols(hm_F_template_eighthrest, input_image, binary_binary_template_eighthrest_sobel_blur, symbols, EIGHTHREST);
	
	write_detection_image("detected5.png", symbols, input_image);
	cout<<"Output written: detected5.png\n\n";
}

int main(int argc, char *argv[]) {
	if (!(argc == 2)) {
		cerr << "usage: " << argv[0] << " input_image" << endl;
		return 1;
	}
	
//******************** Initialize ***************************	
	string input_filename(argv[1]);
	string TEMPLATE_NOTEHEAD = "template1.png";
	string TEMPLATE_QUARTERREST = "template2.png";
	string TEMPLATE_EIGHTHREST = "template3.png";
	int filter_size = 3;

	SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());
	SDoublePlane template1 = SImageIO::read_png_file("template1.png");
													
	// Read the image, template NOTEHEAD, QUARTERREST AND EIGHTHREST
	cout<<"\nReading image and templates...\n";
	SDoublePlane template_notehead = SImageIO::read_png_file(TEMPLATE_NOTEHEAD.c_str());
	SDoublePlane template_quarterrest = SImageIO::read_png_file(TEMPLATE_QUARTERREST.c_str());
	SDoublePlane template_eighthrest = SImageIO::read_png_file(TEMPLATE_EIGHTHREST.c_str());
	
	SDoublePlane row_filter(1, filter_size);
	row_filter = make_row_filter(filter_size);
	
	SDoublePlane col_filter(filter_size, 1);
	col_filter = make_col_filter(filter_size);
	
	SDoublePlane mean_filter(filter_size, filter_size);
	mean_filter = make_mean_filter(filter_size); 
	
//******************** Q2 - Begins ***************************	
	SDoublePlane convoluted_image_2 = convolve_general(input_image, mean_filter);
	//write_detection_image("output_q2.png", convoluted_image_2);
//******************** Q2 - Ends ***************************

	
//******************** Q3 - Begins ***************************
	SDoublePlane convoluted_image_3 = convolve_separable(input_image, row_filter, col_filter);
	//write_detection_image("output_q3.png", convoluted_image_3);
//******************** Q3 - Ends ***************************


//******************** Q4 - Begins ***************************
	// Writing scores4.png //
	// The convolution step is done here to get scores4.png. However, the convolution step is performed inside detectSymbolsHammingDistance	
	SDoublePlane convoluted_image = convolve_general(input_image, mean_filter);
	write_detection_image("scores4.png", convoluted_image);
	cout<<"Output written: scores4.png\n\n";
	detectSymbolsHammingDistance(input_image, col_filter, row_filter, template_notehead, template_quarterrest, template_eighthrest, "detected4");
//******************** Q4 - Ends ***************************	
	
	
//******************** Q5 - Begins ***************************	
	q5_score(input_image, mean_filter, template_notehead, template_quarterrest, template_eighthrest);
//******************** Q5- Ends ***************************	
	
	
//******************** Q6 and Q7 - Begin *****************
	SDoublePlane sobelInputPNG = sobel_filter(input_image);	
	SDoublePlane hough_transform_accu = runHoughTransform(sobelInputPNG);
	vector<Line> linesFromHoughSpace = getLinesFromHoughSpace(hough_transform_accu, sobelInputPNG, 0.75);
	
	write_staves_image("staves.png", sobelInputPNG,linesFromHoughSpace);
	
	// Find Scaling factor
	double avgDistBetweenEveryStaffLine = getAvgDistanceBetweenStaffLines(linesFromHoughSpace);
	double newScaleRatio = template1.rows()/(avgDistBetweenEveryStaffLine - 0.5);
	
	cout<<"Rescaling..\n";	
	SDoublePlane rescaledImage = resize_image(input_image, newScaleRatio);		

	detectSymbolsHammingDistance(input_image, col_filter, row_filter, template_notehead, template_quarterrest, template_eighthrest, "detected7");	
	//******************** Q6 and Q7 - End ****************
	
}

