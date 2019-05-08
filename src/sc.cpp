
#include "sc.h"
#include "math.h"

using namespace cv;
using namespace std;


struct path{
    int x;
    int y;
};

/*calculate the energy for each pixel and get the Mat_energy */
void energy_calculating(int rows, int cols, Mat& Mat_gray, Mat& Mat_energy) {
	//left up---x:right---y:down
	float leftUp_x = Mat_gray.at<float>(0, 1);
	float leftUp_y = Mat_gray.at<float>(1, 0);
	Mat_energy.at<float>(0, 0) = sqrt(leftUp_x*leftUp_x + leftUp_y * leftUp_y);

	//right up---x:left---y:down
	float rightUp_x = Mat_gray.at<float>(0, cols-2);
	float rightUp_y = Mat_gray.at<float>(1, cols-1);
	Mat_energy.at<float>(0, cols - 1) = sqrt(rightUp_x*rightUp_x + rightUp_y * rightUp_y);

	//left down---x:right---y:up
	float leftDown_x = Mat_gray.at<float>(rows - 1, 1);
	float leftDown_y = Mat_gray.at<float>(row - 2, 0);
	Mat_energy.at<float>(rows - 1, 0) = sqrt(leftDown_x*leftDown_x + leftDown_y * leftDown_y);

	//right down---x:left---y:up
	float rightDown_x = Mat_gray.at<float>(rows - 1, cols - 2);
	float rightDown_y = Mat_gray.at<float>(rows - 2, cols - 1);
	Mat_energy.at<float>(rows - 1, cols - 1) = sqrt(rightDown_x*rightDown_x + rightDown_y * rightDown_y);

	//first row without above corner
	for (int i = 1; i < cols - 1; i++) {
		float x = Mat_gray.at<float>(0, i + 1) - Mat_gray.at<float>(0, i - 1);
		float y = Mat_gray.at<float>(1, i);
		Mat_energy.at<float>(0, i) = sqrt(x*x + y * y);
	}

	//last row without above corner
	for (int i = 1; i < cols - 1; i++) {
		float x = Mat_gray.at<float>(rows - 1, i + 1) - Mat_gray.at<float>(rows - 1, i - 1);
		float y = Mat_gray.at<float>(rows - 2, i);
		Mat_energy.at<float>(rows - 1, i) = sqrt(x*x + y * y);
	}

	//first column withouot corner
	for (int i = 1; i < rows - 1; i++) {
		float x = Mat_gray.at<float>(i, 1);
		float y = Mat_gray.at<float>(i + 1, 0) - Mat_gray.at<float>(i - 1, 0);
		Mat_energy.at<float>(i, 0) = sqrt(x*x + y * y);
	}

	//last column without corner
	for (int i = 1; i < rows - 1; i++) {
		float x = Mat_gray.at<float>(i, cols - 2);
		float y = Mat_gray.at<float>(i + 1, cols - 1) - Mat_gray.at<float>(i-1,cols-1);
		Mat_energy.at<float>(i,cols-1) = sqrt(x*x + y * y);
	}

	//middle part
	for (int i = 1; i < rows - 1; i++) {
		for (int j = 1; j < cols - 1; j++) {
			float x = Mat_gray.at<float>(i, j + 1) - Mat_gray.at<float>(i, j - 1);
			float y = Mat_gray.at<float>(i + 1, j) - Mat_gray.at<float>(i - 1, j);
			Mat_energy.at<float>(i,j) = sqrt(x*x + y * y);
		}
	}
}

bool seam_carving(Mat& in_image, int new_width, int new_height, Mat& out_image){

    // some sanity checks
    // Check 1 -> new_width <= in_image.cols
    if(new_width>in_image.cols){
        cout<<"Invalid request!!! new_width has to be smaller than the current size!"<<endl;
        return false;
    }
    if(new_height>in_image.rows){
        cout<<"Invalid request!!! ne_height has to be smaller than the current size!"<<endl;
        return false;
    }
    
    if(new_width<=0){
        cout<<"Invalid request!!! new_width has to be positive!"<<endl;
        return false;

    }
    
    if(new_height<=0){
        cout<<"Invalid request!!! new_height has to be positive!"<<endl;
        return false;
        
    }

    
    return seam_carving_trivial(in_image, new_width, new_height, out_image);
}



// seam carves by removing trivial seams
bool seam_carving_trivial(Mat& in_image, int new_width, int new_height, Mat& out_image){
   
    Mat iimage = in_image.clone();
    Mat oimage = in_image.clone();
	
    //grayscale mat
    Mat Mat_gray = Mat(iimage.rows, iimage.cols, CV_32FC1);
    //cvtColor( in_image, Mat_gray, CV_BGR2GRAY );
   
    for(int i = 0;i < iimage.rows;i++){
	for(int j = 0;j < iimage.cols;j++){
	    Vec3b pixel = in_image.at<Vec3b>(i, j);
            Mat_gray.at<float>(i, j) = pixel[0] * 0.2126 + pixel[1] * 0.7152 + pixel[2] * 0.0722;
	}
    }
    
    //imshow("testGray",Mat_gray);

  /* pass the Mat_gray to the following 2 function to get seams and cut*/

    
    while(iimage.rows!=new_height || iimage.cols!=new_width){
        
	// horizontal seam if needed
        if(iimage.rows > new_height){
	    reduce_horizontal_seam_trivial(iimage, oimage, Mat_gray);
	    //update
            iimage = oimage.clone();
        }
		
        // vertical
        if(iimage.cols > new_width){
	    reduce_vertical_seam_trivial(iimage, oimage, Mat_gray);
            iimage = oimage.clone();
        }
    }
    
    // maybe need to change it later--done
    out_image = oimage.clone();
    return true;
}

// horizontal trivial seam is a seam through the center of the image
bool reduce_horizontal_seam_trivial(Mat& in_image, Mat& out_image, Mat& Mat_gray){

    // retrieve the dimensions of the new image
    int rows = in_image.rows;
    int cols = in_image.cols;

    // the out_image minus one column(Redefine)
    out_image = Mat(rows-1, cols, CV_8UC3);

    //Store the path of the min seam in Mat_S, set the flag: 0--left up,1--left.2--left down.
    Mat Mat_S = Mat(rows,cols,CV_8UC1);

    /*After each cut get the lastest energy for each pixel*/
    Mat Mat_energy = Mat(rows, cols, CV_32FC1);

    //update the gray mat
    /*
    for(int i = 0;i < rows;i++){
	for(int j = 0;j < cols;j++){
	    Vec3b pixel = in_image.at<Vec3b>(i, j);
            Mat_gray.at<float>(i, j) = pixel[0] * 0.2126 + pixel[1] * 0.7152 + pixel[2] * 0.0722;
	}
    }
    */	

    energy_calculating(rows, cols, Mat_gray, Mat_energy);

    //new Mat_energy for store cumulative num
    Mat Mat_energyNew = Mat(rows, cols, CV_32FC1);

    /* get the left-right cumulative minimum energy*/

    // the column 0
    for(int row = 0;row < rows;row++){
	//Mat_energy first column value do not change.
        Mat_energyNew.at<float>(row,0) = Mat_energy.at<float>(row,0);
	Mat_S.at<char>(row,0) = -1;
    }
    // column 1 to column -1
    for(int i = 1;i < cols; i++){

	//row 0--rows -1
    	for(int j = 0;j < rows;j++){
	    if(j == 0){
	   	float min = Mat_energyNew.at<float>(j,i-1);
	    	Mat_S.at<char>(j,i) = 1;

	    	if(min > Mat_energyNew.at<float>(j+1,i-1)){
	        	min = Mat_energyNew.at<float>(j+1,i-1);
			Mat_S.at<char>(j,i) = 2;
	    	}
	    	Mat_energyNew.at<float>(j,i) = Mat_energy.at<float>(j,i) + min;
	    }
	    //last row
    	    else if(j == rows-1){
	    	float min = Mat_energyNew.at<float>(j-1,i-1);
	   	 Mat_S.at<char>(j,i) = 0;
		
	    	if(min > Mat_energyNew.at<float>(j,i-1)){
	        	min = Mat_energyNew.at<float>(j,i-1);
			Mat_S.at<char>(j,i) = 1;
	    	}

	  	Mat_energyNew.at<float>(j,i) = Mat_energy.at<float>(j,i) + min;
            }
	    else{
		float min = Mat_energyNew.at<float>(j-1,i-1);
            	Mat_S.at<char>(j,i) = 0;
		
	    	if(min > Mat_energyNew.at<float>(j,i-1)){
	        	min = Mat_energyNew.at<float>(j,i-1);
	        	Mat_S.at<char>(j,i) = 1;
	    	}
		
            	if(min > Mat_energyNew.at<float>(j+1,i-1)){
			min = Mat_energyNew.at<float>(j+1,i-1);
	        	Mat_S.at<char>(j,i) = 2;
	    	}
 
      	   	Mat_energyNew.at<float>(j,i) = Mat_energy.at<float>(j,i) + min;
	    }	
	    
   	}
	    
    }

 
    /*traverse the last column to get the minimum entry*/

    //min index for row 
    int minIndex = 0;
    float min = Mat_energyNew.at<float>(0,cols-1);

    for(int j = 1;j < rows;j++){
	if(min > Mat_energyNew.at<float>(j,cols-1)){
	    min = Mat_energyNew.at<float>(j,cols-1);
	    minIndex = j;
	}
    }
    //cout << minIndex << endl;

    /*backtrack from this minimum entry*/

    //Store the path
    //vector<path> vecPath;
    //vecPath.resize(rows);

    int j = minIndex;
    int i = cols - 1;

    
    //remove the pixel from each column
    while (i>=0) {
	//each row need to delete 1 celling,so column minus 1.record the delete x,y.
	//vecPath[i].x = i;
	//vecPath[i].y = j;

	for(int imgrow = 0;imgrow < rows; imgrow++){
	    if(imgrow < j){
		out_image.at<Vec3b>(imgrow, i) = in_image.at<Vec3b>(imgrow,i);
	    }
	    else if(imgrow > j){
		out_image.at<Vec3b>(imgrow-1, i) = in_image.at<Vec3b>(imgrow, i);
		//remove pixel from the Mat gray
		Mat_gray.at<float>(imgrow-1, i) = Mat_gray.at<float>(imgrow, i);
	    }
	}
	//left-up row j
	if (Mat_S.at<char>(j,i) == 0) {
		j = j - 1;
	}
	//left
	else if (Mat_S.at<char>(j,i) == 1) {
		//j do not change
		j = j;
	}
	//left-down
	else if (Mat_S.at<char>(j,i) == 2) {
		j = j + 1;
	}
	i--;
    }
       

    /*
    //remove the pixel from each column  j---row
    //1- last
    for(int j = 0;j < rows;j++){
	if(j < minIndex){
	   out_image.at<Vec3b>(j, cols-1) = in_image.at<Vec3b>(j, cols-1); 
	}
	else if(j > minIndex){
	    out_image.at<Vec3b>(j-1, cols-1) = in_image.at<Vec3b>(j, cols-1);
            Mat_gray.at<float>(j-1, cols-1) = Mat_gray.at<float>(j, cols-1);
	}
    }
    //2- rest columns
    for(int i = cols-2; i >= 0;i--){
	if(minIndex == 0){
	    if(Mat_energyNew.at<float>(1,i) < Mat_energyNew.at<float>(0,i)){
		minIndex = 1;
	    }
	}
	//last row
	else if(minIndex == rows-1){
	    if(Mat_energyNew.at<float>(rows-2,i) < Mat_energyNew.at<float>(rows-1, i)){
		minIndex = rows-2;
	    }
	}
	else{
	    if(Mat_energyNew.at<float>(minIndex-1, i) <= Mat_energyNew.at<float>(minIndex, i) && Mat_energyNew.at<float>(minIndex-1, i) <= Mat_energyNew.at<float>(minIndex+1, i)){
		minIndex = minIndex-1;
	    }
	    else if(Mat_energyNew.at<float>(minIndex,i) > Mat_energyNew.at<float>(minIndex+1,i)){
		minIndex = minIndex+1;
	    }
	}

	for(int j = 0;j < rows;j++){
	    if(j < minIndex){
	    	out_image.at<Vec3b>(j, i) = in_image.at<Vec3b>(j, i); 
	    }
	    else if(j > minIndex){
	        out_image.at<Vec3b>(j-1, i) = in_image.at<Vec3b>(j, i);
                Mat_gray.at<float>(j-1, i) = Mat_gray.at<float>(j, i);
	    }
	}
    }
    */

    return true;
}

// vertical trivial seam is a seam through the center of the image
bool reduce_vertical_seam_trivial(Mat& in_image, Mat& out_image, Mat& Mat_gray){
	
    // retrieve the dimensions of the new image
    int rows = in_image.rows;
    int cols = in_image.cols;

    // the out_image minus one column(Redefine)
    out_image = Mat(rows, cols-1, CV_8UC3);

    //Store the path of the min seam in Mat_S, set the flag: 0--left up,1--up.2--right up.
    Mat Mat_S = Mat(rows,cols,CV_8UC1);  
     
    /*After each cut get the lastest energy for each pixel*/
    Mat Mat_energy = Mat(rows, cols, CV_32FC1);

    //update the gray mat
    /*
    for(int i = 0;i < rows;i++){
	for(int j = 0;j < cols;j++){
	    Vec3b pixel = in_image.at<Vec3b>(i, j);
            Mat_gray.at<float>(i, j) = pixel[0] * 0.2126 + pixel[1] * 0.7152 + pixel[2] * 0.0722;
	}
    }	
    */

    energy_calculating(rows, cols, Mat_gray, Mat_energy);

    //new Mat_energy for store cumulative num
    Mat Mat_energyNew = Mat(rows, cols, CV_32FC1);

    /* get the top-down cumulative minimum energy*/

    // the row 0
    for(int col = 0;col < cols;col++){
	Mat_energyNew.at<float>(0,col) = Mat_energy.at<float>(0,col);
	Mat_S.at<char>(0,col) = -1;
    }
    // row 1 to rows-1
    for(int i = 1;i < rows;i++){
	
	//column 0--cols-1
    	for(int j = 0;j < cols;j++){
	  
	    if(j == 0){
	        float min = Mat_energyNew.at<float>(i-1,j);
	        Mat_S.at<char>(i,j) = 1;

	        if(min > Mat_energyNew.at<float>(i-1,j+1)){
	             min = Mat_energyNew.at<float>(i-1,j+1);
	             Mat_S.at<char>(i,j) = 2;
	    	}
	    	Mat_energyNew.at<float>(i,j) = Mat_energy.at<float>(i,j) + min;
	    }

	    //last column
   	    else if(j == cols-1){
	        float min = Mat_energyNew.at<float>(i-1,j-1);
	        Mat_S.at<char>(i,j) = 0;
		
	        if(min > Mat_energyNew.at<float>(i-1,j)){
	            min = Mat_energyNew.at<float>(i-1,j);
		    Mat_S.at<char>(i,j) = 1;
	        }
	        Mat_energyNew.at<float>(i,j) = Mat_energy.at<float>(i,j) + min;
    	    }

	    else{
		float min = Mat_energyNew.at<float>(i-1,j-1);
            	Mat_S.at<char>(i,j) = 0;
		
	    	if(min > Mat_energyNew.at<float>(i-1,j)){
	       	    min = Mat_energyNew.at<float>(i-1,j);
	            Mat_S.at<char>(i,j) = 1;
	        }
		
                if(min > Mat_energyNew.at<float>(i-1,j+1)){
		    min = Mat_energyNew.at<float>(i-1,j+1);
	            Mat_S.at<char>(i,j) = 2;
	        }
      	        Mat_energyNew.at<float>(i,j) = Mat_energy.at<float>(i,j) + min;
	    }
	    
	}
		    
    }

 
    /*traverse the last row to get the minimum entry*/

    int minIndex = 0;
    float min = Mat_energyNew.at<float>(rows-1,0);

    for(int j = 1;j < cols;j++){
	if(min > Mat_energyNew.at<float>(rows-1,j)){
	    min = Mat_energyNew.at<float>(rows-1,j);
	    minIndex = j;
	}

    }
    //cout << minIndex << endl;

    /*backtrack from this minimum entry*/

    //Store the path
    //vector<path> vecPath;
    //vecPath.resize(rows);

    int j = minIndex;
    int i = rows - 1;

    
    while (i>=0) {
	//each row need to delete 1 celling,so column minus 1.record the delete x,y.
	//vecPath[i].x = i;
	//vecPath[i].y = j;

	for(int imgcol = 0;imgcol < cols; imgcol++){
	    if(imgcol < j){
		out_image.at<Vec3b>(i, imgcol) = in_image.at<Vec3b>(i,imgcol);
	    }
	    else if(imgcol > j){
		out_image.at<Vec3b>(i, imgcol-1) = in_image.at<Vec3b>(i,imgcol);
		//remove pixel from the Mat gray
		Mat_gray.at<float>(i, imgcol-1) = Mat_gray.at<float>(i,imgcol);
	    }
	}


	//left-top
	if (Mat_S.at<char>(i,j) == 0) {
		j = j - 1;
	}
	else if (Mat_S.at<char>(i,j) == 1) {
		//j do not change
		j = j;
	}
	else if (Mat_S.at<char>(i,j) == 2) {
		j = j + 1;
	}
	i--;
    }
    
    
    /*
    //remove the pixel for each row
    //1--last row
    for(int j = 0;j < cols; j++){
	if (j < minIndex) {
            out_image.at<Vec3b>(rows-1, j) = in_image.at<Vec3b>(rows-1, j);
        } 
	else if (j > minIndex) {
            out_image.at<Vec3b>(rows-1, j-1) = in_image.at<Vec3b>(rows-1, j);
            Mat_gray.at<float>(rows-1, j-1) = Mat_gray.at<float>(rows-1, j);
        }
    }
    
    //2--rest rows
    for(int i = rows-2 ;i >= 0 ;i--){
	if (minIndex == 0) {
            if (Mat_energyNew.at<float>(i, 0) > Mat_energyNew.at<float>(i, 1)) {
                minIndex = 1;
            }
        }
	else if(minIndex == cols-1){
	    if(Mat_energyNew.at<float>(i, cols-2) < Mat_energyNew.at<float>(i, cols-1)){
		minIndex = cols-2;
	    }
	}
	else{
	    if(Mat_energyNew.at<float>(i, minIndex-1) <= Mat_energyNew.at<float>(i, minIndex) && Mat_energyNew.at<float>(i, minIndex-1) <= Mat_energyNew.at<float>(i, minIndex+1)){
		minIndex = minIndex - 1;
	    }
	    else if(Mat_energyNew.at<float>(i, minIndex) > Mat_energyNew.at<float>(i, minIndex+1)){
		minIndex = minIndex + 1;
	    }
	}
	for(int j = 0;j < cols;j++){
	    if (j < minIndex) {
                out_image.at<Vec3b>(i, j) = in_image.at<Vec3b>(i, j);
            } 
	    else if (j > minIndex) {
                out_image.at<Vec3b>(i, j-1) = in_image.at<Vec3b>(i, j);
                Mat_gray.at<float>(i, j-1) = Mat_gray.at<float>(i, j);
            }
	}
    }
    */
    
    return true;

}

//optimal seam


