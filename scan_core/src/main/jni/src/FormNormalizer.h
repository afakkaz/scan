#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include <iostream>
#include <fstream>
#include <opencv/cv.h>
#include <json/json.h>
using namespace std;
using namespace cv;
class FormNormalizer
{
public:
	int threshold_value = 0;
	int threshold_type = 1;
	int const max_value = 255;
	int const max_type = 4;
	int const max_BINARY_value = 255;

	Mat src, srcCopy, cropCopy, src_gray, dst;
	char* window_name = "Threshold Demo";

	char* trackbar_type = "Type: \n 0: Binary \n 1: Binary Inverted \n 2: Truncate \n 3: To Zero \n 4: To Zero Inverted";
	char* trackbar_value = "Value";

	float height;
	float width;
	float xfactor;
	float yfactor;
	std::vector<Rect> InitialSegments;

	//This method is used to apply filters on source image
	Mat ApplyFilters(Mat bgr_image)
	{

		//CLAHE algorithm
		cv::Mat lab_image;

		cv::cvtColor(bgr_image, lab_image, CV_BGR2Lab);

		// Extract the L channel
		std::vector<cv::Mat> lab_planes(3);
		cv::split(lab_image, lab_planes);  // now we have the L image in lab_planes[0]

		// apply the CLAHE algorithm to the L channel
		cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
		clahe->setClipLimit(3);
		cv::Mat dst;
		clahe->apply(lab_planes[0], dst);

		// Merge the the color planes back into an Lab image
		dst.copyTo(lab_planes[0]);
		cv::merge(lab_planes, lab_image);

		// convert back to RGB
		static cv::Mat image_clahe;
		cv::cvtColor(lab_image, image_clahe, CV_Lab2BGR);
		///////////////////////////////////////////////////////////////////////////////////////////
		cvtColor(image_clahe, src_gray, CV_RGB2GRAY);
		Mat gray = src_gray.clone();
		threshold(gray, dst, 110, max_BINARY_value, threshold_type);

		int erosion_size = 2;
		/// Apply the erosion operation
		Mat element = getStructuringElement(MORPH_ELLIPSE,
			Size(2 * erosion_size + 1, 2 * erosion_size + 1),
			Point(erosion_size, erosion_size));
		//erode(dst, dst, element);
		return dst;
	}
	//Finding all the blobs 
	std::vector<Rect> FindBloobs(Mat imag)
	{
		std::vector<Rect> rects_list;
		std::vector<std::vector<Point>> contours;
		std::vector<Vec4i> hierarchy;
		findContours(imag, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
		for (size_t i = 0; i < contours.size(); i++)
		{
			cv::Point p(boundingRect(contours[i]).x, boundingRect(contours[i]).y);
			rects_list.push_back(boundingRect(contours[i]));
		}
		return rects_list;
	}

	bool valueInRange(int value, int min, int max)
	{
		return (value >= min) && (value <= max);
	}
	//checking to see if two rectangles are overlaping and of almost same in one dimension
	bool AreOverlapingSameSize(Rect A, Rect B)
	{
		int BlobHeight;
		bool IsSameSize;
		if (A.height > A.width){
			BlobHeight = A.width;
			IsSameSize = B.width > (BlobHeight)* 0.1 && B.width < (BlobHeight*1.1);
		}
		else{
			BlobHeight = A.height;
			IsSameSize = B.height >(BlobHeight) * 0.4 && B.height < (BlobHeight*1.1);
		}
		//Checking to see the both rectangles have almose same on dimension

		bool xOverlap = valueInRange(A.x, B.x, B.x + B.width) ||
			valueInRange(B.x, A.x, A.x + A.width);

		bool yOverlap = valueInRange(A.y, B.y, B.y + B.height) ||
			valueInRange(B.y, A.y, A.y + A.height);

		return xOverlap && yOverlap && IsSameSize;
	}
	//Checking from only ovelapping
	bool AreOverlaping(Rect A, Rect B)
	{
		bool xOverlap = valueInRange(A.x, B.x, B.x + B.width) ||
			valueInRange(B.x, A.x, A.x + A.width);

		bool yOverlap = valueInRange(A.y, B.y, B.y + B.height) ||
			valueInRange(B.y, A.y, A.y + A.height);

		return xOverlap && yOverlap;
	}
	//checking to see if a segment needs repositioning
	bool NeedRepositioning(Rect segment, std::vector<Rect> overlapingBlobs)
	{
		int xmin = 10000, xmax = -1;
		for (int i = 1; i < overlapingBlobs.size(); i++)
		{
			if (xmin > overlapingBlobs[i].x)
				xmin = overlapingBlobs[i].x;
			if (xmax < overlapingBlobs[i].x)
				xmax = overlapingBlobs[i].x;
		}

		if (xmax - xmin > segment.width / 2)
			return true;
		else
			false;
	}
	//GETTING all the blobs containing itmes in it. 
	std::vector<Rect> FindElements(std::vector<Rect> rects_listAll)
	{
		//Find All
		std::vector<Rect> rects_listElements;
		for (int i = 0; i < rects_listAll.size(); i++)
		{
			int min = srcCopy.cols * 0.006127;
			int max = srcCopy.cols * 0.00919117;
			if ((rects_listAll[i].width > min && rects_listAll[i].height > min)
				&& rects_listAll[i].width < src.size().width / max
				&& rects_listAll[i].height < src.size().height / max
				&& (rects_listAll[i].width < rects_listAll[i].height*1.5))
			{
				rects_listElements.push_back(Rect(rects_listAll[i].x, rects_listAll[i].y, rects_listAll[i].width, rects_listAll[i].height));
			}
		}

		//Removing Overlaping
		for (int i = 0; i < rects_listElements.size(); i++)
		{
			for (int j = i + 1; j < rects_listElements.size(); j++)
			{
				if (AreOverlaping(rects_listElements[i], rects_listElements[j]))
				{
					int size1 = rects_listElements[i].width * rects_listElements[i].height;
					int size2 = rects_listElements[j].width * rects_listElements[j].height;
					if (size1 > size2)
					{
						rects_listElements.erase(rects_listElements.begin() + j);
						j--;
					}
					else
					{
						rects_listElements.erase(rects_listElements.begin() + i);
						i--;
						break;
					}
				}
			}
		}

		return rects_listElements;
	}
	//Finding from bounding Rectangle
	Rect FindForm(std::vector<Rect> rects_list)
	{
		//Finding the boudary
		int xmin = 9000, xmax = 0, ymin = 9000, ymax = 0;
		for (int i = 0; i < rects_list.size(); i++)
		{
			if ((rects_list[i].width > 200 && rects_list[i].height > 200) && rects_list[i].width < src.size().width / 1.1 && rects_list[i].height < src.size().height / 1.1)
			{
				if (rects_list[i].x < xmin)
					xmin = rects_list[i].x;
				if (rects_list[i].x + rects_list[i].width > xmax)
					xmax = rects_list[i].x + rects_list[i].width;

				if (rects_list[i].y < ymin)
					ymin = rects_list[i].y;
				if (rects_list[i].y + rects_list[i].height > ymax)
					ymax = rects_list[i].y + rects_list[i].height;
			}
		}
		int lineCount = 0;
		bool isblack = false;
		// Traversing top to chack if top boundary is correct if nor reset it
		int y;
		for (y = ymin; y > 0; y--)
		{
			for (int x = xmin; x < xmax; x++)
			{
				Vec3b colour = src_gray.at<uchar>(Point(x, y));
				if (colour.val[0] < 50)
				{
					isblack = true;
					break;
				}
			}
			if (lineCount > 10)
				break;
			else
			if (isblack)
			{
				lineCount = 0;
				isblack = false;
			}
			else
				lineCount++;
		}
		if (lineCount > 10)
			ymin = y;
		// Traversing bottom to chack if bottom boundary is correct if nor reset it
		lineCount = 0;
		for (y = ymax; y < src_gray.size().height; y++)
		{
			for (int x = xmin; x < xmax; x++)
			{
				Vec3b colour = src_gray.at<uchar>(Point(x, y));
				if (colour.val[0] < 50)
				{
					isblack = true;
					break;
				}
			}
			if (lineCount > 10)
				break;
			else
			if (isblack)
			{
				lineCount = 0;
				isblack = false;
			}
			else
				lineCount++;
		}
		if (lineCount > 10)
			ymax = y;

		return Rect(xmin, ymin - 30, xmax - xmin, ymax - ymin);
	}
	//For sorting rectangle using position
	struct byposition {
		bool operator () (const std::vector<Rect> & a, const std::vector<Rect> & b) {
			return ((a[0].y * 1000) + a[0].x) < ((b[0].y * 1000) + b[0].x);
		}
	};
	//for sorting rectangle using x position
	struct byXposition {
		bool operator () (const Rect & a, const Rect & b) {
			return (a.x) < (b.x);
		}
	};
	//for sorting rectangle using y position
	struct byYposition {
		bool operator () (const Rect & a, const Rect & b) {
			return (a.y) < (b.y);
		}
	};
	//for sorting rectangle using area
	struct byHeight {
		bool operator () (const Rect & a, const Rect & b) {
			return (a.height) > (b.height);
		}
	};
	//Used to reagesegmentans item information from json
	std::vector<std::vector<Rect>> readSegments(Rect FormRect, string path){
		std::vector<std::vector<Rect>> FinalSegments;

		Json::Value templateRoot;
		Json::Reader reader;
		//getting jason file 
		std::ifstream  file(path);
		//reading json data from file
		if (!reader.parse(file, templateRoot)){
		}

		//getting height and width of template 
		height = templateRoot.get("height", 0).asInt();
		width = templateRoot.get("width", 0).asInt();
		//creating factor of it properly position it to the detected from
		xfactor = (FormRect.width) / width;
		yfactor = (FormRect.height) / height;

		const Json::Value fields = templateRoot["fields"];

		for (size_t i = 0; i < fields.size(); i++) {
			const Json::Value field = fields[i];
			const Json::Value segments = field["segments"];
			for (size_t j = 0; j < segments.size(); j++) {
				const Json::Value segment = segments[j];

				int x = (segment.get("segment_x", INT_MIN).asInt()) * xfactor + FormRect.x;
				int y = (segment.get("segment_y", INT_MIN).asInt()) * yfactor + FormRect.y;
				int width = (segment.get("segment_width", INT_MIN).asInt()) * xfactor;
				int height = (segment.get("segment_height", INT_MIN).asInt()) * yfactor;
				if (segment.isMember("align_segment"))
				{
					std::vector <Rect> segItems;
					segItems.push_back(cv::Rect(x, y, width, height));
					InitialSegments.push_back(cv::Rect(x, y, width, height));
					if (segment.isMember("items")){
						const Json::Value items = segment["items"];

						for (size_t K = 0; K < items.size(); K++) {
							const Json::Value item = items[K];

							int xx = (item.get("item_x", INT_MIN).asInt());
							int yy = (item.get("item_y", INT_MIN).asInt());
							segItems.push_back(Rect(xx - 1, yy - 1, 2, 2));
						}

					}
					FinalSegments.push_back(segItems);
				}
			}
		}
		return FinalSegments;
	}
	// IT IS UNDER DEVELOPMENT. Used to crop and same repositioned segments and items.
	void cropSaveSegments(std::vector<Rect> SegmentItem, std::vector<Rect> overlapingRects, int index){
		if (overlapingRects.size() >= SegmentItem.size() - 1){
			//Mat seg(SegmentItem[0].height, SegmentItem[0].width, CV_8U, Scalar(255, 255, 255));

			Mat seg = Mat(SegmentItem[0].height + 20, SegmentItem[0].width + 20, CV_8UC1, Scalar(255));
			//Mat seg = Mat(srcCopy, Rect(0, 0, SegmentItem[0].width*3, SegmentItem[0].height*3)); = 
			int startx = InitialSegments[index].x;
			int starty = InitialSegments[index].y;

			std::vector<Rect> overlapingRects1;

			std::sort(overlapingRects.begin(), overlapingRects.end(), byHeight());
			for (int i = 0; i < SegmentItem.size() - 1; i++)
				overlapingRects1.push_back(overlapingRects[i]);

			//cout << overlapingRects.size();
			//cout << "\n";

			int FirstW = overlapingRects[0].width, FirstH = overlapingRects[0].height;

			if (SegmentItem[0].width > SegmentItem[0].height)
				std::sort(overlapingRects1.begin(), overlapingRects1.end(), byXposition());
			else
				std::sort(overlapingRects1.begin(), overlapingRects1.end(), byYposition());

			for (int i = 1; i < SegmentItem.size(); i++){
				//cv::rectangle(seg, Rect(SegmentItem[i].x*xfactor, SegmentItem[i].y*yfactor, SegmentItem[i].width*xfactor, SegmentItem[i].height * yfactor), Scalar(255, 255, 0), 16, 4, 0);

				int dx = FirstW - overlapingRects1[i - 1].width;
				int dy = FirstH - overlapingRects1[i - 1].height;
				int ox, oy, ow, oh;
				ox = overlapingRects1[i - 1].x - dx / 2;
				if (ox < 0)
					ox = overlapingRects1[i - 1].x;
				oy = overlapingRects1[i - 1].y - dy / 2;
				if (oy < 0)
					overlapingRects1[i - 1].y;
				ow = overlapingRects1[i - 1].width + dx;
				if (ox + ow > cropCopy.cols)
					ow = overlapingRects1[i - 1].width;
				oh = overlapingRects1[i - 1].height + dy;
				if (oy + oh > cropCopy.rows)
					oh = overlapingRects1[i - 1].height;

				Rect roiR = Rect(ox, oy, ow, oh);
				Mat roi = Mat(cropCopy, roiR);
				cvtColor(roi, roi, CV_RGB2GRAY);

				int centerX = SegmentItem[i].x *xfactor;
				int centerY = SegmentItem[i].y *yfactor;
				int x = centerX - roiR.width / 2;
				int y = centerY - roiR.height / 2;

				/*	if (x < 0)
				x = 0;
				if (y < 0)
				y = 0;*/

				for (int k = 0; k<roi.rows; k++){
					for (int j = 0; j<roi.cols; j++){
						//Vec3b color = roi.at<Vec3b>(Point(j,k));
						//srcCopy.at<Vec3b>(Point(starty + j + x, startx + k + y)) = color;
						int cx = startx + j + x;
						int cy = starty + k + y;

						if (cx >= 0 && cx < srcCopy.cols &&
							cy >= 0 && cy < srcCopy.rows){
							srcCopy.at<Point3_<char> >(cy, cx).x = roi.at<char>(k, j);
							srcCopy.at<Point3_<char> >(cy, cx).y = roi.at<char>(k, j);
							srcCopy.at<Point3_<char> >(cy, cx).z = roi.at<char>(k, j);
						}
					}

				}
				//cv::rectangle(srcCopy, Rect(startx, starty, 4, 4), Scalar(0, 0, 0), 4, 4, 0);
				//break;
			}
			cout << "/////\n";
		}
	}
	//Used to reposition the segementss
	std::vector<std::vector<Rect>> repositionSegments(std::vector<std::vector<Rect>> FinalSegments, std::vector<Rect> rects_listElements)
	{

		for (int i = 0; i < FinalSegments.size(); i++)
		{

			//copying segment rectangle
			Rect segment = FinalSegments[i][0];
			int blobW = 0;
			int SizeIncs = 0;
			std::vector<Rect> overlapingRects;

			segment.x += 3;
			segment.y += 3;
			segment.width -= 6;
			segment.height -= 6;

			do{
				//increasing its size for better results

				blobW = 0;
				SizeIncs++;

				cv::rectangle(src, segment, Scalar(255, 255, 0), 2, 4, 0);
				//Finding all the overlaping item blobs with the segment
				for (int j = 0; j < rects_listElements.size(); j++)
				{
					if (i == 2)
						cout << blobW;
					if (AreOverlapingSameSize(segment, rects_listElements[j]))
					{
						//cv::rectangle(src, rects_listElements[j], Scalar(0, 0, 0), 6, 4, 0);
						overlapingRects.push_back(Rect(rects_listElements[j].x, rects_listElements[j].y,
							rects_listElements[j].width, rects_listElements[j].height));


						//cv::rectangle(src, overlapingRects[overlapingRects.size()-1], Scalar(255, 255, 0), 6, 4, 0);
						blobW++;
					}
				}
				segment.x -= 3;
				segment.y -= 3;
				segment.width += 6;
				segment.height += 6;
			} while (blobW < FinalSegments[i].size() - 1 && SizeIncs < 10);


			if (overlapingRects.size() > 0)//&& NeedRepositioning(FinalSegments[i][0], overlapingRects))
			{
				//sorting the overlaping items with resoect to the segment's dimensions
				if (segment.width > segment.height)
					std::sort(overlapingRects.begin(), overlapingRects.end(), byXposition());
				else
					std::sort(overlapingRects.begin(), overlapingRects.end(), byYposition());

				Rect LeftTop = Rect(0, 0, 0, 0);
				std::vector<Rect> OverlapingRegion;
				//getting first items overlapping rectangle with segment
				for (int j = 0; j < overlapingRects.size(); j++)
				{
					int x = max(overlapingRects[j].x, segment.x);
					int y = max(overlapingRects[j].y, segment.y);
					int width = min(overlapingRects[j].x + overlapingRects[j].width, segment.x + segment.width) - x;
					int height = min(overlapingRects[j].y + overlapingRects[j].height, segment.y + segment.height) - y;
					if (x > 0 && y > 0 && width > 0 && height > 0)
						OverlapingRegion.push_back(Rect(x, y, width, height));
					else{
						overlapingRects.erase(overlapingRects.begin() + j);
						j--;
					}
				}

				for (int l = 0; l <1; l++){
					if (i == 17)
						cv::rectangle(src, OverlapingRegion[l], Scalar(0, 0, 0), 2, 4, 0);
				}

				/*putText(src, std::to_string(i), cvPoint(segment.x, segment.y),
				FONT_HERSHEY_COMPLEX_SMALL, 5, cvScalar(0, 0, 0), 1, CV_AA);*/

				//FInding the first item to the segment for segment repositioning
				if (overlapingRects.size() > 1)
				{

					if (segment.width > segment.height){

						if (overlapingRects[0].width < overlapingRects[0].height / 1.1 || overlapingRects[1].width < overlapingRects[1].height / 1.1){

							if (OverlapingRegion[0].height >(OverlapingRegion[1].height / 1.2)){
								LeftTop = overlapingRects[0];
							}
							else{
								LeftTop = overlapingRects[1];
							}
						}
						else{

							if (OverlapingRegion[0].width >(OverlapingRegion[1].width / 1.2)){
								LeftTop = overlapingRects[0];
							}
							else{
								LeftTop = overlapingRects[1];

							}
						}
					}
					else{
						//Resolving issue when checkboxes are either height wise or width wise check need
						if (OverlapingRegion[0].width > OverlapingRegion[0].height / 1.1 || OverlapingRegion[1].width > OverlapingRegion[1].height / 1.1){
							if (OverlapingRegion[0].width > (OverlapingRegion[1].width / 1.2))
								LeftTop = overlapingRects[0];
							else{
								LeftTop = overlapingRects[1];
								if (i == 2)
									cv::rectangle(src, LeftTop, Scalar(0, 0, 0), 16, 4, 0);
							}
						}
						else{
							if (OverlapingRegion[0].height > (OverlapingRegion[1].height / 1.2))
								LeftTop = overlapingRects[0];
							else
								LeftTop = overlapingRects[1];
						}
					}
				}
				else
					LeftTop = overlapingRects[0];

				//SHowing the first item found for the segment 
				cv::rectangle(src, LeftTop, Scalar(255, 255, 0), 2, 4, 0);

				//finding the difference while repositioning
				int dx = LeftTop.x - segment.x;
				int dy = LeftTop.y - segment.y;

				//Reposioning the segment
				//if (dx < abs(min(FinalSegments[i][0].height, FinalSegments[i][0].width)) / 1.5)
				FinalSegments[i][0].x = LeftTop.x;
				//if (dy < abs(min(FinalSegments[i][0].height, FinalSegments[i][0].width)) / 2)
				FinalSegments[i][0].y = LeftTop.y;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Finding all the overlaping item blobs with the segment
				OverlapingRegion.clear();
				overlapingRects.clear();

				segment = FinalSegments[i][0];
				do{
					//increasing its size for better results

					blobW = 0;
					SizeIncs++;

					cv::rectangle(src, segment, Scalar(255, 255, 0), 2, 4, 0);
					//Finding all the overlaping item blobs with the segment
					for (int j = 0; j < rects_listElements.size(); j++)
					{
						if (i == 2)
							cout << blobW;
						if (AreOverlapingSameSize(segment, rects_listElements[j]))
						{
							//cv::rectangle(src, rects_listElements[j], Scalar(0, 0, 0), 6, 4, 0);
							overlapingRects.push_back(Rect(rects_listElements[j].x, rects_listElements[j].y,
								rects_listElements[j].width, rects_listElements[j].height));


							//cv::rectangle(src, overlapingRects[overlapingRects.size()-1], Scalar(255, 255, 0), 6, 4, 0);
							blobW++;
						}
					}
					segment.x -= 3;
					segment.y -= 3;
					segment.width += 6;
					segment.height += 6;
				} while (blobW < FinalSegments[i].size() - 1 && SizeIncs < 7);

				//getting first items overlapping rectangle with segment
				for (int j = 0; j < overlapingRects.size(); j++)
				{
					int x = max(overlapingRects[j].x, segment.x);
					int y = max(overlapingRects[j].y, segment.y);
					int width = min(overlapingRects[j].x + overlapingRects[j].width, segment.x + segment.width) - x;
					int height = min(overlapingRects[j].y + overlapingRects[j].height, segment.y + segment.height) - y;
					if (x > 0 && y > 0 && width > 0 && height > 0)
						OverlapingRegion.push_back(Rect(x, y, width, height));
					else{
						overlapingRects.erase(overlapingRects.begin() + j);
						j--;
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////


				//Sorting based upone overlapingregion
				std::vector<Rect> regionindex;
				std::vector<Rect> FinalOverlaping;
				for (int j = 0; j < OverlapingRegion.size(); j++){
					regionindex.push_back(Rect(j, 0, 0, OverlapingRegion[j].width*OverlapingRegion[j].height));
					//intializing all indexes by 0,0,0,0
					FinalOverlaping.push_back(Rect(0, 0, 0, 0));
				}
				std::sort(regionindex.begin(), regionindex.end(), byHeight());

				for (int j = 0; j < regionindex.size(); j++){
					//cout << regionindex[j].x;
					//cout << "\n";
					FinalOverlaping[j] = Rect(overlapingRects[regionindex[j].x].x, overlapingRects[regionindex[j].x].y,
						overlapingRects[regionindex[j].x].width, overlapingRects[regionindex[j].x].height);
				}
				cout << FinalOverlaping.size();
				cout << "-";
				cout << overlapingRects.size();
				cout << "\n";


				if (FinalOverlaping.size() > 1){
					int fx = FinalOverlaping[0].x, fy = FinalOverlaping[0].y;
					for (int l = FinalSegments[i].size() - 1; l < FinalOverlaping.size(); l++)
					{
						/*int fxx = FinalOverlaping[l].x;
						int fyy = FinalOverlaping[l].y;
						int w = FinalOverlaping[l].width;
						int h = FinalOverlaping[l].height;
						if (fxx < fx)
						{*/
						//FinalOverlaping.erase(FinalOverlaping.begin() + l);
						//l--;
						/*}
						else
						fx = fxx;
						fy = fyy;*/

					}
				}

				//Calling to crope and save the current results (Under development)
				if (i == 22)
				for (int l = 0; l < FinalOverlaping.size(); l++)
					cv::rectangle(src, FinalOverlaping[l], Scalar(0, 0, 0), 2, 4, 0);
				/*putText(src, std::to_string(i), cvPoint(segment.x, segment.y),
				FONT_HERSHEY_COMPLEX_SMALL, 5, cvScalar(0, 0, 0), 1, CV_AA);*/

				cropSaveSegments(FinalSegments[i], FinalOverlaping, i);

			}
		}
		imwrite("Segements/FullphotoElements.jpg", srcCopy);
		return FinalSegments;
	}
	Mat PreprocessForm(Mat scr, string TemplatePath)
	{
		srcCopy = src.clone();
		cropCopy = src.clone();
		//converting it to grayscale
		//cvtColor(src, src_gray, CV_RGB2GRAY);

		//namedWindow("Display threshod", WINDOW_NORMAL);// Create a window for display.
		//namedWindow("Display Gray", WINDOW_NORMAL);// Create a window for display.

		//Applying some filter to exrtract blob information
		dst = ApplyFilters(src);

		//dst = dst.inv();
		//showing the grayscale image
		//imshow("Display threshod", dst);

		//Getting all the blobs
		std::vector<Rect> rects_listAll = FindBloobs(dst);

		//Filtring blobs which contains segment itemsin it.
		std::vector<Rect> rects_listElements = FindElements(rects_listAll);

		//Finding the boundary og form later on  htis will be done using your code
		Rect FormRect = Rect(0,0,srcCopy.cols,srcCopy.rows);// FindForm(rects_listAll);

		//Drawing full form bounary in green
		cv::rectangle(src, FormRect, Scalar(0, 255, 0), 6, 4, 0);

		//Drawing elements in Red
		for (int i = 0; i < rects_listElements.size(); i++)
			cv::rectangle(src, rects_listElements[i], Scalar(0, 0, 255), 2, 4, 0);

		//Getting all the segments and items 
		std::vector<std::vector<Rect>> FinalSegments = readSegments(FormRect, TemplatePath);
		//std::vector<std::vector<Rect>> AllItems;

		//sorting all the degments relative to their positions
		//std::sort(FinalSegments.begin(), FinalSegments.end(), byposition());

		//Drawing all the segments in blue on src image
		for (int i = 0; i < FinalSegments.size(); i++)
			cv::rectangle(src, FinalSegments[i][0], Scalar(255, 0, 0), 4, 4, 0);


		//Repositioning Segments to start at the first item they contain and also extracting repositioned items form it
		//Thismethod is under development
		FinalSegments = repositionSegments(FinalSegments, rects_listElements);

		//SHowing all the repositioned segemnts and Items dots.
		for (int i = 0; i < FinalSegments.size(); i++){
			cv::rectangle(src, FinalSegments[i][0], Scalar(0, 255, 0), 6, 4, 0);
			for (int j = 1; j < FinalSegments[i].size(); j++){
				Rect item = Rect(FinalSegments[i][j].x * xfactor + FinalSegments[i][0].x, FinalSegments[i][j].y * yfactor + FinalSegments[i][0].y, 2, 2);
				cv::rectangle(src, item, Scalar(0, 255, 0), 8, 4, 0);
			}
		}

		//imshow("Display Gray", src);
		//imwrite("Segements/FullphotoElements.jpg", srcCopy);
		return srcCopy;
	}
	FormNormalizer::FormNormalizer()
	{
	}
	FormNormalizer::~FormNormalizer()
	{
	}
	//FormNormalizer();
	//Mat ApplyFilters(Mat gray);
	////Finding all the blobs 
	//std::vector<Rect> FindBloobs(Mat imag);

	//bool valueInRange(int value, int min, int max);
	////checking to see if two rectangles are overlaping and of almost same in one dimension
	//bool AreOverlapingSameSize(Rect A, Rect B);
	////Checking from only ovelapping
	//bool AreOverlaping(Rect A, Rect B);
	////checking to see if a segment needs repositioning
	//bool NeedRepositioning(Rect segment, std::vector<Rect> overlapingBlobs);
	////GETTING all the blobs containing itmes in it. 
	//std::vector<Rect> FindElements(std::vector<Rect> rects_listAll);
	////Finding from bounding Rectangle
	//Rect FindForm(std::vector<Rect> rects_list);
	////Used to reagesegmentans item information from json
	//std::vector<std::vector<Rect>> readSegments(Rect FormRect, string path);
	//// IT IS UNDER DEVELOPMENT. Used to crop and same repositioned segments and items.
	//void cropSaveSegments(std::vector<Rect> SegmentItem, std::vector<Rect> overlapingRects, int index);
	////Used to reposition the segementss
	//std::vector<std::vector<Rect>> repositionSegments(std::vector<std::vector<Rect>> FinalSegments, std::vector<Rect> rects_listElements);
	//
	//Mat PreprocessForm(string FormPath, string TemplatePath);
	//~FormNormalizer();
};

